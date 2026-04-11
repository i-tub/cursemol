"""
Functions for editing a molecule, given the cursor position and other inputs
from the UI.

This module does not use curses.
"""

import logging
import math

from rdkit import Chem
from rdkit.Chem import AllChem

from . import canvas
from . import chem
from . import state as state_


def create_or_adjust_bond(state,
                          cursor_x,
                          cursor_y,
                          screen_dims,
                          bond_order,
                          bond_dir=None):
    """
    Create or adjust bond between two nearest atoms at cursor position.
    bond_order: 1 (single), 2 (double), or 3 (triple)
    bond_dir: Optional bond direction (e.g., Chem.BondDir.BEGINWEDGE)
    Returns True if bond was created/modified, False otherwise.
    """
    # Find the two atoms that should be bonded (using screen coordinates)
    atom_pair = canvas.find_bond_atoms(state, cursor_x, cursor_y, screen_dims)

    if atom_pair is not None:
        atom1_idx, atom2_idx = atom_pair
        return chem.modify_bond(state.mol, atom1_idx, atom2_idx, bond_order,
                                bond_dir)

    return False


def adjust_formal_charge(state, cursor_x, cursor_y, screen_dims, delta):
    """
    Adjust formal charge of atom at cursor position by delta.
    Returns True if charge was adjusted, False if no change or no atom found.
    """
    if delta == 0:
        return False  # No change requested

    atom_idx = canvas.find_atom_at_cursor(state, cursor_x, cursor_y,
                                          screen_dims)
    if atom_idx is not None:
        atom = state.mol.GetAtomWithIdx(atom_idx)
        current_charge = atom.GetFormalCharge()
        atom.SetFormalCharge(current_charge + delta)
        return True

    return False


def delete_atoms_in_rect(state, x1, y1, x2, y2, screen_dims):
    """
    Delete all atoms whose screen positions fall within the rectangle.
    Returns True if any atoms were deleted, False otherwise.
    """
    # Normalize rectangle coordinates
    min_x, min_y, max_x, max_y_rect = canvas.normalize_rect(x1, y1, x2, y2)

    # Find atoms within the rectangle
    atoms_to_delete = []
    for atom, screen_x, screen_y in canvas.iter_atom_screen_positions(
            state, screen_dims):
        if min_x <= screen_x <= max_x and min_y <= screen_y <= max_y_rect:
            atoms_to_delete.append(atom.GetIdx())

    # Delete atoms in reverse order to avoid index issues
    if atoms_to_delete:
        for atom_idx in sorted(atoms_to_delete, reverse=True):
            try:
                state.mol.RemoveAtom(atom_idx)
            except Exception:
                logging.exception(f"Error removing atom {atom_idx}")
        return True

    return False


def delete_at_cursor(state, cursor_x, cursor_y, screen_dims):
    """
    Delete atom or bond at cursor position.
    Returns True if something was deleted, False otherwise.
    """
    # First try to find an atom at cursor
    atom_idx = canvas.find_atom_at_cursor(state, cursor_x, cursor_y,
                                          screen_dims)
    if atom_idx is not None:
        try:
            state.mol.RemoveAtom(atom_idx)
            return True
        except Exception:
            logging.exception("Error removing atom (x command)")
            return False
    else:
        # No atom found, try to delete a bond instead
        atom_pair = canvas.find_bond_atoms(state, cursor_x, cursor_y,
                                           screen_dims)
        if atom_pair is not None:
            atom1_idx, atom2_idx = atom_pair
            if chem.modify_bond(state.mol, atom1_idx, atom2_idx, 0):
                return True

    return False


def delete_fragment_at_cursor(state, cursor_x, cursor_y, screen_dims):
    """
    Delete all atoms reachable from the atom at cursor position via BFS.
    Returns True if any atoms were deleted, False otherwise.
    """
    # Find starting atom at cursor
    start_atom_idx = canvas.find_atom_at_cursor(state, cursor_x, cursor_y,
                                                screen_dims)
    if start_atom_idx is None:
        return False

    # BFS to find all reachable atoms
    visited = set()
    queue = [start_atom_idx]
    visited.add(start_atom_idx)

    while queue:
        atom_idx = queue.pop(0)
        atom = state.mol.GetAtomWithIdx(atom_idx)

        # Add all neighbors to queue
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                visited.add(neighbor_idx)
                queue.append(neighbor_idx)

    # Delete atoms in reverse order (highest index first) to avoid index shifting
    if visited:
        try:
            for atom_idx in sorted(visited, reverse=True):
                state.mol.RemoveAtom(atom_idx)
            return True
        except Exception:
            logging.exception("Error deleting fragment (D command)")
            return False

    return False


def connect_sidechain_to_bond(state, bond_atom_pair, start_idx, end_idx,
                              cursor_x, cursor_y, screen_dims):
    """
    Connect a sidechain to a bond by forming a ring.

    The sidechain is inserted between the two bond atoms: the closer atom
    connects to the first sidechain atom, and the farther atom connects
    to the last sidechain atom.

    Args:
        state: Current molecular state
        bond_atom_pair: (atom1_idx, atom2_idx) tuple
        start_idx: Index of first atom in sidechain
        end_idx: Index of last atom in sidechain
        cursor_x, cursor_y: Cursor position (to determine which atom is closer)
        screen_dims: Screen dimensions
    """
    a1_idx, a2_idx = bond_atom_pair

    # Determine which atom is closer to cursor
    conf = state.mol.GetConformer()
    mol_x, mol_y = canvas.screen_to_mol_coords(cursor_x, cursor_y, state.box,
                                               state.scale, screen_dims,
                                               state.y_offset)

    pos1 = conf.GetAtomPosition(a1_idx)
    pos2 = conf.GetAtomPosition(a2_idx)

    dist1 = math.sqrt((pos1.x - mol_x)**2 + (pos1.y - mol_y)**2)
    dist2 = math.sqrt((pos2.x - mol_x)**2 + (pos2.y - mol_y)**2)

    # Order atoms so a1 is closer to cursor
    if dist1 > dist2:
        a1_idx, a2_idx = a2_idx, a1_idx

    # Connect sidechain: a1 (closer) -> first atom, a2 (farther) -> last atom
    state.mol.AddBond(a1_idx, start_idx, Chem.BondType.SINGLE)
    state.mol.AddBond(a2_idx, end_idx, Chem.BondType.SINGLE)


def insert_or_modify_atom(state, cursor_x, cursor_y, screen_dims,
                          element_symbol):
    """
    Handle the 'i' command: insert atom at cursor or modify existing atom.
    If element_symbol is provided, use it; otherwise prompt the user.
    Returns mol if successful, None if no change was made.
    """
    # Check if cursor is on an atom
    atom_idx = canvas.find_atom_at_cursor(state, cursor_x, cursor_y,
                                          screen_dims)

    if element_symbol:
        try:
            if atom_idx is not None:
                # Change existing atom's symbol
                atom = state.mol.GetAtomWithIdx(atom_idx)
                # Check if atom already has this symbol
                if atom.GetSymbol() == element_symbol:
                    return None  # No change needed
                atom.SetAtomicNum(
                    Chem.GetPeriodicTable().GetAtomicNumber(element_symbol))
            else:
                # Add new atom to molecule
                atom_idx = state.mol.AddAtom(Chem.Atom(element_symbol))

                # Convert cursor position to molecule coordinates
                mol_x, mol_y = canvas.screen_to_mol_coords(
                    cursor_x, cursor_y, state.box, state.scale, screen_dims,
                    state.y_offset)

                # Set atom position in conformer
                conf = state.mol.GetConformer()
                conf.SetAtomPosition(atom_idx, [mol_x, mol_y, 0.0])

            return state.mol
        except Exception:
            logging.exception("Error inserting/modifying atom")
            pass

    return None


def cleanup_coordinates(state, screen_dims):
    """
    Regenerate 2D coordinates for the molecule and recenter the view.
    Keeps the current zoom level. Returns True if successful.
    """
    mol = state.mol
    try:
        chem.assign_stereo(mol)

        # Clear bond directions because we'll have to recompute them
        # for the new coordinates.
        for bond in mol.GetBonds():
            bond.SetBondDir(Chem.BondDir.NONE)
        AllChem.Compute2DCoords(mol)
        Chem.WedgeMolBonds(mol, mol.GetConformer())

        if state.scale is not None:
            state.box, state.y_offset = state_.recalculate_box_and_offset(
                mol, state.scale, screen_dims)
        return True
    except Exception:
        logging.exception("Error regenerating coordinates")
        return False
