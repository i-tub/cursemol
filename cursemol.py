#!/usr/bin/env python3
"""
CurseMol - Molecular sketcher for the terminal

Controls:
  h, j, k, l - Move cursor left, down, up, right
  H, J, K, L - Shift molecule left, down, up, down, right
  s          - Enter a SMILES string to replace the current molecule
  S          - Toggle SMILES display
  i          - Insert/modify atom at cursor position
  a          - Append atoms from SMILES to atom or bond under cursor
               (appending to a bond forms a ring by connecting the bond atoms
               to the first and last atoms from the SMILES)
  c, n, o    - Insert/modify carbon/nitrogen/oxygen atom
  x          - Delete atom or bond
  +, -       - Increase/decrease formal charge on atom
  <, >       - Zoom out/in
  1, 2, 3    - Add bond or change bond (order 1/2/3) between nearest atoms
  @          - Clear canvas (reset to blank slate)
  u          - Undo
  r          - Redo
  Ctrl-L     - Clean up (regenerate coordinates)
  ?          - Show this help
  q          - Quit
"""

import argparse
import atexit
import curses
import logging
import math
import os
import sys
from rdkit import Chem
from rdkit import Geometry
from rdkit import RDLogger
from rdkit.Chem import AllChem

MIN_SCALE = 2.0  # columns per angstrom
DEFAULT_SCALE = 8.0  # columns per angstrom
MAX_SCALE = 16.0  # columns per angstrom
ASPECT_RATIO = 0.4  # horizontal / vertical
PADDING = 5

BOND_CHARS = {
    Chem.BondType.SINGLE: '.',
    Chem.BondType.DOUBLE: '=',
    Chem.BondType.TRIPLE: '#'
}


def get_box(conf):
    xyz = conf.GetPositions()
    return (xyz.min(axis=0), xyz.max(axis=0))


def int_coords_for_atom(atom, box, scale, conf, y_offset=0, rows=0):
    pos = conf.GetAtomPosition(atom.GetIdx())
    x = PADDING + int((pos.x - box[0][0]) * scale[0])
    # Flip y coordinate so molecules appear right-side up (not reversed)
    # Higher molecule y -> lower screen y (closer to top)
    y_from_bottom = PADDING + y_offset + int((pos.y - box[0][1]) * scale[1])
    y = rows - 1 - y_from_bottom
    return x, y


def draw_line(screen, char, x1, y1, x2, y2):
    vertical = False
    if abs(x2 - x1) < abs(y2 - y1):
        x1, y1 = y1, x1
        x2, y2 = y2, x2
        vertical = True
    try:
        slope = 1.0 * (y2 - y1) / (x2 - x1)
    except ZeroDivisionError:
        return

    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
    for x in range(x1 + 1, x2):
        y = int(round(y1 + slope * (x - x1)))
        if vertical:
            screen[x][y] = char
        else:
            screen[y][x] = char


def enter_smiles(stdscr, max_y):
    """Prompt user to enter a SMILES string and return it."""
    # Show prompt at the bottom
    stdscr.addstr(max_y - 1, 0, "Enter SMILES: ")
    stdscr.clrtoeol()
    stdscr.refresh()

    # Enable echoing and get string input
    curses.echo()
    try:
        smiles_bytes = stdscr.getstr(max_y - 1, 14)
        smiles = smiles_bytes.decode('utf-8')
    except Exception:
        logging.exception("Error in enter_smiles")
        smiles = ""
    finally:
        curses.noecho()

    return smiles


def enter_element(stdscr, max_y):
    """Prompt user to enter an element symbol and return it."""
    # Show prompt at the bottom
    stdscr.addstr(max_y - 1, 0, "Element symbol: ")
    stdscr.clrtoeol()
    stdscr.refresh()

    # Enable echoing and get string input
    curses.echo()
    try:
        symbol_bytes = stdscr.getstr(max_y - 1, 16)
        symbol = symbol_bytes.decode('utf-8').strip()
    except Exception:
        logging.exception("Error in enter_element")
        symbol = ""
    finally:
        curses.noecho()

    return symbol


def screen_to_mol_coords(cursor_x, cursor_y, box, scale, max_y, y_offset):
    """Convert cursor/terminal coordinates to molecule coordinates."""
    # Screen array has (max_y - 2) rows
    rows = max_y - 2

    # Terminal position is now directly the screen array position
    # Reverse the coordinate transformation from int_coords_for_atom
    mol_x = (cursor_x - PADDING) / scale[0] + box[0][0]

    # Reverse the y flipping: cursor_y -> y_from_bottom -> mol_y
    y_from_bottom = rows - 1 - cursor_y
    mol_y = (y_from_bottom - PADDING - y_offset) / scale[1] + box[0][1]

    return mol_x, mol_y


def find_atom_at_cursor(mol,
                        cursor_x,
                        cursor_y,
                        box,
                        scale,
                        max_y,
                        y_offset,
                        tolerance=1):
    """
    Find an atom at or near the cursor position (within tolerance cells).
    Returns atom index or None if no atom found.
    """
    if mol.GetNumAtoms() == 0:
        return None

    conf = mol.GetConformer()
    rows = max_y - 2

    # Check atoms to find one whose screen position is within tolerance
    for atom in mol.GetAtoms():
        screen_x, screen_y = int_coords_for_atom(atom, box, scale, conf,
                                                 y_offset, rows)

        # Screen position is now directly the terminal position
        # Check if within tolerance
        if (abs(screen_x - cursor_x) <= tolerance and
                abs(screen_y - cursor_y) <= tolerance):
            return atom.GetIdx()

    return None


def find_bond_atoms_screen(mol, cursor_x, cursor_y, box, scale, max_y, y_offset):
    """
    Find the two atoms closest to cursor position such that
    the cursor is roughly between them (angle >= 160 degrees).
    Uses screen coordinates to avoid floating point precision issues.
    Returns (atom1_idx, atom2_idx) or None if no valid pair found.
    """
    if mol.GetNumAtoms() < 2:
        return None

    conf = mol.GetConformer()
    rows = max_y - 2

    # Calculate screen positions and distances for all atoms
    distances = []
    for atom in mol.GetAtoms():
        screen_x, screen_y = int_coords_for_atom(atom, box, scale, conf, y_offset, rows)
        dx = screen_x - cursor_x
        dy = screen_y - cursor_y
        dist = math.sqrt(dx * dx + dy * dy)
        distances.append((dist, atom.GetIdx(), screen_x, screen_y))

    # Sort by distance
    distances.sort()

    # Find all valid pairs and choose the one with smallest combined distance
    best_pair = None
    best_distance_sum = float('inf')

    for i in range(len(distances)):
        dist1, idx1, x1, y1 = distances[i]

        for j in range(i + 1, len(distances)):
            dist2, idx2, x2, y2 = distances[j]

            # Calculate vectors from cursor position to each atom
            v1_x = x1 - cursor_x
            v1_y = y1 - cursor_y
            v2_x = x2 - cursor_x
            v2_y = y2 - cursor_y

            # Calculate lengths
            len1 = math.sqrt(v1_x * v1_x + v1_y * v1_y)
            len2 = math.sqrt(v2_x * v2_x + v2_y * v2_y)

            if len1 == 0 or len2 == 0:
                continue

            # Calculate angle using dot product
            dot = v1_x * v2_x + v1_y * v2_y
            cos_angle = dot / (len1 * len2)
            angle_deg = math.degrees(math.acos(max(-1.0, min(1.0, cos_angle))))

            # If angle is at least 160 degrees, this is a valid pair
            if angle_deg >= 160:
                # Calculate distance between the two atoms on screen
                atom_dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

                # Score: prefer pairs where cursor is close to both atoms
                # AND the atoms themselves are close to each other
                score = dist1 + dist2 + atom_dist

                if score < best_distance_sum:
                    best_distance_sum = score
                    best_pair = (idx1, idx2)

    return best_pair


def modify_bond(mol, atom1_idx, atom2_idx, bond_order):
    """
    Modify or create a bond between two atoms.
    bond_order: 0 (delete), 1 (single), 2 (double), 3 (triple)
    Returns True if successful, False if no change was made or modification would create an invalid molecule.
    """
    bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)

    # Remember the old state in case we need to revert
    old_bond_type = bond.GetBondType() if bond is not None else None

    # Map bond order to BondType
    bond_type_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }

    try:
        if bond_order == 0:
            # Delete bond if it exists
            if bond is not None:
                mol.RemoveBond(atom1_idx, atom2_idx)
            else:
                # Bond doesn't exist, nothing to delete
                return False
        else:
            bond_type = bond_type_map[bond_order]

            if bond is not None:
                # Check if bond already has this type
                if bond.GetBondType() == bond_type:
                    return False  # No change needed
                # Modify existing bond
                bond.SetBondType(bond_type)
            else:
                # Add new bond
                mol.AddBond(atom1_idx, atom2_idx, bond_type)

        # Test if the molecule can be kekulized
        mol_copy = Chem.RWMol(mol)
        Chem.Kekulize(mol_copy, True)
        return True

    except Exception:
        logging.exception("Error in modify_bond, reverting change")
        # Revert the change if kekulization fails
        if bond_order == 0:
            # We deleted the bond, add it back
            if old_bond_type is not None:
                mol.AddBond(atom1_idx, atom2_idx, old_bond_type)
        else:
            bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
            if old_bond_type is not None:
                # We modified an existing bond, restore it
                bond.SetBondType(old_bond_type)
            else:
                # We added a new bond, remove it
                mol.RemoveBond(atom1_idx, atom2_idx)
        return False


def recalculate_box_and_offset(mol, scale, max_x, max_y):
    """
    Recalculate box and y_offset for a molecule at a given scale.
    Centers the view on the molecule's actual bounding box.
    Returns (box, y_offset).
    """
    conf = mol.GetConformer(0)
    actual_box = get_box(conf)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = actual_box

    # Calculate center of molecule
    center_x = (xmin + xmax) / 2
    center_y = (ymin + ymax) / 2

    # Calculate box dimensions that fill the screen at this scale
    screen_width = max_x - 2 * PADDING
    screen_height = max_y - 2 - 2 * PADDING  # Leave room for instructions
    mol_width = screen_width / scale[0]
    mol_height = screen_height / scale[1]

    # Create box centered on molecule center
    box = ((center_x - mol_width / 2, center_y - mol_height / 2, 0.0),
           (center_x + mol_width / 2, center_y + mol_height / 2, 0.0))

    # Calculate vertical offset to center the displayed content
    mol_display_height = int(mol_height * scale[1] + 2 * PADDING)
    available_height = max_y - 2
    y_offset = max(0, (available_height - mol_display_height) // 2)

    return box, y_offset


def calculate_box_and_scale(mol, max_x, max_y):
    """Calculate bounding box and scale for a molecule, centered on screen."""
    conf = mol.GetConformer(0)
    actual_box = get_box(conf)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = actual_box

    # Calculate scale to fit the molecule
    xscale = min((max_x - PADDING * 2) / (xmax - xmin), DEFAULT_SCALE)
    yscale = xscale * ASPECT_RATIO
    scale = (xscale, yscale)

    # Recalculate box and offset at this scale
    box, y_offset = recalculate_box_and_offset(mol, scale, max_x, max_y)

    return box, scale, y_offset


def draw_mol(stdscr, mol, box, scale, max_y, y_offset):
    """Draw the molecule using ASCII art using the given box and scale."""
    # Nothing to draw if molecule is empty
    if mol.GetNumAtoms() == 0:
        return

    try:
        Chem.Kekulize(mol, True)
    except Exception:
        logging.exception("Error kekulizing molecule in draw_mol")
        # If kekulization fails, skip drawing bonds (just show atoms)
        pass

    conf = mol.GetConformer(0)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = box

    # Calculate screen size - make it large enough to handle any atom position
    # Use a generous size to accommodate atoms added outside the original box
    rows = max_y - 2  # Leave room for instructions
    cols = 200  # Generous width

    screen = [[' '] * cols for i in range(rows)]
    screen_colors = [[0] * cols for i in range(rows)]  # 0 = default color

    # Color mapping for elements
    element_colors = {
        'O': 1,  # Red
        'N': 2,  # Blue
        'S': 3,  # Yellow
        'P': 3,  # Yellow
        'F': 4,  # Green
        'Cl': 4,  # Green
        'Br': 4,  # Green
        'I': 4,  # Green
    }

    try:
        for bond in mol.GetBonds():
            x1, y1 = int_coords_for_atom(bond.GetBeginAtom(), box, scale, conf,
                                         y_offset, rows)
            x2, y2 = int_coords_for_atom(bond.GetEndAtom(), box, scale, conf,
                                         y_offset, rows)
            # Only draw if bond type is in our dictionary
            if bond.GetBondType() in BOND_CHARS:
                draw_line(screen, BOND_CHARS[bond.GetBondType()], x1, y1, x2,
                          y2)
    except Exception:
        logging.exception("Error drawing bonds in draw_mol")
        # If there's any issue drawing bonds, continue to draw atoms
        pass

    for atom in mol.GetAtoms():
        x, y = int_coords_for_atom(atom, box, scale, conf, y_offset, rows)
        sym = atom.GetSymbol()
        color = element_colors.get(sym, 0)  # Get color for this element
        # Make all symbols bold except C
        is_bold = sym != 'C'

        for i, c in enumerate(sym):
            # Bounds check to avoid IndexError
            if 0 <= y < rows and 0 <= x + i < cols:
                screen[y][x + i] = c
                # Store color and bold flag (color in lower bits, bold in high bit)
                screen_colors[y][x + i] = color | (0x100 if is_bold else 0)

        # Draw formal charge if non-zero
        charge = atom.GetFormalCharge()
        if charge != 0:
            # Format charge string
            if charge == 1:
                charge_str = "+"
            elif charge == -1:
                charge_str = "-"
            elif charge > 0:
                charge_str = f"{charge}+"
            else:  # charge < 0
                charge_str = f"{abs(charge)}-"

            # Position: one cell above, one cell to the right of the symbol
            # (accounts for symbol length - e.g., "Cl" vs "C")
            charge_x = x + len(sym)
            charge_y = y - 1

            # Draw charge string with same color and bold as atom
            for i, c in enumerate(charge_str):
                if 0 <= charge_y < rows and 0 <= charge_x + i < cols:
                    screen[charge_y][charge_x + i] = c
                    screen_colors[charge_y][charge_x +
                                            i] = color | (0x100
                                                          if is_bold else 0)

    # Display screen without reversing (coordinates are already correct)
    for i in range(len(screen)):
        for j in range(len(screen[i])):
            char = screen[i][j]
            color_data = screen_colors[i][j]
            color = color_data & 0xFF  # Lower 8 bits
            is_bold = (color_data & 0x100) != 0  # Bit 8
            try:
                attr = 0
                if color > 0:
                    attr |= curses.color_pair(color)
                if is_bold:
                    attr |= curses.A_BOLD

                if attr > 0:
                    stdscr.addstr(i, j, char, attr)
                else:
                    stdscr.addstr(i, j, char)
            except curses.error:
                pass


def get_smiles(mol):
    mol_for_smiles = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(mol_for_smiles)
    except:
        mol_for_smiles = mol
    return Chem.MolToSmiles(mol_for_smiles)


def save_state(mol, box, scale, y_offset):
    """Save current state for undo/redo. Returns a deep copy of the state."""
    return {
        'mol': Chem.RWMol(mol) if mol is not None else None,
        'box': box,
        'scale': scale,
        'y_offset': y_offset
    }


def restore_state(state):
    """Restore state from saved state dict. Returns (mol, box, scale, y_offset)."""
    mol = Chem.RWMol(state['mol']) if state['mol'] is not None else None
    return mol, state['box'], state['scale'], state['y_offset']


class UndoHistory:
    """Manages undo/redo history for molecule editing."""

    def __init__(self, mol, box, scale, y_offset):
        self.mol = mol
        self.box = box
        self.scale = scale
        self.y_offset = y_offset
        self._history = [save_state(mol, box, scale, y_offset)]
        self._index = 0

    def __enter__(self):
        """Context manager entry: truncate future history."""
        self._history = self._history[:self._index + 1]
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit: save state if no exception occurred."""
        if exc_type is None:
            self._history.append(
                save_state(self.mol, self.box, self.scale, self.y_offset))
            self._index = len(self._history) - 1
        return False  # Don't suppress exceptions

    def undo(self):
        """Move back in history. Returns True if successful."""
        if self._index > 0:
            self._index -= 1
            self.mol, self.box, self.scale, self.y_offset = restore_state(
                self._history[self._index])
            return True
        return False

    def redo(self):
        """Move forward in history. Returns True if successful."""
        if self._index < len(self._history) - 1:
            self._index += 1
            self.mol, self.box, self.scale, self.y_offset = restore_state(
                self._history[self._index])
            return True
        return False

    def save_to_history(self):
        """Truncate future history and save current state."""
        self._history = self._history[:self._index + 1]
        self._history.append(save_state(self.mol, self.box, self.scale, self.y_offset))
        self._index = len(self._history) - 1


def insert_or_modify_atom(stdscr,
                          mol,
                          box,
                          scale,
                          y_offset,
                          cursor_x,
                          cursor_y,
                          max_y,
                          element_symbol=None):
    """
    Handle the 'i' command: insert atom at cursor or modify existing atom.
    If element_symbol is provided, use it; otherwise prompt the user.
    Returns mol if successful, None if no change was made.
    """
    # Check if cursor is on an atom
    atom_idx = find_atom_at_cursor(mol, cursor_x, cursor_y, box, scale, max_y,
                                   y_offset)

    # Get element symbol if not provided
    if element_symbol is None:
        element_symbol = enter_element(stdscr, max_y)

    if element_symbol:
        try:
            if atom_idx is not None:
                # Change existing atom's symbol
                atom = mol.GetAtomWithIdx(atom_idx)
                # Check if atom already has this symbol
                if atom.GetSymbol() == element_symbol:
                    return None  # No change needed
                atom.SetAtomicNum(
                    Chem.GetPeriodicTable().GetAtomicNumber(element_symbol))
            else:
                # Add new atom to molecule
                atom_idx = mol.AddAtom(Chem.Atom(element_symbol))

                # Convert cursor position to molecule coordinates
                mol_x, mol_y = screen_to_mol_coords(cursor_x, cursor_y, box,
                                                    scale, max_y, y_offset)

                # Set atom position in conformer
                conf = mol.GetConformer()
                conf.SetAtomPosition(atom_idx, [mol_x, mol_y, 0.0])

            return mol
        except Exception:
            logging.exception("Error inserting/modifying atom")
            pass

    return None


def create_or_adjust_bond(history, cursor_x, cursor_y, max_y, bond_order):
    """
    Create or adjust bond between two nearest atoms at cursor position.
    bond_order: 1 (single), 2 (double), or 3 (triple)
    Returns True if bond was created/modified, False otherwise.
    """
    if history.mol is None or history.box is None or history.scale is None:
        return False

    # Find the two atoms that should be bonded (using screen coordinates)
    atom_pair = find_bond_atoms_screen(history.mol, cursor_x, cursor_y,
                                       history.box, history.scale,
                                       max_y, history.y_offset)

    if atom_pair is not None:
        atom1_idx, atom2_idx = atom_pair
        return modify_bond(history.mol, atom1_idx, atom2_idx, bond_order)

    return False


def adjust_formal_charge(history, cursor_x, cursor_y, max_y, delta):
    """
    Adjust formal charge of atom at cursor position by delta.
    Returns True if charge was adjusted, False if no change or no atom found.
    """
    if history.mol is None or history.box is None or history.scale is None:
        return False

    if delta == 0:
        return False  # No change requested

    atom_idx = find_atom_at_cursor(history.mol, cursor_x, cursor_y,
                                   history.box, history.scale,
                                   max_y, history.y_offset)
    if atom_idx is not None:
        atom = history.mol.GetAtomWithIdx(atom_idx)
        current_charge = atom.GetFormalCharge()
        atom.SetFormalCharge(current_charge + delta)
        return True

    return False


def delete_at_cursor(history, cursor_x, cursor_y, max_y):
    """
    Delete atom or bond at cursor position.
    Returns True if something was deleted, False otherwise.
    """
    if history.mol is None or history.box is None or history.scale is None:
        return False

    # First try to find an atom at cursor
    atom_idx = find_atom_at_cursor(history.mol, cursor_x, cursor_y,
                                   history.box, history.scale,
                                   max_y, history.y_offset)
    if atom_idx is not None:
        try:
            history.mol.RemoveAtom(atom_idx)
            return True
        except Exception:
            logging.exception("Error removing atom (x command)")
            return False
    else:
        # No atom found, try to delete a bond instead
        atom_pair = find_bond_atoms_screen(history.mol, cursor_x, cursor_y,
                                           history.box, history.scale,
                                           max_y, history.y_offset)
        if atom_pair is not None:
            atom1_idx, atom2_idx = atom_pair
            if modify_bond(history.mol, atom1_idx, atom2_idx, 0):
                return True

    return False


def clear_canvas(history, max_y):
    """
    Clear the canvas and reset to blank slate with default settings.
    """
    # Create empty molecule
    history.mol = Chem.RWMol()
    conf = Chem.Conformer()
    history.mol.AddConformer(conf)

    # Default box: 20 angstroms centered at origin
    box_size = 10.0
    history.box = ((-box_size, -box_size, 0.0), (box_size, box_size, 0.0))

    # Use default scale
    xscale = DEFAULT_SCALE
    yscale = xscale * ASPECT_RATIO
    history.scale = (xscale, yscale)

    # Center vertically
    mol_height = int(2 * box_size * yscale + 2 * PADDING)
    available_height = max_y - 2
    history.y_offset = max(0, (available_height - mol_height) // 2)


def shift_view(history, dx_sign, dy_sign):
    """
    Shift the view (pan the molecule) by one step in the given direction.
    dx_sign, dy_sign: -1, 0, or +1 indicating direction.
    """
    if history.box is None or history.scale is None:
        return

    (xmin, ymin, zmin), (xmax, ymax, zmax) = history.box

    # Calculate step size based on scale
    dx = dx_sign / history.scale[0] if dx_sign != 0 else 0
    dy = dy_sign / history.scale[1] if dy_sign != 0 else 0

    history.box = ((xmin + dx, ymin + dy, zmin), (xmax + dx, ymax + dy, zmax))


def cleanup_coordinates(history, max_x, max_y):
    """
    Regenerate 2D coordinates for the molecule and recenter the view.
    Keeps the current zoom level. Returns True if successful.
    """
    if history.mol is None or history.mol.GetNumAtoms() == 0:
        return False

    try:
        AllChem.Compute2DCoords(history.mol)
        if history.scale is not None:
            history.box, history.y_offset = recalculate_box_and_offset(
                history.mol, history.scale, max_x, max_y)
        return True
    except Exception:
        logging.exception("Error regenerating coordinates")
        return False


def zoom_view(history, max_x, max_y, zoom_factor):
    """
    Zoom in or out by the given factor.
    zoom_factor > 1 means zoom in, < 1 means zoom out.
    """
    if history.box is None or history.scale is None:
        return

    # Calculate current center of the box
    (xmin, ymin, zmin), (xmax, ymax, zmax) = history.box
    center_x = (xmin + xmax) / 2
    center_y = (ymin + ymax) / 2

    # Adjust scale by zoom factor and clamp to limits
    xscale = history.scale[0] * zoom_factor
    xscale = max(MIN_SCALE, min(MAX_SCALE, xscale))
    yscale = xscale * ASPECT_RATIO
    history.scale = (xscale, yscale)

    # Calculate new box dimensions to show at new scale
    screen_width = max_x - 2 * PADDING
    screen_height = max_y - 2 - 2 * PADDING

    # Molecule coordinate range that fits on screen at new scale
    mol_width = screen_width / xscale
    mol_height = screen_height / yscale

    # New box centered on the same point
    history.box = ((center_x - mol_width / 2, center_y - mol_height / 2, 0.0),
                   (center_x + mol_width / 2, center_y + mol_height / 2, 0.0))

    # Recalculate y_offset for new scale
    mol_display_height = int(mol_height * yscale + 2 * PADDING)
    available_height = max_y - 2
    history.y_offset = max(0, (available_height - mol_display_height) // 2)


def append_smiles_fragment(stdscr, mol, box, scale, y_offset, cursor_x,
                           cursor_y, max_x, max_y):
    """
    Handle the 'a' command: append atoms from SMILES to atom or bond under cursor.
    Returns (mol, box, scale, y_offset) if successful, None if no change was made.
    """
    # Find atom under cursor
    atom_idx = find_atom_at_cursor(mol, cursor_x, cursor_y, box, scale, max_y,
                                   y_offset)

    # Check if on a bond if not on an atom
    bond_atom_pair = None
    if atom_idx is None:
        bond_atom_pair = find_bond_atoms_screen(mol, cursor_x, cursor_y,
                                                box, scale, max_y, y_offset)

    if atom_idx is not None or bond_atom_pair is not None:
        # Prompt for SMILES
        sidechain_smiles = enter_smiles(stdscr, max_y)
        if sidechain_smiles:
            try:
                # Create sidechain molecule
                sidechain = Chem.MolFromSmiles(sidechain_smiles)
                if sidechain is not None:
                    # Get the index where sidechain will start
                    start_idx = mol.GetNumAtoms()
                    end_idx = start_idx + sidechain.GetNumAtoms() - 1

                    # Insert sidechain into molecule
                    mol.InsertMol(sidechain)

                    if atom_idx is not None:
                        # Cursor on atom: connect to first atom of sidechain
                        mol.AddBond(atom_idx, start_idx, Chem.BondType.SINGLE)
                    else:
                        # Cursor on bond: insert sidechain between the two atoms
                        a1_idx, a2_idx = bond_atom_pair

                        # Determine which atom is closer to cursor
                        conf = mol.GetConformer()
                        mol_x, mol_y = screen_to_mol_coords(
                            cursor_x, cursor_y, box, scale, max_y, y_offset)

                        pos1 = conf.GetAtomPosition(a1_idx)
                        pos2 = conf.GetAtomPosition(a2_idx)

                        dist1 = math.sqrt((pos1.x - mol_x)**2 +
                                          (pos1.y - mol_y)**2)
                        dist2 = math.sqrt((pos2.x - mol_x)**2 +
                                          (pos2.y - mol_y)**2)

                        if dist1 > dist2:
                            a1_idx, a2_idx = a2_idx, a1_idx

                        # Connect sidechain: a1 -> first atom, a2 -> last atom
                        mol.AddBond(a1_idx, start_idx, Chem.BondType.SINGLE)
                        mol.AddBond(a2_idx, end_idx, Chem.BondType.SINGLE)

                    # Create coordinate map to keep original atoms fixed
                    coord_map = {}
                    conf = mol.GetConformer()
                    for i in range(start_idx):
                        pos = conf.GetAtomPosition(i)
                        coord_map[i] = Geometry.Point2D(pos.x, pos.y)

                    # Compute 2D coordinates for new atoms only
                    AllChem.Compute2DCoords(mol, coordMap=coord_map)

                    # Update box to show all atoms while keeping same scale
                    box, y_offset = recalculate_box_and_offset(mol, scale, max_x, max_y)

                    return mol, box, scale, y_offset
            except Exception as e:
                logging.exception("Error appending atoms (a command)")
                # Error appending atoms
                pass

    return None


def redraw_screen(stdscr, history, show_smiles, instructions, max_x, max_y):
    """Redraw the entire screen with molecule, SMILES, and instructions."""
    stdscr.clear()

    # Draw molecule if present
    if history.mol is not None and history.box is not None and history.scale is not None:
        draw_mol(stdscr, history.mol, history.box, history.scale, max_y,
                 history.y_offset)

    # Draw SMILES at the top if enabled (after molecule so it's on top)
    if show_smiles and history.mol is not None:
        current_smiles = get_smiles(history.mol)
        # Wrap SMILES to screen width
        row = 0
        for i in range(0, len(current_smiles), max_x - 1):
            chunk = current_smiles[i:i + max_x - 1]
            try:
                stdscr.addstr(row, 0, chunk)
                row += 1
            except curses.error:
                break

    # Draw instructions at the bottom
    for i, line in enumerate(instructions):
        try:
            stdscr.addstr(max_y - len(instructions) + i, 0,
                          line[:max_x - 1])
        except curses.error:
            pass


def main_loop(stdscr, initial_smiles=None):
    # Initialize curses
    curses.curs_set(1)  # Show cursor
    curses.use_default_colors()  # Use terminal's default colors

    # Initialize colors
    curses.start_color()
    curses.init_pair(1, curses.COLOR_RED, -1)  # Oxygen - red
    curses.init_pair(2, curses.COLOR_BLUE, -1)  # Nitrogen - blue
    curses.init_pair(3, curses.COLOR_YELLOW, -1)  # Sulfur - yellow
    curses.init_pair(4, curses.COLOR_GREEN, -1)  # Halogens - green

    stdscr.clear()

    # Get screen dimensions
    max_y, max_x = stdscr.getmaxyx()

    # Starting cursor position (center of screen)
    cursor_y, cursor_x = max_y // 2, max_x // 2

    # SMILES string storage
    smiles = initial_smiles or ""
    mol = Chem.RWMol()
    box = None
    scale = None
    y_offset = 0
    show_smiles = False

    # Load initial molecule if provided
    if initial_smiles:
        m = Chem.MolFromSmiles(initial_smiles)
        Chem.Kekulize(mol)
        if m is not None:
            mol = Chem.RWMol(m)
            AllChem.Compute2DCoords(mol)
            box, scale, y_offset = calculate_box_and_scale(mol, max_x, max_y)
    else:
        # Set up default coordinate system for empty molecule
        # Add an empty conformer
        conf = Chem.Conformer()
        mol.AddConformer(conf)

        # Default box: 20 angstroms centered at origin
        box_size = 10.0
        box = ((-box_size, -box_size, 0.0), (box_size, box_size, 0.0))
        # Use max scale
        xscale = DEFAULT_SCALE
        yscale = xscale * ASPECT_RATIO
        scale = (xscale, yscale)
        # Center vertically
        mol_height = int(2 * box_size * yscale + 2 * PADDING)
        available_height = max_y - 2
        y_offset = max(0, (available_height - mol_height) // 2)

    # Undo/redo history
    history = UndoHistory(mol, box, scale, y_offset)

    # Instructions (try to keep lines under 80 characters and more or less balanced)
    instructions = [
        "hjkl: move | HJKL: shift | s/S: SMILES | i/a/c/n/o: insert | x: del | +/-: chg",
        "<>: zoom | u/r: undo/redo | ^L: clean | 1-3: bond | @: clear | ?: help | q: quit"
    ]

    # Track when we need to redraw the entire screen
    need_redraw = True

    while True:
        # Only redraw everything when necessary
        if need_redraw:
            redraw_screen(stdscr, history, show_smiles, instructions, max_x, max_y)
            need_redraw = False

        # Move cursor to current position
        try:
            stdscr.move(cursor_y, cursor_x)
        except curses.error:
            pass

        stdscr.refresh()

        # Get user input
        key = stdscr.getch()

        # Handle movement (vi-style)
        if key == ord('h'):  # left
            cursor_x = max(0, cursor_x - 1)
        elif key == ord('j'):  # down
            cursor_y = min(max_y - 1, cursor_y + 1)
        elif key == ord('k'):  # up
            cursor_y = max(0, cursor_y - 1)
        elif key == ord('l'):  # right
            cursor_x = min(max_x - 1, cursor_x + 1)

        # Shift molecule (move all atoms)
        elif key == ord('K'):  # shift up
            shift_view(history, 0, -1)
            need_redraw = True

        elif key == ord('H'):  # shift left
            shift_view(history, 1, 0)
            need_redraw = True

        elif key == ord('J'):  # shift down
            shift_view(history, 0, 1)
            need_redraw = True

        elif key == ord('L'):  # shift right
            shift_view(history, -1, 0)
            need_redraw = True

        # Enter SMILES string
        elif key == ord('s'):
            smiles = enter_smiles(stdscr, max_y)
            if smiles:
                m = Chem.MolFromSmiles(smiles)
                if m is not None:
                    with history:
                        history.mol = Chem.RWMol(m)
                        AllChem.Compute2DCoords(history.mol)
                        history.box, history.scale, history.y_offset = calculate_box_and_scale(
                            history.mol, max_x, max_y)

            # Always redraw to clear the prompt
            need_redraw = True

        # Toggle SMILES display
        elif key == ord('S'):
            show_smiles = not show_smiles
            need_redraw = True

        # Insert atom at cursor position or change atom symbol
        elif key == ord('i'):
            if history.mol is not None and history.box is not None and history.scale is not None:
                result = insert_or_modify_atom(stdscr, history.mol, history.box,
                                               history.scale, history.y_offset,
                                               cursor_x, cursor_y, max_y)
                if result is not None:
                    with history:
                        history.mol = result

                # Always redraw to clear the prompt
                need_redraw = True

        # Insert common atoms (c, n, o) - shortcuts, or change atom symbol
        elif key in [ord('c'), ord('n'), ord('o')]:
            if history.mol is not None and history.box is not None and history.scale is not None:
                symbol = chr(key).upper()
                result = insert_or_modify_atom(stdscr, history.mol, history.box,
                                               history.scale, history.y_offset,
                                               cursor_x, cursor_y, max_y,
                                               symbol)
                if result is not None:
                    with history:
                        history.mol = result

                    need_redraw = True

        # Append atoms from SMILES to atom under cursor or bond
        elif key == ord('a'):
            if history.mol is not None and history.box is not None and history.scale is not None:
                result = append_smiles_fragment(stdscr, history.mol,
                                                history.box, history.scale,
                                                history.y_offset, cursor_x,
                                                cursor_y, max_x, max_y)
                if result is not None:
                    with history:
                        history.mol, history.box, history.scale, history.y_offset = result

                # Always redraw to clear the prompt
                need_redraw = True

        # Delete atom or bond at cursor position
        elif key == ord('x'):
            with history:
                if delete_at_cursor(history, cursor_x, cursor_y, max_y):
                    need_redraw = True

        # Increase formal charge
        elif key == ord('+'):
            with history:
                if adjust_formal_charge(history, cursor_x, cursor_y, max_y, 1):
                    need_redraw = True

        # Decrease formal charge
        elif key == ord('-'):
            with history:
                if adjust_formal_charge(history, cursor_x, cursor_y, max_y, -1):
                    need_redraw = True

        # Cleanup/regenerate coordinates (Ctrl-L)
        elif key == 12:  # Ctrl-L
            with history:
                if cleanup_coordinates(history, max_x, max_y):
                    need_redraw = True
                    pass

        # Zoom out
        elif key == ord('<'):
            zoom_view(history, max_x, max_y, 1.0 / 1.2)
            need_redraw = True

        # Zoom in
        elif key == ord('>'):
            zoom_view(history, max_x, max_y, 1.2)
            need_redraw = True

        # Add/modify/delete bond
        elif key in [ord('1'), ord('2'), ord('3')]:
            bond_order = int(chr(key))
            if create_or_adjust_bond(history, cursor_x, cursor_y, max_y, bond_order):
                history.save_to_history()
                need_redraw = True

        # Clear canvas (reset to blank slate)
        elif key == ord('@'):
            with history:
                clear_canvas(history, max_y)
            need_redraw = True

        # Undo
        elif key == ord('u'):
            if history.undo():
                need_redraw = True

        # Redo
        elif key == ord('r'):
            if history.redo():
                need_redraw = True

        # Help
        elif key == ord('?'):
            stdscr.clear()
            # Display the help text from the docstring
            help_text = __doc__.strip()
            help_lines = help_text.split('\n')

            # Display help text
            for i, line in enumerate(help_lines):
                try:
                    stdscr.addstr(i, 0, line[:max_x - 1])
                except curses.error:
                    pass

            # Add "press any key" message
            try:
                stdscr.addstr(
                    len(help_lines) + 1, 0, "Press any key to exit help")
            except curses.error:
                pass

            stdscr.refresh()
            stdscr.getch()  # Wait for any key
            need_redraw = True

        # Quit
        elif key == ord('q'):
            return get_smiles(history.mol)


def main():
    # Set up logging to file (truncate on start)
    logging.basicConfig(
        filename='cursemol.log',
        filemode='w',  # Truncate on open
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s')

    # Silence RDKit warnings
    logger = RDLogger.logger()
    logger.setLevel(RDLogger.CRITICAL)

    parser = argparse.ArgumentParser(
        description='Display molecules in the terminal')
    parser.add_argument(
        'smiles',
        nargs='?',
        help='Initial SMILES string to display (use "-" to read from stdin)')
    args = parser.parse_args()

    # Handle reading from stdin if "-" is provided
    initial_smiles = args.smiles
    if initial_smiles == "-":
        initial_smiles = sys.stdin.readline().strip()

    # If stdin is not a TTY (e.g., piped input), redirect to /dev/tty
    # so curses can read keyboard input
    if not sys.stdin.isatty():
        sys.stdin.close()  # Close the old stdin to avoid resource warning
        tty_fd = os.open('/dev/tty', os.O_RDONLY)
        os.dup2(tty_fd, 0)  # Replace fd 0 (stdin) with /dev/tty
        os.close(tty_fd)
        sys.stdin = os.fdopen(0, 'r')
        # Register cleanup to avoid resource warning on exit
        atexit.register(lambda: sys.stdin.close())

    print(curses.wrapper(main_loop, initial_smiles))


if __name__ == "__main__":
    main()
