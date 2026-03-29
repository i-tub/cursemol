#!/usr/bin/env python3
"""
CurseMol - Molecular sketcher for the terminally insane

Controls:
  h, j, k, l - Move cursor left, down, up, right
  H, J, K, L - Translate molecule left, down, up, down, right
  s          - Enter a SMILES string to replace the current molecule
  S          - Toggle SMILES display
  i          - Insert/modify atom at cursor position
  a          - Append atoms from SMILES to atom or bond under cursor
               (appending to a bond forms a ring by connecting the bond atoms
               to the first and last atoms from the SMILES)
  c, n, o    - Insert/modify carbon/nitrogen/oxygen atom
  x          - Delete atom or bond
  X          - Area delete (select rectangle, Enter to delete, Esc to cancel)
  +, -       - Increase/decrease formal charge on atom
  <, >       - Zoom out/in
  1, 2, 3    - Add bond or change bond (order 1/2/3) between nearest atoms
  w, d       - Add/change to wedge or dash bond (press again to reverse)
  @          - Clear canvas (reset to blank slate)
  u          - Undo
  r          - Redo
  Ctrl-L     - Clean up (regenerate coordinates)
  ?          - Show this help
  q          - Quit and print SMILES to stdout
"""

import argparse
import atexit
import curses
from dataclasses import dataclass
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
    Chem.BondType.SINGLE: '·',
    Chem.BondType.DOUBLE: '=',
    Chem.BondType.TRIPLE: '#',
}

BOND_DIR_CHARS = {
    Chem.BondDir.BEGINWEDGE: '•',
    Chem.BondDir.BEGINDASH: '◦',
}

# Color mapping for elements
ELEMENT_COLORS = {
    'O': 1,  # Red
    'N': 2,  # Blue
    'S': 3,  # Yellow
    'P': 3,  # Yellow
    'F': 4,  # Green
    'Cl': 4,  # Green
    'Br': 4,  # Green
    'I': 4,  # Green
}

# Instructions (try to keep lines under 80 characters and more or less
# balanced)
INSTRUCTIONS = [
    "hjkl: move | HJKL: translate | s/S: SMILES | i/a/c/n/o: insert | x/X: del",
    "+/-: chg | <>: zoom | u/r: undo | ^L: clean | 123/wd: bond | @: clear | ?: help"
]


@dataclass
class State:
    """Molecular drawing state: molecule and its display parameters."""
    mol: Chem.RWMol
    box: tuple  # ((min_x, min_y, min_z), (max_x, max_y, max_z))
    scale: tuple  # (xscale, yscale)
    y_offset: int


class UndoHistory:
    """Manages undo/redo history for molecule editing."""

    def __init__(self, state):
        self.state = state
        self._history = [save_state(state)]
        self._index = 0

    def __enter__(self):
        """Context manager entry: truncate future history."""
        self._history = self._history[:self._index + 1]
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit: save state if no exception occurred."""
        if exc_type is None:
            self._history.append(save_state(self.state))
            self._index = len(self._history) - 1
        return False  # Don't suppress exceptions

    def undo(self):
        """Move back in history. Returns True if successful."""
        if self._index > 0:
            self._index -= 1
            self.state = restore_state(self._history[self._index])
            return True
        return False

    def redo(self):
        """Move forward in history. Returns True if successful."""
        if self._index < len(self._history) - 1:
            self._index += 1
            self.state = restore_state(self._history[self._index])
            return True
        return False

    def push(self):
        """Truncate future history and save current state."""
        self._history = self._history[:self._index + 1]
        self._history.append(save_state(self.state))
        self._index = len(self._history) - 1


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


def draw_line(screen, char, char2, x1, y1, x2, y2):
    vertical = False
    if abs(x2 - x1) < abs(y2 - y1):
        x1, y1 = y1, x1
        x2, y2 = y2, x2
        vertical = True
    try:
        slope = 1.0 * (y2 - y1) / (x2 - x1)
    except ZeroDivisionError:
        return

    rev = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        rev = True
    mid = round((x1 + 1 + x2) / 2)
    for x in range(x1 + 1, x2):
        y = int(round(y1 + slope * (x - x1)))
        if rev:
            c = char if x > mid else char2
        else:
            c = char if x < mid else char2
        if vertical:
            screen[x][y] = c
        else:
            screen[y][x] = c


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


def find_atom_at_cursor(state, cursor_x, cursor_y, max_y, tolerance=1):
    """
    Find an atom at or near the cursor position (within tolerance cells).
    Returns atom index or None if no atom found.
    """
    if state.mol.GetNumAtoms() == 0:
        return None

    conf = state.mol.GetConformer()
    rows = max_y - 2

    # Check atoms to find one whose screen position is within tolerance
    for atom in state.mol.GetAtoms():
        screen_x, screen_y = int_coords_for_atom(atom, state.box, state.scale,
                                                 conf, state.y_offset, rows)

        # Screen position is now directly the terminal position
        # Check if within tolerance
        if (abs(screen_x - cursor_x) <= tolerance and
                abs(screen_y - cursor_y) <= tolerance):
            return atom.GetIdx()

    return None


def find_bond_atoms(state, cursor_x, cursor_y, max_y):
    """
    Find the two atoms closest to cursor position such that
    the cursor is roughly between them (angle >= 160 degrees).
    Uses screen coordinates to avoid floating point precision issues.
    Returns (atom1_idx, atom2_idx) or None if no valid pair found.
    """
    if state.mol.GetNumAtoms() < 2:
        return None

    conf = state.mol.GetConformer()
    rows = max_y - 2

    # Calculate screen positions and distances for all atoms
    distances = []
    for atom in state.mol.GetAtoms():
        screen_x, screen_y = int_coords_for_atom(atom, state.box, state.scale,
                                                 conf, state.y_offset, rows)
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


def reverse_bond(bond):
    """
    Reverse bond by deleting and re-adding with swapped atoms
    we don't know the direction.
    """
    mol = bond.GetOwningMol()
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()
    bond_dir = bond.GetBondDir()
    bond_type = bond.GetBondType()
    mol.RemoveBond(a1, a2)
    bond_idx = mol.AddBond(a2, a1, bond_type) - 1
    bond = mol.GetBondWithIdx(bond_idx)
    bond.SetBondDir(bond_dir)


def modify_bond(mol, atom1_idx, atom2_idx, bond_order, bond_dir=None):
    """
    Modify or create a bond between two atoms.
    bond_order: 0 (delete), 1 (single), 2 (double), 3 (triple)
    bond_dir: Optional bond direction (e.g., Chem.BondDir.BEGINWEDGE)

    If a bond already has the specified direction, it will be reversed
    (atoms swapped).

    Returns True if successful, False if no change was made.
    """
    bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)

    # Map bond order to BondType
    bond_type_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }

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
            current_type = bond.GetBondType()
            current_dir = bond.GetBondDir()

            if (current_type == bond_type and bond_dir is None and
                current_dir == Chem.BondDir.NONE):
                return False

            # Check if bond already has this exact type and direction
            # If so, reverse the bond (swap atoms)
            if (current_type == bond_type and bond_dir is not None and
                    current_dir == bond_dir):
                reverse_bond(bond)
            else:
                # Modify existing bond
                bond.SetBondType(bond_type)
                if bond_dir is not None:
                    bond.SetBondDir(bond_dir)
                else:
                    bond.SetBondDir(Chem.BondDir.NONE)
        else:
            # Add new bond
            mol.AddBond(atom1_idx, atom2_idx, bond_type)
            if bond_dir is not None:
                bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                bond.SetBondDir(bond_dir)

    return True


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


def draw_atom(screen, screen_colors, atom, x, y, rows, cols):
    """
    Draw a single atom with its symbol and charge into the screen buffer.

    Args:
        screen: 2D character array
        screen_colors: 2D color/attribute array
        atom: RDKit atom object
        x, y: screen coordinates for the atom
        rows, cols: screen buffer dimensions
    """
    sym = atom.GetSymbol()
    color = ELEMENT_COLORS.get(sym, 0)  # Get color for this element
    # Make all symbols bold except C
    is_bold = sym != 'C'

    # Draw atom symbol
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
                                        i] = color | (0x100 if is_bold else 0)


def fill_screen_buffer(state, max_y):
    """
    Fill screen buffer with bonds and atoms.
    Returns (screen, screen_colors) tuple of 2D arrays.
    """
    # Calculate screen size
    rows = max_y - 2  # Leave room for instructions
    cols = 200  # Generous width

    screen = [[' '] * cols for i in range(rows)]
    screen_colors = [[0] * cols for i in range(rows)]  # 0 = default color

    conf = state.mol.GetConformer(0)

    # Draw bonds
    try:
        for bond in state.mol.GetBonds():
            x1, y1 = int_coords_for_atom(bond.GetBeginAtom(), state.box,
                                         state.scale, conf, state.y_offset,
                                         rows)
            x2, y2 = int_coords_for_atom(bond.GetEndAtom(), state.box,
                                         state.scale, conf, state.y_offset,
                                         rows)
            # Only draw if bond type is in our dictionary
            if bond_char := BOND_CHARS.get(bond.GetBondType()):
                bond_dir_char = BOND_DIR_CHARS.get(bond.GetBondDir(), bond_char)
                draw_line(screen, bond_char, bond_dir_char, x1, y1, x2, y2)
    except Exception:
        logging.exception("Error drawing bonds")
        # If there's any issue drawing bonds, continue to draw atoms
        pass

    # Draw atoms
    for atom in state.mol.GetAtoms():
        x, y = int_coords_for_atom(atom, state.box, state.scale, conf,
                                   state.y_offset, rows)
        draw_atom(screen, screen_colors, atom, x, y, rows, cols)

    return screen, screen_colors


def render_screen_buffer(stdscr, screen, screen_colors):
    """
    Render a screen buffer to the curses window.

    Args:
        stdscr: curses window
        screen: 2D character array
        screen_colors: 2D color/attribute array
    """
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


def draw_mol(stdscr, state, max_y):
    """Draw the molecule using ASCII art."""
    # Nothing to draw if molecule is empty
    if state.mol.GetNumAtoms() == 0:
        return

    try:
        Chem.Kekulize(state.mol, True)
    except Exception:
        logging.exception("Error kekulizing molecule")
        # If kekulization fails, skip drawing bonds (just show atoms)
        pass

    # Fill screen buffer with molecular structure
    screen, screen_colors = fill_screen_buffer(state, max_y)

    # Render buffer to curses window
    render_screen_buffer(stdscr, screen, screen_colors)


def get_smiles(mol):
    mol_for_smiles = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(mol_for_smiles)
        Chem.DetectBondStereochemistry(mol_for_smiles)
        Chem.AssignChiralTypesFromBondDirs(mol_for_smiles)
        Chem.AssignStereochemistry(mol_for_smiles, force=True)
    except:
        logging.exception("get_smiles error")
        mol_for_smiles = mol
    return Chem.MolToSmiles(mol_for_smiles)


def save_state(state):
    """Save current state for undo/redo. Returns a deep copy of the state."""
    return State(mol=Chem.RWMol(state.mol),
                 box=state.box,
                 scale=state.scale,
                 y_offset=state.y_offset)


def restore_state(saved_state):
    """Restore state from saved state. Returns a new State object."""
    return State(mol=Chem.RWMol(saved_state.mol),
                 box=saved_state.box,
                 scale=saved_state.scale,
                 y_offset=saved_state.y_offset)


def insert_or_modify_atom(stdscr,
                          state,
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
    atom_idx = find_atom_at_cursor(state, cursor_x, cursor_y, max_y)

    # Get element symbol if not provided
    if element_symbol is None:
        element_symbol = enter_element(stdscr, max_y)

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
                mol_x, mol_y = screen_to_mol_coords(cursor_x, cursor_y,
                                                    state.box, state.scale,
                                                    max_y, state.y_offset)

                # Set atom position in conformer
                conf = state.mol.GetConformer()
                conf.SetAtomPosition(atom_idx, [mol_x, mol_y, 0.0])

            return state.mol
        except Exception:
            logging.exception("Error inserting/modifying atom")
            pass

    return None


def create_or_adjust_bond(state,
                          cursor_x,
                          cursor_y,
                          max_y,
                          bond_order,
                          bond_dir=None):
    """
    Create or adjust bond between two nearest atoms at cursor position.
    bond_order: 1 (single), 2 (double), or 3 (triple)
    bond_dir: Optional bond direction (e.g., Chem.BondDir.BEGINWEDGE)
    Returns True if bond was created/modified, False otherwise.
    """
    # Find the two atoms that should be bonded (using screen coordinates)
    atom_pair = find_bond_atoms(state, cursor_x, cursor_y, max_y)

    if atom_pair is not None:
        atom1_idx, atom2_idx = atom_pair
        return modify_bond(state.mol, atom1_idx, atom2_idx, bond_order,
                           bond_dir)

    return False


def adjust_formal_charge(state, cursor_x, cursor_y, max_y, delta):
    """
    Adjust formal charge of atom at cursor position by delta.
    Returns True if charge was adjusted, False if no change or no atom found.
    """
    if delta == 0:
        return False  # No change requested

    atom_idx = find_atom_at_cursor(state, cursor_x, cursor_y, max_y)
    if atom_idx is not None:
        atom = state.mol.GetAtomWithIdx(atom_idx)
        current_charge = atom.GetFormalCharge()
        atom.SetFormalCharge(current_charge + delta)
        return True

    return False


def delete_atoms_in_rect(state, x1, y1, x2, y2, max_y):
    """
    Delete all atoms whose screen positions fall within the rectangle.
    Returns True if any atoms were deleted, False otherwise.
    """
    conf = state.mol.GetConformer()
    rows = max_y - 2

    # Normalize rectangle coordinates
    min_x, max_x = (x1, x2) if x1 <= x2 else (x2, x1)
    min_y, max_y_rect = (y1, y2) if y1 <= y2 else (y2, y1)

    # Find atoms within the rectangle
    atoms_to_delete = []
    for atom in state.mol.GetAtoms():
        screen_x, screen_y = int_coords_for_atom(atom, state.box, state.scale,
                                                 conf, state.y_offset, rows)
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


def delete_at_cursor(state, cursor_x, cursor_y, max_y):
    """
    Delete atom or bond at cursor position.
    Returns True if something was deleted, False otherwise.
    """
    # First try to find an atom at cursor
    atom_idx = find_atom_at_cursor(state, cursor_x, cursor_y, max_y)
    if atom_idx is not None:
        try:
            state.mol.RemoveAtom(atom_idx)
            return True
        except Exception:
            logging.exception("Error removing atom (x command)")
            return False
    else:
        # No atom found, try to delete a bond instead
        atom_pair = find_bond_atoms(state, cursor_x, cursor_y, max_y)
        if atom_pair is not None:
            atom1_idx, atom2_idx = atom_pair
            if modify_bond(state.mol, atom1_idx, atom2_idx, 0):
                return True

    return False


def create_empty_state(max_y):
    """
    Create an empty molecular state with default settings.
    Returns State object.
    """
    # Create empty molecule
    mol = Chem.RWMol()
    conf = Chem.Conformer()
    mol.AddConformer(conf)

    # Default box: 20 angstroms centered at origin
    box_size = 10.0
    box = ((-box_size, -box_size, 0.0), (box_size, box_size, 0.0))

    # Use default scale
    xscale = DEFAULT_SCALE
    yscale = xscale * ASPECT_RATIO
    scale = (xscale, yscale)

    # Center vertically
    mol_height = int(2 * box_size * yscale + 2 * PADDING)
    available_height = max_y - 2
    y_offset = max(0, (available_height - mol_height) // 2)

    return State(mol=mol, box=box, scale=scale, y_offset=y_offset)


def create_molecule_from_smiles(smiles, max_x, max_y):
    """
    Create molecule from SMILES string with 2D coordinates.
    Returns State object if successful, None otherwise.
    """
    m = Chem.MolFromSmiles(smiles)
    if m is not None:
        mol = Chem.RWMol(m)
        AllChem.Compute2DCoords(mol)
        Chem.WedgeMolBonds(mol, mol.GetConformer())
        box, scale, y_offset = calculate_box_and_scale(mol, max_x, max_y)
        return State(mol=mol, box=box, scale=scale, y_offset=y_offset)
    return None


def load_smiles(stdscr, max_x, max_y):
    """
    Prompt user for SMILES string and create molecule from it.
    Returns State object if successful, None otherwise.
    """
    smiles = enter_smiles(stdscr, max_y)
    if smiles:
        return create_molecule_from_smiles(smiles, max_x, max_y)
    return None


def clear_canvas(state, max_y):
    """
    Clear the canvas and reset to blank slate with default settings.
    Modifies state in place.
    """
    # Create empty molecule
    mol = Chem.RWMol()
    conf = Chem.Conformer()
    mol.AddConformer(conf)

    # Default box: 20 angstroms centered at origin
    box_size = 10.0
    box = ((-box_size, -box_size, 0.0), (box_size, box_size, 0.0))

    # Use default scale
    xscale = DEFAULT_SCALE
    yscale = xscale * ASPECT_RATIO
    scale = (xscale, yscale)

    # Center vertically
    mol_height = int(2 * box_size * yscale + 2 * PADDING)
    available_height = max_y - 2
    y_offset = max(0, (available_height - mol_height) // 2)

    # Update state in place
    state.mol = mol
    state.box = box
    state.scale = scale
    state.y_offset = y_offset


def show_help(stdscr, max_x):
    """
    Display help text and wait for user to press a key.
    """
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
        stdscr.addstr(len(help_lines) + 1, 0, "Press any key to exit help")
    except curses.error:
        pass

    stdscr.refresh()
    stdscr.getch()  # Wait for any key


def shift_view(state, dx_sign, dy_sign):
    """
    Shift the view (pan the molecule) by one step in the given direction.
    dx_sign, dy_sign: -1, 0, or +1 indicating direction.
    """
    if state.box is None or state.scale is None:
        return

    (xmin, ymin, zmin), (xmax, ymax, zmax) = state.box

    # Calculate step size based on scale
    dx = dx_sign / state.scale[0] if dx_sign != 0 else 0
    dy = dy_sign / state.scale[1] if dy_sign != 0 else 0

    state.box = ((xmin + dx, ymin + dy, zmin), (xmax + dx, ymax + dy, zmax))


def cleanup_coordinates(state, max_x, max_y):
    """
    Regenerate 2D coordinates for the molecule and recenter the view.
    Keeps the current zoom level. Returns True if successful.
    """
    try:
        AllChem.Compute2DCoords(state.mol)
        if state.scale is not None:
            state.box, state.y_offset = recalculate_box_and_offset(
                state.mol, state.scale, max_x, max_y)
        return True
    except Exception:
        logging.exception("Error regenerating coordinates")
        return False


def zoom_view(state, max_x, max_y, zoom_factor):
    """
    Zoom in or out by the given factor.
    zoom_factor > 1 means zoom in, < 1 means zoom out.
    """
    if state.box is None or state.scale is None:
        return

    # Calculate current center of the box
    (xmin, ymin, zmin), (xmax, ymax, zmax) = state.box
    center_x = (xmin + xmax) / 2
    center_y = (ymin + ymax) / 2

    # Adjust scale by zoom factor and clamp to limits
    xscale = state.scale[0] * zoom_factor
    xscale = max(MIN_SCALE, min(MAX_SCALE, xscale))
    yscale = xscale * ASPECT_RATIO
    state.scale = (xscale, yscale)

    # Calculate new box dimensions to show at new scale
    screen_width = max_x - 2 * PADDING
    screen_height = max_y - 2 - 2 * PADDING

    # Molecule coordinate range that fits on screen at new scale
    mol_width = screen_width / xscale
    mol_height = screen_height / yscale

    # New box centered on the same point
    state.box = ((center_x - mol_width / 2, center_y - mol_height / 2, 0.0),
                 (center_x + mol_width / 2, center_y + mol_height / 2, 0.0))

    # Recalculate y_offset for new scale
    mol_display_height = int(mol_height * yscale + 2 * PADDING)
    available_height = max_y - 2
    state.y_offset = max(0, (available_height - mol_display_height) // 2)


def compute_coords_with_fixed_atoms(mol, num_fixed_atoms):
    """
    Compute 2D coordinates for a molecule, keeping existing atoms fixed.

    This is useful when adding new atoms to a molecule - the original atoms
    maintain their positions while new atoms are positioned around them.

    Args:
        mol: RDKit molecule
        num_fixed_atoms: Number of atoms at the beginning to keep fixed
    """
    # Create coordinate map to keep original atoms fixed
    coord_map = {}
    conf = mol.GetConformer()
    for i in range(num_fixed_atoms):
        pos = conf.GetAtomPosition(i)
        coord_map[i] = Geometry.Point2D(pos.x, pos.y)

    # Compute 2D coordinates for new atoms only
    AllChem.Compute2DCoords(mol, coordMap=coord_map)


def connect_sidechain_to_bond(state, bond_atom_pair, start_idx, end_idx,
                              cursor_x, cursor_y, max_y):
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
        max_y: Screen height
    """
    a1_idx, a2_idx = bond_atom_pair

    # Determine which atom is closer to cursor
    conf = state.mol.GetConformer()
    mol_x, mol_y = screen_to_mol_coords(cursor_x, cursor_y, state.box,
                                        state.scale, max_y, state.y_offset)

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


def append_smiles_fragment(stdscr, state, cursor_x, cursor_y, max_x, max_y):
    """
    Handle the 'a' command: append atoms from SMILES to atom or bond under
    cursor.

    Returns State object if successful, None if no change was made.
    """
    # Find atom under cursor
    atom_idx = find_atom_at_cursor(state, cursor_x, cursor_y, max_y)

    # Check if on a bond if not on an atom
    bond_atom_pair = None
    if atom_idx is None:
        bond_atom_pair = find_bond_atoms(state, cursor_x, cursor_y, max_y)

    if atom_idx is None and bond_atom_pair is None:
        return None  # Nothing to do

    # Prompt for SMILES
    sidechain_smiles = enter_smiles(stdscr, max_y)
    if sidechain_smiles is None:
        return None  # Nothing to do

    try:
        # Create sidechain molecule
        sidechain = Chem.MolFromSmiles(sidechain_smiles)
        if sidechain is None:
            return None  # Bad SMILES

        # Get the index where sidechain will start
        start_idx = state.mol.GetNumAtoms()
        end_idx = start_idx + sidechain.GetNumAtoms() - 1

        # Insert sidechain into molecule
        state.mol.InsertMol(sidechain)

        if atom_idx is not None:
            # Cursor on atom: connect to first atom of sidechain
            state.mol.AddBond(atom_idx, start_idx, Chem.BondType.SINGLE)
        else:
            # Cursor on bond: insert sidechain between the two atoms
            connect_sidechain_to_bond(state, bond_atom_pair, start_idx, end_idx,
                                      cursor_x, cursor_y, max_y)

        # Compute 2D coordinates for new atoms, keeping original atoms fixed
        compute_coords_with_fixed_atoms(state.mol, start_idx)

        # Update box to show all atoms while keeping same scale
        box, y_offset = recalculate_box_and_offset(state.mol, state.scale,
                                                   max_x, max_y)

        return State(mol=state.mol,
                     box=box,
                     scale=state.scale,
                     y_offset=y_offset)
    except Exception as e:
        logging.exception("Error appending atoms (a command)")
        # Error appending atoms
        return None


def draw_selection_rect(stdscr, x1, y1, x2, y2, max_x, max_y):
    """Draw a selection rectangle on the screen."""
    # Normalize coordinates
    min_x, max_x_rect = (x1, x2) if x1 <= x2 else (x2, x1)
    min_y, max_y_rect = (y1, y2) if y1 <= y2 else (y2, y1)

    # Draw rectangle using box drawing characters or simple ASCII
    try:
        # Draw corners
        stdscr.addch(min_y, min_x, '+')
        stdscr.addch(min_y, max_x_rect, '+')
        stdscr.addch(max_y_rect, min_x, '+')
        stdscr.addch(max_y_rect, max_x_rect, '+')

        # Draw horizontal lines
        for x in range(min_x + 1, max_x_rect):
            if x < max_x:
                stdscr.addch(min_y, x, '-')
                stdscr.addch(max_y_rect, x, '-')

        # Draw vertical lines
        for y in range(min_y + 1, max_y_rect):
            if y < max_y:
                stdscr.addch(y, min_x, '|')
                stdscr.addch(y, max_x_rect, '|')
    except curses.error:
        pass


def redraw_screen(stdscr,
                  state,
                  show_smiles,
                  max_x,
                  max_y,
                  selection_mode=False,
                  selection_anchor_x=None,
                  selection_anchor_y=None,
                  cursor_x=None,
                  cursor_y=None):
    """
    Redraw the entire screen with molecule, SMILES, and optional selection.
    """
    stdscr.clear()

    # Draw molecule if present
    draw_mol(stdscr, state, max_y)

    # Draw SMILES at the top if enabled (after molecule so it's on top)
    if show_smiles:
        current_smiles = get_smiles(state.mol)
        # Wrap SMILES to screen width
        row = 0
        for i in range(0, len(current_smiles), max_x - 1):
            chunk = current_smiles[i:i + max_x - 1]
            try:
                stdscr.addstr(row, 0, chunk)
                row += 1
            except curses.error:
                break

    # Draw selection rectangle if in selection mode
    if (selection_mode and selection_anchor_x is not None and
            cursor_x is not None):
        draw_selection_rect(stdscr, selection_anchor_x, selection_anchor_y,
                            cursor_x, cursor_y, max_x, max_y)

    # Draw instructions at the bottom
    for i, line in enumerate(INSTRUCTIONS):
        try:
            stdscr.addstr(max_y - len(INSTRUCTIONS) + i, 0, line[:max_x - 1])
        except curses.error:
            pass


def main_loop(stdscr, initial_smiles=None):
    # Initialize curses
    curses.curs_set(1)  # Show cursor
    curses.use_default_colors()  # Use terminal's default colors

    # Reduce escape key delay (default is 1000ms)
    # This makes Esc key more responsive in selection mode
    try:
        curses.set_escdelay(25)  # 25ms is usually sufficient
    except AttributeError:
        # set_escdelay() not available (Python < 3.9)
        # Can set ESCDELAY environment variable before running instead
        pass

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
    show_smiles = False

    # Load initial molecule if provided
    if initial_smiles:
        state = create_molecule_from_smiles(initial_smiles, max_x, max_y)
        if state is None:
            # Failed to parse SMILES, create empty state
            state = create_empty_state(max_y)
    else:
        # Create empty state
        state = create_empty_state(max_y)

    # Undo/redo history
    history = UndoHistory(state)

    # Track when we need to redraw the entire screen
    need_redraw = True

    # Selection mode state
    selection_mode = False
    selection_anchor_x = None
    selection_anchor_y = None

    while True:
        # Only redraw everything when necessary
        if need_redraw:
            redraw_screen(stdscr, state, show_smiles, max_x, max_y,
                          selection_mode, selection_anchor_x,
                          selection_anchor_y, cursor_x, cursor_y)
            need_redraw = False

        # Move cursor to current position
        try:
            stdscr.move(cursor_y, cursor_x)
        except curses.error:
            pass

        stdscr.refresh()

        # Get user input
        key_code = stdscr.getch()

        # Handle terminal resize
        if key_code == curses.KEY_RESIZE:
            max_y, max_x = stdscr.getmaxyx()
            # Recalculate molecule position for new screen size
            state.box, state.y_offset = recalculate_box_and_offset(
                state.mol, state.scale, max_x, max_y)
            # Clamp cursor to new bounds
            cursor_x = min(cursor_x, max_x - 1)
            cursor_y = min(cursor_y, max_y - 1)
            need_redraw = True
            continue

        # Convert to character (will skip non-char keys)
        try:
            key = chr(key_code)
        except (ValueError, OverflowError):
            # Ignore other special keys we don't handle
            continue

        # Handle movement (vi-style) - works in both normal and selection mode
        if key in 'hjkl':
            if key == 'h':  # left
                cursor_x = max(0, cursor_x - 1)
            elif key == 'j':  # down
                cursor_y = min(max_y - 1, cursor_y + 1)
            elif key == 'k':  # up
                cursor_y = max(0, cursor_y - 1)
            elif key == 'l':  # right
                cursor_x = min(max_x - 1, cursor_x + 1)
            if selection_mode:
                need_redraw = True

        # Special handling for selection mode
        elif selection_mode:
            if key in '\r\n':  # Enter
                # Delete atoms in selection
                if delete_atoms_in_rect(state, selection_anchor_x,
                                        selection_anchor_y, cursor_x, cursor_y,
                                        max_y):
                    history.state = state
                    history.push()
                selection_mode = False
                selection_anchor_x = None
                selection_anchor_y = None
                need_redraw = True
            elif key == '\x1b':  # Escape
                # Cancel selection
                selection_mode = False
                selection_anchor_x = None
                selection_anchor_y = None
                need_redraw = True
            elif key == 'q':
                # Allow quitting from selection mode
                return get_smiles(state.mol)
            # Ignore all other keys in selection mode
            continue

        # Shift molecule (move all atoms)
        elif key in 'HJKL':
            if key == 'H':  # shift left
                shift_view(state, 1, 0)
            elif key == 'J':  # shift down
                shift_view(state, 0, 1)
            elif key == 'K':  # shift up
                shift_view(state, 0, -1)
            elif key == 'L':  # shift right
                shift_view(state, -1, 0)
            need_redraw = True

        # Enter SMILES string
        elif key == 's':
            result = load_smiles(stdscr, max_x, max_y)
            if result is not None:
                state = result
                history.state = state
                history.push()

            # Always redraw to clear the prompt
            need_redraw = True

        # Toggle SMILES display
        elif key == 'S':
            show_smiles = not show_smiles
            need_redraw = True

        # Insert atom at cursor position or change atom symbol
        elif key == 'i':
            result = insert_or_modify_atom(stdscr, state, cursor_x,
                                           cursor_y, max_y)
            if result is not None:
                state.mol = result
                history.state = state
                history.push()

            # Always redraw to clear the prompt
            need_redraw = True

        # Insert common atoms (c, n, o) - shortcuts, or change atom symbol
        elif key in ['c', 'n', 'o']:
            symbol = key.upper()
            result = insert_or_modify_atom(stdscr, state, cursor_x,
                                           cursor_y, max_y, symbol)
            if result is not None:
                state.mol = result
                history.state = state
                history.push()
                need_redraw = True

        # Append atoms from SMILES to atom under cursor or bond
        elif key == 'a':
            result = append_smiles_fragment(stdscr, state, cursor_x,
                                            cursor_y, max_x, max_y)
            if result is not None:
                state = result
                history.state = state
                history.push()

            # Always redraw to clear the prompt
            need_redraw = True

        # Delete atom or bond at cursor position
        elif key == 'x':
            if delete_at_cursor(state, cursor_x, cursor_y, max_y):
                history.state = state
                history.push()
                need_redraw = True

        # Enter area delete (selection) mode
        elif key == 'X':
            selection_mode = True
            selection_anchor_x = cursor_x
            selection_anchor_y = cursor_y
            need_redraw = True

        # Increase formal charge
        elif key == '+':
            if adjust_formal_charge(state, cursor_x, cursor_y, max_y, 1):
                history.state = state
                history.push()
                need_redraw = True

        # Decrease formal charge
        elif key == '-':
            if adjust_formal_charge(state, cursor_x, cursor_y, max_y, -1):
                history.state = state
                history.push()
                need_redraw = True

        # Cleanup/regenerate coordinates (Ctrl-L)
        elif key == '\x0c':  # Ctrl-L
            if cleanup_coordinates(state, max_x, max_y):
                history.state = state
                history.push()
                need_redraw = True

        # Zoom out
        elif key == '<':
            zoom_view(state, max_x, max_y, 1.0 / 1.2)
            need_redraw = True

        # Zoom in
        elif key == '>':
            zoom_view(state, max_x, max_y, 1.2)
            need_redraw = True

        # Add/modify/delete bond
        elif key in ['1', '2', '3']:
            bond_order = int(key)
            if create_or_adjust_bond(state, cursor_x, cursor_y, max_y,
                                     bond_order):
                history.state = state
                history.push()
                need_redraw = True

        # Add/modify wedge bond (single bond with up stereochemistry)
        elif key == 'w':
            if create_or_adjust_bond(state, cursor_x, cursor_y, max_y, 1,
                                     Chem.BondDir.BEGINWEDGE):
                history.state = state
                history.push()
                need_redraw = True

        # Add/modify dash bond (single bond with down stereochemistry)
        elif key == 'd':
            if create_or_adjust_bond(state, cursor_x, cursor_y, max_y, 1,
                                     Chem.BondDir.BEGINDASH):
                history.state = state
                history.push()
                need_redraw = True

        # Clear canvas (reset to blank slate)
        elif key == '@':
            clear_canvas(state, max_y)
            history.state = state
            history.push()
            need_redraw = True

        # Undo
        elif key == 'u':
            if history.undo():
                state = history.state
                need_redraw = True

        # Redo
        elif key == 'r':
            if history.redo():
                state = history.state
                need_redraw = True

        # Help
        elif key == '?':
            show_help(stdscr, max_x)
            need_redraw = True

        # Quit
        elif key == 'q':
            return get_smiles(state.mol)


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
        description='CurseMol - molecular sketcher for the terminally insane')
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
