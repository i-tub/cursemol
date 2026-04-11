"""
Functions that "draw" molecules to a an abstract canvas, that is a 3D array of
characters.

This module does not use curses.
"""

import logging

from rdkit import Chem

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


def mol_y_to_screen_y(mol_y, box_min_y, scale_y, y_offset, rows):
    """
    Convert molecular Y coordinate to screen Y coordinate.

    The Y-axis is flipped: terminal Y-axis points down (0 at top),
    but molecular coordinates point up (higher Y = higher position).
    Therefore higher molecular Y -> lower screen Y (closer to top).
    """
    y_from_bottom = PADDING + y_offset + int((mol_y - box_min_y) * scale_y)
    return rows - 1 - y_from_bottom


def screen_coords_for_atom(atom, state, conf, rows):
    """Calculate screen coordinates for an atom."""
    pos = conf.GetAtomPosition(atom.GetIdx())
    x = PADDING + int((pos.x - state.box[0][0]) * state.scale[0])
    y = mol_y_to_screen_y(pos.y, state.box[0][1], state.scale[1],
                          state.y_offset, rows)
    return x, y


def iter_atom_screen_positions(state, screen_dims):
    """Yield (atom, screen_x, screen_y) for each atom in molecule."""
    if state.mol.GetNumAtoms() == 0:
        return

    conf = state.mol.GetConformer()
    for atom in state.mol.GetAtoms():
        x, y = screen_coords_for_atom(atom, state, conf, screen_dims.rows)
        yield atom, x, y


def draw_line(screen, char, char2, x1, y1, x2, y2):
    """
    Draw a line from (x1, y1) to (x2, y2) using `char` for the first half of the
    line and `char2` for the second half.
    """
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
            c = char if x >= mid else char2
        else:
            c = char if x < mid else char2
        if vertical:
            screen[x][y] = c
        else:
            screen[y][x] = c


def draw_string(screen, screen_colors, s, x, y, rows, cols, color, is_bold):
    """
    Draw a string into the screen buffer with the given color and bold flag.

    Args:
        screen: 2D character array
        screen_colors: 2D color/attribute array
        s: String to draw
        x, y: screen coordinates for the start of the string
        rows, cols: screen buffer dimensions
        color: color code for the string
        is_bold: whether to draw the string in bold
    """
    for i, c in enumerate(s):
        if 0 <= y < rows and 0 <= x + i < cols:
            screen[y][x + i] = c
            screen_colors[y][x + i] = color | (0x100 if is_bold else 0)


def draw_atom(screen, screen_colors, atom, x, y, rows, cols, state, conf):
    """
    Draw a single atom with its symbol and charge into the screen buffer.

    Args:
        screen: 2D character array
        screen_colors: 2D color/attribute array
        atom: RDKit atom object
        x, y: screen coordinates for the atom
        rows, cols: screen buffer dimensions
        state: State object with scale, box, y_offset
        conf: RDKit conformer object
    """
    sym = atom.GetSymbol()
    color = ELEMENT_COLORS.get(sym, 0)  # Get color for this element
    # Make all symbols bold except C
    is_bold = sym != 'C'

    # Draw atom symbol
    draw_string(screen, screen_colors, sym, x, y, rows, cols, color, is_bold)

    # Draw hydrogens if heteroatom
    h_str = ''
    h_on_left = False
    if sym in 'NOPS':
        atom.UpdatePropertyCache()
        h = atom.GetTotalNumHs()
        if h > 0:
            h_str = 'H'
            if h > 1:
                h_str += {2: '₂', 3: '₃', 4: '₄'}.get(h, str(h))

            # Check neighbor positions to determine hydrogen placement
            has_neighbor_on_left = False
            has_neighbor_on_right = False
            neighbors = atom.GetNeighbors()

            for neighbor in neighbors:
                neighbor_pos = conf.GetAtomPosition(neighbor.GetIdx())
                neighbor_x = PADDING + int(
                    (neighbor_pos.x - state.box[0][0]) * state.scale[0])
                if neighbor_x < x:
                    has_neighbor_on_left = True
                if neighbor_x > x:
                    has_neighbor_on_right = True

            # Draw hydrogens on the left if no neighbors on the left AND at
            # least one on the right
            if not has_neighbor_on_left and has_neighbor_on_right:
                h_x = x - len(h_str)
                h_on_left = True
            else:
                h_x = x + len(sym)

            draw_string(screen, screen_colors, h_str, h_x, y, rows, cols, color,
                        is_bold)

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
        # When H is on the left, charge goes to the right of the symbol only
        if h_on_left:
            charge_x = x + len(sym)
        else:
            charge_x = x + len(sym) + len(h_str)
        charge_y = y - 1

        draw_string(screen, screen_colors, charge_str, charge_x, charge_y, rows,
                    cols, color, is_bold)


def fill_screen_buffer(state, screen_dims):
    """
    Fill screen buffer with bonds and atoms.
    Returns (screen, screen_colors) tuple of 2D arrays.
    """
    # Calculate screen size
    rows = screen_dims.rows
    cols = 200  # Generous width

    screen = [[' '] * cols for i in range(rows)]
    screen_colors = [[0] * cols for i in range(rows)]  # 0 = default color

    conf = state.mol.GetConformer(0)

    # Draw bonds
    try:
        for bond in state.mol.GetBonds():
            x1, y1 = screen_coords_for_atom(bond.GetBeginAtom(), state, conf,
                                            rows)
            x2, y2 = screen_coords_for_atom(bond.GetEndAtom(), state, conf,
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
    for atom, x, y in iter_atom_screen_positions(state, screen_dims):
        draw_atom(screen, screen_colors, atom, x, y, rows, cols, state, conf)

    return screen, screen_colors
