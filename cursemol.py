#!/usr/bin/env python3
"""
Cursemol - A simple curses-based program for displaying molecules.

Controls:
  h, j, k, l - Move cursor left, down, up, right (vi-style)
  s          - Enter a SMILES string
  i          - Insert atom at cursor position
  q          - Quit
"""

import argparse
import curses
from rdkit import Chem
from rdkit.Chem import AllChem


MAX_SCALE = 10.0  # columns per angstrom
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


def int_coords_for_atom(atom, box, scale, conf):
    pos = conf.GetAtomPosition(atom.GetIdx())
    x = PADDING + int((pos.x - box[0][0]) * scale[0])
    y = PADDING + int((pos.y - box[0][1]) * scale[1])
    return x, y


def draw_line(screen, char, x1, y1, x2, y2):
    vertical = False
    if abs(x2-x1) < abs(y2-y1):
        x1, y1 = y1, x1
        x2, y2 = y2, x2
        vertical = True
    try:
        slope = 1.0*(y2-y1)/(x2-x1)
    except ZeroDivisionError:
        return

    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
    for x in range(x1+1, x2):
        y = int(round(y1 + slope * (x-x1)))
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
        symbol = ""
    finally:
        curses.noecho()

    return symbol


def screen_to_mol_coords(cursor_x, cursor_y, box, scale, max_y):
    """Convert cursor/terminal coordinates to molecule coordinates."""
    # Screen array has (max_y - 2) rows, displayed reversed
    rows = max_y - 2

    # Convert terminal position to screen array position
    # Terminal row 0 → screen row (rows-1), terminal row 1 → screen row (rows-2), etc.
    screen_y = rows - 1 - cursor_y

    # Reverse the coordinate transformation from int_coords_for_atom
    mol_x = (cursor_x - PADDING) / scale[0] + box[0][0]
    mol_y = (screen_y - PADDING) / scale[1] + box[0][1]

    return mol_x, mol_y

def calculate_box_and_scale(mol, max_x):
    """Calculate bounding box and scale for a molecule."""
    conf = mol.GetConformer(0)
    box = get_box(conf)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = box

    xscale = min((max_x-PADDING*2)/(xmax-xmin), MAX_SCALE)
    yscale = xscale*ASPECT_RATIO
    scale = (xscale, yscale)

    return box, scale


def draw_mol(stdscr, mol, box, scale, max_y):
    """Draw the molecule using ASCII art using the given box and scale."""
    Chem.Kekulize(mol)
    conf = mol.GetConformer(0)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = box

    # Calculate screen size - make it large enough to handle any atom position
    # Use a generous size to accommodate atoms added outside the original box
    rows = max_y - 2  # Leave room for instructions
    cols = 200  # Generous width

    screen = [[' '] * cols for i in range(rows)]

    for bond in mol.GetBonds():
        x1, y1 = int_coords_for_atom(bond.GetBeginAtom(), box, scale, conf)
        x2, y2 = int_coords_for_atom(bond.GetEndAtom(), box, scale, conf)
        draw_line(screen, BOND_CHARS[bond.GetBondType()], x1, y1, x2, y2)

    for atom in mol.GetAtoms():
        x, y = int_coords_for_atom(atom, box, scale, conf)
        sym = atom.GetSymbol()
        for i, c in enumerate(sym):
            # Bounds check to avoid IndexError
            if 0 <= y < rows and 0 <= x+i < cols:
                screen[y][x+i] = c

    for i, line in enumerate(reversed(screen)):
        try:
            stdscr.addstr(i, 0, ''.join(line))
        except curses.error:
            pass

def main_loop(stdscr, initial_smiles=None):
    # Initialize curses
    curses.curs_set(1)  # Show cursor
    stdscr.clear()

    # Get screen dimensions
    max_y, max_x = stdscr.getmaxyx()

    # Starting cursor position
    cursor_y, cursor_x = 0, 0

    # SMILES string storage
    smiles = initial_smiles or ""
    mol = None
    box = None
    scale = None

    # Load initial molecule if provided
    if initial_smiles:
        m = Chem.MolFromSmiles(initial_smiles)
        if m is not None:
            mol = Chem.RWMol(m)
            AllChem.Compute2DCoords(mol)
            box, scale = calculate_box_and_scale(mol, max_x)

    # Instructions
    instructions = [
        "Cursemol - Display molecules",
        "h/j/k/l: move | s: SMILES | i: insert atom | q: quit"
    ]

    # Track when we need to redraw the entire screen
    need_redraw = True

    while True:
        # Only redraw everything when necessary
        if need_redraw:
            stdscr.clear()

            # Draw molecule if present
            if mol is not None and box is not None and scale is not None:
                draw_mol(stdscr, mol, box, scale, max_y)

            # Draw instructions at the bottom
            for i, line in enumerate(instructions):
                try:
                    stdscr.addstr(max_y - len(instructions) + i, 0, line[:max_x-1])
                except curses.error:
                    pass

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

        # Enter SMILES string
        elif key == ord('s'):
            smiles = enter_smiles(stdscr, max_y)
            if smiles:
                m = Chem.MolFromSmiles(smiles)
                if m is not None:
                    mol = Chem.RWMol(m)
                    AllChem.Compute2DCoords(mol)
                    box, scale = calculate_box_and_scale(mol, max_x)
                    need_redraw = True

        # Insert atom at cursor position
        elif key == ord('i'):
            if mol is not None and box is not None and scale is not None:
                symbol = enter_element(stdscr, max_y)
                if symbol:
                    try:
                        # Add atom to molecule
                        atom_idx = mol.AddAtom(Chem.Atom(symbol))

                        # Convert cursor position to molecule coordinates
                        mol_x, mol_y = screen_to_mol_coords(cursor_x, cursor_y, box, scale, max_y)

                        # Set atom position in conformer
                        conf = mol.GetConformer()
                        conf.SetAtomPosition(atom_idx, [mol_x, mol_y, 0.0])

                        need_redraw = True
                    except Exception:
                        # Invalid element symbol or other error
                        pass

        # Quit
        elif key == ord('q'):
            break


def main():
    parser = argparse.ArgumentParser(description='Display molecules in the terminal')
    parser.add_argument('smiles', nargs='?', help='Initial SMILES string to display')
    args = parser.parse_args()

    curses.wrapper(main_loop, args.smiles)


if __name__ == "__main__":
    main()
