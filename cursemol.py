#!/usr/bin/env python3
"""
Cursemol - A simple curses-based program for placing markers on screen.

Controls:
  h, j, k, l - Move cursor left, down, up, right (vi-style)
  x          - Place an 'x' at current cursor position
  c          - Show coordinates of all x'es
  s          - Enter a SMILES string
  q          - Quit
"""

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


def show_coordinates(stdscr, x_positions, max_y):
    """Display coordinates of all X markers."""
    stdscr.clear()
    stdscr.addstr(0, 0, "Coordinates of all X markers:")
    stdscr.addstr(1, 0, "-" * 40)

    if x_positions:
        for i, (y, x) in enumerate(sorted(x_positions.keys()), start=2):
            coord_str = f"X #{i-1}: (y={y}, x={x})"
            try:
                stdscr.addstr(i, 0, coord_str)
            except curses.error:
                break
    else:
        stdscr.addstr(2, 0, "No X markers placed yet.")

    stdscr.addstr(max_y - 1, 0, "Press any key to continue...")
    stdscr.refresh()
    stdscr.getch()


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

def draw_mol(stdscr, mol, max_x):
    """Draw the molecule using ASCII art"""
    Chem.Kekulize(mol)
    conf = mol.GetConformer(0)
    box = get_box(conf)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = box

    xscale = min((max_x-PADDING*2)/(xmax-xmin), MAX_SCALE)
    yscale = xscale*ASPECT_RATIO
    scale = (xscale, yscale)
    rows = int(PADDING*2 + (ymax-ymin) * yscale)
    cols = int(round(xscale * (xmax-xmin) + 2 * PADDING))
    screen = [[' '] * cols for i in range(rows)]

    for bond in mol.GetBonds():
        x1, y1 = int_coords_for_atom(bond.GetBeginAtom(), box, scale, conf)
        x2, y2 = int_coords_for_atom(bond.GetEndAtom(), box, scale, conf)
        draw_line(screen, BOND_CHARS[bond.GetBondType()], x1, y1, x2, y2)

    for atom in mol.GetAtoms():
        x, y = int_coords_for_atom(atom, box, scale, conf)
        sym = atom.GetSymbol()
        for i, c in enumerate(sym):
            screen[y][x+i] = c

    for i, line in enumerate(reversed(screen)):
        try:
            stdscr.addstr(i, 0, ''.join(line))
        except curses.error:
            pass

def main(stdscr):
    # Initialize curses
    curses.curs_set(1)  # Show cursor
    stdscr.clear()

    # Get screen dimensions
    max_y, max_x = stdscr.getmaxyx()

    # Starting cursor position
    cursor_y, cursor_x = 0, 0

    # Dictionary to store positions of x'es: {(y, x): 'x'}
    x_positions = {}

    # SMILES string storage
    smiles = ""
    mol = None

    # Instructions
    instructions = [
        "Cursemol - Place markers on screen",
        "h/j/k/l: move | x: place X | c: show coords | s: SMILES | q: quit"
    ]

    # Track when we need to redraw the entire screen
    need_redraw = True

    while True:
        # Only redraw everything when necessary
        if need_redraw:
            stdscr.clear()

            # Draw molecule if present
            if mol is not None:
                draw_mol(stdscr, mol, max_x)

            # Draw instructions at the bottom
            for i, line in enumerate(instructions):
                try:
                    stdscr.addstr(max_y - len(instructions) + i, 0, line[:max_x-1])
                except curses.error:
                    pass

            # Draw all placed x'es
            for (y, x), char in x_positions.items():
                try:
                    stdscr.addstr(y, x, char)
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

        # Place an 'x' at current position
        elif key == ord('x'):
            x_positions[(cursor_y, cursor_x)] = 'X'
            try:
                stdscr.addstr(cursor_y, cursor_x, 'X')
            except curses.error:
                pass

        # Show coordinates of all x'es
        elif key == ord('c'):
            show_coordinates(stdscr, x_positions, max_y)
            need_redraw = True

        # Enter SMILES string
        elif key == ord('s'):
            smiles = enter_smiles(stdscr, max_y)
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    AllChem.Compute2DCoords(mol)
                    need_redraw = True

        # Quit
        elif key == ord('q'):
            break


if __name__ == "__main__":
    curses.wrapper(main)
