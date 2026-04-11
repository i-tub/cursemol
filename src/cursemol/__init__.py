"""
CurseMol - Molecular sketcher for the terminally committed

Controls:
  h, j, k, l       - Move cursor left/down/up/right (arrow keys supported too)
  H, J, K, L       - Move cursor faster (10 cells horizontal, 4 cells vertical)
  Space            - Snap cursor to nearest atom
  m                - Enter move mode (hjkl moves molecule, Esc to exit)
  s                - Enter a SMILES string to replace the current molecule
  S                - Toggle SMILES display
  i                - Insert/modify atom at cursor position
  a                - Append atoms from SMILES to atom or bond under cursor
                     (appending to a bond forms a ring by connecting the bond
                     atoms to the first and last atoms from the SMILES)
  c, n, o          - Insert/modify carbon/nitrogen/oxygen atom
  x                - Delete atom or bond
  D                - Delete fragment (all atoms connected to cursor atom)
  X                - Area delete (select rectangle)
  +, -             - Increase/decrease formal charge on atom
  <, >             - Zoom out/in
  b                - Add bond mode (add atom and move it, Enter to accept)
  1, 2, 3          - Add bond or change bond (order 1/2/3) between nearest atoms
  w, d             - Add/change to wedge or dash bond (press again to reverse)
  @                - Clear canvas (reset to blank slate)
  u, r             - Undo/redo
  Ctrl-L           - Clean up (regenerate coordinates)
  ?                - Show this help
  q                - Quit and print SMILES to stdout
"""

__version__ = "4.1.0"

import argparse
import atexit
import curses
import io
import logging
import os
import sys

from rdkit import Chem

from . import canvas
from . import chem
from . import edit
from . import ui
from .state import Mode
from .state import ScreenDimensions
from .state import State
from .state import UndoHistory

rdkit_logger = logging.getLogger('rdkit')

MIN_SCALE = 2.0  # columns per angstrom
DEFAULT_SCALE = 8.0  # columns per angstrom
MAX_SCALE = 16.0  # columns per angstrom
ASPECT_RATIO = 0.4  # horizontal / vertical
PADDING = 5
ZOOM_STEP = 1.2


# Training wheels for those who haven't converted to vi. :-)
ARROW_KEY_MAP = {
    curses.KEY_LEFT: 'h',
    curses.KEY_DOWN: 'j',
    curses.KEY_UP: 'k',
    curses.KEY_RIGHT: 'l',
}


def load_smiles(stdscr, screen_dims):
    """
    Prompt user for SMILES string and create molecule from it.
    Returns tuple: (State object or None, error_message string)
    Error message is empty string if no error occurred.
    """
    smiles = ui.enter_smiles(stdscr, screen_dims.max_y)
    if smiles:
        state, msg = State.fromSmiles(smiles, screen_dims)
        # Return (state, error_message) - error if state is None
        if state is None:
            return (None, msg)
        return (state, "")
    # User cancelled - no error
    return (None, "")


def clear_canvas(state, screen_dims):
    """
    Clear the canvas and reset to blank slate with default settings.
    Modifies state in place.
    """
    empty = State.createEmpty(screen_dims)
    state.mol = empty.mol
    state.box = empty.box
    state.scale = empty.scale
    state.y_offset = empty.y_offset


def show_help(stdscr, screen_dims):
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
            stdscr.addstr(i, 0, line[:screen_dims.max_x - 1])
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


def zoom_view(state, screen_dims, zoom_factor):
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
    screen_width = screen_dims.max_x - 2 * PADDING
    screen_height = screen_dims.rows - 2 * PADDING

    # Molecule coordinate range that fits on screen at new scale
    mol_width = screen_width / xscale
    mol_height = screen_height / yscale

    # New box centered on the same point
    state.box = ((center_x - mol_width / 2, center_y - mol_height / 2, 0.0),
                 (center_x + mol_width / 2, center_y + mol_height / 2, 0.0))

    # Recalculate y_offset for new scale
    mol_display_height = int(mol_height * yscale + 2 * PADDING)
    state.y_offset = max(0, (screen_dims.rows - mol_display_height) // 2)


def append_smiles_fragment(stdscr, state, cursor_x, cursor_y, screen_dims):
    """
    Handle the 'a' command: append atoms from SMILES to atom or bond under
    cursor.

    Returns State object if successful, None if no change was made.
    """
    # Find atom under cursor
    atom_idx = canvas.find_atom_at_cursor(state, cursor_x, cursor_y,
                                          screen_dims)

    # Check if on a bond if not on an atom
    bond_atom_pair = None
    if atom_idx is None:
        bond_atom_pair = canvas.find_bond_atoms(state, cursor_x, cursor_y,
                                                screen_dims)

    if atom_idx is None and bond_atom_pair is None:
        return None, ""  # Nothing to do

    # Prompt for SMILES
    sidechain_smiles = ui.enter_smiles(stdscr, screen_dims.max_y)
    if sidechain_smiles is None:
        return None, ""  # Nothing to do

    try:
        # Create sidechain molecule
        with chem.CaptureRDKitLog() as log:
            sidechain = Chem.MolFromSmiles(sidechain_smiles)
        if sidechain is None:
            return None, log.getMessage()  # Bad SMILES

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
            edit.connect_sidechain_to_bond(state, bond_atom_pair, start_idx,
                                           end_idx, cursor_x, cursor_y,
                                           screen_dims)

        # Compute 2D coordinates for new atoms, keeping original atoms fixed
        chem.compute_coords_with_fixed_atoms(state.mol, start_idx)

        # Update box to show all atoms while keeping same scale
        box, y_offset = canvas.recalculate_box_and_offset(
            state.mol, state.scale, screen_dims)

        return State(mol=state.mol,
                     box=box,
                     scale=state.scale,
                     y_offset=y_offset), ""
    except Exception as e:
        logging.exception("Error appending atoms (a command)")
        # Error appending atoms
        return None, str(e)


def main_loop(stdscr, initial_smiles=None):
    ui.init_curses(stdscr)

    # Get screen dimensions
    max_y, max_x = stdscr.getmaxyx()
    screen_dims = ScreenDimensions(max_x=max_x, max_y=max_y)

    # Starting cursor position (center of screen)
    cursor_y, cursor_x = screen_dims.rows // 2, screen_dims.max_x // 2

    # SMILES string storage
    show_smiles = False

    # Error message to display (empty string means no error)
    error_message = ""

    # Load initial molecule if provided
    if initial_smiles:
        state, error_message = State.fromSmiles(initial_smiles, screen_dims)
        if state is None:
            # Failed to parse SMILES, create empty state
            state = State.createEmpty(screen_dims)
    else:
        # Create empty state
        state = State.createEmpty(screen_dims)

    # Undo/redo history
    history = UndoHistory(state)

    # Track when we need to redraw the entire screen
    need_redraw = True

    # UI mode and selection state
    mode = Mode.NORMAL
    selection_anchor_x = None
    selection_anchor_y = None

    # Bond mode state
    bond_start_atom_idx = None
    bond_new_atom_idx = None

    while True:
        # Only redraw everything when necessary
        if need_redraw:
            ui.redraw_screen(stdscr, state, show_smiles, screen_dims, mode,
                             selection_anchor_x, selection_anchor_y, cursor_x,
                             cursor_y, error_message)
            need_redraw = False

        # Move cursor to current position
        try:
            stdscr.move(cursor_y, cursor_x)
        except curses.error:
            pass

        stdscr.refresh()

        # Get user input
        key_code = stdscr.getch()

        # Clear error message on any key press
        if error_message:
            error_message = ""
            need_redraw = True
            continue

        # Handle terminal resize
        if key_code == curses.KEY_RESIZE:
            max_y, max_x = stdscr.getmaxyx()
            screen_dims = ScreenDimensions(max_x=max_x, max_y=max_y)
            # Recalculate molecule position for new screen size
            state.box, state.y_offset = canvas.recalculate_box_and_offset(
                state.mol, state.scale, screen_dims)
            # Clamp cursor to new bounds
            cursor_x = min(cursor_x, screen_dims.max_x - 1)
            cursor_y = min(cursor_y, screen_dims.rows - 1)
            need_redraw = True
            continue

        # Convert to character (will skip non-char keys)
        try:
            key = ARROW_KEY_MAP.get(key_code) or chr(key_code)
        except (ValueError, OverflowError):
            # Ignore other special keys we don't handle
            continue

        # Handle movement (vi-style)
        if key in 'hjklHJKL':
            # Determine delta (1 for hjkl, 10 for HJKL)
            delta = 10 if key.isupper() else 1
            delta_y = int(delta * ASPECT_RATIO) if key.isupper() else delta

            key_lower = key.lower()

            if mode == Mode.MOVE:
                # In move mode, hjkl translates the molecule
                if key_lower == 'h':  # shift left
                    shift_view(state, delta, 0)
                elif key_lower == 'j':  # shift down
                    shift_view(state, 0, delta_y)
                elif key_lower == 'k':  # shift up
                    shift_view(state, 0, -delta_y)
                elif key_lower == 'l':  # shift right
                    shift_view(state, -delta, 0)
                need_redraw = True
            else:
                # Normal/select/bond mode: move cursor
                if key_lower == 'h':  # left
                    cursor_x = max(0, cursor_x - delta)
                elif key_lower == 'j':  # down
                    cursor_y = min(screen_dims.rows - 1, cursor_y + delta_y)
                elif key_lower == 'k':  # up
                    cursor_y = max(0, cursor_y - delta_y)
                elif key_lower == 'l':  # right
                    cursor_x = min(screen_dims.max_x - 1, cursor_x + delta)

                # Update new atom position in bond mode
                if mode == Mode.BOND and bond_new_atom_idx is not None:
                    mol_x, mol_y = canvas.screen_to_mol_coords(
                        cursor_x, cursor_y, state.box, state.scale, screen_dims,
                        state.y_offset)
                    conf = state.mol.GetConformer()
                    conf.SetAtomPosition(bond_new_atom_idx, [mol_x, mol_y, 0.0])

                if mode in (Mode.SELECT, Mode.BOND):
                    need_redraw = True

        # Special handling for move mode
        elif mode == Mode.MOVE:
            if key in '\x1b\n':  # Escape or Enter
                # Exit move mode
                mode = Mode.NORMAL
                need_redraw = True
            elif key == 'q':
                # Allow quitting from move mode
                return chem.get_smiles(state.mol)
            # Ignore all other keys in move mode
            continue

        # Special handling for selection mode
        elif mode == Mode.SELECT:
            if key in '\nx':  # Enter or x commits delete
                # Delete atoms in selection
                if edit.delete_atoms_in_rect(state, selection_anchor_x,
                                             selection_anchor_y, cursor_x,
                                             cursor_y, screen_dims):
                    history.push(state)
                mode = Mode.NORMAL
                selection_anchor_x = None
                selection_anchor_y = None
                need_redraw = True
            elif key == '\x1b':  # Escape
                # Cancel selection
                mode = Mode.NORMAL
                selection_anchor_x = None
                selection_anchor_y = None
                need_redraw = True
            elif key == 'q':
                # Allow quitting from selection mode
                return chem.get_smiles(state.mol)
            else:
                # Ignore all other keys in selection mode
                continue

        # Snap to nearest atom
        elif key == ' ':
            # In bond mode, exclude the new atom from snap search
            exclude_idx = bond_new_atom_idx if mode == Mode.BOND else None
            result = canvas.find_nearest_atom(state, cursor_x, cursor_y,
                                              screen_dims, exclude_idx)
            if result is not None:
                _, screen_x, screen_y = result
                cursor_x = screen_x
                cursor_y = screen_y

                # Update new atom position in bond mode
                if mode == Mode.BOND and bond_new_atom_idx is not None:
                    mol_x, mol_y = canvas.screen_to_mol_coords(
                        cursor_x, cursor_y, state.box, state.scale, screen_dims,
                        state.y_offset)
                    conf = state.mol.GetConformer()
                    conf.SetAtomPosition(bond_new_atom_idx, [mol_x, mol_y, 0.0])
                    need_redraw = True

        # Special handling for bond mode
        elif mode == Mode.BOND:
            if key == '\n':  # Enter - accept bond
                # Check if cursor is on another atom (excluding the new atom)
                target_atom_idx = canvas.find_atom_at_cursor(
                    state, cursor_x, cursor_y, screen_dims)
                if target_atom_idx is not None and target_atom_idx != bond_new_atom_idx:
                    # Delete new atom and form bond between start and target
                    state.mol.RemoveAtom(bond_new_atom_idx)
                    # Adjust indices if needed (if new atom was before target)
                    if bond_new_atom_idx < target_atom_idx:
                        target_atom_idx -= 1
                    adjusted_start_idx = bond_start_atom_idx
                    if bond_new_atom_idx < bond_start_atom_idx:
                        adjusted_start_idx -= 1
                    # Create bond (only if not self-loop)
                    if adjusted_start_idx != target_atom_idx:
                        chem.modify_bond(state.mol, adjusted_start_idx,
                                         target_atom_idx, 1)
                # else: leave new atom where it is

                # Push undo and exit bond mode
                history.push(state)
                mode = Mode.NORMAL
                bond_start_atom_idx = None
                bond_new_atom_idx = None
                need_redraw = True
            elif key == '\x1b':  # Escape - cancel bond
                # Delete new atom and exit bond mode (no undo entry)
                if bond_new_atom_idx is not None:
                    state.mol.RemoveAtom(bond_new_atom_idx)
                mode = Mode.NORMAL
                bond_start_atom_idx = None
                bond_new_atom_idx = None
                need_redraw = True
            elif key == 'q':
                # Allow quitting from bond mode (clean up first)
                if bond_new_atom_idx is not None:
                    state.mol.RemoveAtom(bond_new_atom_idx)
                return chem.get_smiles(state.mol)
            else:
                # Ignore all other keys in bond mode
                continue

        # Enter move mode
        elif key == 'm':
            mode = Mode.MOVE
            need_redraw = True

        # Enter bond mode
        elif key == 'b':
            # Find atom at cursor to start bond from
            start_idx = canvas.find_atom_at_cursor(state, cursor_x, cursor_y,
                                                   screen_dims)
            if start_idx is not None:
                # Add new C atom at cursor position
                bond_start_atom_idx = start_idx
                bond_new_atom_idx = state.mol.AddAtom(Chem.Atom('C'))

                # Set position of new atom
                mol_x, mol_y = canvas.screen_to_mol_coords(
                    cursor_x, cursor_y, state.box, state.scale, screen_dims,
                    state.y_offset)
                conf = state.mol.GetConformer()
                conf.SetAtomPosition(bond_new_atom_idx, [mol_x, mol_y, 0.0])

                # Add bond between start atom and new atom
                state.mol.AddBond(bond_start_atom_idx, bond_new_atom_idx,
                                  Chem.BondType.SINGLE)

                # Enter bond mode
                mode = Mode.BOND
                need_redraw = True

        # Enter SMILES string
        elif key == 's':
            result, error_message = load_smiles(stdscr, screen_dims)
            if result is not None:
                state = result
                history.push(state)

            # Always redraw to clear the prompt
            need_redraw = True

        # Toggle SMILES display
        elif key == 'S':
            show_smiles = not show_smiles
            need_redraw = True

        # Insert atom at cursor position or change atom symbol
        elif key == 'i':
            symbol = ui.enter_element(stdscr, screen_dims.max_y)
            result = edit.insert_or_modify_atom(state, cursor_x, cursor_y,
                                                screen_dims, symbol)
            if result is not None:
                state.mol = result
                history.push(state)

            # Always redraw to clear the prompt
            need_redraw = True

        # Insert common atoms (c, n, o) - shortcuts, or change atom symbol
        elif key in ['c', 'n', 'o']:
            symbol = key.upper()
            result = edit.insert_or_modify_atom(state, cursor_x, cursor_y,
                                                screen_dims, symbol)
            if result is not None:
                state.mol = result
                history.push(state)
                need_redraw = True

        # Append atoms from SMILES to atom under cursor or bond
        elif key == 'a':
            result, error_message = append_smiles_fragment(
                stdscr, state, cursor_x, cursor_y, screen_dims)
            if result is not None:
                state = result
                history.push(state)

            # Always redraw to clear the prompt
            need_redraw = True

        # Delete atom or bond at cursor position
        elif key == 'x':
            if edit.delete_at_cursor(state, cursor_x, cursor_y, screen_dims):
                history.push(state)
                need_redraw = True

        # Delete fragment (all atoms connected to cursor atom)
        elif key == 'D':
            if edit.delete_fragment_at_cursor(state, cursor_x, cursor_y,
                                              screen_dims):
                history.push(state)
                need_redraw = True

        # Enter area delete (selection) mode
        elif key == 'X':
            mode = Mode.SELECT
            selection_anchor_x = cursor_x
            selection_anchor_y = cursor_y
            need_redraw = True

        # Adjust formal charge
        elif key in '+-':
            dq = 1 if key == '+' else -1
            if edit.adjust_formal_charge(state, cursor_x, cursor_y, screen_dims,
                                         dq):
                history.push(state)
                need_redraw = True

        # Cleanup/regenerate coordinates (Ctrl-L)
        elif key == '\x0c':  # Ctrl-L
            if edit.cleanup_coordinates(state, screen_dims):
                history.push(state)
                need_redraw = True

        # Zoom
        elif key in '<>':
            zoom = ZOOM_STEP if key == '>' else 1.0 / ZOOM_STEP
            zoom_view(state, screen_dims, zoom)
            need_redraw = True

        # Add/modify/delete bond
        elif key in ['1', '2', '3']:
            bond_order = int(key)
            if edit.create_or_adjust_bond(state, cursor_x, cursor_y,
                                          screen_dims, bond_order):
                history.push(state)
                need_redraw = True

        # Add/modify wedge bond (single bond with up stereochemistry)
        elif key == 'w':
            if edit.create_or_adjust_bond(state, cursor_x, cursor_y,
                                          screen_dims, 1,
                                          Chem.BondDir.BEGINWEDGE):
                history.push(state)
                need_redraw = True

        # Add/modify dash bond (single bond with down stereochemistry)
        elif key == 'd':
            if edit.create_or_adjust_bond(state, cursor_x, cursor_y,
                                          screen_dims, 1,
                                          Chem.BondDir.BEGINDASH):
                history.push(state)
                need_redraw = True

        # Clear canvas (reset to blank slate)
        elif key == '@':
            clear_canvas(state, screen_dims)
            history.push(state)
            need_redraw = True

        # Undo
        elif key == 'u':
            if history.undo():
                state = history.state
                need_redraw = True

        # Redo
        elif key in 'r\x12':  # Also support Ctrl-R for vim muscle memory
            if history.redo():
                state = history.state
                need_redraw = True

        # Help
        elif key == '?':
            show_help(stdscr, screen_dims)
            need_redraw = True

        # Quit
        elif key == 'q':
            return chem.get_smiles(state.mol)


def parse_args():
    parser = argparse.ArgumentParser(
        description='CurseMol - molecular sketcher for the terminally committed'
    )
    parser.add_argument(
        'smiles',
        nargs='?',
        help='Initial SMILES string to display (use "-" to read from stdin)')
    return parser.parse_args()


def setup_tty():
    """
    If stdin/stdout are not TTYs (e.g., piped input/output), redirect to /dev/tty
    so curses can read keyboard input and display to the terminal.
    Returns the original stdout fd if it was redirected (for final SMILES output).
    """
    original_stdout_fd = None

    # Handle stdin
    if not sys.stdin.isatty():
        sys.stdin.close()  # Close the old stdin to avoid resource warning
        tty_fd = os.open('/dev/tty', os.O_RDONLY)
        os.dup2(tty_fd, 0)  # Replace fd 0 (stdin) with /dev/tty
        os.close(tty_fd)
        sys.stdin = os.fdopen(0, 'r')
        # Register cleanup to avoid resource warning on exit
        atexit.register(lambda: sys.stdin.close())

    # Handle stdout
    if not sys.stdout.isatty():
        original_stdout_fd = os.dup(1)  # Duplicate fd 1 before redirecting
        tty_fd = os.open('/dev/tty', os.O_WRONLY)
        os.dup2(tty_fd, 1)  # Replace fd 1 (stdout) with /dev/tty
        os.close(tty_fd)
        sys.stdout = sys.__stdout__ = os.fdopen(
            1, 'w')  # Update both stdout and __stdout__

    return original_stdout_fd


def main():
    # Set up logging to file (truncate on start)
    logging.basicConfig(
        filename='cursemol.log',
        filemode='w',  # Truncate on open
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s')

    # Capture RDKit warnings
    rdkit_logger.setLevel(logging.ERROR)
    rdkit_logger.handlers[0].setStream(io.StringIO())
    Chem.rdBase.LogToPythonLogger()

    args = parse_args()

    # Handle reading from stdin if "-" is provided
    initial_smiles = args.smiles
    if initial_smiles == "-":
        initial_smiles = sys.stdin.readline().strip()

    original_stdout_fd = setup_tty()

    smiles = curses.wrapper(main_loop, initial_smiles)

    # Print to original stdout if it was redirected, otherwise to current stdout
    if original_stdout_fd is not None:
        output = os.fdopen(original_stdout_fd, 'w')
        print(smiles, file=output)
        output.close()
    else:
        print(smiles)


if __name__ == "__main__":
    main()
