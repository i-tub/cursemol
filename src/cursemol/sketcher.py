"""
Implementation of the sketcher widget. Main loop and functions implementing
various commands.
"""

import curses
import logging

from rdkit import Chem

from . import canvas
from . import config
from . import chem
from . import edit
from . import ui
from .state import Mode
from .state import ScreenDimensions
from .state import State
from .state import UndoHistory
from .state import recalculate_box_and_offset

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
        box, y_offset = recalculate_box_and_offset(state.mol, state.scale,
                                                   screen_dims)

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

        # Most commands need redraw; the few that don't will set it to False.
        need_redraw = True

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
            continue

        # Handle terminal resize
        if key_code == curses.KEY_RESIZE:
            max_y, max_x = stdscr.getmaxyx()
            screen_dims = ScreenDimensions(max_x=max_x, max_y=max_y)
            # Recalculate molecule position for new screen size
            state.box, state.y_offset = recalculate_box_and_offset(
                state.mol, state.scale, screen_dims)
            # Clamp cursor to new bounds
            cursor_x = min(cursor_x, screen_dims.max_x - 1)
            cursor_y = min(cursor_y, screen_dims.rows - 1)
            continue

        # Convert to character (will skip non-char keys)
        try:
            key = ARROW_KEY_MAP.get(key_code) or chr(key_code)
        except (ValueError, OverflowError):
            # Ignore other special keys we don't handle
            need_redraw = False
            continue

        # Handle movement (vi-style)
        if key in 'hjklHJKL':
            # Determine delta (1 for hjkl, 10 for HJKL)
            delta = 10 if key.isupper() else 1
            delta_y = int(delta *
                          config.ASPECT_RATIO) if key.isupper() else delta

            key_lower = key.lower()

            if mode == Mode.MOVE:
                # In move mode, hjkl translates the molecule
                if key_lower == 'h':  # shift left
                    canvas.shift_view(state, delta, 0)
                elif key_lower == 'j':  # shift down
                    canvas.shift_view(state, 0, delta_y)
                elif key_lower == 'k':  # shift up
                    canvas.shift_view(state, 0, -delta_y)
                elif key_lower == 'l':  # shift right
                    canvas.shift_view(state, -delta, 0)
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

                if mode == Mode.NORMAL:
                    need_redraw = False

        # Special handling for move mode
        elif mode == Mode.MOVE:
            if key in '\x1b\n':  # Escape or Enter
                # Exit move mode
                mode = Mode.NORMAL
            elif key == 'q':
                # Allow quitting from move mode
                return chem.get_smiles(state.mol)
            # Ignore all other keys in move mode
            need_redraw = False

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
            elif key == '\x1b':  # Escape
                # Cancel selection
                mode = Mode.NORMAL
                selection_anchor_x = None
                selection_anchor_y = None
            elif key == 'q':
                # Allow quitting from selection mode
                return chem.get_smiles(state.mol)
            else:
                # Ignore all other keys in selection mode
                need_redraw = False

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
            elif key == '\x1b':  # Escape - cancel bond
                # Delete new atom and exit bond mode (no undo entry)
                if bond_new_atom_idx is not None:
                    state.mol.RemoveAtom(bond_new_atom_idx)
                mode = Mode.NORMAL
                bond_start_atom_idx = None
                bond_new_atom_idx = None
            elif key == 'q':
                # Allow quitting from bond mode (clean up first)
                if bond_new_atom_idx is not None:
                    state.mol.RemoveAtom(bond_new_atom_idx)
                return chem.get_smiles(state.mol)
            else:
                # Ignore all other keys in bond mode
                need_redraw = False

        # Enter move mode
        elif key == 'm':
            mode = Mode.MOVE

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

        # Enter SMILES string
        elif key == 's':
            result, error_message = load_smiles(stdscr, screen_dims)
            if result is not None:
                state = result
                history.push(state)

        # Toggle SMILES display
        elif key == 'S':
            show_smiles = not show_smiles

        # Insert atom at cursor position or change atom symbol
        elif key == 'i':
            symbol = ui.enter_element(stdscr, screen_dims.max_y)
            result = edit.insert_or_modify_atom(state, cursor_x, cursor_y,
                                                screen_dims, symbol)
            if result is not None:
                state.mol = result
                history.push(state)

        # Insert common atoms (c, n, o) - shortcuts, or change atom symbol
        elif key in ['c', 'n', 'o']:
            symbol = key.upper()
            result = edit.insert_or_modify_atom(state, cursor_x, cursor_y,
                                                screen_dims, symbol)
            if result is not None:
                state.mol = result
                history.push(state)

        # Append atoms from SMILES to atom under cursor or bond
        elif key == 'a':
            result, error_message = append_smiles_fragment(
                stdscr, state, cursor_x, cursor_y, screen_dims)
            if result is not None:
                state = result
                history.push(state)

        # Delete atom or bond at cursor position
        elif key == 'x':
            if edit.delete_at_cursor(state, cursor_x, cursor_y, screen_dims):
                history.push(state)

        # Delete fragment (all atoms connected to cursor atom)
        elif key == 'D':
            if edit.delete_fragment_at_cursor(state, cursor_x, cursor_y,
                                              screen_dims):
                history.push(state)

        # Enter area delete (selection) mode
        elif key == 'X':
            mode = Mode.SELECT
            selection_anchor_x = cursor_x
            selection_anchor_y = cursor_y

        # Adjust formal charge
        elif key in '+-':
            dq = 1 if key == '+' else -1
            if edit.adjust_formal_charge(state, cursor_x, cursor_y, screen_dims,
                                         dq):
                history.push(state)

        # Cleanup/regenerate coordinates (Ctrl-L)
        elif key == '\x0c':  # Ctrl-L
            if edit.cleanup_coordinates(state, screen_dims):
                history.push(state)

        # Zoom
        elif key in '<>':
            zoom = config.ZOOM_STEP if key == '>' else 1.0 / config.ZOOM_STEP
            canvas.zoom_view(state, screen_dims, zoom)

        # Add/modify/delete bond
        elif key in ['1', '2', '3']:
            bond_order = int(key)
            if edit.create_or_adjust_bond(state, cursor_x, cursor_y,
                                          screen_dims, bond_order):
                history.push(state)

        # Add/modify wedge bond (single bond with up stereochemistry)
        elif key == 'w':
            if edit.create_or_adjust_bond(state, cursor_x, cursor_y,
                                          screen_dims, 1,
                                          Chem.BondDir.BEGINWEDGE):
                history.push(state)

        # Add/modify dash bond (single bond with down stereochemistry)
        elif key == 'd':
            if edit.create_or_adjust_bond(state, cursor_x, cursor_y,
                                          screen_dims, 1,
                                          Chem.BondDir.BEGINDASH):
                history.push(state)

        # Clear canvas (reset to blank slate)
        elif key == '@':
            edit.clear_canvas(state, screen_dims)
            history.push(state)

        # Undo
        elif key == 'u':
            if history.undo():
                state = history.state

        # Redo
        elif key in 'r\x12':  # Also support Ctrl-R for vim muscle memory
            if history.redo():
                state = history.state

        # Help
        elif key == '?':
            ui.show_help(stdscr, screen_dims)

        # Quit
        elif key == 'q':
            return chem.get_smiles(state.mol)

        else:
            # Ignore unknown keys
            need_redraw = False


def run(initial_smiles):
    return curses.wrapper(main_loop, initial_smiles)
