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
from dataclasses import dataclass
from enum import Enum
import io
import logging
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

from . import canvas
from . import chem
from . import edit

rdkit_logger = logging.getLogger('rdkit')

MIN_SCALE = 2.0  # columns per angstrom
DEFAULT_SCALE = 8.0  # columns per angstrom
MAX_SCALE = 16.0  # columns per angstrom
ASPECT_RATIO = 0.4  # horizontal / vertical
PADDING = 5
ZOOM_STEP = 1.2


class Mode(Enum):
    """UI mode for the main loop."""
    NORMAL = "normal"
    MOVE = "move"
    SELECT = "select"
    BOND = "bond"


# Instructions (try to keep lines under 80 characters and more or less balanced)
INSTRUCTIONS = {
    Mode.NORMAL: [
        "hjkl: move | HJKL: fast | SPC: snap | m: move mol | s/S: SMILES | i/a/c/n/o: ins",
        "x/X/D: del | +/-: chg | <>: zoom | u/r: undo | ^L: clean | b/123/wd: bond | ?: help"
    ],
    Mode.MOVE: [
        "[Move molecule mode]",
        "hjkl/HJKL: move molecule | Esc/Enter: leave move mode | q: quit"
    ],
    Mode.SELECT: [
        "[Area delete mode]",
        "hjkl: move | HJKL: fast | Enter/x: delete | Esc: cancel | q: quit"
    ],
    Mode.BOND: [
        "[Add bond mode]", "hjkl/HJKL/SPC: move | Enter: accept | Esc: cancel"
    ]
}

# Number of instruction lines to reserve at the bottom of the screen.
INSTRUCTION_LINES = max(len(v) for v in INSTRUCTIONS.values())

# Training wheels for those who haven't converted to vi. :-)
ARROW_KEY_MAP = {
    curses.KEY_LEFT: 'h',
    curses.KEY_DOWN: 'j',
    curses.KEY_UP: 'k',
    curses.KEY_RIGHT: 'l',
}


@dataclass
class State:
    """Molecular drawing state: molecule and its display parameters."""
    mol: Chem.RWMol
    box: tuple  # ((min_x, min_y, min_z), (max_x, max_y, max_z))
    scale: tuple  # (xscale, yscale)
    y_offset: int

    def copy(self):
        """Create deep copy for undo/redo."""
        return State(mol=Chem.RWMol(self.mol),
                     box=self.box,
                     scale=self.scale,
                     y_offset=self.y_offset)


@dataclass
class ScreenDimensions:
    """Terminal screen dimensions and derived values."""
    max_x: int
    max_y: int

    @property
    def rows(self):
        """Drawable rows (excluding instruction lines)."""
        return self.max_y - INSTRUCTION_LINES


class UndoHistory:
    """Manages undo/redo history for molecule editing."""

    def __init__(self, state):
        self.state = state
        self._history = [state.copy()]
        self._index = 0

    def undo(self):
        """Move back in history. Returns True if successful."""
        if self._index > 0:
            self._index -= 1
            self.state = self._history[self._index].copy()
            return True
        return False

    def redo(self):
        """Move forward in history. Returns True if successful."""
        if self._index < len(self._history) - 1:
            self._index += 1
            self.state = self._history[self._index].copy()
            return True
        return False

    def push(self, state):
        """Truncate future history and save current state."""
        self.state = state
        self._history = self._history[:self._index + 1]
        self._history.append(self.state.copy())
        self._index = len(self._history) - 1


def prompt_user_input(stdscr, max_y, prompt_text):
    """Display prompt and get user input. Returns empty string on error."""
    stdscr.addstr(max_y - 1, 0, prompt_text)
    stdscr.clrtoeol()
    stdscr.refresh()

    curses.echo()
    try:
        input_bytes = stdscr.getstr(max_y - 1, len(prompt_text))
        return input_bytes.decode('utf-8').strip()
    except Exception:
        logging.exception(f"Error in prompt: {prompt_text}")
        return ""
    finally:
        curses.noecho()


def enter_smiles(stdscr, max_y):
    """Prompt user to enter a SMILES string and return it."""
    return prompt_user_input(stdscr, max_y, "Enter SMILES: ")


def enter_element(stdscr, max_y):
    """Prompt user to enter an element symbol and return it."""
    return prompt_user_input(stdscr, max_y, "Element symbol: ")


def recalculate_box_and_offset(mol, scale, screen_dims):
    """
    Recalculate box and y_offset for a molecule at a given scale.
    Centers the view on the molecule's actual bounding box.
    Returns (box, y_offset).
    """
    conf = mol.GetConformer(0)
    actual_box = chem.get_box(conf)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = actual_box

    # Calculate center of molecule
    center_x = (xmin + xmax) / 2
    center_y = (ymin + ymax) / 2

    # Calculate box dimensions that fill the screen at this scale
    screen_width = screen_dims.max_x - 2 * PADDING
    screen_height = screen_dims.rows - 2 * PADDING
    mol_width = screen_width / scale[0]
    mol_height = screen_height / scale[1]

    # Create box centered on molecule center
    box = ((center_x - mol_width / 2, center_y - mol_height / 2, 0.0),
           (center_x + mol_width / 2, center_y + mol_height / 2, 0.0))

    # Calculate vertical offset to center the displayed content
    mol_display_height = int(mol_height * scale[1] + 2 * PADDING)
    y_offset = max(0, (screen_dims.rows - mol_display_height) // 2)

    return box, y_offset


def calculate_box_and_scale(mol, max_x, max_y):
    """Calculate bounding box and scale for a molecule, centered on screen."""
    conf = mol.GetConformer(0)
    actual_box = chem.get_box(conf)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = actual_box

    # Calculate scale to fit the molecule
    xscale = min((max_x - PADDING * 2) / (xmax - xmin), DEFAULT_SCALE)
    yscale = xscale * ASPECT_RATIO
    scale = (xscale, yscale)

    # Recalculate box and offset at this scale
    screen_dims = ScreenDimensions(max_x=max_x, max_y=max_y)
    box, y_offset = recalculate_box_and_offset(mol, scale, screen_dims)

    return box, scale, y_offset


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


def draw_mol(stdscr, state, screen_dims):
    """Draw the molecule using ASCII art."""
    # Nothing to draw if molecule is empty
    if state.mol.GetNumAtoms() == 0:
        return

    # Fill screen buffer with molecular structure
    screen, screen_colors = canvas.fill_screen_buffer(state, screen_dims)

    # Render buffer to curses window
    render_screen_buffer(stdscr, screen, screen_colors)


def create_empty_state(screen_dims):
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
    y_offset = max(0, (screen_dims.rows - mol_height) // 2)

    return State(mol=mol, box=box, scale=scale, y_offset=y_offset)


def create_molecule_from_smiles(smiles, screen_dims):
    """
    Create molecule from SMILES string with 2D coordinates.
    Returns (state, error_message); in case of error, state will be None.
    """
    with chem.CaptureRDKitLog() as log:
        m = Chem.MolFromSmiles(smiles)

    if m is None:
        return None, 'Invalid SMILES\n' + log.getMessage()

    Chem.Kekulize(m, True)
    mol = Chem.RWMol(m)
    AllChem.Compute2DCoords(mol)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    box, scale, y_offset = calculate_box_and_scale(mol, screen_dims.max_x,
                                                   screen_dims.max_y)
    return State(mol=mol, box=box, scale=scale, y_offset=y_offset), ""


def load_smiles(stdscr, screen_dims):
    """
    Prompt user for SMILES string and create molecule from it.
    Returns tuple: (State object or None, error_message string)
    Error message is empty string if no error occurred.
    """
    smiles = enter_smiles(stdscr, screen_dims.max_y)
    if smiles:
        state, msg = create_molecule_from_smiles(smiles, screen_dims)
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
    empty = create_empty_state(screen_dims)
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
    sidechain_smiles = enter_smiles(stdscr, screen_dims.max_y)
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


def draw_selection_rect(stdscr, x1, y1, x2, y2, screen_dims):
    """Draw a selection rectangle on the screen."""
    # Normalize coordinates
    min_x, min_y, max_x_rect, max_y_rect = canvas.normalize_rect(x1, y1, x2, y2)

    # Draw rectangle using box drawing characters or simple ASCII
    try:
        # Draw corners
        stdscr.addch(min_y, min_x, '+')
        stdscr.addch(min_y, max_x_rect, '+')
        stdscr.addch(max_y_rect, min_x, '+')
        stdscr.addch(max_y_rect, max_x_rect, '+')

        # Draw horizontal lines
        for x in range(min_x + 1, max_x_rect):
            if x < screen_dims.max_x:
                stdscr.addch(min_y, x, '-')
                stdscr.addch(max_y_rect, x, '-')

        # Draw vertical lines
        for y in range(min_y + 1, max_y_rect):
            if y < screen_dims.max_y:
                stdscr.addch(y, min_x, '|')
                stdscr.addch(y, max_x_rect, '|')
    except curses.error:
        pass


def draw_instructions(stdscr, screen_dims, mode):
    """Draw instructions at the bottom of the screen."""
    instructions = INSTRUCTIONS[mode]

    for i, line in enumerate(instructions):
        try:
            stdscr.addstr(screen_dims.max_y - INSTRUCTION_LINES + i, 0,
                          line[:screen_dims.max_x - 1])
        except curses.error:
            pass


def draw_error_message(stdscr, screen_dims, error_message):
    """
    Draw error message at the bottom of the screen.

    Args:
        stdscr: curses screen object
        screen_dims: Screen dimensions
        error_message: Error message string (will be split into lines)
    """
    # Split message into lines
    error_lines = [
        *error_message.strip().split('\n'), '[Press any key to clear]'
    ]

    for i, line in enumerate(error_lines):
        # Calculate row position (last line of error is on last screen row)
        row = screen_dims.max_y - len(error_lines) + i
        display_line = line[:screen_dims.max_x - 1].ljust(screen_dims.max_x - 1)
        try:
            stdscr.addstr(row, 0, display_line, curses.color_pair(1))
        except curses.error:
            pass


def redraw_screen(stdscr,
                  state,
                  show_smiles,
                  screen_dims,
                  mode,
                  selection_anchor_x=None,
                  selection_anchor_y=None,
                  cursor_x=None,
                  cursor_y=None,
                  error_message=""):
    """
    Redraw the entire screen with molecule, SMILES, and optional selection.

    Args:
        error_message: Optional error message string to display at bottom
    """
    stdscr.clear()

    # Draw molecule if present
    draw_mol(stdscr, state, screen_dims)

    # Draw SMILES at the top if enabled (after molecule so it's on top)
    if show_smiles:
        current_smiles = chem.get_smiles(state.mol)
        # Wrap SMILES to screen width
        row = 0
        for i in range(0, len(current_smiles), screen_dims.max_x - 1):
            chunk = current_smiles[i:i + screen_dims.max_x - 1]
            try:
                stdscr.addstr(row, 0, chunk)
                row += 1
            except curses.error:
                break

    # Draw selection rectangle if in selection mode
    if mode == Mode.SELECT and selection_anchor_x is not None and cursor_x is not None:
        draw_selection_rect(stdscr, selection_anchor_x, selection_anchor_y,
                            cursor_x, cursor_y, screen_dims)

    # Draw error message or instructions at the bottom
    if error_message:
        draw_error_message(stdscr, screen_dims, error_message)
    else:
        draw_instructions(stdscr, screen_dims, mode)


def init_curses(stdscr):
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


def main_loop(stdscr, initial_smiles=None):
    init_curses(stdscr)

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
        state, error_message = create_molecule_from_smiles(
            initial_smiles, screen_dims)
        if state is None:
            # Failed to parse SMILES, create empty state
            state = create_empty_state(screen_dims)
    else:
        # Create empty state
        state = create_empty_state(screen_dims)

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
            redraw_screen(stdscr, state, show_smiles, screen_dims, mode,
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
            state.box, state.y_offset = recalculate_box_and_offset(
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
            symbol = enter_element(stdscr, screen_dims.max_y)
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
