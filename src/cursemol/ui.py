"""
Module implementing the curses-specific parts of the UI.
"""

import curses
import logging

from . import canvas
from . import chem
from .state import Mode

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
