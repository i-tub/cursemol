"""
Module implementing the curses-specific parts of the UI.
"""

from __future__ import annotations

import curses
import logging

from . import canvas
from . import config
from . import chem
from .canvas import Coords
from .state import Mode, State, ScreenDimensions

HELP_TEXT = """
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


def init_curses(stdscr: curses.window) -> None:
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


def render_screen_buffer(stdscr: curses.window, screen: list[list[str]],
                         screen_colors: list[list[int]]) -> None:
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


def draw_mol(stdscr: curses.window, state: State,
             screen_dims: ScreenDimensions) -> None:
    """Draw the molecule using ASCII art."""
    # Nothing to draw if molecule is empty
    if state.mol.GetNumAtoms() == 0:
        return

    # Fill screen buffer with molecular structure
    screen, screen_colors = canvas.fill_screen_buffer(state, screen_dims)

    # Render buffer to curses window
    render_screen_buffer(stdscr, screen, screen_colors)


def draw_selection_rect(stdscr: curses.window, corner1: Coords, corner2: Coords,
                        screen_dims: ScreenDimensions) -> None:
    """Draw a selection rectangle on the screen."""
    # Normalize coordinates
    min_x, min_y, max_x_rect, max_y_rect = canvas.normalize_rect(
        corner1, corner2)

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


def draw_instructions(stdscr: curses.window, screen_dims: ScreenDimensions,
                      mode: Mode) -> None:
    """Draw instructions at the bottom of the screen."""
    instructions = config.INSTRUCTIONS[mode.value]

    for i, line in enumerate(instructions):
        try:
            stdscr.addstr(screen_dims.max_y - config.INSTRUCTION_LINES + i, 0,
                          line[:screen_dims.max_x - 1])
        except curses.error:
            pass


def draw_error_message(stdscr: curses.window, screen_dims: ScreenDimensions,
                       error_message: str) -> None:
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


def redraw_screen(stdscr: curses.window,
                  state: State,
                  show_smiles: bool,
                  screen_dims: ScreenDimensions,
                  mode: Mode,
                  selection_anchor: Coords | None = None,
                  cursor: Coords | None = None,
                  error_message: str = "") -> None:
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
    if mode == Mode.SELECT and selection_anchor is not None and cursor is not None:
        draw_selection_rect(stdscr, selection_anchor, cursor, screen_dims)

    # Draw error message or instructions at the bottom
    if error_message:
        draw_error_message(stdscr, screen_dims, error_message)
    else:
        draw_instructions(stdscr, screen_dims, mode)


def prompt_user_input(stdscr: curses.window, max_y: int,
                      prompt_text: str) -> str:
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


def enter_smiles(stdscr: curses.window, max_y: int) -> str:
    """Prompt user to enter a SMILES string and return it."""
    return prompt_user_input(stdscr, max_y, "Enter SMILES: ")


def enter_element(stdscr: curses.window, max_y: int) -> str:
    """Prompt user to enter an element symbol and return it."""
    return prompt_user_input(stdscr, max_y, "Element symbol: ")


def show_help(stdscr: curses.window, screen_dims: ScreenDimensions) -> None:
    """
    Display help text and wait for user to press a key.
    """
    stdscr.clear()
    # Display the help text from the docstring
    help_lines = HELP_TEXT.strip().split('\n')

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
