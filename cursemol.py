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

    # Instructions
    instructions = [
        "Cursemol - Place markers on screen",
        "h/j/k/l: move | x: place X | c: show coords | s: SMILES | q: quit"
    ]

    while True:
        stdscr.clear()

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

        # Show coordinates of all x'es
        elif key == ord('c'):
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

        # Enter SMILES string
        elif key == ord('s'):
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

        # Quit
        elif key == ord('q'):
            break


if __name__ == "__main__":
    curses.wrapper(main)
