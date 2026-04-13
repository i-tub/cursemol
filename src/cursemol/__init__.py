"""
CurseMol - Molecular sketcher for the terminally committed

This module has the main function and helpers for starting up the program; the
guts of the application are in the `sketcher` module.
"""

from __future__ import annotations

__version__ = "4.1.1"

import argparse
import atexit
import logging
import os
import sys

from rdkit import Chem

from . import chem
from . import sketcher

rdkit_logger = logging.getLogger('rdkit')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='CurseMol - molecular sketcher for the terminally committed'
    )
    parser.add_argument(
        'smiles',
        nargs='?',
        help='Initial SMILES string to display (use "-" to read from stdin)')
    return parser.parse_args()


def setup_tty() -> int | None:
    """
    If stdin/stdout are not TTYs (e.g., piped input/output), redirect to /dev/tty
    so curses can read keyboard input and display to the terminal.
    Returns the original stdout fd if it was redirected (for final SMILES output).
    """
    original_stdout_fd: int | None = None

    # Handle stdin
    if not sys.stdin.isatty():
        sys.stdin.close()  # Close the old stdin to avoid resource warning
        tty_fd = os.open('/dev/tty', os.O_RDONLY)
        os.dup2(tty_fd, 0)  # Replace fd 0 (stdin) with /dev/tty
        os.close(tty_fd)
        sys.stdin = os.fdopen(0, 'r')
        # Register cleanup to avoid resource warning on exit
        atexit.register(sys.stdin.close)

    # Handle stdout
    if not sys.stdout.isatty():
        original_stdout_fd = os.dup(1)  # Duplicate fd 1 before redirecting
        tty_fd = os.open('/dev/tty', os.O_WRONLY)
        os.dup2(tty_fd, 1)  # Replace fd 1 (stdout) with /dev/tty
        os.close(tty_fd)
        sys.stdout = os.fdopen(1, 'w')
        atexit.register(sys.stdout.close)

    return original_stdout_fd


def main() -> None:
    # Set up logging to file (truncate on start)
    logging.basicConfig(
        filename='cursemol.log',
        filemode='w',  # Truncate on open
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s')

    # Capture RDKit warnings
    rdkit_logger.setLevel(logging.ERROR)
    Chem.rdBase.LogToPythonLogger()

    args = parse_args()

    # Handle reading from stdin if "-" is provided
    initial_smiles = args.smiles
    if initial_smiles == "-":
        initial_smiles = sys.stdin.readline().strip()

    with chem.CaptureRDKitLog():
        original_stdout_fd = setup_tty()

    smiles = sketcher.run(initial_smiles)

    # Print to original stdout if it was redirected, otherwise to current stdout
    if original_stdout_fd is not None:
        output = os.fdopen(original_stdout_fd, 'w')
        print(smiles, file=output)
        output.close()
    else:
        print(smiles)
