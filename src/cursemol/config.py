"""
Global configuration constants.
"""
# TODO: allow the user to override with a config file.

from rdkit import Chem

MIN_SCALE = 2.0  # columns per length unit
DEFAULT_SCALE = 8.0  # columns per length unit
MAX_SCALE = 16.0  # columns per length unit
ASPECT_RATIO = 0.4  # horizontal / vertical
PADDING = 5
ZOOM_STEP = 1.2

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

# Instructions (try to keep lines under 80 characters and more or less balanced)
INSTRUCTIONS = {
    "normal": [
        "hjkl: move | HJKL: fast | SPC: snap | m: move mol | s/S: SMILES | i/a/c/n/o: ins",
        "x/X/D: del | +/-: chg | <>: zoom | u/r: undo | ^L: clean | b/123/wd: bond | ?: help",
    ],
    "move": [
        "[Move molecule mode]",
        "hjkl/HJKL: move molecule | Esc/Enter: leave move mode | q: quit",
    ],
    "select": [
        "[Area delete mode]",
        "hjkl: move | HJKL: fast | Enter/x: delete | Esc: cancel | q: quit",
    ],
    "bond": [
        "[Add bond mode]",
        "hjkl/HJKL/SPC: move | Enter: accept | Esc: cancel",
    ]
}

# Number of instruction lines to reserve at the bottom of the screen.
INSTRUCTION_LINES = max(len(v) for v in INSTRUCTIONS.values())
