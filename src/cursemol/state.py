"""
Classes that represent the current state of the sketcher, as well as the
undo/redo history.
"""

from dataclasses import dataclass
from enum import Enum

from rdkit import Chem

from . import chem
from . import config


class Mode(Enum):
    """UI mode for the main loop."""
    NORMAL = "normal"
    MOVE = "move"
    SELECT = "select"
    BOND = "bond"


@dataclass
class State:
    """Molecular drawing state: molecule and its display parameters."""
    mol: Chem.RWMol
    box: tuple  # ((min_x, min_y, min_z), (max_x, max_y, max_z))
    scale: tuple  # (xscale, yscale)

    @classmethod
    def createEmpty(cls, screen_dims):
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
        xscale = config.DEFAULT_SCALE
        yscale = xscale * config.ASPECT_RATIO
        scale = (xscale, yscale)

        return cls(mol=mol, box=box, scale=scale)

    @classmethod
    def fromSmiles(cls, smiles, screen_dims):
        """
        Create molecule from SMILES string with 2D coordinates.
        Returns (state, error_message); in case of error, state will be None.
        """
        with chem.CaptureRDKitLog() as log:
            mol = chem.get_mol(smiles)

        if mol is None:
            return None, 'Invalid SMILES\n' + log.getMessage()

        box, scale = calculate_box_and_scale(mol, screen_dims.max_x,
                                             screen_dims.max_y)
        return cls(mol=mol, box=box, scale=scale), ""

    def copy(self):
        """Create deep copy for undo/redo."""
        return State(mol=Chem.RWMol(self.mol), box=self.box, scale=self.scale)


@dataclass
class ScreenDimensions:
    """Terminal screen dimensions and derived values."""
    max_x: int
    max_y: int

    @property
    def rows(self):
        """Drawable rows (excluding instruction lines)."""
        return self.max_y - config.INSTRUCTION_LINES


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


def recalculate_box_and_offset(mol, scale, screen_dims):
    """
    Recalculate box for a molecule at a given scale.
    Centers the view on the molecule's actual bounding box.
    Returns box.
    """
    actual_box = chem.get_box(mol)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = actual_box

    # Calculate center of molecule
    center_x = (xmin + xmax) / 2
    center_y = (ymin + ymax) / 2

    # Calculate box dimensions that fill the screen at this scale
    screen_width = screen_dims.max_x - 2 * config.PADDING
    screen_height = screen_dims.rows - 2 * config.PADDING
    mol_width = screen_width / scale[0]
    mol_height = screen_height / scale[1]

    # Create box centered on molecule center
    box = ((center_x - mol_width / 2, center_y - mol_height / 2, 0.0),
           (center_x + mol_width / 2, center_y + mol_height / 2, 0.0))

    return box


def calculate_box_and_scale(mol, max_x, max_y):
    """Calculate bounding box and scale for a molecule, centered on screen."""
    actual_box = chem.get_box(mol)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = actual_box

    # Calculate scale to fit the molecule
    xscale = min((max_x - config.PADDING * 2) / (xmax - xmin),
                 config.DEFAULT_SCALE)
    yscale = xscale * config.ASPECT_RATIO
    scale = (xscale, yscale)

    # Recalculate box at this scale
    screen_dims = ScreenDimensions(max_x=max_x, max_y=max_y)
    box = recalculate_box_and_offset(mol, scale, screen_dims)

    return box, scale
