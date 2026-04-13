# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

CurseMol is a terminal-based molecular sketcher built with Python curses and
RDKit. It's a Python package (~1941 lines) that allows interactive drawing
and editing of chemical structures in the terminal.

## Running the Application

```bash
# Install the package
uv sync

# Run the installed command
cursemol

# Load a SMILES string on startup
cursemol "CCO"

# See all options
cursemol --help
```

The application outputs the final SMILES string to stdout when you quit (q).

## Development Setup

```bash
# Install dependencies and package in editable mode
uv sync

# Format code with yapf (Google style as configured in .style.yapf)
yapf -i src/cursemol/*.py

# Run tests with pytest
pytest

# Run tests with coverage
pytest --cov=cursemol

# Lint code with flake8
flake8 src/cursemol
```

## Architecture

### Project Structure

The package is organized into the following modules:

**Core modules:**
- `src/cursemol/__init__.py` (95 lines): Entry point with `main()` function,
  argument parsing, TTY setup, and logging configuration
- `src/cursemol/sketcher.py` (485 lines): Main event loop and command
  handlers (load, append, help, etc.)
- `src/cursemol/state.py` (170 lines): State management classes (State,
  UndoHistory, ScreenDimensions, Mode enum)
- `src/cursemol/canvas.py` (411 lines): Screen rendering, coordinate
  conversion, zoom/pan operations
- `src/cursemol/edit.py` (272 lines): Molecule editing operations (insert,
  delete, modify atoms/bonds)
- `src/cursemol/chem.py` (179 lines): RDKit integration and chemistry
  utilities
- `src/cursemol/ui.py` (270 lines): User interface helpers (prompts, help
  screen, status line)
- `src/cursemol/config.py` (59 lines): Configuration constants (scales, bond
  characters, element colors)

**Tests:**
- `tests/test_chem.py` (153 lines): Unit tests for chemistry functions
- `tests/test_canvas.py` (67 lines): Unit tests for canvas rendering

**Build configuration:**
- `pyproject.toml`: Package configuration using hatchling build backend
- Entry point: `cursemol` command (configured in pyproject.toml)

### Module Dependencies

The modules are organized with clear separation of concerns:

```
__init__.py (entry point)
    └── sketcher.py (main loop & commands)
        ├── canvas.py (rendering & coordinates)
        │   ├── state.py (data structures)
        │   └── config.py (constants)
        ├── edit.py (molecule editing)
        │   ├── chem.py (RDKit integration)
        │   ├── state.py
        │   └── config.py
        ├── ui.py (user interface helpers)
        │   └── config.py
        └── chem.py
```

Key dependencies:
- `sketcher.py` orchestrates all other modules
- `canvas.py` and `edit.py` are independent of each other
- `state.py`, `config.py`, and `chem.py` are leaf modules with no internal
  dependencies
- All modules depend on RDKit for chemistry operations

### Core Data Structures

Defined in `state.py`:
- **State**: Main application state containing the RDKit molecule, scale,
  and box (bounding coordinates)
- **UndoHistory**: Manages undo/redo stack by saving/restoring serialized
  molecule states
- **ScreenDimensions**: Holds screen dimensions (rows, cols, max_y, max_x)
- **Mode**: Enum for UI modes (NORMAL, MOVE, SELECT, BOND)

Defined in `chem.py`:
- **CaptureRDKitLog**: Context manager for capturing RDKit logging output

### Coordinate System

The application maps between three coordinate spaces:
1. **Molecular coordinates**: Actual 3D coordinates from RDKit (Angstroms)
2. **Screen coordinates**: Terminal cursor position (columns/rows)
3. **Scaled coordinates**: Adjusted for zoom level and aspect ratio

Key constants (defined in `config.py`):
- `ASPECT_RATIO = 0.4`: Compensates for terminal character width/height ratio
- `DEFAULT_SCALE = 8.0`: Columns per angstrom
- `MIN_SCALE/MAX_SCALE`: Zoom limits
- `BOND_CHARS`: Mapping of bond types to display characters
- `BOND_DIR_CHARS`: Characters for wedge/dash stereochemistry bonds
- `ELEMENT_COLORS`: Color codes for different elements

Functions to convert between coordinate systems (in `canvas.py`):
- `screen_to_mol_coords()`: cursor position → molecular coordinates
- `screen_coords_for_atom()`: atom → screen position

Functions in `state.py`:
- `recalculate_box_and_offset()`: Centers molecule after modifications

### RDKit Integration

The `chem.py` module wraps RDKit's molecular editing capabilities:
- Uses `RWMol` (read-write molecule) for editing operations
- `rdDepictor.Compute2DCoords()` for 2D coordinate generation
- `Chem.WedgeMolBonds()` to assign wedge/dash stereochemistry from 3D coords
- Stereochemistry detection via `DetectBondStereochemistry()`,
  `AssignChiralTypesFromBondDirs()`, and `AssignStereochemistry()`
- `Chem.MolToSmiles()` / `Chem.MolFromSmiles()` for serialization
- `compute_coords_with_fixed_atoms()`: Custom coordinate generation that
  preserves existing atom positions when appending fragments
- `get_box()`: Calculate bounding box from molecule coordinates
- `CaptureRDKitLog`: Context manager for capturing RDKit warnings/errors

### Bond Representation

Bonds are rendered as ASCII characters between atoms (defined in
`config.py`, rendering in `canvas.py`):
- Single bonds: dots (·) or middle dots (·)
- Double bonds: equals signs (=)
- Triple bonds: hash marks (#)
- Wedge bonds: solid bullets (●) pointing from narrow to wide end
- Dash bonds: hollow bullets (○) pointing from narrow to wide end

The `reverse_bond()` function in `edit.py` handles stereochemistry direction
changes.

### State Management

The `state.py` module manages application state with the `UndoHistory` class.
All operations that modify the molecule use save/restore pattern:
1. Save current state to undo history
2. Perform operation
3. Recalculate box and offset to keep molecule centered

Key editing operations (in `edit.py`):
- `insert_or_modify_atom()`: Add/change atoms at cursor
- `modify_bond()`: Add/change bonds between atoms
- `delete_at_cursor()`: Remove atoms or bonds
- `connect_sidechain_to_bond()`: Special logic for ring formation when
  appending to bonds

Sketcher-level operations (in `sketcher.py`):
- `append_smiles_fragment()`: Attach SMILES fragments to atoms or bonds
- `load_smiles()`: Load molecule from SMILES string
- `cleanup_coordinates()`: Regenerate 2D coordinates

### Screen Rendering

The rendering pipeline (in `canvas.py`):
1. `fill_screen_buffer()`: Creates arrays of characters and colors
2. Renders atoms (with element symbols, charges)
3. Renders bonds using Bresenham line algorithm in `draw_line()`
4. `render_screen_buffer()`: Writes buffer to screen with colors

UI elements (in `ui.py`):
- `draw_instructions()`: Mode-specific status line
- `draw_smiles()`: Optional SMILES display
- `draw_error_message()`: Error message display

### Application Flow

Entry point (in `__init__.py`):
- `main()`: Parses arguments and initializes the application
- `setup_tty()`: Ensures proper TTY handling for curses
- `parse_args()`: Command-line argument parsing

Main sketcher loop (in `sketcher.py`):
- `run()`: Top-level function that sets up curses and runs the main loop
- `init_curses()`: Configures curses environment (colors, cursor visibility)
- `main_loop()`: Main event loop handling keyboard input and screen updates

Rendering (in `canvas.py`):
- `redraw_screen()`: Coordinates rendering of molecule, UI elements, and
  messages
- `shift_view()`: Pan the view in screen coordinates
- `zoom_view()`: Adjust zoom level

## Testing

The project uses pytest for unit testing. Tests are organized by module:

- `tests/test_chem.py`: Tests for chemistry and RDKit integration functions
  - Coordinate bounding box calculation
  - RDKit log capturing
  - SMILES parsing and validation
  - Stereochemistry handling
  - Fragment attachment logic

- `tests/test_canvas.py`: Tests for canvas rendering functions
  - Screen buffer creation
  - Atom and bond rendering
  - Coordinate transformations
  - Zoom and pan operations

Run tests with:
```bash
pytest                      # Run all tests
pytest tests/test_chem.py  # Run specific test file
pytest -v                  # Verbose output
pytest --cov=cursemol      # With coverage report
```

The tests use fixtures to create mock molecules and screen dimensions, and
include helper functions like `compare_chars()` for validating rendered
output.

## Key Implementation Details

### Cursor Navigation

Functions in `edit.py`:
- The cursor can be on empty space, atoms, or bonds
- `find_atom_at_cursor()`: Locates atom within tolerance
- `find_nearest_atom()`: Finds closest atom for snapping or bond creation
- `find_bond_atoms()`: Identifies which bond the cursor is over by checking
  proximity to line segments

### Molecule Cleanup

Functions in `sketcher.py`:
- Ctrl-L triggers `cleanup_coordinates()` which regenerates all 2D
  coordinates using `Compute2DCoords()`
- This is useful after manual edits that distort the structure

### Area Delete (X key)

Functions in `edit.py`:
- Implements rectangle selection in screen coordinates
- `delete_atoms_in_rect()` converts screen rect to molecular coordinates
  and deletes atoms within the rectangle

### SMILES Fragment Appending

Functions in `sketcher.py` and `edit.py`:
- When appending to a bond (not an atom), the fragment forms a ring by
  connecting its first and last atoms to the two bond atoms. This allows
  easy ring fusion operations.
- `append_smiles_fragment()` in `sketcher.py` handles user interaction
- `connect_sidechain_to_bond()` in `edit.py` implements ring formation logic

### UI Helpers

Functions in `ui.py`:
- `prompt_user_input()`: Generic prompt for text input
- `enter_smiles()`: Specialized prompt for SMILES strings
- `enter_element()`: Specialized prompt for element symbols
- `show_help()`: Displays help screen with keybindings
- `draw_instructions()`: Shows mode-specific status line
- `draw_error_message()`: Displays error messages to user
- `draw_smiles()`: Displays current SMILES string
