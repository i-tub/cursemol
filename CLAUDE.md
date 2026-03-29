# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

CurseMol is a terminal-based molecular sketcher built with Python curses and
RDKit. It's a single-file application (~1460 lines) that allows interactive
drawing and editing of chemical structures in the terminal.

## Running the Application

```bash
# Run directly
./cursemol.py

# Or with Python
python cursemol.py

# Load a SMILES string on startup
./cursemol.py "CCO"

# See all options
./cursemol.py --help
```

The application outputs the final SMILES string to stdout when you quit (q).

## Development Setup

```bash
# Install dependencies (uses uv)
uv sync

# Activate virtual environment
source .venv/bin/activate

# Format code with yapf (Google style as configured in .style.yapf)
yapf -i cursemol.py
```

## Architecture

### Core Data Structures

- **State**: Main application state containing the RDKit molecule, scale,
  box (bounding coordinates), y_offset, and SMILES display toggle
- **UndoHistory**: Manages undo/redo stack by saving/restoring serialized
  molecule states

### Coordinate System

The application maps between three coordinate spaces:
1. **Molecular coordinates**: Actual 3D coordinates from RDKit (Angstroms)
2. **Screen coordinates**: Terminal cursor position (columns/rows)
3. **Scaled coordinates**: Adjusted for zoom level and aspect ratio

Key constants:
- `ASPECT_RATIO = 0.4`: Compensates for terminal character width/height ratio
- `DEFAULT_SCALE = 8.0`: Columns per angstrom
- `MIN_SCALE/MAX_SCALE`: Zoom limits

Functions to convert between coordinate systems:
- `screen_to_mol_coords()`: cursor position → molecular coordinates
- `int_coords_for_atom()`: atom → screen position
- `recalculate_box_and_offset()`: Centers molecule after modifications

### RDKit Integration

The application wraps RDKit's molecular editing capabilities:
- Uses `RWMol` (read-write molecule) for editing operations
- `AllChem.Compute2DCoords()` for 2D coordinate generation
- `Chem.WedgeMolBonds()` to assign wedge/dash stereochemistry from 3D coords
- Stereochemistry detection via `DetectBondStereochemistry()`,
  `AssignChiralTypesFromBondDirs()`, and `AssignStereochemistry()`
- `Chem.MolToSmiles()` / `Chem.MolFromSmiles()` for serialization
- `compute_coords_with_fixed_atoms()`: Custom coordinate generation that
  preserves existing atom positions when appending fragments

### Bond Representation

Bonds are rendered as ASCII characters between atoms:
- Single bonds: dots (·) or middle dots (·)
- Double bonds: equals signs (=)
- Triple bonds: hash marks (#)
- Wedge bonds: solid bullets (●) pointing from narrow to wide end
- Dash bonds: hollow bullets (○) pointing from narrow to wide end

The `reverse_bond()` function handles stereochemistry direction changes.

### State Management

All operations that modify the molecule use save/restore pattern:
1. Save current state to undo history
2. Perform operation
3. Recalculate box and offset to keep molecule centered

Key editing operations:
- `insert_or_modify_atom()`: Add/change atoms at cursor
- `modify_bond()`: Add/change bonds between atoms
- `delete_at_cursor()`: Remove atoms or bonds
- `append_smiles_fragment()`: Attach SMILES fragments to atoms or bonds
- `connect_sidechain_to_bond()`: Special logic for ring formation when
  appending to bonds

### Screen Rendering

The rendering pipeline:
1. `fill_screen_buffer()`: Creates arrays of characters and colors
2. Renders atoms (with element symbols, charges)
3. Renders bonds using Bresenham line algorithm in `draw_line()`
4. `render_screen_buffer()`: Writes buffer to screen with colors
5. Adds status line and optional SMILES display

## Key Implementation Details

### Cursor Navigation
- The cursor can be on empty space, atoms, or bonds
- `find_atom_at_cursor()`: Locates atom within tolerance
- `find_bond_atoms()`: Identifies which bond the cursor is over by checking
  proximity to line segments

### Molecule Cleanup
- Ctrl-L triggers `cleanup_coordinates()` which regenerates all 2D
  coordinates using `Compute2DCoords()`
- This is useful after manual edits that distort the structure

### Area Delete (X key)
- Implements rectangle selection in screen coordinates
- `delete_atoms_in_rect()` converts screen rect to molecular coordinates
  and deletes atoms within the rectangle

### SMILES Fragment Appending
When appending to a bond (not an atom), the fragment forms a ring by
connecting its first and last atoms to the two bond atoms. This allows
easy ring fusion operations.
