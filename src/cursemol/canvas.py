"""
Functions that "draw" molecules to an abstract canvas, that is, a 2D array of
characters. Also, helper functions for converting between real atom coordinates
and canvas coordinates.

This module does not use curses.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Iterator

import numpy as np
from rdkit import Chem

from . import config
from .state import State, ScreenDimensions


@dataclass
class Coords:
    """2D screen coordinates (integer cursor/terminal position)."""
    x: int
    y: int


def screen_to_mol_coords(
    cursor: Coords,
    box: tuple[np.ndarray, np.ndarray],
    scale: tuple[float, float],
    screen_dims: ScreenDimensions,
) -> tuple[float, float]:
    """Convert cursor/terminal coordinates to molecule coordinates."""
    # Reverse the coordinate transformation from screen_coords_for_atom
    mol_x = (cursor.x - config.PADDING) / scale[0] + box[0][0]
    mol_y = screen_y_to_mol_y(cursor.y, box[0][1], scale[1], screen_dims.rows)

    return mol_x, mol_y


def screen_y_to_mol_y(screen_y: int, box_min_y: float, scale_y: float,
                      rows: int) -> float:
    """
    Convert screen Y coordinate to molecular Y coordinate.

    Reverses the Y-axis flip: terminal Y points down but molecular Y points up.
    """
    y_from_bottom = rows - 1 - screen_y
    return (y_from_bottom - config.PADDING) / scale_y + box_min_y


def mol_y_to_screen_y(mol_y: float, box_min_y: float, scale_y: float,
                      rows: int) -> int:
    """
    Convert molecular Y coordinate to screen Y coordinate.

    The Y-axis is flipped: terminal Y-axis points down (0 at top),
    but molecular coordinates point up (higher Y = higher position).
    Therefore higher molecular Y -> lower screen Y (closer to top).
    """
    y_from_bottom = config.PADDING + int((mol_y - box_min_y) * scale_y)
    return rows - 1 - y_from_bottom


def shift_view(state: State, dx: float, dy: float) -> None:
    """
    Shift the view (pan the molecule) by the given (dx, dy), in screen units.
    """
    (xmin, ymin, zmin), (xmax, ymax, zmax) = state.box

    # Calculate step size based on scale
    dx = dx / state.scale[0]
    dy = dy / state.scale[1]

    state.box = (
        np.array([xmin + dx, ymin + dy, zmin]),
        np.array([xmax + dx, ymax + dy, zmax]),
    )


def zoom_view(state: State, screen_dims: ScreenDimensions,
              zoom_factor: float) -> None:
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
    xscale = max(config.MIN_SCALE, min(config.MAX_SCALE, xscale))
    yscale = xscale * config.ASPECT_RATIO
    state.scale = (xscale, yscale)

    # Calculate new box dimensions to show at new scale
    screen_width = screen_dims.max_x - 2 * config.PADDING
    screen_height = screen_dims.rows - 2 * config.PADDING

    # Molecule coordinate range that fits on screen at new scale
    mol_width = screen_width / xscale
    mol_height = screen_height / yscale

    # New box centered on the same point
    state.box = (
        np.array([center_x - mol_width / 2, center_y - mol_height / 2, 0.0]),
        np.array([center_x + mol_width / 2, center_y + mol_height / 2, 0.0]),
    )


def screen_coords_for_atom(atom: Chem.Atom, state: State, conf: Chem.Conformer,
                           rows: int) -> tuple[int, int]:
    """Calculate screen coordinates for an atom."""
    pos = conf.GetAtomPosition(atom.GetIdx())
    x = config.PADDING + int((pos.x - state.box[0][0]) * state.scale[0])
    y = mol_y_to_screen_y(pos.y, state.box[0][1], state.scale[1], rows)
    return x, y


def normalize_rect(corner1: Coords,
                   corner2: Coords) -> tuple[int, int, int, int]:
    """Normalize rectangle coordinates to (min_x, min_y, max_x, max_y)."""
    return (min(corner1.x, corner2.x), min(corner1.y, corner2.y),
            max(corner1.x, corner2.x), max(corner1.y, corner2.y))


def iter_atom_screen_positions(
        state: State,
        screen_dims: ScreenDimensions) -> Iterator[tuple[Chem.Atom, int, int]]:
    """Yield (atom, screen_x, screen_y) for each atom in molecule."""
    if state.mol.GetNumAtoms() == 0:
        return

    conf = state.mol.GetConformer()
    for atom in state.mol.GetAtoms():
        x, y = screen_coords_for_atom(atom, state, conf, screen_dims.rows)
        yield atom, x, y


def find_nearest_atom(
        state: State,
        cursor: Coords,
        screen_dims: ScreenDimensions,
        exclude_atom_idx: int | None = None) -> tuple[int, Coords] | None:
    """
    Find the atom nearest to the cursor position.
    Returns (atom_index, coords) or None if no atoms.
    """
    min_dist_sq = float('inf')
    nearest_atom = None
    nearest_pos = None

    for atom, screen_x, screen_y in iter_atom_screen_positions(
            state, screen_dims):
        atom_idx = atom.GetIdx()

        # Skip excluded atom
        if exclude_atom_idx is not None and atom_idx == exclude_atom_idx:
            continue

        # Calculate squared distance (avoid sqrt for performance)
        dist_sq = (screen_x - cursor.x)**2 + (screen_y - cursor.y)**2

        if dist_sq < min_dist_sq:
            min_dist_sq = dist_sq
            nearest_atom = atom_idx
            nearest_pos = Coords(screen_x, screen_y)

    if nearest_atom is not None and nearest_pos is not None:
        return (nearest_atom, nearest_pos)
    return None


def find_atom_at_cursor(state: State,
                        cursor: Coords,
                        screen_dims: ScreenDimensions,
                        tolerance: int = 1) -> int | None:
    """
    Find an atom at or near the cursor position (within tolerance cells).
    Returns atom index or None if no atom found.
    """
    for atom, screen_x, screen_y in iter_atom_screen_positions(
            state, screen_dims):
        if (abs(screen_x - cursor.x) <= tolerance and
                abs(screen_y - cursor.y) <= tolerance):
            return atom.GetIdx()

    return None


def find_bond_atoms(state: State, cursor: Coords,
                    screen_dims: ScreenDimensions) -> tuple[int, int] | None:
    """
    Find the two atoms closest to cursor position such that
    the cursor is roughly between them (angle >= 160 degrees).
    Uses screen coordinates to avoid floating point precision issues.
    Returns (atom1_idx, atom2_idx) or None if no valid pair found.
    """
    if state.mol.GetNumAtoms() < 2:
        return None

    # Calculate screen positions and distances for all atoms
    distances = []
    for atom, screen_x, screen_y in iter_atom_screen_positions(
            state, screen_dims):
        dx = screen_x - cursor.x
        dy = screen_y - cursor.y
        dist = math.sqrt(dx * dx + dy * dy)
        distances.append((dist, atom.GetIdx(), screen_x, screen_y))

    # Sort by distance
    distances.sort()

    # Find all valid pairs and choose the one with smallest combined distance
    best_pair = None
    best_distance_sum = float('inf')

    for i in range(len(distances)):
        dist1, idx1, x1, y1 = distances[i]

        for j in range(i + 1, len(distances)):
            dist2, idx2, x2, y2 = distances[j]

            # Calculate vectors from cursor position to each atom
            v1_x = x1 - cursor.x
            v1_y = y1 - cursor.y
            v2_x = x2 - cursor.x
            v2_y = y2 - cursor.y

            # Calculate lengths
            len1 = math.sqrt(v1_x * v1_x + v1_y * v1_y)
            len2 = math.sqrt(v2_x * v2_x + v2_y * v2_y)

            if len1 == 0 or len2 == 0:
                continue

            # Calculate angle using dot product
            dot = v1_x * v2_x + v1_y * v2_y
            cos_angle = dot / (len1 * len2)
            angle_deg = math.degrees(math.acos(max(-1.0, min(1.0, cos_angle))))

            # If angle is at least 160 degrees, this is a valid pair
            if angle_deg >= 160:
                # Calculate distance between the two atoms on screen
                atom_dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

                # Score: prefer pairs where cursor is close to both atoms
                # AND the atoms themselves are close to each other
                score = dist1 + dist2 + atom_dist

                if score < best_distance_sum:
                    best_distance_sum = score
                    best_pair = (idx1, idx2)

    return best_pair


def draw_line(screen: list[list[str]], char: str, char2: str, x1: int, y1: int,
              x2: int, y2: int) -> None:
    """
    Draw a line from (x1, y1) to (x2, y2) using `char` for the first half of the
    line and `char2` for the second half.
    """
    vertical = False
    if abs(x2 - x1) < abs(y2 - y1):
        x1, y1 = y1, x1
        x2, y2 = y2, x2
        vertical = True
    try:
        slope = 1.0 * (y2 - y1) / (x2 - x1)
    except ZeroDivisionError:
        return

    rev = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        rev = True
    mid = round((x1 + 1 + x2) / 2)
    for x in range(max(x1 + 1, 0), x2):
        y = int(round(y1 + slope * (x - x1)))
        if y < 0:
            continue
        if rev:
            c = char if x >= mid else char2
        else:
            c = char if x < mid else char2
        try:
            if vertical:
                screen[x][y] = c
            else:
                screen[y][x] = c
        except IndexError:
            pass  # Ignore out-of-bounds.


def draw_string(screen: list[list[str]], screen_colors: list[list[int]], s: str,
                x: int, y: int, rows: int, cols: int, color: int,
                is_bold: bool) -> None:
    """
    Draw a string into the screen buffer with the given color and bold flag.

    Args:
        screen: 2D character array
        screen_colors: 2D color/attribute array
        s: String to draw
        x, y: screen coordinates for the start of the string
        rows, cols: screen buffer dimensions
        color: color code for the string
        is_bold: whether to draw the string in bold
    """
    for i, c in enumerate(s):
        if 0 <= y < rows and 0 <= x + i < cols:
            screen[y][x + i] = c
            screen_colors[y][x + i] = color | (0x100 if is_bold else 0)


def draw_atom(screen: list[list[str]], screen_colors: list[list[int]],
              atom: Chem.Atom, x: int, y: int, rows: int, cols: int,
              state: State, conf: Chem.Conformer) -> None:
    """
    Draw a single atom with its symbol and charge into the screen buffer.

    Args:
        screen: 2D character array
        screen_colors: 2D color/attribute array
        atom: RDKit atom object
        x, y: screen coordinates for the atom
        rows, cols: screen buffer dimensions
        state: State object with scale and box
        conf: RDKit conformer object
    """
    sym = atom.GetSymbol()
    color = config.ELEMENT_COLORS.get(sym, 0)  # Get color for this element
    # Make all symbols bold except C
    is_bold = sym != 'C'

    # Draw atom symbol
    draw_string(screen, screen_colors, sym, x, y, rows, cols, color, is_bold)

    # Draw hydrogens if heteroatom
    h_str = ''
    h_on_left = False
    if sym in 'NOPS':
        try:
            atom.UpdatePropertyCache()
        except Exception:
            # Ignore exception; H count may be wrong but it's up to the user to
            # fix the charge if necessary.
            pass
        h = atom.GetTotalNumHs()
        if h > 0:
            h_str = 'H'
            if h > 1:
                h_str += {2: '₂', 3: '₃', 4: '₄'}.get(h, str(h))

            # Check neighbor positions to determine hydrogen placement
            has_neighbor_on_left = False
            has_neighbor_on_right = False
            neighbors = atom.GetNeighbors()

            for neighbor in neighbors:
                neighbor_pos = conf.GetAtomPosition(neighbor.GetIdx())
                neighbor_x = config.PADDING + int(
                    (neighbor_pos.x - state.box[0][0]) * state.scale[0])
                if neighbor_x < x:
                    has_neighbor_on_left = True
                if neighbor_x > x:
                    has_neighbor_on_right = True

            # Draw hydrogens on the left if no neighbors on the left AND at
            # least one on the right
            if not has_neighbor_on_left and has_neighbor_on_right:
                h_x = x - len(h_str)
                h_on_left = True
            else:
                h_x = x + len(sym)

            draw_string(screen, screen_colors, h_str, h_x, y, rows, cols, color,
                        is_bold)

    # Draw formal charge if non-zero
    charge = atom.GetFormalCharge()
    if charge != 0:
        # Format charge string
        if charge == 1:
            charge_str = "+"
        elif charge == -1:
            charge_str = "-"
        elif charge > 0:
            charge_str = f"{charge}+"
        else:  # charge < 0
            charge_str = f"{abs(charge)}-"

        # Position: one cell above, one cell to the right of the symbol
        # (accounts for symbol length - e.g., "Cl" vs "C")
        # When H is on the left, charge goes to the right of the symbol only
        if h_on_left:
            charge_x = x + len(sym)
        else:
            charge_x = x + len(sym) + len(h_str)
        charge_y = y - 1

        draw_string(screen, screen_colors, charge_str, charge_x, charge_y, rows,
                    cols, color, is_bold)


def fill_screen_buffer(
        state: State, screen_dims: ScreenDimensions
) -> tuple[list[list[str]], list[list[int]]]:
    """
    Fill screen buffer with bonds and atoms.
    Returns (screen, screen_colors) tuple of 2D arrays.
    """
    # Calculate screen size
    rows = screen_dims.rows
    cols = screen_dims.max_x

    screen = [[' '] * cols for i in range(rows)]
    screen_colors = [[0] * cols for i in range(rows)]  # 0 = default color

    conf = state.mol.GetConformer(0)

    # Draw bonds
    try:
        for bond in state.mol.GetBonds():
            x1, y1 = screen_coords_for_atom(bond.GetBeginAtom(), state, conf,
                                            rows)
            x2, y2 = screen_coords_for_atom(bond.GetEndAtom(), state, conf,
                                            rows)
            # Only draw if bond type is in our dictionary
            if bond_char := config.BOND_CHARS.get(bond.GetBondType()):
                bond_dir_char = config.BOND_DIR_CHARS.get(
                    bond.GetBondDir(), bond_char)
                draw_line(screen, bond_char, bond_dir_char, x1, y1, x2, y2)
    except Exception:
        logging.exception("Error drawing bonds")
        # If there's any issue drawing bonds, continue to draw atoms
        pass

    # Draw atoms
    for atom, x, y in iter_atom_screen_positions(state, screen_dims):
        draw_atom(screen, screen_colors, atom, x, y, rows, cols, state, conf)

    return screen, screen_colors
