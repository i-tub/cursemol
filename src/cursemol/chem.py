"""
Functions related to the underlying chemistry model (RDKit).
"""

from __future__ import annotations

import io
import re
import logging

from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import AllChem
import numpy as np

rdkit_logger = logging.getLogger('rdkit')


class CaptureRDKitLog:
    """
    Context manager to capture RDKit log messages.
    """

    def __enter__(self) -> 'CaptureRDKitLog':
        self._stream = io.StringIO()
        self._old_stream = rdkit_logger.handlers[0].setStream(self._stream)
        return self

    def __exit__(self, *a) -> None:
        rdkit_logger.handlers[0].setStream(self._old_stream)

    def getMessage(self) -> str:
        """Return log messages after stripping them of timestamps"""
        return re.sub(r'\[..:..:..] ', '', self._stream.getvalue())


def get_box(mol: Chem.Mol) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the bounding box for the molecule.
    """
    conf = mol.GetConformer(0)
    xyz = conf.GetPositions()
    return (xyz.min(axis=0), xyz.max(axis=0))


def reverse_bond(bond: Chem.Bond) -> None:
    """
    Reverse bond by deleting and re-adding with swapped atoms.
    """
    mol = bond.GetOwningMol()
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()
    bond_dir = bond.GetBondDir()
    bond_type = bond.GetBondType()
    mol.RemoveBond(a1, a2)
    bond_idx = mol.AddBond(a2, a1, bond_type) - 1
    bond = mol.GetBondWithIdx(bond_idx)
    bond.SetBondDir(bond_dir)


def modify_bond(mol: Chem.RWMol,
                atom1_idx: int,
                atom2_idx: int,
                bond_order: int,
                bond_dir: Chem.BondDir | None = None) -> bool:
    """
    Modify or create a bond between two atoms.
    bond_order: 0 (delete), 1 (single), 2 (double), 3 (triple)
    bond_dir: Optional bond direction (e.g., Chem.BondDir.BEGINWEDGE)

    If a bond already has the specified direction, it will be reversed
    (atoms swapped).

    Returns True if successful, False if no change was made.
    """
    bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)

    # Map bond order to BondType
    bond_type_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }

    if bond_order == 0:
        # Delete bond if it exists
        if bond is not None:
            mol.RemoveBond(atom1_idx, atom2_idx)
        else:
            # Bond doesn't exist, nothing to delete
            return False
    else:
        bond_type = bond_type_map[bond_order]

        if bond is not None:
            current_type = bond.GetBondType()
            current_dir = bond.GetBondDir()

            if (current_type == bond_type and bond_dir is None and
                    current_dir == Chem.BondDir.NONE):
                return False

            # Check if bond already has this exact type and direction
            # If so, reverse the bond (swap atoms)
            if (current_type == bond_type and bond_dir is not None and
                    current_dir == bond_dir):
                reverse_bond(bond)
            else:
                # Modify existing bond
                bond.SetBondType(bond_type)
                bond.SetBondDir(Chem.BondDir.NONE if bond_dir is
                                None else bond_dir)
        else:
            # Add new bond
            mol.AddBond(atom1_idx, atom2_idx, bond_type)
            if bond_dir is not None:
                bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                bond.SetBondDir(bond_dir)

    return True


def assign_stereo(mol: Chem.RWMol) -> None:
    """
    Assign stereochemical state to the molecule based on its geometry and bond
    directions (wedges/dashes).
    """
    # Don't set aromaticity because we want to keep the Kekule representation
    # for sketching.
    flags = (Chem.SanitizeFlags.SANITIZE_ALL &
             ~Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    Chem.SanitizeMol(mol, flags)
    Chem.DetectBondStereochemistry(mol)
    Chem.AssignChiralTypesFromBondDirs(mol)
    Chem.AssignStereochemistry(mol, force=True)


def get_smiles(mol: Chem.Mol) -> str:
    """
    Generate a SMILES deriving the stereochemical configuration from atomic
    coordinates and bond directions.
    """
    mol_for_smiles = Chem.Mol(mol)
    try:
        assign_stereo(mol_for_smiles)
        mol_for_smiles = Chem.RemoveHs(mol_for_smiles)
    except Exception:
        logging.exception("get_smiles error")
        mol_for_smiles = mol
    return Chem.MolToSmiles(mol_for_smiles)


def get_mol(smiles: str) -> Chem.RWMol | None:
    """
    Create a Mol from a SMILES, in the state expected by the sketcher. This
    includes kekulization, 2D coordinates, and wedge bonds.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        Chem.Kekulize(mol, True)
        mol = Chem.RWMol(mol)
        AllChem.Compute2DCoords(mol)
        Chem.WedgeMolBonds(mol, mol.GetConformer())
    return mol


def compute_coords_with_fixed_atoms(mol: Chem.RWMol,
                                    num_fixed_atoms: int) -> None:
    """
    Compute 2D coordinates for a molecule, keeping existing atoms fixed.

    This is useful when adding new atoms to a molecule - the original atoms
    maintain their positions while new atoms are positioned around them.

    Args:
        mol: RDKit molecule
        num_fixed_atoms: Number of atoms at the beginning to keep fixed
    """
    # Create coordinate map to keep original atoms fixed
    coord_map = {}
    conf = mol.GetConformer()
    for i in range(num_fixed_atoms):
        pos = conf.GetAtomPosition(i)
        coord_map[i] = Geometry.Point2D(pos.x, pos.y)

    # Compute 2D coordinates for new atoms only
    AllChem.Compute2DCoords(mol, coordMap=coord_map)
