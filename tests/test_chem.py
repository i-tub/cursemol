import re

import numpy as np
import pytest
from rdkit import Chem
from rdkit import rdBase

from cursemol import chem


@pytest.fixture
def rdkit_log_to_python():
    rdBase.LogToPythonLogger()
    try:
        yield
    finally:
        rdBase.LogToCppStreams()


def test_get_box():
    mol = Chem.MolFromSmiles('CCC')
    conf = Chem.Conformer(3)
    conf.SetPositions(np.array([[1., 2.], [3., 4.], [5., 6.]]))
    mol.AddConformer(conf)
    box = chem.get_box(mol)
    assert (box == np.array([[1., 2., 0.], [5., 6., 0.]])).all()


def test_CaptureRDKitLog(rdkit_log_to_python):
    with chem.CaptureRDKitLog() as log:
        mol = Chem.MolFromSmiles('x')
        assert mol is None
    msg = log.getMessage()
    assert 'SMILES Parse Error' in msg
    assert not re.search(r'\[..:..:..\]', msg)


def test_reverse_bond():
    mol = chem.get_mol('C[C@H](O)CC')

    # One bond should be dash/wedge, but which/what are up to RDKit.
    chiral_bonds = [
        b for b in mol.GetBonds()
        if b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
    ]
    assert len(chiral_bonds) == 1

    bond = chiral_bonds[0]
    bond_dir = bond.GetBondDir()
    begin_atom_idx = bond.GetBeginAtomIdx()
    end_atom_idx = bond.GetEndAtomIdx()

    chem.reverse_bond(bond)

    # Direction shouldn't change, but atoms should be swapped.
    assert bond.GetBondDir() == bond_dir
    assert bond.GetBeginAtomIdx() == end_atom_idx
    assert bond.GetEndAtomIdx() == begin_atom_idx


@pytest.mark.parametrize('bond_order,expected_return,expected_smiles', [
    (1, False, 'CC'),
    (2, True, 'C=C'),
    (3, True, 'C#C'),
])
def test_modify_bond_order(bond_order, expected_return, expected_smiles):
    mol = chem.get_mol('CC')
    assert chem.modify_bond(mol, 0, 1, bond_order) == expected_return
    assert Chem.MolToSmiles(mol) == expected_smiles


def test_modify_bond_dir():
    mol = chem.get_mol('CC')

    # Wedges in CC are  meaningless in terms of stereo, but the skecher allows
    # it because wedges are not interpreted until final conversion to SMILES.
    assert chem.modify_bond(mol,
                            0,
                            1,
                            bond_order=1,
                            bond_dir=Chem.BondDir.BEGINWEDGE)
    assert mol.GetBondWithIdx(0).GetBondDir() == Chem.BondDir.BEGINWEDGE
    assert mol.GetBondWithIdx(0).GetBeginAtomIdx() == 0
    assert mol.GetBondWithIdx(0).GetEndAtomIdx() == 1

    # Calling with the same args should swap the begin/end end atoms
    assert chem.modify_bond(mol,
                            0,
                            1,
                            bond_order=1,
                            bond_dir=Chem.BondDir.BEGINWEDGE)
    assert mol.GetBondWithIdx(0).GetBondDir() == Chem.BondDir.BEGINWEDGE
    assert mol.GetBondWithIdx(0).GetBeginAtomIdx() == 1
    assert mol.GetBondWithIdx(0).GetEndAtomIdx() == 0

    # Calling with a different direction should only change the direction
    assert chem.modify_bond(mol,
                            0,
                            1,
                            bond_order=1,
                            bond_dir=Chem.BondDir.BEGINDASH)
    assert mol.GetBondWithIdx(0).GetBondDir() == Chem.BondDir.BEGINDASH
    assert mol.GetBondWithIdx(0).GetBeginAtomIdx() == 1
    assert mol.GetBondWithIdx(0).GetEndAtomIdx() == 0

    # Calling with no direction should set direction to NONE
    assert chem.modify_bond(mol, 0, 1, bond_order=1)
    assert mol.GetBondWithIdx(0).GetBondDir() == Chem.BondDir.NONE
    assert mol.GetBondWithIdx(0).GetBeginAtomIdx() == 1
    assert mol.GetBondWithIdx(0).GetEndAtomIdx() == 0
