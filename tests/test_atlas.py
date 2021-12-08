import pickle
import unittest
from pathlib import Path

import numpy as np
from Bio.PDB import Vector, Residue, Atom, PDBParser
from numpy.testing import assert_allclose

from atligator.atlas import AtlasDatapoint, generate_atlas

dp_l_restype = "GLU"
dp_b_restype = "THR"
dp_l_origin = '1o6o/0/A/GLU2'
dp_b_origin = '1o6o/0/A/THR5'
dp_l_atoms = {'N': Vector(-0.54, -0.70, -1.20), 'CA': Vector(0.00, 0.00, 0.00), 'C': Vector(-0.49, 1.43, 0.00),
              'O': Vector(-0.93, 1.95, -1.03), 'CB': Vector(1.53, 0.00, 0.00), 'CG': Vector(2.15, -1.37, -0.04),
              'CD': Vector(1.78, -2.13, -1.29), 'OE1': Vector(2.08, -1.63, -2.39), 'OE2': Vector(1.17, -3.21, -1.16)}

dp_b_atoms = {'N': Vector(1.96, 4.13, -1.97), 'CA': Vector(1.62, 4.21, -3.38), 'C': Vector(0.58, 5.28, -3.66),
              'O': Vector(0.79, 6.14, -4.52), 'CB': Vector(1.09, 2.85, -3.86), 'OG1': Vector(2.06, 1.83, -3.58),
              'CG2': Vector(0.82, 2.87, -5.36)}

dp_l_atoms_g = {'N': Vector(-0.54, -0.70, -1.20), 'CA': Vector(0.00, 0.00, 0.00), 'C': Vector(-0.49, 1.43, 0.00),
                'O': Vector(-0.93, 1.95, -1.03), 'CB': Vector(1.53, 0.00, 0.00)}

dp_b_atoms_g = {'N': Vector(1.96, 4.13, -1.97), 'CA': Vector(1.62, 4.21, -3.38), 'C': Vector(0.58, 5.28, -3.66),
                'O': Vector(0.79, 6.14, -4.52), 'CB': Vector(1.09, 2.85, -3.86)}


class TestAtlasDatapoint(unittest.TestCase):

    def setUp(self) -> None:
        self.dp = AtlasDatapoint(
            ligand_restype=dp_l_restype,
            binder_restype=dp_b_restype,
            ligand_origin=dp_l_origin,
            binder_origin=dp_b_origin,
            ligand_atoms=dp_l_atoms,
            binder_atoms=dp_b_atoms
        )

    def test_atoms(self):
        self.assertEqual(self.dp.binder_calpha(), dp_b_atoms['CA'])
        self.assertEqual(self.dp.binder_cbeta(), dp_b_atoms['CB'])
        self.assertEqual(self.dp.binder_co(), dp_b_atoms['C'])

    def test_orientations(self):
        cb_ca = dp_b_atoms['CB'] - dp_b_atoms['CA']
        co_ca = dp_b_atoms['C'] - dp_b_atoms['CA']
        assert_allclose(self.dp.binder_orientation().get_array(), cb_ca.get_array())
        assert_allclose(self.dp.binder_secondary_orientation().get_array(), co_ca.get_array())

    def test_to_pdb(self):
        def assert_def(coords):
            self.assertIsInstance(residue, Residue.Residue)
            self.assertEqual(residue.get_id()[1], 15)
            self.assertIsInstance(residue["CA"], Atom.Atom)
            assert_allclose(residue["CA"].get_coord().get_array(), coords)

        residue = self.dp.binder_to_pdb_residue(15)
        assert_def(dp_b_atoms["CA"].get_array())
        residue = self.dp.ligand_to_pdb_residue(15)
        assert_def(dp_l_atoms["CA"].get_array())

    def test_to_pdb_gly(self):
        def assert_gly():
            self.assertIsInstance(residue, Residue.Residue)
            self.assertEqual(residue.get_id()[1], 25)
            with self.assertRaises(KeyError):
                return residue["CB"]

        self.dp_g = AtlasDatapoint(
            ligand_restype="GLY",
            binder_restype="GLY",
            ligand_origin='1o6o/0/A/GLY2',
            binder_origin='1o6o/0/A/THR5',
            ligand_atoms=dp_l_atoms_g,
            binder_atoms=dp_b_atoms_g
        )
        residue = self.dp_g.binder_to_pdb_residue(25)
        assert_gly()
        residue = self.dp_g.ligand_to_pdb_residue(25)
        assert_gly()


class TestAtlas(unittest.TestCase):

    def setUp(self) -> None:
        self.filename = "5AEI_pair.pdb"
        self.filename_big = "5AEI_short.pdb"
        self.filepath = Path(__file__).parent.joinpath(self.filename).__str__()
        self.filepath_big = Path(__file__).parent.joinpath(self.filename_big).__str__()
        self.parameters = {"ir_default": 4.0, "ir_hbond": 6.0, "ir_aromatic": 6.0, "ir_ionic": 8.0,
                           "minlen": 1, "maxlen": 20, "include_hydrogens": False,
                           "include_alternative_models": False, "alternative_lig_aa": None,
                           "restrict_to_alternative": True, "allow_self_interactions": False, "skip_bb_atoms": False}
        self.atlas = generate_atlas([self.filepath], **self.parameters, n_workers=1, observer=None)
        big_atlas_filename = Path(__file__).parent.joinpath("big.atlas")
        if False:#big_atlas_filename.exists():
            with open(big_atlas_filename, "rb") as rf:
                self.big_atlas = pickle.load(rf)
        else:
            self.big_atlas = generate_atlas([self.filepath_big], **self.parameters, n_workers=1, observer=None)
            with open(big_atlas_filename, "wb") as wf:
                pickle.dump(self.big_atlas, wf)

    def test_generation(self):
        # Assert all parameters are safed into Atlas features
        self.assertDictEqual(self.atlas.features, self.parameters)

    def test_find_datapoints(self):
        def check_ligand_chirality(dp: AtlasDatapoint):
            res: Residue = dp.ligand_to_pdb_residue(1)
            # TODO remove print(res["CB"].coord.get_array(), np.full([3], 0.0))
            assert_allclose(res["CA"].coord.get_array(), np.full([3], 0.0))
            assert_allclose(res["CB"].coord.get_array(), np.array([1.5, 0.0, 0.0]), rtol=0.5, atol=1e-10)
            assert_allclose(res["C"].coord.get_array(), np.array([-0.5, 1.5, 0.0]), rtol=0.5, atol=1e-10)
            assert_allclose(res["N"].coord.get_array(), np.array([-0.5, -0.7, -1.0]), rtol=0.5)

        # Assert the number of datapoints is as expected for pair
        # TODO do we want two datapoints from a pair of residues or one? This happens if all residues are found
        #  as ligand residues (length of polypeptide chain). Current answer: Yes, because both can be seen as
        #  ligand and binder and do not influence each other at Atlas maps/pages or Pockets
        self.assertEqual(len(self.atlas.datapoints), 2)

        structure = PDBParser(QUIET=True).get_structure("5AEI", self.filepath)

        for s in self.atlas.datapoints:
            bori = s.binder_origin.split("/")
            lori = s.ligand_origin.split("/")
            model = int(bori[1])
            chain = bori[2]
            resid = int(bori[3][3:])
            bresname = bori[3][:3]
            lresname = lori[3][:3]

            # Assert the residue name in origin is the same as in restype
            self.assertEqual(lresname, s.ligand_restype)
            self.assertEqual(bresname, s.binder_restype)

            # Tests if binder atoms are the same as all atoms at this position
            self.assertListEqual(list(s.binder_atoms.keys()),
                                 [a.name for a in structure[model][chain][resid].get_atoms()
                                  if a.name[0] != "H"])
            del model, chain, resid, bori, lresname, bresname

            # Tests if chirality of ligand residue is still including L-amino acids only
            check_ligand_chirality(s)

    def test_group_by_ligand_origin(self):
        grouped = self.atlas.group_by_ligand_origin()
        for group in grouped.keys():
            self.assertEqual(group, grouped[group][0].ligand_origin)

    def test_filter(self):
        arg_asn = self.big_atlas.filter("ARG", "ASN")
        arg_glu = self.big_atlas.filter("ARG", "GLU")
        arg = self.big_atlas.filter("ARG")
        self.assertGreater(arg.datapoints.__len__(), arg_glu.datapoints.__len__())
        self.assertEqual(len(arg_asn.datapoints), 5)
        self.assertEqual(len(arg_glu.datapoints), 5)


if __name__ == "__main__":
    unittest.main()
