import unittest
from itertools import combinations_with_replacement, product

from Bio import PDB
from Bio.PDB.Chain import Chain
from Bio.SeqUtils import seq1

from atligator.pdb_util import canonical_amino_acids, canonical_amino_acids_short, aa_3to1_conv, ABC, \
    get_interaction_types, get_chain_by_id, find_ligands, is_canonical_amino_acid, is_amino_acid_residue, \
    has_amino_acid_atoms, get_path, get_nonbinder_residues_within_radius


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:
        # Loads reference structure and copies model to have an alternative model
        self.binder_chain = "A"
        self.ligand_chain = "D"
        self.struc = PDB.PDBParser(QUIET=True).get_structure("5AEI_short", "5AEI_short.pdb")
        new_model = self.struc[0].copy()
        new_model.id = "1"
        self.struc.add(new_model)
        self.residue = self.struc[0][self.ligand_chain][2]
        self.residue_gly = self.struc[0][self.binder_chain][137]

    def test_static(self):
        self.assertEqual(20, len(canonical_amino_acids_short))
        self.assertEqual(20, len(canonical_amino_acids))
        self.assertEqual([seq1(aa) for aa in canonical_amino_acids], [aa_3to1_conv[aa] for aa in canonical_amino_acids])
        self.assertListEqual([chr(sign) for sign in range(ord("A"), ord("Z") + 1)] +
                             [chr(sign) for sign in range(ord("a"), ord("z") + 1)] +
                             [chr(sign) for sign in range(ord("0"), ord("9") + 1)],
                             list(ABC))

    def test_get_interaction_type(self):
        ats = {"ALA": ["C", "CA", "CB", "N", "O"],
               "ARG": ["C", "CA", "CB", "CG", "CD", "CZ", "N", "NE", "NH1", "NH2", "O"],
               "ASP": ["C", "CA", "CB", "CG", "N", "O", "OD1", "OD2"],
               "ASN": ["C", "CA", "CB", "CG", "N", "ND2", "O", "OD1"],
               "CYS": ["C", "CA", "CB", "N", "O", "SG"],
               "GLU": ["C", "CA", "CB", "CG", "CD", "N", "O", "OE1", "OE2"],
               "GLN": ["C", "CA", "CB", "CG", "CD", "N", "NE2", "O", "OE1"],
               "GLY": ["C", "CA", "N", "O"],
               "HIS": ["C", "CA", "CB", "CG", "CD2", "CE1", "N", "ND1", "NE2", "O"],
               "ILE": ["C", "CA", "CB", "CG1", "CG2", "CD1", "N", "O"],
               "LEU": ["C", "CA", "CB", "CG", "CD1", "CD2", "N", "O"],
               "LYS": ["C", "CA", "CB", "CG", "CD", "CE", "N", "NZ", "O"],
               "MET": ["C", "CA", "CB", "CG", "CE", "N", "O", "SD"],
               "PHE": ["C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O"],
               "PRO": ["C", "CA", "CB", "CG", "CD", "N", "O"],
               "SER": ["C", "CA", "CB", "N", "O", "OG"],
               "THR": ["C", "CA", "CB", "CG2", "N", "O", "OG1"],
               "TRP": ["C", "CA", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "N", "NE1", "O"],
               "TYR": ["C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "OH"],
               "VAL": ["C", "CA", "CB", "CG1", "CG2", "N", "O"]
               }

        prs = {'ALA': {'ALA': {None, 'hbond'}, 'ARG': {None, 'hbond'}, 'ASP': {None, 'hbond'}, 'ASN': {None, 'hbond'},
                       'CYS': {None, 'hbond'}, 'GLU': {None, 'hbond'}, 'GLN': {None, 'hbond'}, 'GLY': {None, 'hbond'},
                       'HIS': {None, 'hbond'}, 'ILE': {None, 'hbond'}, 'LEU': {None, 'hbond'}, 'LYS': {None, 'hbond'},
                       'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'}, 'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'},
                       'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'ARG': {'ARG': {None, 'hbond'}, 'ASP': {None, 'hbond', 'ionic'}, 'ASN': {None, 'hbond'},
                       'CYS': {None, 'hbond'}, 'GLU': {None, 'hbond', 'ionic'}, 'GLN': {None, 'hbond'},
                       'GLY': {None, 'hbond'}, 'HIS': {None, 'hbond'}, 'ILE': {None, 'hbond'}, 'LEU': {None, 'hbond'},
                       'LYS': {None, 'hbond'}, 'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'}, 'PRO': {None, 'hbond'},
                       'SER': {None, 'hbond'}, 'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'},
                       'VAL': {None, 'hbond'}},
               'ASP': {'ASP': {None, 'hbond'}, 'ASN': {None, 'hbond'}, 'CYS': {None, 'hbond'}, 'GLU': {None, 'hbond'},
                       'GLN': {None, 'hbond'}, 'GLY': {None, 'hbond'}, 'HIS': {None, 'hbond', 'ionic'},
                       'ILE': {None, 'hbond'}, 'LEU': {None, 'hbond'}, 'LYS': {None, 'hbond', 'ionic'},
                       'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'}, 'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'},
                       'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'ASN': {'ASN': {None, 'hbond'}, 'CYS': {None, 'hbond'}, 'GLU': {None, 'hbond'}, 'GLN': {None, 'hbond'},
                       'GLY': {None, 'hbond'}, 'HIS': {None, 'hbond'}, 'ILE': {None, 'hbond'}, 'LEU': {None, 'hbond'},
                       'LYS': {None, 'hbond'}, 'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'}, 'PRO': {None, 'hbond'},
                       'SER': {None, 'hbond'}, 'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'},
                       'VAL': {None, 'hbond'}},
               'CYS': {'CYS': {None, 'hbond'}, 'GLU': {None, 'hbond'}, 'GLN': {None, 'hbond'}, 'GLY': {None, 'hbond'},
                       'HIS': {None, 'hbond'}, 'ILE': {None, 'hbond'}, 'LEU': {None, 'hbond'}, 'LYS': {None, 'hbond'},
                       'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'}, 'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'},
                       'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'GLU': {'GLU': {None, 'hbond'}, 'GLN': {None, 'hbond'}, 'GLY': {None, 'hbond'},
                       'HIS': {None, 'hbond', 'ionic'}, 'ILE': {None, 'hbond'}, 'LEU': {None, 'hbond'},
                       'LYS': {None, 'hbond', 'ionic'}, 'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'},
                       'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'}, 'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'},
                       'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'GLN': {'GLN': {None, 'hbond'}, 'GLY': {None, 'hbond'}, 'HIS': {None, 'hbond'}, 'ILE': {None, 'hbond'},
                       'LEU': {None, 'hbond'}, 'LYS': {None, 'hbond'}, 'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'},
                       'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'}, 'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'},
                       'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'GLY': {'GLY': {None, 'hbond'}, 'HIS': {None, 'hbond'}, 'ILE': {None, 'hbond'}, 'LEU': {None, 'hbond'},
                       'LYS': {None, 'hbond'}, 'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'}, 'PRO': {None, 'hbond'},
                       'SER': {None, 'hbond'}, 'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'},
                       'VAL': {None, 'hbond'}},
               'HIS': {'HIS': {None, 'hbond', 'aromatic'}, 'ILE': {None, 'hbond'}, 'LEU': {None, 'hbond'},
                       'LYS': {None, 'hbond'}, 'MET': {None, 'hbond'}, 'PHE': {None, 'hbond', 'aromatic'},
                       'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'}, 'THR': {None, 'hbond'},
                       'TRP': {None, 'hbond', 'aromatic'}, 'TYR': {None, 'hbond', 'aromatic'}, 'VAL': {None, 'hbond'}},
               'ILE': {'ILE': {None, 'hbond'}, 'LEU': {None, 'hbond'}, 'LYS': {None, 'hbond'}, 'MET': {None, 'hbond'},
                       'PHE': {None, 'hbond'}, 'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'}, 'THR': {None, 'hbond'},
                       'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'LEU': {'LEU': {None, 'hbond'}, 'LYS': {None, 'hbond'}, 'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'},
                       'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'}, 'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'},
                       'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'LYS': {'LYS': {None, 'hbond'}, 'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'}, 'PRO': {None, 'hbond'},
                       'SER': {None, 'hbond'}, 'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'},
                       'VAL': {None, 'hbond'}},
               'MET': {'MET': {None, 'hbond'}, 'PHE': {None, 'hbond'}, 'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'},
                       'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'PHE': {'PHE': {None, 'hbond', 'aromatic'}, 'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'},
                       'THR': {None, 'hbond'}, 'TRP': {None, 'hbond', 'aromatic'}, 'TYR': {None, 'hbond', 'aromatic'},
                       'VAL': {None, 'hbond'}},
               'PRO': {'PRO': {None, 'hbond'}, 'SER': {None, 'hbond'}, 'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'},
                       'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'SER': {'SER': {None, 'hbond'}, 'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'},
                       'VAL': {None, 'hbond'}},
               'THR': {'THR': {None, 'hbond'}, 'TRP': {None, 'hbond'}, 'TYR': {None, 'hbond'}, 'VAL': {None, 'hbond'}},
               'TRP': {'TRP': {None, 'hbond', 'aromatic'}, 'TYR': {None, 'hbond', 'aromatic'}, 'VAL': {None, 'hbond'}},
               'TYR': {'TYR': {None, 'hbond', 'aromatic'}, 'VAL': {None, 'hbond'}}, 'VAL': {'VAL': {None, 'hbond'}}}

        def get_interactions(aa1, aa2):
            """
            For a pair of amino acids return the given interaction types.
            :param aa1: First aa
            :param aa2: Second aa
            :return: Set of interaction types
            """
            feedback = set()
            for atp in product(ats[aa1], ats[aa2]):
                for result in get_interaction_types(aa1, atp[0], aa2, atp[1]):
                    feedback.add(result)
            return feedback

        # This is testing, if all expected interactions of two amino acid types are found correctly.
        for aas in combinations_with_replacement(ats.keys(), 2):
            self.assertSetEqual(get_interactions(*aas), prs[aas[0]][aas[1]])

        # TODO test for specific atom-to-atom interactions?

    def test_get_chain_by_id(self):
        # Test the number of chains
        self.assertEqual(2, len(self.struc.child_list))
        # Test if the chain is returned correctly
        self.assertEqual(get_chain_by_id(self.struc[0], self.binder_chain), self.struc[0].child_dict[self.binder_chain])
        self.assertEqual(get_chain_by_id(self.struc[0], self.binder_chain).id, self.binder_chain)
        # Test if other entities (Structure or Chain instead of Model) can not be processed
        with self.assertRaises(TypeError):
            get_chain_by_id(self.struc, self.binder_chain)
        with self.assertRaises(TypeError):
            get_chain_by_id(self.struc[0][self.binder_chain], self.binder_chain)

    def test_get_path(self):
        path_res = get_path(self.residue)
        self.assertEqual(["5AEI_short", "0", self.ligand_chain, "ARG2"], path_res.split("/"))
        path_res_gly = get_path(self.residue_gly)
        self.assertEqual(["5AEI_short", "0", self.binder_chain, "GLY137"], path_res_gly.split("/"))

    def test_get_matching_sequences(self):
        # TODO implement? low priority since use is low.
        pass

    def test_is_canonical_amino_acid(self):
        self.assertTrue(is_canonical_amino_acid("THR"))
        self.assertFalse(is_canonical_amino_acid("THE"))

    def test_has_amino_acid_atoms(self):
        self.assertTrue(has_amino_acid_atoms(self.residue))
        self.assertTrue(has_amino_acid_atoms(self.residue, glycine=True))
        self.assertFalse(has_amino_acid_atoms(self.residue_gly))
        self.assertTrue(has_amino_acid_atoms(self.residue_gly, glycine=True))

    def test_is_amino_acid_residue(self):
        residue = self.residue.copy()
        residue_gly = self.residue_gly.copy()
        self.assertTrue(is_amino_acid_residue(residue))
        self.assertTrue(is_amino_acid_residue(residue_gly))
        residue.resname = "THREST"
        self.assertFalse(is_amino_acid_residue(residue))
        self.assertTrue(is_amino_acid_residue(residue, extended="THREST"))
        residue.resname = "Gly"
        self.assertFalse(is_amino_acid_residue(residue))
        self.assertTrue(is_amino_acid_residue(residue, extended="Gly"))

    def test_find_ligands(self):
        for min_len, max_len in [(1, 10), (10, 1000), (1, 1000)]:
            chain_lens = [len(ch.child_list) for ch in self.struc[0].child_list if
                          min_len <= len(ch.child_list) <= max_len]
            ligs = list(find_ligands(self.struc[0], min_len, max_len))
            self.assertIsInstance(ligs[0], Chain)
            self.assertEqual(len(ligs), len(chain_lens))
        # TODO does not yet test lengths of peptides with alternative AAs.

    def test_get_nonbinder_residues_within_radius(self):
        # Testing edge cases for 5AEI_short interaction.
        # first setup parameter is segmentation (only nonbinder residues in a line of at least 3 res are returned)
        # the second one is if hydrogen atoms are also included (interactions are found in short distances)
        # The third parameter is the maximum radius around the binder chain to look for residues
        setups = {(0, 1, 1.8): 0,  # dict values are the expected number of residues to find
                  (0, 1, 1.9): 3,
                  (0, 1, 2.0): 4,
                  (0, 1, 2.3): 5,
                  (1, 1, 1.9): 0,
                  (1, 1, 2.0): 3,
                  (1, 1, 2.3): 5,
                  (0, 0, 2.7): 0,
                  (0, 0, 2.8): 3,
                  (0, 0, 2.9): 5,
                  (0, 0, 5.0): 5,
                  (1, 0, 2.8): 0,
                  (1, 0, 2.9): 5,
                  (1, 0, 5.0): 5,
                  }
        for setup, exp_result in setups.items():
            segmentation, incl_hydrogens, radius = setup
            actual_result = get_nonbinder_residues_within_radius(self.struc[0],
                                                                 binder_chain=self.struc[0][self.binder_chain],
                                                                 ir_max=radius,
                                                                 ligand_min_len=3,
                                                                 segmentation=bool(segmentation),
                                                                 include_hydrogens=bool(incl_hydrogens))
            self.assertEqual(exp_result, len(list(actual_result)))


if __name__ == '__main__':
    unittest.main()
