"""This module enables selection of pdb structures based on their secondary structure content.
The secondary structure prediction is based on a very simple detection of ramachandran space.
I might update this at some point or remove it in favour of alternative algorithms (stride, dssp, etc.).


:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-04-11
"""
import logging
from typing import List, Tuple, Dict, Callable, Union

import matplotlib.path as mpl_path
import numpy as np
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Chain import Chain
from scipy import spatial

from atligator.pdb_util import is_amino_acid_residue, degrees

logger = logging.getLogger(__name__)

class GetRowInPolygon:
    def __init__(self):
        self.ret = []
        self.inner = []
        self.outer = []
        self.outlayers: List[int] = []
        self.current_row: List[int] = []
        self.last_wrong = False
        self.before_last_wrong = False

    def clear(self):
        self.inner.clear()
        self.outer.clear()
        self.outlayers.clear()
        self.current_row.clear()
        self.last_wrong = False

    def g_inner(self):
        return len(self.inner)

    def g_outer(self):
        return len(self.outer)

    def g_both(self):
        return self.g_inner() * 2 + self.g_outer()

    def discard(self, residue: int):
        try:
            self.inner.remove(residue)
        except ValueError:
            pass
        try:
            self.outer.remove(residue)
        except ValueError:
            pass

    def get_r_in_polys(self, inside1, inside2, nums, no_mercy=False, loose=False, helix=False, extended=()) -> List:
        t = 5
        if loose:
            t = 5  # TODO ?

        # TODO get rid of loose and helix and extended

        # TODO combine neighbouring helices

        def check_row():
            if self.g_both() > t:
                if self.last_wrong:
                    self.current_row.pop()
                    self.outlayers.pop()
                    if self.before_last_wrong:
                        self.current_row.pop()
                        self.outlayers.pop()

                # double outlayer: out - in - out
                if len(self.outlayers) > 1:
                    single = set()
                    double = []
                    double_i = []
                    triple_or_higher = []
                    triple_or_higher_i = []

                    for n, r in enumerate(self.current_row):
                        if r not in self.inner + self.outer:
                            single.add(r)
                            if r - 2 in single:
                                double.append(r)
                                double_i.append(n)
                                if r - 2 in double or r - 2 in triple_or_higher:
                                    triple_or_higher.append(r)
                                    triple_or_higher_i.append(n)
                                    try:
                                        double.remove(r - 2)
                                        double_i.remove(n - 2)
                                    except ValueError:
                                        pass
                                    double.remove(r)
                                    double_i.remove(n)
                    split_rows = []
                    done = 0
                    for i in triple_or_higher_i:
                        split_rows.append(self.current_row[done:i - 2])
                        done = i - 1
                    split_rows.append(self.current_row[done:])
                    second_split = []
                    for i in split_rows:
                        done = 0
                        for j in double:
                            if j in i:
                                indx = i.index(j)
                                second_split.append(i[done:indx - 2])
                                done = indx + 1
                        second_split.append(i[done:])

                    def evaluate(res_list: List[int]):
                        eval_points = 0
                        for res in res_list:
                            if res in self.inner:
                                eval_points += 2
                            elif res in self.outer:
                                eval_points += 1
                        return eval_points

                    self.current_row.clear()
                    for i in second_split:
                        if evaluate(i) > t:
                            for j in i:
                                self.current_row.append(j)

                for r in self.current_row:
                    self.ret.append(r)
            self.clear()

        def check_gap(number):
            return not (len(self.current_row) == 0 or number == self.current_row[-1] + 1)

        def move(in_inner, in_outer, number, is_last: bool, in_extended):
            if check_gap(number):
                check_row()
            if (in_inner or in_outer or in_extended) and not is_last:
                self.last_wrong = False
                self.before_last_wrong = False
                self.current_row.append(number)
                if in_inner:
                    self.inner.append(number)
                else:
                    self.outer.append(number)
            else:
                if not (no_mercy or is_last or self.last_wrong) and self.g_both() > 0:
                    self.outlayers.append(number)
                    self.last_wrong = True
                    self.current_row.append(number)
                elif helix and self.g_both() > 0 and not (is_last or self.before_last_wrong):
                    self.outlayers.append(number)
                    self.before_last_wrong = True
                    self.current_row.append(number)
                else:
                    check_row()


        if len(extended) == 0:
            extended = (False for _ in inside1)
        for n, (i, i2, nu, ex) in enumerate(zip(inside1, inside2, nums, extended)):
            is_last_residue = True if n == len(inside1) - 1 else False
            move(i, i2, nu, is_last_residue, ex)
        return self.ret


def find_sheets(residues: List[int], chain: Chain):
    def get_coor(c):
        return chain[int(c)]["CA"].get_coord()

    """
    TODO solve:
      File exiAxatligator/src/alarms/construct_selection.py", line 206, in get_coor
        return struc[0][chain][int(c)]["CA"].get_coord()
      File lib/python3.6/site-packages/Bio/PDB/Chain.py", line 94, in __getitem__
        return Entity.__getitem__(self, id)
      File lib/python3.6/site-packages/Bio/PDB/Entity.py", line 41, in __getitem__
        return self.child_dict[id]
    KeyError: (' ', 1, ' ')
    The current chain is:  B in object  1jt7B -
    """

    residues = np.asarray(residues)
    coor_list = []
    for res_id in residues:
        coor_list.append(get_coor(res_id))
    all_coor = np.asarray(coor_list)

    connected_strand_res_ids = []
    mytree = spatial.cKDTree(all_coor)
    for res_index, res_id_ in enumerate(residues):
        # get residue indices in 6 A distance of current residue
        res_in_dist = mytree.query_ball_point(all_coor[res_index], 6)
        for res_in_dist_index in res_in_dist:
            # prevents neighbours from being found
            if not (-2 <= residues[res_index] - residues[res_in_dist_index] <= 2):
                connected_strand_res_ids.append(res_id_)
                break
    return connected_strand_res_ids


def get_main_chain(pdb_path, pdb_struc, check_if_amino_acid: bool = False) -> Tuple[Chain, Callable]:
    """
    Returns the chain which corresponds the last letter of the pdb name if possible.
    Else returns the chain with the most residues and the function it used to calculate the length.
    :param pdb_path: Path of the file resulting in pdb_struc, if None binder chain is only selected by size.
    :param pdb_struc: Pdb structure to inspect for binder chain
    :param check_if_amino_acid: If every counted residue should be checked for being a valid amino acid
    """

    if check_if_amino_acid:
        def get_len(prot_chain):
            return len(list(x for x in prot_chain if is_amino_acid_residue(x)))
    else:
        get_len = len

    if pdb_path is not None:
        chain = pdb_path.split(".")[0].split("/")[-1]
        # find last character of pdb file in chain names
        if chain in (ch.get_id() for ch in pdb_struc.get_chains()):
            for model in pdb_struc:
                try:
                    return model[chain], get_len
                except KeyError:
                    pass
    # If not possible take the chain with the most residues
    return max((ch for ch in pdb_struc.get_chains()), key=get_len), get_len


def get_helix(degs: np.ndarray, res_ids: np.ndarray, fancy_ramachandran_space: bool = False):
    """
    Counts residues that seem to be part of an alpha helix based on a list of phi/psi angles.
    :param degs: An (x, 2) array with phi and psi angles (degrees)
    :param res_ids: A list of residue ids of the corresponding residues
    :param fancy_ramachandran_space: Mode of describing Ramachandran space: If True a more refined space is taken
    :return: returns number of counted residues
    """

    # define internal and external polygons that represent ramachandran space of alpha helices
    if fancy_ramachandran_space:
        alpha_path_inner = np.array([-127.5, 52.5, -127.5, 47.5, -112.5, 47.5, -112.5, 42.5, -102.5, 42.5, -102.5, 37.5,
                                     -92.5, 37.5, -92.5, 32.5, -87.5, 32.5, -87.5, 22.5, -82.5, 22.5, -82.5, 17.5,
                                     -77.5,
                                     17.5, -77.5, 12.5, -67.5, 12.5, -67.5, 7.5, -62.5, 7.5, -62.5, 2.5, -57.5, 2.5,
                                     -57.5,
                                     -7.5, -52.5, -7.5, -52.5, -12.5, -47.5, -12.5, -47.5, -22.5, -42.5, -22.5, -42.5,
                                     -32.5, -37.5, -32.5, -37.5, -62.5, -42.5, -62.5, -42.5, -67.5, -77.5, -67.5, -77.5,
                                     -62.5, -117.5, -62.5, -117.5, -57.5, -122.5, -57.5, -122.5, -47.5, -127.5, -47.5,
                                     -127.5, -37.5, -132.5, -37.5, -132.5, -17.5, -137.5, -17.5, -137.5, 2.5, -142.5,
                                     2.5,
                                     -142.5, 32.5, -137.5, 32.5, -137.5, 52.5, -127.5, 52.5]).reshape((-1, 2))

        alpha_path_outer = np.array(
            [-152.5, -97.5, -152.5, -92.5, -157.5, -92.5, -157.5, -82.5, -162.5, -82.5, -162.5, -52.5,
             -157.5, -52.5, -157.5, -37.5, -162.5, -37.5, -162.5, -7.5, -167.5, -7.5, -167.5, 32.5,
             -172.5, 32.5, -172.5, 52.5, -62.5, 52.5, -62.5, 42.5, 62.5, 27.5, -57.5, 27.5, -57.5,
             22.5, -52.5, 22.5, -52.5, 12.5, -47.5, 12.5, -47.5, 7.5, -42.5, 7.5, -42.5, 2.5, -37.5,
             2.5, -37.5, -7.5, -32.5, -7.5, -32.5, -12.5, -27.5, -12.5, -27.5, -27.5, -22.5, -27.5,
             -22.5, -47.5, -17.5, -47.5, -17.5, -67.5, -22.5, -67.5, -22.5, -77.5, -27.5, -77.5,
             -27.5, -82.5, -47.5, -82.5, -47.5, -87.5, -77.5, -87.5, -77.5, -92.5, -87.5, -92.5,
             -87.5, -97.5, -152.5, -97.5]).reshape((-1, 2))

    else:
        alpha_path_inner = np.array([-156.5, -60.4, -54.7, -60.4, -54.7, -40.2, -100.4, -40.2,
                                     -123.9, -51.0, -156.5, -51.0, -156.5, -60.4]).reshape((-1, 2))
        alpha_path_outer = np.array([-180.0, 30.0, 0.00, 30.0, 0.00, -100.0, -180.0, -100.0]).reshape((-1, 2))
    path_inner = mpl_path.Path(alpha_path_inner)
    inside_inner = path_inner.contains_points(degs)
    path_outer = mpl_path.Path(alpha_path_outer)
    inside_outer = path_outer.contains_points(degs)
    helix = GetRowInPolygon().get_r_in_polys(inside_inner, inside_outer, res_ids, helix=False, no_mercy=True)
    return len(helix)


def get_sheet(degs: np.ndarray, res_ids: np.ndarray, chain: Chain, fancy_ramachandran_space: bool = False):
    """
    Counts residues that seem to be part of an beta sheet based on a list of phi/psi angles.
    :param degs: an (x, 2) array with phi and psi angles (degrees)
    :param res_ids: A list of residue ids of the corresponding residues
    :param chain: the binder chain
    :param fancy_ramachandran_space: Mode of describing Ramachandran space: If True a more refined space is taken
    :return: returns number of counted residues
    """

    # define internal and external polygons that represent ramachandran space of beta sheets
    if fancy_ramachandran_space:
        beta_poly = np.array([-62.5, 180.0, -62.5, 172.5, -57.5, 172.5, -57.5, 167.5, -52.5, 167.5, -52.5, 157.5,
                              -47.5, 157.5, -47.5, 147.5, -42.5, 147.5, -42.5, 137.5, -37.5, 137.5, -37.5, 122.5,
                              -42.5, 122.5, -42.5, 117.5, -47.5, 117.5, -47.5, 112.5, -57.5, 112.5, -57.5, 107.5,
                              -62.5, 107.5, -62.5, 102.5, -67.5, 102.5, -67.5, 97.5, -72.5, 97.5, -72.5, 62.5,
                              -77.5, 62.5, -77.5, 52.5, -87.5, 52.5, -87.5, 47.5, -92.5, 47.5, -92.5, 52.5,
                              -97.5, 52.5, -97.5, 67.5, -102.5, 67.5, -102.5, 77.5, -107.5, 77.5, -107.5, 82.5,
                              -112.5, 82.5, -112.5, 72.5, -117.5, 72.5, -117.5, 62.5, -122.5, 62.5, -122.5, 52.5,
                              -142.5, 52.5, -142.5, 57.5, -147.5, 57.5, -147.5, 67.5, -152.5, 67.5, -152.5, 77.5,
                              -147.5, 77.5, -147.5, 87.5, -152.5, 87.5, -152.5, 97.5, -157.5, 97.5, -157.5, 112.5,
                              -162.5, 112.5, -162.5, 122.5, -167.5, 122.5, -167.5, 132.5, -172.5, 132.5,
                              -172.5, 142.5, -180.0, 142.5, -180.0, 180.0, -62.5, 180.0]).reshape((-1, 2))

        beta_ex = np.array([-172.5, 32.5, -172.5, 52.5, -177.5, 52.5, -177.5, 77.5, -180.0, 77.5, -180.0, 180.0, -42.5,
                            180.0, -42.5, 172.5, -42.5, 172.5, -37.5, 172.5, -37.5, 167.5, -32.5, 167.5, 32.5, 157.5,
                            -27.5, 157.5, -27.5, 147.5, -22.5, 147.5, -22.5, 127.5, -17.5, 127.5, -17.5, 112.5, -22.5,
                            112.5, -22.5, 107.5, -27.5, 107.5, -27.5, 102.5, -32.5, 102.5, -32.5, 97.5, -47.5, 97.5,
                            -47.5, 92.5, -52.5, 92.5, -52.5, 72.5, -57.5, 72.5, -57.5, 42.5, -62.5, 42.5, -62.5, 32.5,
                            -172.5, 32.5]).reshape((-1, 2))

        beta_sw = np.array([-177.5, -180.0, -177.5, -177.5, -172.5, -177.5, -172.5, -172.5, -167.5, -172.5, -167.5,
                            -167.5, -127.5, -167.5, -127.5, -172.5, -97.5, -172.5, -97.5, -167.5, -77.5, -167.5, -77.5,
                            -172.5, -72.5, -172.5, -72.5, -177.5, -67.5, -177.5, -67.5, -180.0, -177.5, -180.0]
                           ).reshape((-1, 2))

        beta_ne = np.array([162.5, 180.0, 162.5, 147.5, 167.5, 147.5, 167.5, 132.5, 172.5, 132.5, 172.5, 117.5,
                            177.5, 117.5, 177.5, 77.5, 180.0, 77.5, 180.0, 180.0, 162.5, 180.0]).reshape((-1, 2))

        beta_se = np.array([162.5, -180.0, 162.5, -177.5, 167.5, -177.5, 167.5, -167.5, 172.5, -167.5, 172.5, -157.5,
                            177.5, -157.5, 177.5, -147.5, 180.0, -147.5, 180.0, -180.0, 162.5, -180.0]).reshape((-1, 2))

        path4 = mpl_path.Path(beta_ne)
        insideb_ne = path4.contains_points(degs)
        path5 = mpl_path.Path(beta_se)
        insideb_se = path5.contains_points(degs)

    else:
        beta_poly = np.array([-156.5, 91.3, -70.4, 91.3, -54.7, 112.8, -54.7, 173.2, -136.9, 173.2, -136.9, 155.8,
                              -156.5, 135.6, -156.5, 91.3]).reshape((-1, 2))
        beta_ex = np.array([-180.0, 42.9, -140.8, 16.1, -86.0, 16.1, -74.3, 45.6, -74.3, 72.5, -44.3, 102.0, -44.3,
                            161.1, -46.9, 179.9, -180.0, 180.0, -180.0, 42.9]).reshape((-1, 2))
        beta_sw = np.array([-180.0, -163.8, -75.6, -163.8, -46.9, -180.0, -180.0, -180.0, -180.0, -163.8]
                           ).reshape((-1, 2))

    path = mpl_path.Path(beta_poly)
    insideb = path.contains_points(degs)
    path2 = mpl_path.Path(beta_ex)
    insideb_ex = path2.contains_points(degs)
    path3 = mpl_path.Path(beta_sw)
    insideb_sw = path3.contains_points(degs)

    insideb2 = []
    if fancy_ramachandran_space:
        for i, w in enumerate(insideb_ex):
            if w:
                insideb2.append(w)
            elif insideb_sw[i]:
                insideb2.append(insideb_sw[i])
            elif insideb_se[i]:
                print("SE is TRUE!")
                insideb2.append(insideb_se[i])
            elif insideb_ne[i]:
                print("NE is TRUE!")
                insideb2.append(insideb_ne[i])
            else:
                insideb2.append(False)
    else:
        for i, w in enumerate(insideb_ex):
            if w:
                insideb2.append(w)
            else:
                insideb2.append(insideb_sw[i])

    sth_changes = True
    iteration = 0
    loose_strands = []

    while sth_changes:
        start = loose_strands.copy()
        # get residue ids of potential beta strands
        loose_strands = GetRowInPolygon().get_r_in_polys(insideb, insideb2, res_ids, loose=True)
        if len(loose_strands) == 0:
            break
        # find connected strands
        loose_sheet = find_sheets(loose_strands, chain=chain)
        if loose_strands == start:
            sth_changes = False
        else:
            for n, number in enumerate(res_ids):
                if number not in loose_sheet:
                    insideb[n] = False
                    insideb2[n] = False
        iteration += 1
        if iteration > 10:
            print("Too many iterations...")
            import sys
            sys.exit()

    sheet = loose_strands
    return len(sheet)


def read_phi_psi(pdb_entity, check_if_amino_acid: bool = False):
    """
    returns the phi/psi angles of the backbone of a certain pdb_entity.
    :param pdb_entity: Bio.PDB chain/structure/model object to read residue features (number, dihedral angles) from
    :param check_if_amino_acid: verbose testing for canonical amino acids (slower) before counting
    :return: returns two array: an (x, 2) array with phi and psi angles (degrees) and one containing the
    corresponding residue ids
    """
    degs = []
    res_ids = []
    polypeptides = PPBuilder().build_peptides(pdb_entity)
    for poly in polypeptides:
        phipsi_list: List = poly.get_phi_psi_list()
        for res, phipsi_res_rad in zip(poly, phipsi_list):
            if not check_if_amino_acid or is_amino_acid_residue(res):
                res_ids.append(res.get_id()[1])
                phipsi_res_deg: List[float] = []
                for phi_or_psi in phipsi_res_rad:
                    if phi_or_psi is not None:
                        phipsi_res_deg.append(degrees(phi_or_psi))
                    else:
                        res_ids.pop()
                        break
                if len(phipsi_res_deg) == 2:
                    degs.extend(phipsi_res_deg)
    degs = np.asarray(degs).reshape((-1, 2))
    res_ids = np.asarray(res_ids)
    return degs, res_ids


def get_alpha_beta_total(pdb, quiet: bool = True, check_if_amino_acid: bool = False,
                         take_last_character_as_main_chain: bool = True):
    if not quiet:
        print("The current pdb file is ", pdb, end=" - ")
    pdb_struc = PDBParser(QUIET=True).get_structure(pdb, pdb)
    binder_chain, get_len = get_main_chain(pdb_path=pdb if take_last_character_as_main_chain else None,
                                           pdb_struc=pdb_struc,
                                           check_if_amino_acid=check_if_amino_acid)
    phi_psi, residue_ids = read_phi_psi(binder_chain, check_if_amino_acid=check_if_amino_acid)

    if not len(phi_psi):
        return 0, 0, 0

    alpha = get_helix(degs=phi_psi, res_ids=residue_ids)
    beta = get_sheet(degs=phi_psi, res_ids=residue_ids, chain=binder_chain)
    total = get_len(binder_chain)
    return alpha, beta, total


def check_pdbs_for_2s_content(pdbs: List, quiet: bool = True, check_if_amino_acid: bool = False) -> Dict[str, Dict]:
    """
    Objects in self.pdbs are checked for secondary structure content.
    Returns a dictionary with alpha/beta content for each input pdb.

    :param pdbs: List of pdb file paths
    :param quiet: If true, prints are suppressed
    :param check_if_amino_acid: verbose testing for canonical amino acids (slower) before counting
    :return: Returns a dict with alpha beta content per pdb file name
    """

    alpha_beta_content = {pdb: {} for pdb in pdbs}
    for pdb in pdbs:
        alpha, beta, total = get_alpha_beta_total(pdb, quiet=quiet, check_if_amino_acid=check_if_amino_acid)
        alpha_beta_content[pdb]['alpha'] = alpha
        alpha_beta_content[pdb]['beta'] = beta
        alpha_beta_content[pdb]['total'] = total
    return alpha_beta_content


class SelectionProcess:
    """
    Performs the ALARMS selection step for a set of (processed) pdb structures.
    """

    def __init__(self, pdbs: List[str], quiet: bool = True):
        self.pdbs = pdbs
        self.quiet = quiet

    def filter_pdbs(self, alphabeta: Union[None, float] = None, betaalpha: Union[None, float] = None,
                    alphatotal: Union[None, float] = None, betatotal: Union[None, float] = None,
                    take_last_character_as_main_chain: bool = True):
        """
        Return a list of filtered pdb paths which agree to the parameters (greater than or equals!).
        The pdb paths are given by self.pdbs (List[str])
        all parameters are ratios between first part and second part:
         (alphabeta = residues of alpha helices/ residues of beta strands)

        Only the main chain is taken into 2s consideration! The main chain is either given by the last character of the
        structure name or by the longest chain.

        :param alphabeta: desired ratio for alpha/beta residues
        :param betaalpha: desired ratio for beta/alpha residues
        :param alphatotal: desired ratio for alpha/total residues
        :param betatotal: desired ratio for beta/total residues
        :param take_last_character_as_main_chain: If False it does not try to use the last character of the input file
        as the main chain.
        :return: List of pdb files which agree to the parameters
        """

        filtered_list = []
        for pdb in self.pdbs:
            alpha, beta, total = get_alpha_beta_total(
                pdb=pdb, take_last_character_as_main_chain=take_last_character_as_main_chain)
            if total:
                if alphatotal is not None and alphatotal > alpha / total:
                    continue
                if betatotal is not None and betatotal > beta / total:
                    continue
                if alphabeta is not None:
                    if beta and alpha / beta > alphabeta:
                        continue
                    elif not beta and not alpha:
                        continue
                if betaalpha is not None:
                    if alpha and beta / alpha > betaalpha:
                        continue
                    elif not alpha and not beta:
                        continue

                filtered_list.append(pdb)
        return filtered_list

    def check_pdbs_for_2s_content(self) -> Dict[str, Dict]:
        return check_pdbs_for_2s_content(pdbs=self.pdbs, quiet=self.quiet)
