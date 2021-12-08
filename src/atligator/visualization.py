"""Utility classes, functions, and data for the visualization of atlas data using matplotlib.

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-03-05
"""
from collections import defaultdict
from itertools import combinations
from typing import Dict, List, Union, Tuple, Iterator

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from Bio.PDB.vectors import Vector
# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import axes3d, Axes3D  # do not remove this import!

from atligator.acomplex import LigandBinderComplex
from atligator.atlas import Atlas, AtlasDatapoint
from atligator.pdb_util import canonical_amino_acids
from atligator.pocket_miner import Pocket
from atligator.prediction_matrix import proximity_score, BinderPredictionMatrix
from atligator.structure import get_icoor_from_aresidue, DesignedLigandResidue

pio.renderers.default = "browser"

# coloring scheme for protein residues, see doc/color-scheme.html
color_map: Dict[str, str] = {'ALA': '#82B2B8', 'ARG': '#035A12', 'ASN': '#A24689',
                             'ASP': '#DA675E', 'CYS': '#F7E748', 'GLN': '#570340',
                             'GLU': '#750C04', 'GLY': '#A6B9C1', 'HIS': '#FFD7B2',
                             'ILE': '#043F48', 'LEU': '#1B5F68', 'LYS': '#4A8658',
                             'MET': '#DDCF3E', 'PHE': '#DA995E', 'PRO': '#74D73E',
                             'SER': '#5B5BB3', 'THR': '#343477', 'TRP': '#753A04',
                             'TYR': '#AA672A', 'VAL': '#3B7C85'}

# coloring scheme for protein residues, see doc/color-scheme.html
color_map_light: Dict[str, str] = {'ALA': '#94cad1', 'ARG': '#047317', 'ASN': '#bd52a0',
                                   'ASP': '#F27268', 'CYS': '#FFF05C', 'GLN': '#700453',
                                   'GLU': '#8F0F05', 'GLY': '#BDD2DB', 'HIS': '#FFE5CC',
                                   'ILE': '#055561', 'LEU': '#227782', 'LYS': '#59A16A',
                                   'MET': '#F7E845', 'PHE': '#DA995E', 'PRO': '#F2AA68',
                                   'SER': '#6868CC', 'THR': '#404091', 'TRP': '#8F4705',
                                   'TYR': '#C47731', 'VAL': '#46939E'}

color_map_rasmol: Dict[str, str] = {'ALA': '#C8C8C8', 'ARG': '#145AFF', 'ASN': '#00DCDC',
                                    'ASP': '#E60A0A', 'CYS': '#E6E600', 'GLN': '#00DCDC',
                                    'GLU': '#E60A0A', 'GLY': '#EBEBEB', 'HIS': '#8282D2',
                                    'ILE': '#0F820F', 'LEU': '#0F820F', 'LYS': '#145AFF',
                                    'MET': '#E6E600', 'PHE': '#3232AA', 'PRO': '#DC9682',
                                    'SER': '#FA9600', 'THR': '#FA9600', 'TRP': '#B45AB4',
                                    'TYR': '#3232AA', 'VAL': '#0F820F'}

color_map_rasmol_light: Dict[str, str] = {'ALA': '#D7D7D7', 'ARG': '#3070FF', 'ASN': '#00FAF7', 'ASP': '#F41919',
                                          'CYS': '#FFFC04', 'GLN': '#00FAF7', 'GLU': '#F41919', 'GLY': '#FAFAFA',
                                          'HIS': '#9898D9', 'ILE': '#129B12', 'LEU': '#129B12', 'LYS': '#3070FF',
                                          'MET': '#FFFC04', 'PHE': '#3838BF', 'PRO': '#E2A899', 'SER': '#FFA019',
                                          'THR': '#FFA019', 'TRP': '#BC6EBD', 'TYR': '#3838BF', 'VAL': '#129B12'}

color_map_shapely: Dict[str, str] = {'ALA': '#8CFF8C', 'ARG': '#00007C', 'ASN': '#FF7C70',
                                     'ASP': '#A00042', 'CYS': '#FFFF70', 'GLN': '#FF4C4C',
                                     'GLU': '#660000', 'GLY': '#DFDFDF', 'HIS': '#7070FF',
                                     'ILE': '#004C00', 'LEU': '#455E45', 'LYS': '#4747B8',
                                     'MET': '#B8A042', 'PHE': '#534C42', 'PRO': '#525252',
                                     'SER': '#FF7042', 'THR': '#B84C00', 'TRP': '#4F4600',
                                     'TYR': '#8C704C', 'VAL': '#FF8CFF'}

color_map_shapely_light: Dict[str, str] = {'ALA': '#A9FFA9', 'ARG': '#00009A', 'ASN': '#FE958D', 'ASP': '#BE0050',
                                           'CYS': '#FEFD8D', 'GLN': '#FF6969', 'GLU': '#840000', 'GLY': '#EDEDED',
                                           'HIS': '#8D8DFE', 'ILE': '#006A00', 'LEU': '#516E51', 'LYS': '#5B5BC0',
                                           'MET': '#C2A955', 'PHE': '#635A4E', 'PRO': '#616161', 'SER': '#FF845F',
                                           'THR': '#D65500', 'TRP': '#6C5E00', 'TYR': '#9F7D56', 'VAL': '#FDA9FF'}

connections = {'ALA': ('CA', 'CB'),
               'ARG': ('CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'CZ', 'NH2'),
               'ASN': ('CA', 'CB', 'CG', 'OD1', 'CG', 'ND2'),
               'ASP': ('CA', 'CB', 'CG', 'OD1', 'CG', 'OD2'),
               'CYS': ('CA', 'CB', 'SG'),
               'GLN': ('CA', 'CB', 'CG', 'CD', 'OE1', 'CD', 'NE2'),
               'GLU': ('CA', 'CB', 'CG', 'CD', 'OE1', 'CD', 'OE2'),
               'GLY': (),
               'HIS': ('CA', 'CB', 'CG', 'CD2', 'NE2', 'CE1', 'ND1', 'CG'),
               'ILE': ('CA', 'CB', 'CG2', 'CB', 'CG1', 'CD1'),
               'LEU': ('CA', 'CB', 'CG', 'CD1', 'CG', 'CD2'),
               'LYS': ('CA', 'CB', 'CG', 'CD', 'CE', 'NZ'),
               'MET': ('CA', 'CB', 'CG', 'SD', 'CE'),
               'PHE': ('CA', 'CB', 'CG', 'CD2', 'CE2', 'CZ', 'CE1', 'CD1', 'CG'),
               'PRO': ('N', 'CD', 'CG', 'CB', 'CA'),
               'SER': ('CA', 'CB', 'OG'),
               'THR': ('CA', 'CB', 'OG1', 'CB', 'CG2'),
               'TRP': ('CA', 'CB', 'CG', 'CD1', 'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'CE2', 'CD2', 'CG'),
               'TYR': ('CA', 'CB', 'CG', 'CD1', 'CE1', 'CZ', 'OH', 'CZ', 'CE2', 'CD2', 'CG'),
               'VAL': ('CA', 'CB', 'CG2', 'CB', 'CG1')}

cmaps = {"0": color_map,
         "1": color_map_rasmol,
         "2": color_map_shapely}

cmaps_light = {"0": color_map_light,
               "1": color_map_rasmol_light,
               "2": color_map_shapely_light}


def get_color_map(color_mode: str = "0"):
    """
    :param color_mode: If specified this changes the color mapper to a different color scheme.
    :return: The associated color as hexcode, or #808080 if no corresponding entry found in the color map
    """
    try:
        return cmaps[color_mode]
    except KeyError:
        return cmaps["0"]


def get_color(restype: str, color_mode: str = "0") -> str:
    """Using color_map, this function looks up a color for the given amino acid residue.
    :param restype: the residue type in three-letter code to look up
    :param color_mode: If specified this changes the color mapper to a different color scheme.
    :return: The associated color as hexcode, or #808080 if no corresponding entry found in the color map
    """
    _color_map = get_color_map(color_mode)
    try:
        return _color_map[restype.upper()]
    except KeyError:
        return "#808080"


def get_light_color(restype: str, color_mode: str = "0") -> str:
    """Using color_map, this function looks up a color for the given amino acid residue.
    :param restype: the residue type in three-letter code to look up
    :param color_mode: If specified this changes the color mapper to a different color scheme.
    :return: The associated color as hexcode, or #808080 if no corresponding entry found in the color map
    """
    try:
        _color_map_light = cmaps_light[color_mode]
    except KeyError:
        _color_map_light = cmaps_light["0"]
    try:
        return _color_map_light[restype.upper()]
    except KeyError:
        return "#808080"


def visualize_atlas(atlas: Atlas, method: str = "plotly", **kwargs) -> None:
    """Shows the atlas datapoints with a plotly plot (within web broswer) or matplotlib (separate window).
    For additional parameters such as 'ligand_restype' compare the subfunctions used here.

    :param atlas: the atlas to visualize
    :param method: Either "plotly" or "matplotlib". If not plotly matplotlib is used.
    """
    if method == "plotly":
        return visualize_atlas_plotly(atlas=atlas, **kwargs)
    else:
        return visualize_atlas_matplotlib(atlas=atlas, **kwargs)


def visualize_atlas_plotly(atlas: Atlas = None, ligand_restype: Union[str, None] = None,
                           binder_restype: Union[str, None] = None, **_) -> None:
    """Shows a 3d plotly plot with in the web broswer that contains both the specified residue ligand
    and the atlas datapoints associated with the corresponding residue types. All residues are shown as
    Ca Cb markers with the Cb being smaller than Ca.

    :param atlas: the atlas to visualize
    :param ligand_restype: binder residues interacting with this ligand residue types will be displayed. If unset, no
    filtering is applied.
    :param binder_restype: binder residues of this type will be displayed. If unset, no filtering is applied.
    :param _: to prevent failure if matplotlib parameters are given.
    :return: None
    """

    if atlas is None:
        import pickle
        with open("/home/josef/Documents/alarms_test/L_a.118.1_.55a_.1b.atlas", "rb") as rf:
            atlas = pickle.load(rf)

    # Filter atlas, if filtering for ligand or binder residue is desired.
    if ligand_restype is not None or binder_restype is not None:
        atlas = atlas.filter(ligand_restype=ligand_restype, binder_restype=binder_restype)
    if ligand_restype is None:
        ligand_restype = "LIG"

    # collect all coordinates, residue types and the size of the marker datapoints in lists.
    x = []
    y = []
    z = []
    restypes = []
    size = []

    # Two datapoints representing ligand Calpha (x,y,z = 0,0,0) and Cbeta (1,0,0) are placed in the center.
    x.append(0.0)
    y.append(0.0)
    z.append(0.0)
    restypes.append("ligand " + ligand_restype)
    size.append(80)
    x.append(1.0)
    y.append(0.0)
    z.append(0.0)
    restypes.append("ligand " + ligand_restype + " cb")
    size.append(25)

    # For every datapoint save coords, restype and size (Ca and Cb)
    for dp in atlas.datapoints:
        x.append(dp.binder_calpha()[0])
        y.append(dp.binder_calpha()[1])
        z.append(dp.binder_calpha()[2])
        restypes.append(dp.binder_restype)
        size.append(20)

        x.append(dp.binder_cbeta()[0])
        y.append(dp.binder_cbeta()[1])
        z.append(dp.binder_cbeta()[2])
        restypes.append(dp.binder_restype + " cb")
        size.append(7)

    df = pd.DataFrame.from_dict({
        "x": x,
        "y": y,
        "z": z,
        "restype": restypes,
        "size": size,
    })

    # concatenate color dictionaries to address Ca (full color) and Cb (light color)
    color_map_light_cb = {k + " cb": v for k, v in color_map_light.items()}
    color_map_full = {**color_map, **color_map_light_cb}
    color_map_full["ligand " + ligand_restype] = color_map_full[ligand_restype] \
        if ligand_restype != "LIG" else "#808080"
    color_map_full["ligand " + ligand_restype + " cb"] = color_map_full[ligand_restype + " cb"] \
        if ligand_restype != "LIG" else "#808080"

    fig = px.scatter_3d(df, x='x', y='y', z='z',
                        color='restype', color_discrete_map=color_map_full, size="size")
    # size also introduces opacity and borders...

    # removing border, changing opacity.
    fig.update_traces(marker=dict(opacity=1.0,
                                  line=dict(width=0)),
                      selector=dict(name="ligand " + ligand_restype + " cb"))
    fig.update_traces(marker=dict(opacity=1.0,
                                  line=dict(width=0)),
                      selector=dict(name=ligand_restype + " cb"))
    for aa in canonical_amino_acids:
        fig.update_traces(marker=dict(opacity=1.0,
                                      line=dict(width=0)),
                          selector=dict(name=aa + " cb"))
        fig.update_traces(marker=dict(opacity=1.0,
                                      line=dict(width=0)),
                          selector=dict(name=aa))

    fig.show()


def visualize_atlas_matplotlib(atlas: Atlas, ligand_restype: str or None, binder_restype: str or None,
                               title: str = None, dimensions: float = 10.0, draw_secor: bool = False) -> None:
    """Shows a 3d quiver plot that contains both the specified residue ligand (including backbone) and the atlas
    datapoints associated with the corresponding ligand residue types.
    :param atlas: the atlas to visualize
    :param ligand_restype: binder residues interacting with this ligand residue types will be displayed. If unset, no
    filtering is applied.
    :param binder_restype: binder residues of this type will be displayed. If unset, no filtering is applied.
    :param title: this will be shown as the title of the plot
    :param dimensions: max negative and positive coordinates to display, in Angstrom
    :param draw_secor: if True, secondary orientation vectors are drawn (without arrowheads)
    """
    if title is None:
        if ligand_restype is None:
            title = 'Atlas'
        else:
            title = ligand_restype
    filtered_atlas = atlas.filter(ligand_restype, binder_restype)
    # plot the data as matplotlib 3d quiver
    fig = plt.figure()
    fig.suptitle(title)
    plt.suptitle(title)
    ax: Axes3D = fig.gca(projection='3d')
    ax.set_xlim(-dimensions, dimensions)
    ax.set_ylim(-dimensions, dimensions)
    ax.set_zlim(-dimensions, dimensions)
    # draw ligand residue
    ax.quiver([-1.0], [0.0], [0.0], [3.0], [0.0], [0.0],
              color=get_color(ligand_restype), arrow_length_ratio=0.4, linewidth=8)
    # process binder residues
    for binrestype in color_map:
        # for each potential interacting binder residue, paint the quivers in individual colors.
        relevant_datapoints = [d for d in filtered_atlas.datapoints if d.binder_restype == binrestype]
        col = color_map[binrestype]
        x: List[float] = []
        y: List[float] = []
        z: List[float] = []
        u: List[float] = []
        v: List[float] = []
        w: List[float] = []
        u2: List[float] = []
        v2: List[float] = []
        w2: List[float] = []
        for d in relevant_datapoints:
            calpha = d.binder_calpha()
            x.append(calpha[0])
            y.append(calpha[1])
            z.append(calpha[2])
            orient = d.binder_orientation()
            u.append(orient[0])
            v.append(orient[1])
            w.append(orient[2])
            if draw_secor:
                secor = d.binder_secondary_orientation()
                u2.append(secor[0])
                v2.append(secor[1])
                w2.append(secor[2])
        # draw the C_alpha -> C_beta vectors of the residues
        ax.quiver(x, y, z, u, v, w, color=col, arrow_length_ratio=0.5, linewidth=1.0)
        if draw_secor:
            # draw the C_alpha -> C_O vectors of the residues
            ax.quiver(x, y, z, u2, v2, w2, color=col, arrow_length_ratio=0.0, linewidth=1.0)
    plt.show()


def visualize_acomplex(acomplex: LigandBinderComplex, dimensions: float = 10.0, draw_secor: bool = True,
                       filter_score: bool = True, **score_args) -> None:
    """Shows a 3d quiver plot that contains the designed ligand (including backbone), the designed binder residues (in
    black) and the atlas datapoints associated with the corresponding ligand residue types.
    :param acomplex: the complex to visualize
    :param dimensions: sets the display limit of the plot. This spacing will be applied to every dimension of the
            outermost residue.
    :param draw_secor: if True, secondary orientation vectors are drawn (without arrowheads)
    :param filter_score: if True, only those aligned residues whose score is below 0.0 are displayed
    :param score_args: additional keyword arguments that will be passed on to prediction_matrix.proximity_score
    """
    # plot the data as matplotlib 3d quiver
    fig = plt.figure()
    plt.suptitle(acomplex.ligand.get_title())
    ax: Axes3D = fig.gca(projection='3d')
    min_dim = acomplex.ligand.get_min_dim()
    max_dim = acomplex.ligand.get_max_dim()
    ax.set_xlim(min_dim[0] - dimensions, max_dim[0] + dimensions)
    ax.set_ylim(min_dim[1] - dimensions, max_dim[1] + dimensions)
    ax.set_zlim(min_dim[2] - dimensions, max_dim[2] + dimensions)
    bbx: List[float] = []
    bby: List[float] = []
    bbz: List[float] = []
    bbu: List[float] = []
    bbv: List[float] = []
    bbw: List[float] = []
    bbprev: DesignedLigandResidue = None
    # process ligand's backbone-residue structure
    for ligres in acomplex.ligand.residues:
        # create an internal coordinate system for the considered residue
        licoor = get_icoor_from_aresidue(ligres)
        # show C_alpha->C_beta vector of each ligand residue
        lca = licoor.internal_to_external(Vector(0.0, 0.0, 0.0), False)
        lor = licoor.internal_to_external(Vector(4.5, 0.0, 0.0), True)
        # draw ligand residue
        ax.quiver([lca[0]], [lca[1]], [lca[2]], [lor[0]], [lor[1]], [lor[2]],
                  color=get_color(ligres.restype), arrow_length_ratio=0.4, linewidth=9.0)
        # fill backbone info
        if bbprev is not None:
            bbx.append(ligres.calpha()[0])
            bby.append(ligres.calpha()[1])
            bbz.append(ligres.calpha()[2])
            bbu.append(bbprev.calpha()[0] - ligres.calpha()[0])
            bbv.append(bbprev.calpha()[1] - ligres.calpha()[1])
            bbw.append(bbprev.calpha()[2] - ligres.calpha()[2])
        bbprev = ligres
    # draw backbone
    ax.quiver(bbx, bby, bbz, bbu, bbv, bbw, color="#000000", arrow_length_ratio=0.0, linewidth=9.0)
    # process binder's designed residue positions
    for binres in acomplex.binder.residues:
        # show C_alpha->C_beta vector of each designed binder residue
        bca = binres.calpha()
        bor = binres.orientation() ** 3.0
        # draw binder residue
        ax.quiver([bca[0]], [bca[1]], [bca[2]], [bor[0]], [bor[1]], [bor[2]],
                  color="#606060", arrow_length_ratio=0.7, linewidth=9.0)
        if draw_secor:
            bso = binres.secondary_orientation() ** 3.0
            ax.quiver([bca[0]], [bca[1]], [bca[2]], [bso[0]], [bso[1]], [bso[2]],
                      color="#606060", arrow_length_ratio=0.0, linewidth=9.0)
    # process the residues of the atlas
    for binrestype in color_map:
        # for each potential interacting binder residue, paint the quivers in
        # individual colors.
        relevant_datapoints = [d for d in acomplex.aligned_residues if d.binder_restype == binrestype]
        col = color_map[binrestype]
        x: List[float] = []
        y: List[float] = []
        z: List[float] = []
        u: List[float] = []
        v: List[float] = []
        w: List[float] = []
        u2: List[float] = []
        v2: List[float] = []
        w2: List[float] = []
        for d in relevant_datapoints:
            for binder_res in acomplex.binder.residues:
                score = proximity_score(binder_res, d, **score_args)
                if not filter_score or score > 0.0:
                    orient = d.orientation()
                    secor = d.secondary_orientation()
                    x.append(d.calpha()[0])
                    y.append(d.calpha()[1])
                    z.append(d.calpha()[2])
                    u.append(orient[0])
                    v.append(orient[1])
                    w.append(orient[2])
                    if draw_secor:
                        u2.append(secor[0])
                        v2.append(secor[1])
                        w2.append(secor[2])
        # draw the C_alpha -> C_beta vectors of the residues
        ax.quiver(x, y, z, u, v, w, color=col, arrow_length_ratio=0.5, linewidth=1.0)
        if draw_secor:
            # draw the C_alpha -> C_O vectors of the residues
            ax.quiver(x, y, z, u2, v2, w2, color=col, arrow_length_ratio=0.0, linewidth=1.0)
    plt.show()


def visualize_pmatrix(pmatrix: BinderPredictionMatrix) -> None:
    """Visualizes a given binder predictrion matrix as a stacked bar plot. Each stacked bar represents a residue of the
    designed binder. The sub-bars represent potential target residues each, where the type is coded by the bar color,
    the height encodes the score, and the width encodes the certainty.
    :param pmatrix: the prediction matrix to visualize
    """
    maxwidth = max(1, max([max([col.certainty for col in pline.columns], default=1) for pline in pmatrix.lines]))
    for ind, pline in enumerate(pmatrix.lines):
        bot = 0.0
        for pcolumn in pline.columns:
            col = get_color(pcolumn.residue_type)
            # height of the bar reflects score, width n_origins.
            plt.bar([ind], [pcolumn.score*100.0], pcolumn.certainty*(0.9/float(maxwidth)), color=col, bottom=[bot])
            bot += pcolumn.score*100.0
        # plot bars for reference amino acids
        reference_col = get_color(pline.original_restype)
        plt.bar([ind], 10.0, 0.9, color=reference_col, bottom=[-20])
    plt.xticks(range(0, len(pmatrix.lines)),
               [pline.original_restype + str(pline.residue_id) for pline in pmatrix.lines], rotation=60)
    plt.yticks([-15] + list(range(0, 100, 10)), ["Reference"] + list([f"{i}%" for i in range(0, 100, 10)]))
    plt.show()


def visualize_atlas_stats(atlas: Atlas, method: str = "plotly", **kwargs) -> None:
    """Visualizes the datapoints of an atlas as a stacked bar plot. Every stack represents a ligand residue type, and
    the individual bars represent the interacting binder residue types.

    :param atlas: the atlas to visualize
    :param method: Either "plotly" or "matplotlib". If not given plotly matplotlib is used.
    """
    if method == "plotly":
        return visualize_atlas_stats_plotly(atlas=atlas, **kwargs)
    else:
        return visualize_atlas_stats_matplotlib(atlas=atlas, **kwargs)


def visualize_atlas_stats_plotly(atlas: Atlas, ligand_per_binder: bool = False, color_mode: str = "0") -> None:
    """Visualizes the datapoints of an atlas as a stacked bar plot. Every stack represents a ligand residue type, and
    the individual bars represent the interacting binder residue types.
    param atlas: The atlas to take the data from
    param invert: if True, the plot maps ligand (x axis) to binder (y axis), else vice versa.
    """
    stats = atlas.get_stats(ligand_per_binder=ligand_per_binder)
    df = pd.DataFrame.from_dict(stats)
    df = df.transpose()
    cmap = get_color_map(color_mode)
    fig = px.bar(df,
                 x=df.index.name,
                 y=canonical_amino_acids,
                 color_discrete_map=cmap,
                 title="Wide-Form Input")

    fig.update_layout(
        title="Atlas statistics",
        xaxis_title="ligand amino acid" if not ligand_per_binder else "binder amino acid",
        yaxis_title="ligand amino acid" if ligand_per_binder else "binder amino acid")
    fig.show()


def visualize_atlas_stats_matplotlib(atlas: Atlas, ligand_per_binder: bool = False, color_mode: str = "0") -> None:
    """Visualizes the datapoints of an atlas as a stacked bar plot. Every stack represents a ligand residue type, and
    the individual bars represent the interacting binder residue types.
    param atlas: The atlas to take the data from
    param ligand_per_binder: if True, inverts the plot: ligand residues (y axis) per binder residues (x axis),
     else vice versa.
    """
    # TODO viualize alternative aas too...
    stats = atlas.get_stats(ligand_per_binder)
    max_height = max([sum([val for val in line.values()]) for line in stats.values()])
    for ind, ligres in enumerate(canonical_amino_acids):
        if ligres not in stats:
            continue
        bot = 0.0
        for binres in sorted(canonical_amino_acids, key=lambda aa: stats[ligres].get(aa, 0.0), reverse=True):
            if binres not in stats[ligres]:
                continue
            count = stats[ligres][binres]
            height = float(count) / float(max_height) * 100.0
            bincol = get_color(binres, color_mode)
            plt.bar([ind], [height], 0.8, color=bincol, bottom=[bot])
            bot += height
        ligand_col = get_color(ligres, color_mode)
        plt.bar([ind], 10.0, 0.9, color=ligand_col, bottom=[-20.0])
    plt.xticks(range(0, 20), canonical_amino_acids, rotation=60)
    plt.yticks([-15] + list(range(0, 101, 10)),
               ["Binder" if ligand_per_binder else "Ligand"] +
               list([f"{i}" for i in range(0, max_height, int(max_height / 10))]))
    plt.show()


def visualize_pocket_atlas(pocket: Pocket, method: str = "plotly", **kwargs):
    """Shows the datapoints of a pocket atlas with a plotly plot (within web broswer) or matplotlib (separate window).
        For additional parameters such as 'ligand_restype' compare the subfunctions used here.

        :param pocket: the pocket to visualize
        :param method: Either "plotly" or "matplotlib". If not plotly matplotlib is used.
    """
    return visualize_atlas(atlas=pocket.to_atlas(), method=method, **kwargs)


class DatapointNotFound(Exception):
    """raised if datapoint is not found in atlas."""

    def __init__(self, message="The datapoint you entered cannot be found within this atlas!"):
        self.message = message
        super().__init__(self.message)


def get_best_connections_for_dp(atoms_: List[Tuple[str, Vector, str]], restype_: str, trace_name: str, color_: str) \
       -> Iterator[Dict]:
    """
    Returns an Iterator including a dict with a trace for plotting (plotly?) bonds between the given atoms.
    Only works for canonical amino acids! # TODO non-canonical typical ones ?
    The backbone atoms will be grey, the sidechain atoms will be the assigned color.

    :param atoms_: The atoms that should be connected
    :param restype_: The residue type the atoms belong to
    :param trace_name: The name of the returned trace
    :param color_: The color of sidechain atoms
    :return: Iterator of Dicts including information about the trace to plot
    """

    def get_lines(canonical: bool):
        if not canonical:
            for bond in get_connected_atom_pairs(atoms_ + [("CA", Vector(0, 0, 0))]):
                yield {'type': "scatter3d", 'mode': 'lines',
                       'name': trace_name,
                       'x': [bond[0][0], bond[1][0]], 'y': [bond[0][1], bond[1][1]],
                       'z': [bond[0][2], bond[1][2]],
                       'line': {'width': 20, 'color': color_}}
        else:
            # For every connected line safed in best_connection define a bond of vectors
            for connected_atom_names in best_connection:
                bond_vectors = []
                if not len(connected_atom_names):
                    continue

                # Look for right atom names to match the predefined atoms in a line
                for connected_atom in connected_atom_names:
                    for (atom_name, atom_vector, origin) in atoms_:
                        # If name is found extend the connected line with vector
                        if atom_name == connected_atom:
                            bond_vectors.append(atom_vector)
                            break

                # TODO make warning? Assertion fails it atoms are missing in the input structure
                #  (e.g. Ser w/o OG - 3jb9 i/S9)
                # assert len(connected_atom_names) == len(bond_vectors)

                current_color = "#BBBBBB" if connected_atom_names[0] == "N" else color_

                yield {'type': "scatter3d", 'mode': 'lines', 'hoverinfo': 'none', 'name': trace_name,
                       'x': [b[0] for b in bond_vectors],
                       'y': [b[1] for b in bond_vectors], 'z': [b[2] for b in bond_vectors],
                       'line': {'width': 20, 'color': current_color}}

    # Set residue types and get atoms of residue
    try:
        # For canonical amino acids, the lines are predefined to use less traces (faster)
        best_connection = [("N", "CA", "C", "O"), connections[restype_]]
    except KeyError:
        return get_lines(canonical=False)

    return get_lines(canonical=True)


def get_connected_atom_pairs(all_atoms: Union[List[Tuple[str, Vector, str]], List[Tuple[str, Vector]]]):
    for a1, a2 in combinations(all_atoms, 2):
        if 0.5 < (a1[1] - a2[1]).norm() < 2:
            yield a1[1], a2[1]


def visualize_single_pocket(pocket: Pocket, datapoint: AtlasDatapoint, color_by_element: bool = False):

    # Filter pocket atlas for those who share the ligand origin with the datapoint parameter
    filtered_atlas = Atlas([dp for dp in pocket.to_atlas().datapoints if dp.ligand_origin == datapoint.ligand_origin])
    if not len(filtered_atlas.datapoints):
        raise DatapointNotFound
    dps = filtered_atlas.datapoints

    # define special (non-C) heavy atoms and their color
    hatom_names = ("N", "S", "P", "O")
    hatom_atoms = dict(zip(hatom_names, (defaultdict(list) for _ in hatom_names)))
    hatoms_colours = dict(zip(hatom_names, ("#0000DD", "#DDDD00", "#AB0000", "#DD0000")))

    # collect all coordinates, residue types and the size of the marker datapoints in lists.
    x = []
    y = []
    z = []
    restypes = []
    size = []
    bonds_list = []

    # For every datapoint save coords, restype and size (Ca and Cb)
    for dp in dps:
        bin_atoms = []

        # Define atom bubbles
        for atom_name, atom_vector in dp.binder_atoms.items():
            bin_atoms.append((atom_name, atom_vector, dp.binder_origin))
            # if color_by_element, save heavy atoms differently to plot different atom colours.
            if not (color_by_element and atom_name[0] in hatom_names):
                x.append(atom_vector[0])
                y.append(atom_vector[1])
                z.append(atom_vector[2])
                restypes.append(dp.binder_restype)
                size.append(8)
            else:
                hatom_atoms[atom_name[0]]['x'].append(atom_vector[0])
                hatom_atoms[atom_name[0]]['y'].append(atom_vector[1])
                hatom_atoms[atom_name[0]]['z'].append(atom_vector[2])

        # define connections between residue atoms for plotting of covalent bonds.
        for binder_line in get_best_connections_for_dp(bin_atoms, restype_=dp.binder_restype,
                                                       trace_name=dp.binder_restype + " bond",
                                                       color_=get_color(dp.binder_restype)):
            bonds_list.append(binder_line)

    # Define dataframe for markers of binder atoms
    df = pd.DataFrame.from_dict({
        "x": x,
        "y": y,
        "z": z,
        "restype": restypes,
        "size": size,
    })

    # concatenate color dictionaries to address Ca (full color) and Cb (light color)
    color_map_light_cb = {k + " cb": v for k, v in color_map_light.items()}
    color_map_full = {**color_map, **color_map_light_cb}

    # Define figure and include all atoms (if color_by_element: only carbon atoms) as markers
    fig = px.scatter_3d(df, x='x', y='y', z='z',
                        color='restype', color_discrete_map=color_map_full)

    # Special heavy atom coloring: O, N, S, P as markers
    if color_by_element:
        for atoms_name, atoms in hatom_atoms.items():
            fig.add_trace(go.Scatter3d(x=atoms['x'], y=atoms['y'], z=atoms['z'], name=atoms_name,
                                       mode='markers', showlegend=False,
                                       marker={'line': {'width': 20, 'color': hatoms_colours[atoms_name]}}))

    # Define bond lines
    for bond in bonds_list:
        fig.add_trace(go.Scatter3d(x=bond['x'], y=bond['y'], z=bond['z'], name=bond['name'], mode='lines',
                                   showlegend=False, line=bond['line'], hoverinfo='skip'))
    fig.show()
