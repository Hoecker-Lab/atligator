"""This module contains functions providing old-school static HTML output for the integrated visualization of binding
pockets obtained by the pocket miner.

:Deprecated: This module has been deprecated by the atligator web interface. See project atligator_web

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-06-06
"""
import datetime
from pathlib import Path
from typing import List, Dict

from Bio.SeqUtils import seq1

from atligator.atlas import Atlas
from atligator.pocket_miner import Pocket


def generate_html_pocket_index(output_folder: str, atlas: Atlas, pockets: Dict[str, List[Pocket]]) -> None:
    """Generates a HTML index page showing summarizing information on the atlas and the mined pockets. The page also
    contains links to individual pocket descriptions which are supposed to be generated by the function
    generate_html_pocket_description.
    :param output_folder: where the index.html file is generated
    :param atlas: the atlas the pockets were mined from
    :param pockets: the mined pockets
    """
    html: List[str] = list()
    html.append('<!DOCTYPE html>')
    html.append('<html lang="en">')
    html.append('<head>')
    html.append('    <meta charset="utf-8">')
    html.append('</head>')
    html.append('</body>')
    html.append('    <h1>Atlas Binding Pockets</h1>')
    html.append(f"    <p>generated by ATLIGATOR's pocket miner script on {datetime.datetime.now()}.</p>")
    html.append('    <h2>General Information</h2>')
    html.append(f"    <p>number of datapoints in the atlas: {len(atlas.datapoints)}.</p>")
    html.append('    <h2>Binding Pockets</h2>')
    html.append('    <p>select one of the ligand residue types below:</p>')
    html.append('    <ul>')
    for lig_restype, pocketlist in sorted(pockets.items(), key=lambda e: e[0]):
        html.append(f'        <li><a href={lig_restype}.html>{lig_restype}</a> ({len(pocketlist)} pockets '
                    f'from {len([d for d in atlas.datapoints if d.ligand_restype == lig_restype])} datapoints)</li>')
    html.append('    </ul>')
    html.append('</body>')
    html.append('</html>')
    with open(Path(output_folder) / 'index.html', 'w') as fo:
        for line in html:
            fo.write(line + '\n')


def generate_html_pocket_description(output_folder: str, ligand_restype: str, pockets: List[Pocket],
                                     visualize_raw_pdb: bool, visualize_clustered_pdb: bool) -> None:
    """Generates a detailed HTML page showing pockets connected to a specfic ligand residue type. Optionally, the
    superimposed and clustered PDB structures are presented in a NGL.js viewer. The name of the generated file depends
    on the ligand residue type, e.g. LYS.html
    :param output_folder: will contain the generated file.
    :param ligand_restype: three-letter code of the residue type
    :param pockets: list of pockets belonging to the residue type
    :param visualize_raw_pdb: whether to visualize the superimposed PDB structure
    :param visualize_clustered_pdb: whether to visualize the clustered PDB structure
    """
    html: List[str] = list()
    html.append('<!DOCTYPE html>')
    html.append('<html lang="en">')
    html.append('<head>')
    html.append('    <meta charset="utf-8">')
    html.append('</head>')
    html.append('</body>')
    html.append(f'    <h3>{ligand_restype} Pockets</h3>')
    if visualize_raw_pdb or visualize_clustered_pdb:
        html.append('    <script src="https://unpkg.com/ngl"></script>')
        html.append(f'    <p>in the visualizations below, C atoms of the ligand are presented in white,'
                    f' of the binder in gray, respectively. CA atoms are in green.</p>')
        html.append(f'    <p>hint: double-click onto the structures for full-screen representation.</p>')
    for i, pocket in enumerate(pockets):
        pocket_code = f"{seq1(ligand_restype)}_{pocket.pocket_id()}_{int(pocket.support * 100.0)}"
        html.append(f'    <h4>Pocket {i + 1}: {pocket.pocket_id()}</h4>')
        html.append(f'    <p>from {pocket.n_ligand_origins()} datapoints; support: {pocket.support:.3f}</p>')
        if visualize_raw_pdb:
            html.append('    <p>superimposed pocket representation:</p>')
            html.append('    <script>')
            html.append('        document.addEventListener("DOMContentLoaded", function () {')
            html.append(f'            var stage = new NGL.Stage("viewport-{pocket_code}");')
            html.append('            stage.viewer.container.addEventListener("dblclick", function () {')
            html.append('                stage.toggleFullscreen();')
            html.append('            });')
            html.append('            var colorScheme = NGL.ColormakerRegistry.addSelectionScheme([')
            html.append('                ["darkgreen", ".CA"],')
            html.append('                ["white", ":A and _C"],')
            html.append('                ["gray", "(not :A) and _C"],')
            html.append('                ["red", "_O"],')
            html.append('                ["blue", "_N"],')
            html.append('                ["yellow", "_S"],')
            html.append('                ["green", "*"]')
            html.append(f'            ], "{pocket_code}");')
            html.append(f'            stage.loadFile("{pocket_code}.pdb").then(function (o) ' + '{')
            html.append('                o.addRepresentation("ball+stick", {color: colorScheme});')
            html.append('                o.autoView();')
            html.append('            });')
            html.append('        });')
            html.append('    </script>')
            html.append(f'    <div id="viewport-{pocket_code}" style="width:400px; height:300px;"></div>')
        if visualize_clustered_pdb:
            html.append('    <p>pocket centroid representation:</p>')
            html.append('    <script>')
            html.append('        document.addEventListener("DOMContentLoaded", function () {')
            html.append(f'            var stage = new NGL.Stage("viewport-{pocket_code}_clustered");')
            html.append('            stage.viewer.container.addEventListener("dblclick", function () {')
            html.append('                stage.toggleFullscreen();')
            html.append('            });')
            html.append('            var colorScheme = NGL.ColormakerRegistry.addSelectionScheme([')
            html.append('                ["darkgreen", ".CA"],')
            html.append('                ["white", ":A and _C"],')
            html.append('                ["gray", "(not :A) and _C"],')
            html.append('                ["red", "_O"],')
            html.append('                ["blue", "_N"],')
            html.append('                ["yellow", "_S"],')
            html.append('                ["green", "*"]')
            html.append(f'            ], "{pocket_code}_clustered");')
            html.append(f'            stage.loadFile("{pocket_code}_clustered.pdb").then(function (o) ' + '{')
            html.append('                o.addRepresentation("ball+stick", {color: colorScheme});')
            html.append('                o.autoView();')
            html.append('            });')
            html.append('        });')
            html.append('    </script>')
            html.append(f'    <div id="viewport-{pocket_code}_clustered" style="width:400px; height:300px;"></div>')
    html.append('    </ul>')
    html.append('    <p><a href="index.html">back to index</a></p>')
    html.append('</body>')
    html.append('</html>')
    with open(Path(output_folder) / f"{ligand_restype}.html", 'w') as fo:
        for line in html:
            fo.write(line + '\n')
