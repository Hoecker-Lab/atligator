# ATLIGATOR

ATLIGATOR is created to analyse protein-protein or protein-peptide interactions.

# Installation

You can find the atligator in the Python Packaging Index.
Thus, to install atligator in your python 3.7+ environment just type:

```pip install atligator```

#### Installation Notes
If installation fails because of the build process of dependencies use the '--only-binary=:all:' flag:
```
python -m pip install --only-binary=:all: atligator
```

For the compilation of scipy and numpy you might need to install the compilers first 
(or try to use only wheel files, see above! This might depend on your OS environment):
```
sudo apt-get install gfortran libopenblas-dev liblapack-dev
```

#### What is included?

The package from PyPI delivers all options to work with ATLIGATOR as a python API.
The CLI scripts can be received from github as an aditional entry point. 

# Usage

You can use ATLIGATOR either with the predefined scripts or the API. 
Most use cases can be accessed via both options.
Note: Most scripts will also give you a help option using the `-h` flag.

### Citation

If you use ATLIGATOR please cite it. The article is currently being published, so stay tuned
for the correct citation.

### Why?

Why you should use ATLIGATOR?

#### Atlases & Alarms
With ATLIGATOR you can build up interaction databases (called atlases) which store the interactions between pairs of
amino acid residues. By combining many of these pairwise interactions you can get a better understanding of how amino 
acids tend to interact with each other.

All the input structures can be selected manually or by using fold-specific SCOPe-identifiers to only use protein 
structures with similar architecture. Additionally, input structures can be preprocessed to only contain relevant 
protein chains (a central binder chain and ligand chains in close distance) and filtered for a certain secondary
structure content.

#### Pockets
Atlas information can be exploited further by using the pocket mining functionality. Pockets are frequent groups of
interactions which can be detected within an atlas. With Pockets you can directly see potential binding pockets for
designing amino acid binders in a structural context.

#### Design
Pockets can also be grafted onto scaffold proteins to get starting points for new protein or peptide binding abilities 
in existing proteins.

## Structures and Preprocessing

To be able to build a new atlas database prepare a directory with input pdb structures.

#### Downloading Structures
If you want to use ATLIGATORS SCOPe-matching functionality, either use the predefined script:

*to download all files matched with SCOPe id a.118.11 to the directory ./structures*
```shell script
python3 download_pdbs_by_scop_query.py -o ./structures a.118.11
```
or write a python script/use the interpreter:

```python
from atligator.pdb_util import download_pdbs_by_scop_query

download_pdbs_by_scop_query(query="a.118.11",
                            dest_dir="./structures")
# If you have a local copy of the pdb database you can try setting the parameter: pdb_db_path
# If you have a local copy of the SCOPe database you can try setting the parameter: scop_dir_cla_file
# Both will prevent downloading those files from the web.
```

or for single PDB files:

```python
from atligator.pdb_util import get_pdb_structure

get_pdb_structure("1z5s", dest_dir="./strucures/")
```

#### Preprocessing Structures
To preprocess the input structure files you can use:

```shell script
python3 process_chains.py -o ./structures/processed -pro -d 8.0 -mb 30 -ml 3 structures/*.pdb
```

*The -d option defines the maximum distance two atoms of a residue pair can have to still be included. It's recommended
to use at least the maximum distance you use for atlas generation lateron.*

*The options -mb and -ml define the length of the polypeptide chains to be considered for binder (central chain) or 
ligand (any other chain it interacts with).*

or write a python script/use the interpreter:

```python
from alarms.chain_processing import MultiProcessChainProcessing
from glob import iglob

pdbs = [x for x in iglob("./structures/*.pdb")]
proc = MultiProcessChainProcessing(pdbs=pdbs, output_path="./structures/processed/", min_binder_len=30, 
                                   min_ligand_len=3, max_distance=8.0, n_workers=2, progress=True, quiet=False)
proc.main()
```

To process single files: 

```python
from alarms.chain_processing import process_pdb

process_pdb(pdb="./structures/7ev7.pdb", output_path="./structures/processed/", min_binder_len=30, 
            min_ligand_len=3, max_distance=8.0, quiet=False)
```

#### Filtering Structures

You can filter structures by secondary structure content:
1. total alpha helix
2. total beta sheet
3. alpha helix/beta sheet ratio
4. beta sheet/alpha-helix ratio

*E.g. at least 60% alpha helical residues*

```shell script
python3 select_by_sec_struc.py -o ./structures/processed/selected -at 0.6 -cp structures/processed/*.pdb
```

*the `-cp` flag enables copying the filtered files into the output directory. Otherwise there will be 
only a text file 'filtered_objects.txt' including all selected residues.*

Within python:

```python
from alarms.construct_selection import SelectionProcess
from glob import iglob

pdbs = list(x for x in iglob("./structures/processed/*.pdb"))
select = SelectionProcess(pdbs=pdbs)

# You can retrieve the secondary structure content:
print(select.check_pdbs_for_2s_content())

# Or you can directly filter the files
filtered = select.filter_pdbs(alphatotal=0.8)
print("These fall under the criteria:", filtered)
```

Or to only retrieve the secondary structure content:

```python
from alarms.construct_selection import check_pdbs_for_2s_content
pdbs = list(x for x in iglob("./structures/processed/*.pdb"))
check_pdbs_for_2s_content(pdbs=pdbs, quiet=False)
```


## Atlas

### Atlas Generation

To generate an atlas from pdb input files, just:

```shell script
python3 generate_atlas.py ./atlas.atlas ./structures/processed/*.pdb -b
```

*The `-b` flag excludes interactions of ligand backbone atoms and is recommended to strenghen the influence of the 
corresponding ligand residue side chain.*

Within python:

```python
from atligator.atlas import generate_atlas
from glob import iglob

pdbs = list(x for x in iglob("./structures/processed/*.pdb"))
atlas = generate_atlas(filenames=pdbs, skip_bb_atoms=True)
```

### Atlas Objects

Apart from visualisation or further processing atlas object can be already inspected in python:

```python
# Find a list of AtlasDatapoints in:
atlas.datapoints
# These datapoints contain all atoms of the residue pair (ligand and binder) as well as their residue type and their
# origin

# You can also filter an atlas to only contain specific residue type in ligand and binder position.
filtered_atlas = atlas.filter(ligand_restype = None, binder_restype = None)
```

### Atlas Visualisation

Visualisation of an atlas can be of two types: structural or statistical. 

##### Statistics

The statistical visualization contains a plot of the atlas content, where the ligand residue type is plotted against
the binder residue type and vice versa. 

To generate an atlas from pdb input files, just:


```shell script
python3 visualize_atlas_stats.py ./atlas.atlas -m plotly
```

Within python:
```python
from atligator.visualization import visualize_atlas_stats

visualize_atlas_stats(atlas, method="plotly", ligand_per_binder=False)
```

*There is two ways of plotting atlas statistics: plotly and matplotlib. Plotly is using the browser,
matplotlib creates a new window.*
*To transform the plot to see which ligand residue types interact with the binder residue types, set ligand_per_binder 
to True or use the `-i` flag with the script*.


##### Structural Visualization: Atlas or Atlas subset

The structural representation delivers a good overview about what an atlas is composed of.
The default representation will always draw a 3D plot where only the Calpha and Cbeta atoms of the ligand residue
and the binder residues are present as bigger (Calpha) and smaller (Cbeta) bubbles. 

The ligand residue is typically positioned in the center of the plot to demonstrate where the binder residues are
relative to the ligand residue.


```shell script
python3 visualize_atlas.py ./atlas.atlas -m plotly -l ALA -b TYR
```

Within python:

```python
from atligator.visualization import visualize_atlas

# You can define the ligand residue type
visualize_atlas(atlas, method="plotly", ligand_restype="ALA")

# Or both residue types
visualize_atlas(atlas, method="plotly", ligand_restype="ALA", binder_restype="TYR")
```

### File transfer

#### Not recommended for data exchange: pickle
Atlases can be stored with pickle (https://docs.python.org/3/library/pickle.html) by deserialisation of the atlas
python object: 

```python
import pickle
# for import
with open("./atlas.atlas", "rb") as rf:
    atlas = pickle.load(rf)

# for export
with open("./atlas.atlas", "wb") as wf:
    pickle.dump(atlas, wf)
```

However, in any environment where files are exchanged with others it's not recommended to use pickle. Pickle files
can contain any python code and are executed after reading in the file. Thus, untrusted sources could hide an unpleasant
surprise in their files and we need an alternative.

#### json Atlases 

We have the option to import and export Atlases as json files which are parsed internally. By doing so, we have clear 
text files for immediate inspection and prevent undesired code execution.

```python
from atligator import atlas
# for import
with open("./atlas.json", "r") as rf:
    atlas = atlas.Atlas.from_json(rf)

# for export
with open("./atlas.json", "w") as wf:
    atlas.to_json_file(wf)
```

## Pockets

### Pocket Mining

To mine pockets from an existing atlas:

```shell script
python3 apply_pocket_miner.py ./atlas.atlas -o ./
```

Within python:
```python
from atligator.pocket_miner import mine_pockets

pockets = mine_pockets(atlas)
```
*Additional parameters can be used for example to define how important the size (cardinality) of the pocket is, or 
how many pockets per ligand residue type are found at maximum.*

### Pocket Handling

Apart from visualization or grafting Pockets can be already inspected in python:

```python
# pockets is a dictionary of ligand residue types. To access the underlying list of pockets just:
pockets['TYR']
# or print all content:
for residue_type, pocket_list in pockets.items():
    print(residue_type, pocket_list)
```

### Pocket Visualization

To visualize pockets you can choose between plotting the pocket atlas or a single pocket.
- The pocket atlas is a filtered version of the original atlas. It only includes datapoints that 
can be found as a part of the pocket.   
For example: You mined a Ile to REST pocket - including Ile as the ligand and Arg, Glu, Ser and Thr as binder residues. 
The corresponding pocket atlas will only contain atlas datapoints that include Ile as a ligand and one of the four 
binder residues as a binder residue. Additionally, only those datapoints that are found in combination with the three 
other datapoints (the other three binder residue types) are taken into account (and filtered based on the clustering 
parameters during pocket mining).

- A single pocket is a multi-residue structure including one ligand residue and all binder residues surrounding it. 
It will only contain those residues that natively interact in one of the input structures. Thus, it is just a strictly
filtered representation of one input structure.

**A pocket atlas** is plotted as an atlas, including Calpha and Cbeta bubbles for the central ligand residue
as well as binder residues around.

```shell script
python3 visualize_pocket.py ./pockets.pockets -l TYR -p RLD
```
*The `-l` flag defines the  ligand residue type and the `-p` flag defines the itemset of the desired pocket 
in one-letter-code (In this case Arg, Leu and Asp).*

Within python:

```python
from atligator.visualization import visualize_pocket_atlas

# First we could check which Tyr pockets we have
print(pockets["TYR"])

tyr_pocket = pockets["TYR"][0]
# You can define the ligand residue type
visualize_pocket_atlas(tyr_pocket, method="plotly")
```



**A single pocket** is plotted with all atoms to enable viewing the interaction in all details.

Visualizing a single pocket is not yet available as a script.

Within python:

```python
from atligator.visualization import visualize_single_pocket

tyr_pocket = pockets["TYR"][0]
# First we have to pick one datapoint of this pocket. We could for example determine the most representative member of 
# one cluster (in this case the first cluster):
# Let's first see which binder residue type is present in this cluster
print(tyr_pocket.clusters[0].binder_restype)
# Let's pick the most representative datapoint
most_representative_member = tyr_pocket.clusters[0].get_most_representative_member(2.0, 1.0, 1.0)

# we could also just pick any member of any cluster
any_member = tyr_pocket.clusters[0].members[0]

# Let's visualize this single pocket with plotly
visualize_single_pocket(pocket=tyr_pocket, datapoint=most_representative_member)
```

### File Transfer

#### Not recommended for data exchange: pickle
Pocket collections can be stored with pickle (https://docs.python.org/3/library/pickle.html) by deserialisation of the 
python object: 

```python
import pickle
# for import
with open("./pockets.pockets", "rb") as rf:
    pockets = pickle.load(rf)

# for export
with open("./pockets.pockets", "wb") as wf:
    pickle.dump(pockets, wf)
```

However, in any environment where files are exchanged with others it's not recommended to use pickle. Pickle files
can contain any python code and are executed after reading in the file. Thus, untrusted sources could hide an unpleasant
surprise in their files and we need an alternative.

#### json Pocket Collections

We have the option to import and export Pocket Collections as json files which are parsed internally. By doing so, we 
have clear text files for immediate inspection and prevent undesired code execution.

```python
from atligator import pocket_miner
# for import
with open("./pockets.json", "r") as rf:
    pockets = pocket_miner.json_to_pockets(rf)

# for export
with open("./pockets.json", "w") as wf:
    pocket_miner.pockets_to_json_file(pockets, wf)
```


## Pocket Grafting

To apply the knowledge gained by ATLIGATOR atlas and pockets directly to your design pockets
can be grafted onto new scaffolds.

**Note**: All grafted pockets will be based on natural pockets from the input structures and the side chain
conformations (rotamers) will remain as they were in the original structure. Thus, the grafted rotamers will **not**
perfectly fit into the new environment. It is highly recommended to minimize these rotamers with a common protocol
(e.g. Rosetta fixbb).

### Grafting One Binding Pocket

Individual pockets can be grafted onto scaffold proteins directly. Just define a protein of choice, mutable binder 
residue positions, a ligand residue (the basis for the grafted pocket) and the pocket you want to graft.
ATLIGATOR will automatically pick the best matching members of this pocket and the best fitting positions to mutate.

```python
from atligator.grafting import graft_pocket_onto_scaffold

# First define which positions in which chain are mutable
design_positions = {"C": ["510", "513", "514", "515", "517", "518", "523", "524", "557", "558", "559", "561", "562"]}
# Pick a pocket of the desired residue type and graft it onto design positions of scaffold - based on ligand_residue
# If you define an output_name, the new pdb file will be stored at this path.
design = graft_pocket_onto_scaffold(pocket=pockets["THR"][1], scaffold_name="./1z5s.pdb", ligand_residue=("A", "131"), 
                                    design_positions=design_positions, output_name="./design.pdb")
```

### Quickgraft
#### Quickgraft - single pocket

Quickgraft allows to graft the pocket in a pocket collection which fits the best on the available design positions.
In this case you don't have to select the pocket yourself, but the residue type the ligand 
residue should be mutated to.

```python
from atligator.grafting import graft_best_pocket_onto_scaffold

# First define which positions in which chain are mutable
design_positions = {"C": ["510", "513", "514", "517", "523", "558", "559"]}
# Define which residue (from which chain) should act as a ligand residue and which mutation should be applied.
ligand_residues = {"A": {"131": "THR"}}
# Graft the best matching pocket onto design_positions of scaffold - based on ligand_residues and the new ligand restype
# If you define an output_name, the new pdb file will be stored at this path.
design, graft = graft_best_pocket_onto_scaffold(pockets=pockets, scaffold_name="./1z5s.pdb", 
                                                ligand_residues=ligand_residues, 
                                                design_positions=design_positions, output_name="./design.pdb")
print(graft)
```

#### Quickgraft - multi pocket

If you want to mutate more than one ligand residue and want to get matching pocket grafts for all ligand mutations 
ATLIGATOR offers a multi pocket grafting option.

```python
from atligator.grafting import graft_best_pocket_onto_scaffold

# First define which positions in which chain are mutable
design_positions = {"C": ["510", "513", "514", "517", "523", "558", "559"]}
# Define which residues (from which chain) should act as ligand residues and which mutations should be applied.
ligand_residues = {"A": {"131": "THR", "132": "ASN"}}
# Graft the best matching pocket onto design_positions of scaffold - based on ligand_residue and the new ligand_restype
# If you define an output_name, the new pdb file will be stored at this path.
design, graft = graft_best_pocket_onto_scaffold(pockets=pockets, scaffold_name="./1z5s.pdb", 
                                                ligand_residues=ligand_residues, 
                                                design_positions=design_positions, output_name="./design.pdb")
print(graft)
```

#### Quickgraft - more grafting solutions

If you want multiple grafts based on your settings use:

```python
from atligator.grafting import multi_graft_best_pocket_onto_scaffold

# First define ligand residues and design positions as before
design_positions = {"C": ["510", "513", "514", "517", "523", "558", "559"]}
ligand_residues = {"A": {"131": "THR", "132": "ASN"}}

# In this case the ouput_name will not be taken as is, but extended with an index (starting from 1) for each 
# grafting result: ./design1.pdb ./design2.pdb ./design3.pdb ...
result = multi_graft_best_pocket_onto_scaffold(pockets=pockets, scaffold_name="./1z5s.pdb", 
                                               ligand_residues=ligand_residues, design_positions=design_positions, 
                                               output_name="./design.pdb", n_solutions=5)
# You can inspect the results:
for res_i, (penalty, graft, design) in enumerate(result):
    print(f"Result {res_i} with penalty of {penalty} includes the following grafted mutations:\n{graft}")
```

# Acknowledgments

This tool is the work of the Protein Design group of Prof. Dr. Birte Höcker at University of Bayreuth.
We thankfully use software packages to enable ATLIGATOR functionality.
This is a non-exhaustive list of requirements:
- biopython
- numpy
- plotly
- matplotlib
- scipy
- pandas
- tqdm

This work was enabled by the European Research Council (H2020-FETOpen-RIA grant 764434 'Pre-ART').
