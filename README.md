# stapler
A Hash Based Method for Matching Crosslinkers into Peptides and Proteins for PyRosetta

### Installation
`pip install pystapler`

### Theory
This protocol is developed to quickly identify pairs of residues that can accomidate disulfides given the 3D structure of a protein. 

To achieve this goal, examples of 30,000 native disulfide structures were procured from the PDB, and the relative positions of the backbone atoms (N, CA, C) were calculated, hashed, and stored into a database.

Therefore, when you have a candidate protein structure and you would like to know if any residue pairs from the structure can accomidate a disulfide, you can quickly calculate all of the relative positions of backbone atoms in the protein and ask which transformations are the same as an example in the database of native examples.

The protocol is scalable/adaptable to other crosslinkers given a set of conformations of the new crosslinker. For disulfides, this was given by the PDB, but for crosslinks that do not have a native examples, the examples must be generated from scratch. Also, this protocol is not limited to sidechain-to-sidechain crosslinking. One could implement loop closure with some tweaks.

![Theory Example Image](/image2.png)

### Example
```python
import pyrosetta
pyrosetta.init()

from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector

from stapler import NativeDisulfideStapler

default_residue_selectors = [TrueResidueSelector(), TrueResidueSelector()]
ab_residue_selectors = [ChainSelector('A'), ChainSelector('B')]
aa_or_ab_residue_selectors = [ChainSelector('A'), ChainSelector('A,B')]

native_disulfide_stapler = NativeDisulfideStapler(
    residue_selectors=default_residue_selectors,
    minimum_sequence_distance=4
)

pose = pyrosetta.pose_from_file('input.pdb')

for i, stapled_pose in enumerate(native_disulfide_stapler.apply(pose)):
    stapled_pose.dump_pdb(f'output_{i}.pdb')
```

![Protein/Disulfide Example Image](/image1.png)
