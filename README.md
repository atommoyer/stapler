### stapler
A Motif Hash Based Method for Matching Crosslinkers into Peptides and Proteins

# Installation
`pip install pystapler`

# Example
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
    residue_selectors=ab_residue_selectors,
    minimum_sequence_distance=4
)

pose = pyrosetta.pose_from_file('input.pdb')

for i, stapled_pose in enumerate(native_disulfide_stapler.apply(pose)):
    stapled_pose.dump_pdb(f'output_{i}.pdb')
```
