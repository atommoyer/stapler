import glob
import os

import pyrosetta
pyrosetta.init('--beta --extra_res_fa CHL.params')

from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.select.residue_selector import SecondaryStructureSelector

from pyrosetta.rosetta.protocols.symmetry import DetectSymmetry
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover

from stapler import ChlorophyllStapler

# Preset ResidueSelectors
default_residue_selectors = [TrueResidueSelector(), TrueResidueSelector()]
interface_residue_selectors = [ChainSelector('A'), ChainSelector('B')]
interface_or_internal_residue_selectors = [ChainSelector('A'), ChainSelector('A,B')]
only_binder_residue_selectors = [ChainSelector('B'), ChainSelector('B')]
not_on_loops = [SecondaryStructureSelector('HE'), SecondaryStructureSelector('HE')]
not_on_loops_across_interface = [AndResidueSelector(SecondaryStructureSelector('HE'),ChainSelector('A')),
                                 AndResidueSelector(SecondaryStructureSelector('HE'),ChainSelector('B'))]

# Initialize the native disulfide stapler with defaults.
chlorophyll_stapler = ChlorophyllStapler(
    residue_selectors=interface_residue_selectors
)

for pdb in glob.glob('_chl_inputs/*.pdb'):
    pdb_filename = os.path.splitext(os.path.basename(pdb))[0]
    pose = pyrosetta.pose_from_file(pdb)

    for i, crosslinked_pose in enumerate(chlorophyll_stapler.apply(pose)):
        crosslinked_pose.dump_pdb(f'_chl_outputs/{pdb_filename}_staple_{i}.pdb')
