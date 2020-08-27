import glob
import os
import shutil

def filename(path):
    return

import pyrosetta
pyrosetta.init()

from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.select.residue_selector import SecondaryStructureSelector

from pyrosetta.rosetta.protocols.symmetry import DetectSymmetry
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover

from stapler import NativeDisulfideStapler

# Preset ResidueSelectors
default_residue_selectors = [TrueResidueSelector(), TrueResidueSelector()]
interface_residue_selectors = [ChainSelector('A'), ChainSelector('B')]
interface_or_internal_residue_selectors = [ChainSelector('A'), ChainSelector('A,B')]
only_binder_residue_selectors = [ChainSelector('B'), ChainSelector('B')]
not_on_loops = [SecondaryStructureSelector('HE'), SecondaryStructureSelector('HE')]
not_on_loops_across_interface = [AndResidueSelector(SecondaryStructureSelector('HE'),ChainSelector('A')),
                                 AndResidueSelector(SecondaryStructureSelector('HE'),ChainSelector('B'))]

# Initialize the native disulfide stapler with defaults.
native_disulfide_stapler = NativeDisulfideStapler(
    residue_selectors=default_residue_selectors,
    minimum_sequence_distance=4
)

for pdb in glob.glob('_inputs/*.pdb'):
    pdb_filename = os.path.splitext(os.path.basename(pdb))[0]

    pose = pyrosetta.pose_from_file(pdb)

    # Preset Movers for Symmetry
    # DetectSymmetry().apply(pose)
    # SetupForSymmetryMover('C2.sym').apply(pose)

    for i, crosslinked_pose in enumerate(native_disulfide_stapler.apply(pose)):
        crosslinked_pose.dump_pdb(f'_outputs/{pdb_filename}_staple_{i}.pdb')
