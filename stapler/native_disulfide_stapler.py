import pkg_resources

import numpy as np

import pyrosetta
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.conformation import form_disulfide

from .stapler import Stapler

class NativeDisulfideStapler(Stapler):
    def __init__(
        self,
        residue_selectors=(TrueResidueSelector(), TrueResidueSelector()),
        atom_selectors=(('N', 'CA', 'C'), ('N', 'CA', 'C')),
        minimum_sequence_distance=1,
        maximum_neighborhood_distance=20.0,
        hash_function_kwargs={'cart_resl' : 1.0, 'ori_resl' : 15.0, 'cart_bound' : 512.0},
        hash_table_kwargs={'key_type' : np.dtype('i8'), 'value_type': np.dtype('i4'), 'filename': pkg_resources.resource_filename('stapler', 'hash_tables/default_native_disulfide.bin')}
    ):

        super(NativeDisulfideStapler, self).__init__(
            residue_selectors=residue_selectors,
            atom_selectors=atom_selectors,
            minimum_sequence_distance=minimum_sequence_distance,
            maximum_neighborhood_distance=maximum_neighborhood_distance,
            hash_function_kwargs=hash_function_kwargs,
            hash_table_kwargs=hash_table_kwargs
        )


    def decode(self, i, j, data):
        return (
            (i, 'CYS', (np.array([data]).view(np.dtype('i2'))[0]*30.0,)),
            (j, 'CYS', (np.array([data]).view(np.dtype('i2'))[1]*30.0,))
        )


    def place(self, pose, staple):
        form_disulfide(pose.conformation(), staple[0][0], staple[1][0])

        for i, chi in enumerate(staple[0][2]):
            pose.set_chi(i+1, staple[0][0], chi)

        for i, chi in enumerate(staple[1][2]):
            pose.set_chi(i+1, staple[1][0], chi)

        return pose
