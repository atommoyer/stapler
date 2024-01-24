import pkg_resources

import numpy as np

import pyrosetta
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

from .stapler import Stapler

class ChlorophyllStapler(Stapler):
    def __init__(
        self,
        residue_selectors=(TrueResidueSelector(), TrueResidueSelector()),
        atom_selectors=(('N', 'CA', 'C'), ('N', 'CA', 'C')),
        minimum_sequence_distance=1,
        maximum_neighborhood_distance=30.0,
        hash_function_kwargs={'cart_resl' : 0.5, 'ori_resl' : 10.0, 'cart_bound' : 512.0},
        hash_table_kwargs={'key_type' : np.dtype('i8'), 'value_type': np.dtype('i8'), 'filename': pkg_resources.resource_filename('stapler', 'hash_tables/chl_hashtable.bin')}
    ):

        super(ChlorophyllStapler, self).__init__(
            residue_selectors=residue_selectors,
            atom_selectors=atom_selectors,
            minimum_sequence_distance=minimum_sequence_distance,
            maximum_neighborhood_distance=maximum_neighborhood_distance,
            hash_function_kwargs=hash_function_kwargs,
            hash_table_kwargs=hash_table_kwargs
        )


    def decode(self, i, j, data):
        chis = tuple(np.array([data]).view(np.dtype('i2'))*10)

        return (
            (i, 'CHL', chis),
            (j, 'CHL', chis)
        )


    def place(self, pose, staple):
        MutateResidue(staple[0][0], staple[0][1]).apply(pose)

        for i, chi in enumerate(staple[0][2][:3], start=1):
            print(chi)
            pose.set_chi(i, staple[0][0], chi)

        MutateResidue(staple[1][0], staple[1][1]).apply(pose)

        for i, chi in enumerate(staple[1][2][:3], start=1):
            pose.set_chi(i, staple[1][0], chi)

        return pose


    def select(self, pose):
        assert len(np.argwhere(self.residue_selectors[0].apply(pose))) == len(np.argwhere(self.residue_selectors[1].apply(pose)))

        sele = np.array([
            [i, j] for i, j in zip(
                np.squeeze(np.argwhere(self.residue_selectors[0].apply(pose))),
                np.squeeze(np.argwhere(self.residue_selectors[1].apply(pose)))
            )
        ])

        return sele
