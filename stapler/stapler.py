import itertools
import pkg_resources

import numpy as np
import getpy as gp

import xbin

from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector

from .homog import hstub

class Stapler(object):
    def __init__(
        self,
        residue_selectors=(TrueResidueSelector(), TrueResidueSelector()),
        atom_selectors=(('N', 'CA', 'C'), ('N', 'CA', 'C')),
        minimum_sequence_distance=1,
        maximum_neighborhood_distance=20.0,
        hash_function_kwargs={'cart_resl' : 1.0, 'ori_resl' : 15.0, 'cart_bound' : 512.0},
        hash_table_kwargs={'key_type' : np.dtype('i8'), 'value_type': np.dtype('i8'), 'filename': None}
    ):
        print(hash_table_kwargs['filename'])
        self.hash_function = xbin.XformBinner(**hash_function_kwargs)
        self.hash_table = gp.MultiDict(**hash_table_kwargs)

        self.residue_selectors = residue_selectors
        self.atom_selectors = atom_selectors

        self.minimum_sequence_distance = minimum_sequence_distance
        self.maximum_neighborhood_distance = maximum_neighborhood_distance


    def decode(self, i, j, data):
        raise NotImplementedError


    def place(self, pose, staple):
        raise NotImplementedError


    def apply(self, pose):
        xyzs_a = np.array([[residue.atom(atom).xyz() for atom in self.atom_selectors[0]] for residue in pose.residues if not residue.is_virtual_residue()])
        xyzs_b = np.array([[residue.atom(atom).xyz() for atom in self.atom_selectors[0]] for residue in pose.residues if not residue.is_virtual_residue()])

        sele = np.array([
            [i, j] for i, j in itertools.product(np.argwhere(self.residue_selectors[0].apply(pose)), np.argwhere(self.residue_selectors[1].apply(pose)))
            if np.abs(i - j) >= self.minimum_sequence_distance and np.linalg.norm(xyzs_a[i,1,:] - xyzs_b[j,1,:]) <= self.maximum_neighborhood_distance
        ])

        sele_a = np.squeeze(sele[:,0])
        sele_b = np.squeeze(sele[:,1])

        stubs_a = hstub(xyzs_a[sele_a,0,:], xyzs_a[sele_a,1,:], xyzs_a[sele_a,2,:])
        stubs_b = hstub(xyzs_b[sele_b,0,:], xyzs_b[sele_b,1,:], xyzs_b[sele_b,2,:])

        xforms_ab = np.linalg.inv(stubs_a) @ stubs_b
        xforms_ba = np.linalg.inv(stubs_b) @ stubs_a

        keys_ab = self.hash_function.get_bin_index(xforms_ab)
        keys_ba = self.hash_function.get_bin_index(xforms_ba)

        staples = set(
            [frozenset(self.decode(int(i+1), int(j+1), staple_ab)) for (i, j), staples_ab in zip(sele, self.hash_table[keys_ab]) for staple_ab in staples_ab] +
            [frozenset(self.decode(int(i+1), int(j+1), staple_ba)) for (j, i), staples_ba in zip(sele, self.hash_table[keys_ba]) for staple_ba in staples_ba]
        )

        for staple in staples:
            print(tuple(staple))
            yield self.place(pose.clone(), tuple(staple))
