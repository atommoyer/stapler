import glob

import tqdm

import numpy as np
import getpy as gp

import nerf

import homog
import xbin

from rdkit import Chem

chi_resolution = 30.0
assert 360.0 / chi_resolution < 256.0

hash_function_kwargs = {'cart_resl' : 1.0, 'ori_resl' : 15.0, 'cart_bound' : 512.0}
hash_function = xbin.XformBinner(**hash_function_kwargs)

hash_table_kwargs = {'key_type' : np.dtype('i8'), 'value_type': np.dtype('i4')}
hash_table = gp.MultiDict(**hash_table_kwargs)

for pdb in tqdm.tqdm(glob.glob('/home/apmoyer/stapler/data/native_disulfide/*.pdb')):
    mol = Chem.MolFromPDBFile(pdb)

    xyz = np.array([
        mol.GetConformer(0).GetAtomPosition(0), # N1
        mol.GetConformer(0).GetAtomPosition(1), # CA1
        mol.GetConformer(0).GetAtomPosition(2), # C1

        mol.GetConformer(0).GetAtomPosition(4), # CB1
        mol.GetConformer(0).GetAtomPosition(5), # SG1
        mol.GetConformer(0).GetAtomPosition(11), # SG2
        mol.GetConformer(0).GetAtomPosition(10), # CB2

        mol.GetConformer(0).GetAtomPosition(7), # CA2
        mol.GetConformer(0).GetAtomPosition(6), # N2
        mol.GetConformer(0).GetAtomPosition(8), # C2
    ])

    dep = np.array([
        [0, 0, 0], #
        [0, 0, 0], #
        [1, 0, 0], #
        [1, 0, 2], #
        [3, 1, 0], # <- Perturb SG1-CB1-CA1-N1
        [4, 3, 1], # <- Perturb SG2-SG1-CB1-CA1
        [5, 4, 3], # <- Perturb CB2-SG2-SG1-CB1
        [6, 5, 4], # <- Perturb CA2-CB2-SG2-SG1
        [7, 6, 5], # <- Perturb N2-CA2-CB2-SG2
        [7, 6, 8], #
    ])

    dof = nerf.iNeRF(xyz, dependency=dep)

    repeats = 100
    dofs = np.repeat(dof[np.newaxis], repeats, axis=0)
    for i in [4, 5, 6, 7, 8]: # [CA1-CB1, CB1-SG1, SG1-SG2, SG2-CB2, CB2-CA2]
        dofs[:,i,2] += np.random.uniform(-5, 5, repeats)
        dofs[:,i,2] = (dofs[:,i,2] + (180.0 + chi_resolution/2.0)) % 360.0 - (180.0 + chi_resolution/2.0)

    xyzs = nerf.NeRF(dofs, dependency=dep)

    stubs_a = homog.hstub(xyzs[:,0,:], xyzs[:,1,:], xyzs[:,2,:]) # N1 CA1 C1
    stubs_b = homog.hstub(xyzs[:,-2,:], xyzs[:,-3,:], xyzs[:,-1,:]) # N2 CA2 C2

    xforms_ab = np.linalg.inv(stubs_a) @ stubs_b
    xforms_ba = np.linalg.inv(stubs_b) @ stubs_a

    keys_ab = hash_function.get_bin_index(xforms_ab)
    keys_ba = hash_function.get_bin_index(xforms_ba)

    chis_ab = np.around(dofs[:,[4,8],2]/chi_resolution).astype('i2').flatten().view(np.dtype('i4'))
    chis_ba = np.around(dofs[:,[8,4],2]/chi_resolution).astype('i2').flatten().view(np.dtype('i4'))

    hash_table[keys_ab] = chis_ab
    hash_table[keys_ba] = chis_ba

    tqdm.tqdm.write(str(len(hash_table)))

hash_table.dump('../stapler/hash_tables/default_native_disulfide.bin')
