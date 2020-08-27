import glob

import tqdm

import numpy as np
import getpy as gp

import nerf

import homog
import xbin

from rdkit import Chem

chi_resolution = 10.0
assert 360.0 / chi_resolution < 256.0

hash_function_kwargs = {'cart_resl' : 1.0, 'ori_resl' : 15.0, 'cart_bound' : 512.0}
hash_function = xbin.XformBinner(**hash_function_kwargs)

hash_table_kwargs = {'key_type' : np.dtype('i8'), 'value_type': np.dtype('i4')}
hash_table = gp.MultiDict(**hash_table_kwargs)

for pdb in tqdm.tqdm(glob.glob('/home/apmoyer/stapler/data/native_disulfide/*.pdb')):
    mol = Chem.MolFromPDBFile(pdb)

    abc = np.array([[
        mol.GetConformer(0).GetAtomPosition(2), # C1
        mol.GetConformer(0).GetAtomPosition(0), # N1
        mol.GetConformer(0).GetAtomPosition(1), # CA1
    ]])

    xyz = np.array([[
        mol.GetConformer(0).GetAtomPosition(4), # CB1
        mol.GetConformer(0).GetAtomPosition(5), # SG1
        mol.GetConformer(0).GetAtomPosition(11), # SG2
        mol.GetConformer(0).GetAtomPosition(10), # CB2
        mol.GetConformer(0).GetAtomPosition(7), # CA2
        mol.GetConformer(0).GetAtomPosition(6), # N2
        mol.GetConformer(0).GetAtomPosition(8), # C2
    ]])

    dof = nerf.iNeRF(abc, xyz)

    abcs = np.repeat(abc, 1000, axis=0)
    dofs = np.repeat(dof, 1000, axis=0)

    for i in [1, 2, 3, 4, 5]: # [CA1-CB1, CB1-SG1, SG1-SG2, SG2-CB2, CB2-CA2]
        dofs[:,i,2] += np.random.uniform(-5, 5, dofs.shape[0])
        dofs[:,i,2] = (dofs[:,i,2] + (180.0 + chi_resolution/2.0)) % 360.0 - (180.0 + chi_resolution/2.0)

    xyzs = nerf.NeRF(abcs, dofs)

    stubs_a = homog.hstub(abcs[:,1,:], abcs[:,2,:], abcs[:,0,:]) # N1 CA1 C1
    stubs_b = homog.hstub(xyzs[:,-2,:], xyzs[:,-3,:], xyzs[:,-1,:]) # N2 CA2 C2

    xforms_ab = homog.hinv(stubs_a) @ stubs_b
    xforms_ba = homog.hinv(stubs_b) @ stubs_a

    keys_ab = hash_function.get_bin_index(xforms_ab)
    keys_ba = hash_function.get_bin_index(xforms_ba)

    chis_ab = np.around(dofs[:,[1,5],2]/chi_resolution).astype('i2').flatten().view(np.dtype('i4'))
    chis_ba = np.around(dofs[:,[5,1],2]/chi_resolution).astype('i2').flatten().view(np.dtype('i4'))

    hash_table[keys_ab] = chis_ab
    hash_table[keys_ba] = chis_ba

hash_table.dump('default_native_disulfide.bin')
