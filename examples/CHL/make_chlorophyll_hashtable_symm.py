# Adam Moyer and Nate Novy 7/27/21
# to run: /home/apmoyer/miniconda3/envs/default/bin/python make_chlorophyll_hashtable_symm.py

import os
import sys
import glob

import tqdm

import numpy as np
import getpy as gp

import nerf
import align

import homog
import xbin

from rdkit import Chem

import numpy as np


np.set_printoptions(suppress=True, threshold=sys.maxsize, linewidth=250)

# chi_resolution condenses the dataset later by preventing the addition of conformers that are almost identical
# A high chi_resolution will cause more conformers to be collapsed to the same angles of rotation
# This was originally 30, but we changed it to 10.
chi_resolution = 10.0
assert 360.0 / chi_resolution < 256.0

##### This section was not modified for the purposes of building a chlorophyll hash table
hash_function_kwargs = {'cart_resl' : 0.5, 'ori_resl' : 10.0, 'cart_bound' : 512.0}
hash_function = xbin.XformBinner(**hash_function_kwargs)

hash_table_kwargs = {'key_type' : np.dtype('i8'), 'value_type': np.dtype('i8')}
hash_table = gp.MultiDict(**hash_table_kwargs)
#####

# This loads in the pdb and converts it to a format where each atom is defined partially by its index
# chlorophyll_reordered_symm.pdb is one chlorophyll of the QA0 dimer with three floating atoms of the second QA0
# which will be used later to generate and align a symmetry mate for the first chlorophyll.
# These three floating atoms will be removed when the chloropylls are rebuilt into the docked proteins.
mol = Chem.MolFromPDBFile('chlorophyll_reordered_symm.pdb')
conf = mol.GetConformer(0)

xyz = conf.GetPositions()

# adj_matrix is an NxN matrix where N is the number of atoms in the pose. The Matrix
# describes which atoms are connected by placing a 1 to indicate a bond and a 0 to indicate no bond.
# This section shouldnt need to be modified if you are using chlorophyll_reordered_symm.pdb
adj_matrix = Chem.rdmolops.GetAdjacencyMatrix(mol)

adj_matrix[9,10] = 1
adj_matrix[10,9] = 1

adj_matrix[11,10] = 1
adj_matrix[10,11] = 1

adj_matrix[54,55] = 1
adj_matrix[55,54] = 1

adj_matrix[55,56] = 1
adj_matrix[56,55] = 1

adj_matrix[55,57] = 1
adj_matrix[57,55] = 1


# The depenency table (dep) is a tree descibing which atoms are bonded together based
# on the indices. dep[5] will equal [4,1,0], meaning atom 5 is bonded to 4 which bonds to 1, which bonds to 0.
# You can edit dep if you want to adjust an angle that is not desribed by the table. If you wanted to adjust
# the dihedral at atom index 5, then you would be adjusting the dihedral made from atoms 5, 4, 1, and 0.
dep = nerf.build_deps(adj_matrix)
## print(dep)

# dof is another way to describe the posisstions of atoms using their [distances, angle,  dihedral], stored  in that format.
dof = nerf.iNeRF(xyz, deps=dep)

# Here, we are setting the angle made by atoms indexed as 10, 9, 8 to 125 degrees.
dof[10,1] = 125.5
dof[11,1] = 100.0075
# Here we are setting the dihedral made by atoms 12, 11, 10, 9 to 85 degrees.
dof[12,2] = 85.0

## print(dof)

# repeats is the maximum number of entries you will find in one run. The final number of entries is
# less than repeats because duplicates are not stored in the table
repeats = 250000
# This expands dof to dofs, where each value in dofs is a dof for a different hash_table entry
dofs = np.repeat(dof[np.newaxis], repeats, axis=0)

# The following three lines cause chis 1 and 2 of histidine to change (5,6) and Cause
# the histidine to spin about the NE2-Mg bond (11).
# Each time you run this script, only one of the  angles in random.choice() will be selected,
# so I recommended doing many smaller batches instead of one large batch so that all rotomers will
# be sampled.
dofs[:,11,2] = np.random.uniform(-180, 180, repeats)
dofs[:,5,2] = np.random.normal(0, 10, repeats) + np.random.choice([60,180,-60], size=repeats)
dofs[:,6,2] = np.random.normal(0, 15, repeats) + np.random.choice([100,-100,0], size=repeats)


# The following three lines cause the Chl-His angle to change. You can comment these out if you
# wish to generate a hash table where all entries are set at 90 degrees, but this significantly reduces
# the size of the hash table.
#dofs[:,10,1] += (np.random.normal(0, 8, repeats))
#dofs[:,11,1] += (np.random.normal(0, 8, repeats))
#dofs[:,10,2] += (np.random.normal(0, 8, repeats))

# Gets the xyz values of the atoms in each entry
xyzs1 = nerf.NeRF(dofs, deps=dep)
# Makes a copy of xyzs1 and aligns it onto the three floating atoms. This is how the symmetry
# mate is generated so that the histidines have the same conformation.
xyzs2 = align.align(xyzs1, xyzs1, mapping=[[55,56,57,10,11,15],[10,11,15,55,56,57]])

#### You can use this for debuging or making figures for outputing several states
# of monoA or monoB to a pdb

for k in range(repeats):
   for i, j in enumerate(xyzs1[k]):
       conf.SetAtomPosition(i, j)
   with open('monoA.xyz', 'a') as text_file:
       text_file.write(Chem.MolToXYZBlock(mol))

   for i, j in enumerate(xyzs2[k]):
       conf.SetAtomPosition(i, j)
   with open('monoB.xyz', 'a') as text_file:
       text_file.write(Chem.MolToXYZBlock(mol))

####

##### This section was not modified for the chlorophyll hash table
# this takes the xyz positions of just the backbone atoms of the two histidinees
stubs_a = homog.hstub(xyzs1[:,0,:], xyzs1[:,1,:], xyzs1[:,2,:]) # N1 CA1 C1
stubs_b = homog.hstub(xyzs2[:,0,:], xyzs2[:,1,:], xyzs2[:,2,:]) # N2 CA2 C2

# Calculate transformation between backbone atoms of histidines
xforms_ab = np.linalg.inv(stubs_a) @ stubs_b
xforms_ba = np.linalg.inv(stubs_b) @ stubs_a

# The key for each entry in the hash table is calculated based on the transformation
keys_ab = hash_function.get_bin_index(xforms_ab)
keys_ba = hash_function.get_bin_index(xforms_ba)

# This likely wont need to be modified, but it is lowering the resolution of chi 1 and 2
# to reduce the number of similar entries in the hash table.
# This is the value component of the Key-Value dictionary pair
chis = np.concatenate([dofs[:,[5,6,11],2], np.zeros((repeats, 1))], axis=-1)
encoded_chis = np.around(chis/chi_resolution).astype('i2').flatten().view(np.dtype('i8'))

# Load the hash table from bin file
filename = 'chl_hashtable_test.bin'

if os.path.exists(filename):
    hash_table.load(filename)

# Add the two histidine conformations as entries to the hash table.
hash_table[keys_ab] = encoded_chis
hash_table[keys_ba] = encoded_chis

# Dump hash table into bin file
hash_table.dump(filename)

# Print the number of entries in the hash table
print(len(hash_table))
