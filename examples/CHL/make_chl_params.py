import os
import sys
import glob
from collections import defaultdict

import tqdm

import numpy as np
import getpy as gp

import nerf
import align

import homog
import xbin

from rdkit import Chem

import numpy as np
import random

np.set_printoptions(suppress=True, threshold=sys.maxsize, linewidth=250, precision=4)

mol = Chem.MolFromPDBFile('chlorophyll_reordered_symm_H.pdb', removeHs=False)
conf = mol.GetConformer(0)
xyz = conf.GetPositions()

adj_matrix = Chem.rdmolops.GetAdjacencyMatrix(mol)

adj_matrix[9,10] = 1
adj_matrix[10,9] = 1

adj_matrix[11,10] = 1
adj_matrix[10,11] = 1

print(adj_matrix)


deps = nerf.build_deps(adj_matrix)
print(deps)

dofs = nerf.iNeRF(xyz, deps=deps)
print(dofs)

counter = defaultdict(int)
atomic_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

names = []
for i, atomic_symbol in enumerate(atomic_symbols):
    if i == 10:
        names.append('Mg1')
    elif i == 59:
        names.append('HN1')
    elif atomic_symbol == 'H':
        names.append(atomic_symbol + 'C' + str(1 + counter[atomic_symbol]))
        counter[atomic_symbol] += 1
    else:
        names.append(atomic_symbol + str(1 + counter[atomic_symbol]))
        counter[atomic_symbol] += 1
names = np.array(names)
print(names)

bonds = np.argwhere(adj_matrix == 1)
bonds_names = np.empty_like(bonds, dtype=np.dtype(('U', 4)))
for i, name in enumerate(names):
    bonds_names[np.where(bonds == i)] = name
print(bonds_names)

for bonds_name in bonds_names:
    print(f'BOND {bonds_name[0]} {bonds_name[1]}')

icoors = []
for i, dep in enumerate(deps):
    icoors.append([i, dep[0], dep[1], dep[2]])
icoors = np.array(icoors)

icoors_names = np.empty_like(icoors, dtype=np.dtype(('U', 4)))
for i, name in enumerate(names):
    icoors_names[np.where(icoors == i)] = name
print(icoors_names)

for icoors_name, dof in zip(icoors_names, dofs):
    print(f'ICOOR_INTERNAL    {icoors_name[0]}  {dof[2]:.6f}     {180.0 - dof[1]:.6f}     {dof[0]:.6f}     {icoors_name[1]}   {icoors_name[2]}   {icoors_name[3]}')
