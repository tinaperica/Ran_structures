import numpy as np
import mdtraj as md

trajectory = md.load('~/Documents/Gsp1_bioinformatics/Ran_structures/pdbs/clean/1i2m.pdb')
sasa = md.shrake_rupley(trajectory)

print(trajectory)
print('sasa data shape', sasa.shape)
total_sasa = sasa.sum(axis = 1)
print(total_sasa.shape)
