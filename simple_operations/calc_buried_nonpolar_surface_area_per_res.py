#!/usr/bin/env python3
'''Calculate buried nonpolar surface area normalized by
the number of residues.
Usage:
    ./calc_buried_nonpolar_surface_area_per_res.py pdb_file
'''

import sys

import pyrosetta
from pyrosetta import rosetta

def calc_buried_np_SASA(pose, residue_types=None):
    '''Calculate the buried nonpolar surface area in the designed
    structure on for given residue types'''
    # Calculate the hydrophobic SASA
    
    rsd_sasa = pyrosetta.rosetta.utility.vector1_double()
    rsd_hydrophobic_sasa = pyrosetta.rosetta.utility.vector1_double()
    rosetta.core.scoring.calc_per_res_hydrophobic_sasa(pose, rsd_sasa, rsd_hydrophobic_sasa, 1.4) #The last arguement is the probe radius

    # Calculate the hydrophobic buried SASA

    for i in range(1, pose.size() + 1):
        total_sasa = rosetta.core.scoring.normalizing_area_total_hydrophobic_atoms_only(pose.residue(i).name1())
        rsd_hydrophobic_sasa[i] = total_sasa - rsd_hydrophobic_sasa[i]

    # Sum the SASA of given residues

    if residue_types:
        residues = [i for i in range(1, pose.size() + 1) if pose.residue(i).name1() in residue_types]
        return sum(rsd_hydrophobic_sasa[i] for i in residues)
    
    return sum(rsd_hydrophobic_sasa)

def calc_buried_nonpolar_surface_area_per_res(pose):
    '''Calculate the buried nonpolar surface area normalized by the number of residues.'''
    return calc_buried_np_SASA(pose) / pose.size()

if __name__ == '__main__':
    pyrosetta.init()

    pdb_file = sys.argv[1]

    pose = rosetta.core.import_pose.pose_from_file(pdb_file)

    print('The buried nonpolar surface area per residue is {0}.'.format(calc_buried_nonpolar_surface_area_per_res(pose)))

