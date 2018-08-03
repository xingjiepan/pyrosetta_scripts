#!/usr/bin/env python3
'''Get the amino acid composition of a set of structures
Usage:
    ./get_amino_acid_composition.py pdb_file1 [pdb_file2 ...]
'''

import sys

import pyrosetta
from pyrosetta import rosetta


def get_amino_acid_composition_for_pdb_file(pdb_file):
    '''Get a dictionary of amino acid counts.'''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    aa_d = {}

    for i in range(1, pose.size() + 1):
        name = pose.residue(i).name3()

        if name in aa_d.keys():
            aa_d[name] += 1
        else:
            aa_d[name] = 1

    return aa_d


if __name__ == '__main__':
    pyrosetta.init()

    pdb_files = sys.argv[1:]

    # Get the dictionary for amino acid composition

    aa_d = {}

    for pdb_file in pdb_files:
        if not (pdb_file.endswith('.pdb') or pdb_file.endswith('.pdb.gz')):
            continue
        
        aa_d_single = get_amino_acid_composition_for_pdb_file(pdb_file)

        for k in aa_d_single.keys():
            if k in aa_d.keys():
                aa_d[k] += aa_d_single[k]
            else:
                aa_d[k] = aa_d_single[k]

    # Print the result 

    total_aas = sum(aa_d[k] for k in aa_d.keys())

    for k in sorted(aa_d.keys()):
        print('{0}\t{1}\t{2:.2f}%'.format(k, aa_d[k], 100 * aa_d[k] / total_aas))
