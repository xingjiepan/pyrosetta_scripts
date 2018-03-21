#!/usr/bin/env python2
'''Clean a fuzzball pose and dump a cleaned fuzzball
Usage:
    ./clean_fuzzball.py fuzzball_pdb ligand_params_file ligand_id
'''

import sys

import pyrosetta
from pyrosetta import rosetta

import preprocessing


if __name__ == '__main__':
    fuzzball_pdb = sys.argv[1]
    ligand_params_file = sys.argv[2]
    ligand_id = int(sys.argv[3])

    pyrosetta.init(options='-extra_res_fa {0}'.format(ligand_params_file))

    fuzz_pose = preprocessing.load_cleaned_filtered_fuzz_pose(
            fuzzball_pdb, ligand_id)

    print "Number of motif residues =", fuzz_pose.size() - 1

    fuzz_pose.dump_pdb(fuzzball_pdb[:-4] + '_cleaned.pdb')
