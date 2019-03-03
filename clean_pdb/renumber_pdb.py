#!/usr/bin/env python2
'''Renumber a pdb file. Run as:
    ./renumber_pdb.py input.pdb output.pdb [params_file]
'''

import sys

import pyrosetta
from pyrosetta import rosetta


if __name__ == '__main__':
    ipdb = sys.argv[1]
    opdb = sys.argv[2]

    if len(sys.argv) > 3:
        pyrosetta.init(options='-extra_res_fa {0} -ignore_unrecognized_res true'.format(sys.argv[3]))
    else:
        pyrosetta.init(options='-ignore_unrecognized_res true')

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, ipdb)

    pose.pdb_info().obsolete(True)

    pose.dump_pdb(opdb)
