#!/usr/bin/env python2
'''Score a structure.
Usage:
    ./score_structure.py pdb_file.pdb
'''

import sys
import os

import pyrosetta
from pyrosetta import rosetta


if __name__ == '__main__':
    pyrosetta.init()

    # Load pose

    pdb_file = sys.argv[1]
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    # Score

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)

    pose.dump_pdb(os.path.basename(pdb_file) + '.scored.pdb')
