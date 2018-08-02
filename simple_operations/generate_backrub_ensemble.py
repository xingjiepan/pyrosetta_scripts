#!/usr/bin/env python3
'''Generate a backbrub ensemble for a given protein.
Usage:
    ./generate_backrub_ensemble.py input_pdb
'''

import sys

import pyrosetta
from pyrosetta import rosetta


if __name__ == '__main__':
    pyrosetta.init()
    
    input_pdb = sys.argv[1]

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, input_pdb)


    br_mover = rosetta.protocols.backrub.BackrubMover()
    
    pose.dump_pdb('before_br.pdb')
   
    # Dump 20 structures

    for i in range(20):
        tmp_pose = pose.clone()
        br_mover.apply(tmp_pose)
        tmp_pose.dump_pdb('after_br_{0}.pdb'.format(i))
