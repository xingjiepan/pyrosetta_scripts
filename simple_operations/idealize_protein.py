#!/usr/bin/env python3
'''Idealize structure.
Usage:
    ./idealize_protein.py pdb_file
'''

import sys
import os

import pyrosetta
from pyrosetta import rosetta


def idealize_protein(pose):
    '''Idealize a pose.'''
    # Idealize the pose

    im = rosetta.protocols.idealize.IdealizeMover()
    im.apply(pose)

    # Renumber the protein

    pose.pdb_info().obsolete()
    

if __name__ == '__main__':
    pdb_file = sys.argv[1]

    pyrosetta.init()

    pose = rosetta.core.import_pose.pose_from_file(pdb_file)

    idealize_protein(pose)

    pose.dump_file(os.path.join(os.path.dirname(pdb_file),
        'idealized_' + os.path.basename(pdb_file)))
