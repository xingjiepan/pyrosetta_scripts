#!/usr/bin/env python3
'''Symmetrize a monomer structure.
Usage:
    ./symmetrize_monomer.py pdb_file repeat_unit_start_point repeat_unit_stop_point
'''

import sys
import os

import pyrosetta
from pyrosetta import rosetta


def symmetrize_monomer(pose, repeat_unit_start_point, repeat_unit_stop_point):
    '''Symmetrize a pose.'''
    # Idealize the pose

    for i in range(1, pose.size() + 1):
        rosetta.core.conformation.idealize_position(i, pose.conformation())

    # Symmetrize the monomer

    phis = []
    psis = []
    omegas = []

    for i in range(repeat_unit_start_point, repeat_unit_stop_point + 1):
        phis.append(pose.phi(i))
        psis.append(pose.psi(i))
        omegas.append(pose.omega(i))

    repeat_length = repeat_unit_stop_point - repeat_unit_start_point + 1

    for seqpos in range(1, pose.size()):
        i = (seqpos - repeat_unit_start_point) % repeat_length

        pose.set_phi(seqpos, phis[i])
        pose.set_psi(seqpos, psis[i])
        pose.set_omega(seqpos, omegas[i])
    

if __name__ == '__main__':
    pdb_file = sys.argv[1]
    repeat_unit_start_point = int(sys.argv[2])
    repeat_unit_stop_point = int(sys.argv[3])

    pyrosetta.init()

    pose = rosetta.core.import_pose.pose_from_file(pdb_file)

    symmetrize_monomer(pose, repeat_unit_start_point, repeat_unit_stop_point)

    pose.dump_file(os.path.join(os.path.dirname(pdb_file),
        'symmetrized_' + os.path.basename(pdb_file)))
