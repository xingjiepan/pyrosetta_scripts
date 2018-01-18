#!/usr/bin/env python2

import pyrosetta
from pyrosetta import rosetta


def get_rosetta_loop_from_pdb_indices(pdb_file, start_id, stop_id, output_pdb_file, output_loop_file):
    '''Get the rosetta loop from the pdb indices. The start
    and stop indices are represented as pairs of (chain, atom_id).
    Write the structure and loop definition into new files.
    '''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    start_seqpos = pose.pdb_info().pdb2pose(start_id[0], start_id[1])
    stop_seqpos = pose.pdb_info().pdb2pose(stop_id[0], stop_id[1])

    with open(output_loop_file, 'w') as f:
        f.write('Loop {0} {1} {1} 0 1'.format(start_seqpos, stop_seqpos))

    pose.dump_pdb(output_pdb_file)
