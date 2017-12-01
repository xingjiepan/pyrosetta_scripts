#!/usr/bin/env python2.7

import pyrosetta
from pyrosetta import rosetta

import lh_basic


if __name__ == '__main__':
    pyrosetta.init()

    loop_size = 14
    lh_db_path = "lh_db_helix"
    #lh_db_path = "/home/xingjie/DataBases/loophash_db"
    
    lh_library = lh_basic.load_lh_library([loop_size], lh_db_path) 

    pose = rosetta.core.pose.Pose()
    rosetta.core.pose.make_pose_from_sequence(pose, 'A' * loop_size, 'fa_standard')

    hashmap = lh_library.gethash(loop_size)

    #for i in range(hashmap.n_loops()):
    for i in range(min(100, hashmap.n_loops())):
        bb_seg, sequence = lh_basic.extract_lh_fragment(lh_library, loop_size, i)
        
        lh_basic.apply_bb_seg_to_pose(pose, bb_seg, 1)

        pose.dump_pdb('outputs/out_{0}.pdb'.format(i))

