#!/usr/bin/env python2.7

import pyrosetta
from pyrosetta import rosetta


def remove_terminal_variants(pose):
    '''Remove terminal residue variants of a pose.'''
    for i in range(1, pose.size() + 1):
        if rosetta.core.pose.is_upper_terminus(pose, i):
            rosetta.core.pose.remove_upper_terminus_type_from_pose_residue(pose, i)
        if rosetta.core.pose.is_lower_terminus(pose, i):
            rosetta.core.pose.remove_lower_terminus_type_from_pose_residue(pose, i)


if __name__ == '__main__':
    pyrosetta.init()

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, 'inputs/ke07_active_site_no_substrate.pdb')
    remove_terminal_variants(pose)
    #rosetta.core.pose.correctly_add_cutpoint_variants(pose)
    
    rotamers = rosetta.core.pack.rotamer_set.bb_independent_rotamers( pose.residue(1).type(), True ) 

    print len(rotamers)

    for i in range(1, len(rotamers) + 1):
        pose.replace_residue(1, rotamers[i], False) 
        #pose.replace_residue(1, rotamers[i], True) 

        pose.dump_pdb('test_{0}.pdb'.format(i))
