#!/usr/bin/env python2.7

import numpy as np

import pyrosetta
from pyrosetta import rosetta


def set_up_fold_tree_for_root_residue(pose, root_seqpos):
    '''Set up a fold tree for a pose such that one residue
    is the root and all the other residues are its direct children.
    '''
    ft = rosetta.core.kinematics.FoldTree()
   
    other_residues = [i for i in range(1, pose.size() + 1) if i != root_seqpos] 

    for i, res in enumerate(other_residues):
        ft.add_edge(root_seqpos, res, i + 1)

    pose.fold_tree(ft)

def replace_intra_residue_torsions(pose, seqpos, ref_residue):
    '''Replace the intra residue torsions at seqpos
    by the torsions from the reference residue.
    '''
    ref_chis = ref_residue.chi()
    
    for i in range(1, pose.residue(seqpos).nchi() + 1):
        pose.set_chi(i, seqpos, ref_chis[i])

def set_inverse_rotamer(pose, seqpos, ref_residue):
    '''Set inverse rotamer at a position according by
    applying the torsions from the reference residue.
    '''
    # Update the rotamer
    
    last_chi = pose.residue(seqpos).chi_atoms(pose.residue(seqpos).nchi()) 
    ref_atoms = [last_chi[i] for i in range(2, 5)]
    ref_xyzs = [rosetta.numeric.xyzVector_double_t(pose.residue(seqpos).xyz(a)) for a in ref_atoms]
    ref_stub = rosetta.core.kinematics.Stub(ref_xyzs[0], ref_xyzs[1], ref_xyzs[2])

    replace_intra_residue_torsions(pose, seqpos, ref_residue)

    # Find the jump id for this residue

    jump_id = 0

    for i in range(1, pose.fold_tree().num_jump() + 1):
        if pose.fold_tree().downstream_jump_residue(i) == seqpos:
            jump_id = i

    assert(jump_id > 0)

    # Align the end points of the rotamers
    
    current_xyzs = [rosetta.numeric.xyzVector_double_t(pose.residue(seqpos).xyz(a)) for a in ref_atoms]
    current_stub = rosetta.core.kinematics.Stub(current_xyzs[0], current_xyzs[1], current_xyzs[2])

    current_inv_M = current_stub.M.inverse()
    current_inv_v = (current_inv_M * (-1)) * current_stub.v

    transform_M = ref_stub.M * current_inv_M
    transform_v = ref_stub.M * current_inv_v + ref_stub.v
   
    stub_up = pose.conformation().upstream_jump_stub(jump_id)
    stub_down = pose.conformation().downstream_jump_stub(jump_id)

    M = transform_M * stub_down.M
    v = transform_M * stub_down.v + transform_v

    new_stub_down = rosetta.core.kinematics.Stub(M, v)
    new_jump = rosetta.core.kinematics.Jump(stub_up, new_stub_down)

    pose.set_jump(jump_id, new_jump)

def match_anchor_position(target_pose, target_anchor_seqpos, movable_pose, movable_anchor_seqpos):
    '''Match the movable pose to the anchor
    position of the target pose.
    '''
    res_t = target_pose.residue(target_anchor_seqpos)
    res_m = movable_pose.residue(movable_anchor_seqpos)

    current_frame = rosetta.numeric.xyzTransform_double_t(res_m.xyz('N'), res_m.xyz('CA'), res_m.xyz('C'))
    inv_cf = current_frame.inverse()
    target_frame = rosetta.numeric.xyzTransform_double_t(res_t.xyz('N'), res_t.xyz('CA'), res_t.xyz('C'))

    R = target_frame.R * inv_cf.R
    v = target_frame.R * inv_cf.t + target_frame.t

    movable_pose.apply_transform_Rx_plus_v(R, v) 

def find_matched_rotamers_for_anchored_fuzz_ball(target_pose, target_matching_seqposes,
        fuzz_pose, ligand_residue, motif_residues, fuzz_anchor):
    '''Find all rotamers that could match the anchored fuzz ball to the target pose.
    Args:
        target_pose: The target pose.
        target_matching_seqposes: The residues on the target pose that are available for
            motif matching.
        fuzz_pose: The fuzz ball pose.
        ligand_residue: The id of the ligand residue.
        motif_residues: A list of motif residues.
        fuzz_anchor: The residue on the fuzz ball that is anchored to the target

    Return:
        A list of matches defined as (fuzz_ball_anchor_residue, fuzz_ball_matched_residue, 
            target_matched_residue, rotamer_id)
    '''
    pass

def find_matched_rotamers_for_fuzz_ball(target_pose, target_anchor_seqpos, target_matching_seqposes,
        fuzz_pose, ligand_residue, motif_residues):
    '''Find all rotamers that could match the fuzz ball to the target pose.
    Args:
        target_pose: The target pose.
        target_anchor_seqpos: The achor point on the target pose.
        target_matching_seqposes: The residues on the target pose that are available for
            motif matching.
        fuzz_pose: The fuzz ball pose.
        ligand_residue: The id of the ligand residue.
        motif_residues: A list of motif residues.
        
    Return:
        A list of matches defined as (fuzz_ball_anchor_residue, fuzz_ball_matched_residue, 
            target_matched_residue, rotamer_id)
    '''
   
    # Set up the fold tree that roots on the ligand residue

    set_up_fold_tree_for_root_residue(fuzz_pose, ligand_residue)

    # Find the matches for each of the motif residues achored

    for fuzz_anchor in motif_residues:

        rotamer_set = rosetta.core.pack.rotamer_set.bb_independent_rotamers( fuzz_pose.residue(fuzz_anchor).type(), True )

        for i in range(1, len(rotamer_set) + 1):
            set_inverse_rotamer(fuzz_pose, fuzz_anchor, rotamer_set[i])
            match_anchor_position(target_pose, target_anchor_seqpos, fuzz_pose, fuzz_anchor)
            
            find_matched_rotamers_for_anchored_fuzz_ball(target_pose, target_matching_seqposes,
                fuzz_pose, ligand_residue, motif_residues, fuzz_anchor)
            
            fuzz_pose.dump_pdb('debug/test_fuzz_{0}_{1}.pdb'.format(fuzz_anchor, i)) ###DEBUG
            
    target_pose.dump_pdb('debug/target.pdb') ###DEBUG


if __name__ == '__main__':
    pyrosetta.init(options='-extra_res_fa inputs/REN_no_charge_from_mol2.params')
    
    fuzz_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(fuzz_pose, 'inputs/binding_site_from_james_renumbered.pdb')

    target_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(target_pose, 'inputs/2lvb_no_terms.pdb')

    #target_residues = list(range(3, 9)) + list(range(37, 50)) + list(range(63, 77))

    find_matched_rotamers_for_fuzz_ball(target_pose, 6, list(range(37, 50)) + list(range(63, 77)),
            fuzz_pose, 1, list(range(2, fuzz_pose.size() + 1)))

