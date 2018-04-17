#!/usr/bin/env python2.7

import os

import numpy as np

import pyrosetta
from pyrosetta import rosetta

from fuzzball_matching_basic import *


def residue_contact_degree(target_residue, pose_to_compare_with, cutoff_distance=10):
    '''Calculate the contact degrees of a target residue against
    residues of a pose to compare with. For one residue, the contact degree
    is the number of residues to which the NBR-NBR distance is within a 
    given cutoff_distance.
    '''
    contact_degree = 0

    for r2 in range(1, pose_to_compare_with.size() + 1):
        if target_residue.nbr_atom_xyz().distance(pose_to_compare_with.residue(r2).nbr_atom_xyz()) < cutoff_distance:
            contact_degree += 1

    return contact_degree

def anchor_is_good(target_pose, fuzz_pose, target_anchor, fuzz_anchor, ligand_residue, ligand_contact_degree_cutoff=10):
    '''Return true if an anchor is good.'''
    if intra_residue_heavy_atom_clash(fuzz_pose.residue(fuzz_anchor)):
        return False
    
    for res in range(1, target_pose.size()):
        if res == target_anchor: continue

        if residue_heavy_atom_clashes(fuzz_pose.residue(ligand_residue), target_pose.residue(res), res2_sc=False):
            return False

        if -1 <= target_anchor - res <=1:
            if residue_heavy_atom_clashes(fuzz_pose.residue(fuzz_anchor), target_pose.residue(res), res1_bb=False, res2_sc=False):
                return False

        elif residue_heavy_atom_clashes(fuzz_pose.residue(fuzz_anchor), target_pose.residue(res), res2_sc=False):
            return False

    if residue_contact_degree(fuzz_pose.residue(ligand_residue), target_pose) < ligand_contact_degree_cutoff:
        return False

    return True

def find_matched_rotamers_for_residue_pairs(target_pose, target_seqpos, fuzz_pose, motif_res, allowed_rotamers, rmsd_cutoff=2):
    '''Find all rotamers that matches the residue on the target pose to the motif residue.
    
    Return:
        The rotamer id that has lowest RMSD and the RMSD.
    '''
    # Muatet the target residue to the motif residue

    mutater = rosetta.protocols.simple_moves.MutateResidue()
    mutater.set_res_name(fuzz_pose.residue(motif_res).name3())
    mutater.set_target(target_seqpos)
    mutater.apply(target_pose)

    # Test the rotamers

    matched_rotamers = []

    rotamer_set = rosetta.core.pack.rotamer_set.bb_independent_rotamers( fuzz_pose.residue(motif_res).type(), True )

    for i in allowed_rotamers:
        replace_intra_residue_torsions(target_pose, target_seqpos, rotamer_set[i])
        rmsd = sc_heavy_atom_rmsd(target_pose.residue(target_seqpos), fuzz_pose.residue(motif_res))
        if rmsd < rmsd_cutoff:
            matched_rotamers.append((i, rmsd))
            #target_pose.dump_pdb('debug/target_matched_{0}_{1}_{2}.pdb'.format(target_seqpos, motif_res, i))###DEBUG
            #fuzz_pose.dump_pdb('debug/fuzz_matched_{0}_{1}_{2}.pdb'.format(target_seqpos, motif_res, i))#DEBUG

    if len(matched_rotamers) == 0: return None

    return min(matched_rotamers, key=lambda x : x[1])


def find_matched_rotamers_for_anchored_fuzz_ball(target_pose, target_matching_seqposes,
        fuzz_pose, ligand_residue, motif_residues, fuzz_anchor, bb_compatible_rotamers):
    '''Find all rotamers that could match the anchored fuzz ball to the target pose.
    Args:
        target_pose: The target pose.
        target_matching_seqposes: The residues on the target pose that are available for
            motif matching.
        fuzz_pose: The fuzz ball pose.
        ligand_residue: The id of the ligand residue.
        motif_residues: A list of motif residues.
        fuzz_anchor: The residue on the fuzz ball that is anchored to the target
        bb_compatible_rotamers: The backbone compatible rotamer ids on the target pose
            at each position

    Return:
        A list of Matches
    '''
    def ca_nbr_distance(residue1, residue2):
        '''Calculate the distance between the CA atom of residue1
        and the NBR atom of residue2.
        '''
        return residue1.xyz('CA').distance(residue2.nbr_atom_xyz())

    # Find matches for each motif residue

    matches = []

    for motif_res in motif_residues:
        if motif_res == fuzz_anchor: continue

        # Find all matchable positions with in 5 angstrom

        positions_to_test = [i for i in target_matching_seqposes 
                if ca_nbr_distance(target_pose.residue(i), fuzz_pose.residue(motif_res)) < 5]

        # Find all matches

        for target_seqpos in positions_to_test:
            matched_rotamer = find_matched_rotamers_for_residue_pairs(target_pose, target_seqpos, fuzz_pose,
                    motif_res, bb_compatible_rotamers[target_seqpos][fuzz_pose.residue(motif_res).name3()])
            if matched_rotamer:
                matches.append(Match())
                
                matches[-1].fuzz_ball_anchor_residue = fuzz_anchor
                matches[-1].target_matched_residue = target_seqpos
                matches[-1].fuzz_ball_matched_residue = motif_res
                matches[-1].target_matched_rotamer = matched_rotamer[0]
                matches[-1].match_rmsd = matched_rotamer[1]

    return matches

def find_matched_rotamers_for_fuzz_ball(target_pose_original, target_matching_seqposes,
        fuzz_pose, ligand_residue, motif_residues, bb_compatible_rotamers):
    '''Find all rotamers that could match the fuzz ball to the target pose.
    Args:
        target_pose_original: The target pose.
        target_matching_seqposes: The residues on the target pose that are available for
            motif matching.
        fuzz_pose: The fuzz ball pose.
        ligand_residue: The id of the ligand residue.
        motif_residues: A list of motif residues.
        bb_compatible_rotamers: The backbone compatible rotamer ids on the target pose
            at each position

    Return:
        A list of lists for matches of each anchor.
    '''
    # Make a local clone of the target pose

    target_pose = target_pose_original.clone()

    # Set up the fold tree that roots on the ligand residue

    set_up_fold_tree_for_root_residue(fuzz_pose, ligand_residue)

    # Find the matches for each of the motif residues anchored
    
    matches = []

    for fuzz_anchor in motif_residues:

        rotamer_set = rosetta.core.pack.rotamer_set.bb_independent_rotamers( fuzz_pose.residue(fuzz_anchor).type(), True )

        for target_anchor_seqpos in target_matching_seqposes:
            
            for i in bb_compatible_rotamers[target_anchor_seqpos][fuzz_pose.residue(fuzz_anchor).name3()]:
                
                set_inverse_rotamer(fuzz_pose, fuzz_anchor, rotamer_set[i])
            
                match_anchor_position(target_pose, target_anchor_seqpos, fuzz_pose, fuzz_anchor)
          
                if not anchor_is_good(target_pose, fuzz_pose, target_anchor_seqpos, fuzz_anchor, ligand_residue):
                    continue

                matches_for_anchor = find_matched_rotamers_for_anchored_fuzz_ball(target_pose, target_matching_seqposes,
                    fuzz_pose, ligand_residue, motif_residues, fuzz_anchor, bb_compatible_rotamers)
                for m in matches_for_anchor:
                    m.target_anchor_residue = target_anchor_seqpos 
                    m.fuzz_ball_anchor_rotamer = i
              
                if len(matches_for_anchor) > 0:
                    matches.append(matches_for_anchor)
        
                #print len(matches_for_anchor)###DEBUG

    return matches


if __name__ == '__main__':
    import preprocessing
    
    pyrosetta.init(options='-extra_res_fa test/inputs/REN_no_charge_from_mol2.params -mute all')
   
    print 'Start preprocessing.'

    fuzz_pose = preprocessing.load_cleaned_filtered_fuzz_pose('test/inputs/fuzz_balls_small/binding_site_from_james_renumbered.pdb', 1)
    
    target_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(target_pose, 'test/inputs/target_pdbs/1svx.pdb')
    matchable_positions = preprocessing.find_interface_seqposes_noGP(target_pose, 'A', 'B', cutoff_distance=9)
    bb_compatible_rotamers = preprocessing.get_bb_compatible_rotamers_for_pose(target_pose, matchable_positions)

    print 'Start matching.'

    find_matched_rotamers_for_fuzz_ball(target_pose, matchable_positions, fuzz_pose, 1, list(range(2, fuzz_pose.size() + 1)), bb_compatible_rotamers)
