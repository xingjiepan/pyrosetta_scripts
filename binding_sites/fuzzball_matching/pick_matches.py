#!/usr/bin/env python2.7

import pyrosetta
from pyrosetta import rosetta

from fuzzball_matching_basic import *


def pick_lowest_rmsd_matches(matches):
    '''Find the subset of matches that don't share same
    motif/position and have lowest RMSDs.
    '''
    if 0 == len(matches): return []

    picked_matches = []

    matched_target_residues = [matches[0].target_anchor_residue]
    matched_fuzz_residues = [matches[0].fuzz_ball_anchor_residue]

    sorted_matches = sorted(matches, key=lambda m : m.match_rmsd)

    for m in sorted_matches:
        if m.target_matched_residue in matched_target_residues:
            continue
        if m.fuzz_ball_matched_residue in matched_fuzz_residues:
            continue

        picked_matches.append(m)
        matched_target_residues.append(m.target_matched_residue)
        matched_fuzz_residues.append(m.fuzz_ball_matched_residue)

    return picked_matches 

def pick_non_clashing_lowest_rmsd_matches(target_pose_original, fuzz_pose_original, matches, ligand_residue):
    '''Find the subset of matches that don't share same
    motif/position, don't clash and have lowest RMSDs.
    '''
    if 0 == len(matches): return []
    picked_matches = []
    
    # Make copies of the poses
  
    target_pose = target_pose_original.clone()
    fuzz_pose = fuzz_pose_original.clone()
    
    # Mutate all residues on the target to ALA
    
    mutate_residues(target_pose, range(1, target_pose.size() + 1), 'ALA')

    # Align the anchor

    set_rotamer_and_match_anchor(target_pose, matches[0].target_anchor_residue, fuzz_pose,
            matches[0].fuzz_ball_anchor_residue, matches[0].fuzz_ball_anchor_rotamer)
    
    # Find the compatible matched residues

    matches_with_anchor = [create_a_pseudo_match_for_anchor(matches[0])] + matches

    matched_target_residues = []
    matched_fuzz_residues = []
    sorted_matches = sorted(matches_with_anchor, key=lambda m : m.match_rmsd)

    for m in sorted_matches:
        fuzz_res = m.fuzz_ball_matched_residue 
        target_seqpos = m.target_matched_residue
        
        if target_seqpos in matched_target_residues:
            continue
        if fuzz_res in matched_fuzz_residues:
            continue

        # Mutate the residue and apply the rotamer

        apply_match(target_pose, fuzz_pose, m)

        # Check clashes
    
        clash = intra_residue_heavy_atom_clash(target_pose.residue(target_seqpos))
        
        if residue_heavy_atom_clashes(fuzz_pose.residue(ligand_residue), target_pose.residue(target_seqpos)):
            clash = True
        
        for res in range(1, target_pose.size()):
            if res == target_seqpos: continue

            if -1 <= target_seqpos - res <=1:
                if residue_heavy_atom_clashes(target_pose.residue(target_seqpos), target_pose.residue(res), res1_bb=False):
                    clash = True
                    break

            elif residue_heavy_atom_clashes(target_pose.residue(target_seqpos), target_pose.residue(res)):
                clash = True
                break

        # Accept the non-clashing match

        if clash:
            mutate_residues(target_pose, [target_seqpos], 'ALA') 
        else:
            picked_matches.append(m)
            matched_target_residues.append(target_seqpos)
            matched_fuzz_residues.append(fuzz_res)

    return picked_matches 

def h_bond_maintained(target_pose, target_res, target_ligand, fuzz_pose, fuzz_res, fuzz_ligand, cutoff_ratio=0.5):
    '''Return true if the hydrogen bond is maintained.
    If the H-bond energy is raised above the cutoff ratio,
    the H-bond between two residues is considered as not maintained.
    '''
    # Score the poses
    
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(target_pose)
    sfxn(fuzz_pose)

    # Check if the H-bond is maintained

    e_hb_target = residue_residue_hbond_energy(target_pose, target_res, target_ligand)
    e_hb_fuzz = residue_residue_hbond_energy(fuzz_pose, fuzz_res, fuzz_ligand)

    return e_hb_target <= cutoff_ratio * e_hb_fuzz 

def get_scores_for_matches(target_pose_original, fuzz_pose_original, matches, ligand_residue):
    '''Get the scores for the matched residues.
    Return:
        M_scores: a matrix of of scores. The diagonal is the one body scores
            of matches and the element M_scores[i][j] (j != i) is the interaction energy between
            two matched residues.
        ligand_interaction_scores: a list of the total residue-residue scores between a
            residue and the ligand
        matches_with_anchor: the maches with the anchor prepended
        match_ids_problematic: the list of matches that breaks H-bonds or have ligand
            interaction energy above -1
    '''
    # Make copies of the poses
  
    target_pose = target_pose_original.clone()
    fuzz_pose = fuzz_pose_original.clone()
    
    # Mutate all residues on the target to ALA
    
    mutate_residues(target_pose, range(1, target_pose.size() + 1), 'ALA')

    # Align the anchor and insert the ligand to the target pose

    set_rotamer_and_match_anchor(target_pose, matches[0].target_anchor_residue, fuzz_pose,
            matches[0].fuzz_ball_anchor_residue, matches[0].fuzz_ball_anchor_rotamer)
    target_pose.append_residue_by_jump(fuzz_pose.residue(ligand_residue), 1)
    
    # Find the compatible matched residues

    matches_with_anchor = [create_a_pseudo_match_for_anchor(matches[0])] + matches

    # Set up scoring

    sfxn = rosetta.core.scoring.get_score_function()
    
    sfxn(target_pose)
    ref_score = target_pose.energies().total_energy()

    M_scores = [[0 for j in range(len(matches_with_anchor))] for i in range(len(matches_with_anchor))] 
    ligand_interaction_scores = []
    match_ids_problematic = []

    # Calculate all single residue scores
   
    for i in range(len(matches_with_anchor)):
        apply_match(target_pose, fuzz_pose, matches_with_anchor[i])
        
        sfxn(target_pose)
        M_scores[i][i] = target_pose.energies().total_energy() - ref_score
        
        # Check if the match has good ligand interaction energy
        
        ligand_interaction_scores.append(residue_residue_total_energy(target_pose,
            matches_with_anchor[i].target_matched_residue, target_pose.size()))

        if ligand_interaction_scores[-1] > -1:
            match_ids_problematic.append(i)

        # Check if the match breaks hbonds

        elif not h_bond_maintained(target_pose, matches_with_anchor[i].target_matched_residue, target_pose.size(),
                fuzz_pose, matches_with_anchor[i].fuzz_ball_matched_residue, ligand_residue, cutoff_ratio=0.5):
            match_ids_problematic.append(i)

        # Set check if the match has a hbond interaction

        if residue_residue_hbond_energy(target_pose, matches_with_anchor[i].target_matched_residue,
                target_pose.size()) < 0:
            matches_with_anchor[i].hbond_match = True

        mutate_residues(target_pose, [matches_with_anchor[i].target_matched_residue], 'ALA') 

    # Calculate all pairwise residue scores

    for i in range(len(matches_with_anchor)):
        if i in match_ids_problematic: continue

        for j in range(i + 1, len(matches_with_anchor)):
            if j in match_ids_problematic: continue
            
            if matches_with_anchor[i].target_matched_residue == matches_with_anchor[j].target_matched_residue:
                M_scores[i][j] = 99999 # Set the interaction energy to 'infinity' if two matches have the same position
                M_scores[j][i] = M_scores[i][j]
                continue
            
            apply_match(target_pose, fuzz_pose, matches_with_anchor[i])
            apply_match(target_pose, fuzz_pose, matches_with_anchor[j])

            sfxn(target_pose)
            M_scores[i][j] = target_pose.energies().total_energy() - (ref_score + M_scores[i][i] + M_scores[j][j])
            M_scores[j][i] = M_scores[i][j]

            mutate_residues(target_pose, [matches_with_anchor[i].target_matched_residue], 'ALA') 
            mutate_residues(target_pose, [matches_with_anchor[j].target_matched_residue], 'ALA') 
    
    return M_scores, ligand_interaction_scores, matches_with_anchor, match_ids_problematic

def pick_lowest_score_matches_greedy(target_pose_original, fuzz_pose_original, matches, ligand_residue,
        cutoff_score=10, min_hbond_match=1, ligand_interaction_score_cutoff=0):
    '''Pick lowest score matches using a greedy algorithm.
    Reject the binding site if the number of hbond matches
    is smaller than the min_hbond_match value.
    Reject the binding site if the total interaction energy
    betwen the ligand and motif residues are above the
    ligand_interaction_score_cutoff.

    Return the picked matches and the total_ligand_interaction_score
    '''
    def score_change_after_adding_match(M_scores, accepted_match_ids, new_id, match_ids_problematic):
        '''Return the score change after adding a match.
        If a match breaks hbonds, reject it.
        '''
        if new_id in accepted_match_ids + match_ids_problematic:
            return float('inf')

        change = M_scores[new_id][new_id] 
        for i in accepted_match_ids:
            change += M_scores[i][new_id]

        return change

    if 0 == len(matches): return [], 0
   
    M_scores, ligand_interaction_scores, matches_with_anchor, match_ids_problematic = get_scores_for_matches(
            target_pose_original, fuzz_pose_original, matches, ligand_residue)
    total_ligand_interaction_score = 0

    accepted_match_ids = []
    for i in range(len(matches_with_anchor)):

        # Find the best match to add

        score_changes = [(j, score_change_after_adding_match(M_scores, accepted_match_ids, j, match_ids_problematic))
                for j in range(len(matches_with_anchor))]

        best_match = min(score_changes, key=lambda x : x[1])

        if best_match[1] > cutoff_score:
            break

        accepted_match_ids.append(best_match[0])
        total_ligand_interaction_score += ligand_interaction_scores[best_match[0]]

    # Reject the binding site if not enough hbond interaction is found

    if sum(1 if matches_with_anchor[i].hbond_match else 0 for i in accepted_match_ids) < min_hbond_match:
        return [], 0

    if total_ligand_interaction_score > ligand_interaction_score_cutoff:
        return [], 0

    return [matches_with_anchor[i] for i in accepted_match_ids], total_ligand_interaction_score

def dump_matches_for_an_anchor(target_pose_original, fuzz_pose_original, ligand_residue, matches,
        target_output_file, matched_fuzz_output_file):
    '''Dump matches for an anchor.'''
    
    # Make copies of the poses
  
    target_pose = target_pose_original.clone()
    fuzz_pose = fuzz_pose_original.clone()

    # Mutate all residues on the target to ALA
    
    mutate_residues(target_pose, range(1, target_pose.size() + 1), 'ALA')

    # Align the anchor

    set_rotamer_and_match_anchor(target_pose, matches[0].target_anchor_residue, fuzz_pose,
            matches[0].fuzz_ball_anchor_residue, matches[0].fuzz_ball_anchor_rotamer)
    
    # Set the rotamers on the matched fuzz pose and target pose

    for match in matches:
        apply_match(target_pose, fuzz_pose, match)

    # Extract the matched part of the fuzz pose

    matched_fuzz_pose = rosetta.core.pose.Pose(fuzz_pose, ligand_residue, ligand_residue)
    for m in matches:
        matched_fuzz_pose.append_residue_by_jump(fuzz_pose.residue(m.fuzz_ball_matched_residue).clone(), 1)

    target_pose.dump_pdb(target_output_file)
    matched_fuzz_pose.dump_pdb(matched_fuzz_output_file)

