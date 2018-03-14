#!/usr/bin/env python2.7

import os

import numpy as np

import pyrosetta
from pyrosetta import rosetta


class Match:
    '''A record of a match.'''
    def __init__(self):
        self.target_anchor_residue = 0
        self.fuzz_ball_anchor_residue = 0
        self.fuzz_ball_anchor_rotamer = 0
        
        self.target_matched_residue = 0
        self.target_matched_rotamer = 0
        self.fuzz_ball_matched_residue = 0

        self.match_rmsd = 999

    def __repr__(self):
        return str(self.__dict__)

def clean_fuzz_pose(fuzz_pose, ligand_residue):
    '''The fuzz pose might have a lot of residues that are
    nearly identical, remove the redundancy.
    After cleaning, the ligand will be the first residue.
    '''
    cleaned_pose = rosetta.core.pose.Pose(fuzz_pose, ligand_residue, ligand_residue)

    for i in range(1, fuzz_pose.size() + 1):    
        if i == ligand_residue: continue

        redundant = False
        for j in range(2, cleaned_pose.size() + 1):
            if fuzz_pose.residue(i).name3() == cleaned_pose.residue(j).name3():
                if sc_heavy_atom_rmsd(fuzz_pose.residue(i), cleaned_pose.residue(j)) < 0.2:
                    redundant = True

        if not redundant:
            cleaned_pose.append_residue_by_jump(fuzz_pose.residue(i).clone(), 1)

    return cleaned_pose

def residue_residue_total_energy(pose, res1, res2):
    '''Get the total interaction energy between two residues.'''
    e_edge = pose.energies().energy_graph().find_energy_edge(res1, res2)

    if e_edge:
        return e_edge.dot(pose.energies().weights())
    else:
        return 0

def filter_motif_residues(fuzz_pose, ligand_residue, energy_cutoff=-2):
    '''Filter the motif residues by their interaction energy with
    the ligand residue.
    Return a new pose with the bad motifs filtered out.
    '''
    new_pose = rosetta.core.pose.Pose(fuzz_pose, ligand_residue, ligand_residue)
    
    sfxn = rosetta.core.scoring.get_score_function()
   
    sfxn(fuzz_pose)

    for i in range(1, fuzz_pose.size() + 1):
        if i == ligand_residue: continue
            
        if residue_residue_total_energy(fuzz_pose, ligand_residue, i) < energy_cutoff:
            new_pose.append_residue_by_jump(fuzz_pose.residue(i).clone(), 1)

    return new_pose

def intra_residue_heavy_atom_clash(residue, cutoff_distance=2.5):
    '''Return true if there are clashes inside one residue.'''
    N = residue.nheavyatoms()
   
    for i in range(1, N):
        for j in range(i + 1, N + 1):
            if residue.path_distance(i, j) > 3:
                if residue.xyz(i).distance(residue.xyz(j)) < cutoff_distance:
                    return True
    
    return False

def residue_heavy_atom_clashes(residue1, residue2, cutoff_distance=2.5,
        res1_bb=True, res1_sc=True, res2_bb=True, res2_sc=True):
    '''Return true if two residues have heavy atoms that clash.
    NOTE: CB is also considered in BB clash checking.
    '''
    # No clash if two residues are far apart

    if residue1.nbr_atom_xyz().distance(residue2.nbr_atom_xyz()) > 15:
        return False
   
    # Check atom clashes

    res1_start = 1 if res1_bb else residue1.first_sidechain_atom()
    res1_stop = residue1.nheavyatoms() if res1_sc else (residue1.nheavyatoms() if residue1.name3() == 'GLY' else residue1.atom_index('CB'))
    res2_start = 1 if res2_bb else residue2.first_sidechain_atom()
    res2_stop = residue2.nheavyatoms() if res2_sc else (residue2.nheavyatoms() if residue2.name3() == 'GLY' else residue2.atom_index('CB'))

    for i in range(res1_start, res1_stop + 1):
        p1 = residue1.xyz(i)
        for j in range(res2_start, res2_stop + 1):
            p2 = residue2.xyz(j)

            if p1.distance(p2) < cutoff_distance:
                return True

    return False

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
    if pose.residue(seqpos).name3() == 'ALA' or pose.residue(seqpos).name3() == 'GLY':
        return
    
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

def mutate_residues(pose, residues, res_name):
    '''Mutate residues in a pose to a given type.'''
    mutater = rosetta.protocols.simple_moves.MutateResidue()
    mutater.set_res_name(res_name)

    for res in residues:
        mutater.set_target(res)
        mutater.apply(pose)

def create_a_pseudo_match_for_anchor(match):
    '''Given a match, create a pseudo match for its anchor.'''
    m = Match()
    m.target_anchor_residue = match.target_anchor_residue
    m.fuzz_ball_anchor_residue = match.fuzz_ball_anchor_residue
    m.fuzz_ball_anchor_rotamer = match.fuzz_ball_anchor_rotamer
    m.target_matched_residue = match.target_anchor_residue
    m.target_matched_rotamer = match.fuzz_ball_anchor_rotamer
    m.fuzz_ball_matched_residue = match.fuzz_ball_anchor_residue
    m.match_rmsd = 0

    return m

def set_rotamer_and_match_anchor(target_pose, target_anchor_seqpos, fuzz_pose, fuzz_anchor, rotamer_id):
    '''Set rotamer on the fuzz ball anchor and match the fuzz pose anchor
    to the target_pose.
    '''
    rotamer_set = rosetta.core.pack.rotamer_set.bb_independent_rotamers( fuzz_pose.residue(fuzz_anchor).type(), True )
    set_inverse_rotamer(fuzz_pose, fuzz_anchor, rotamer_set[rotamer_id])

    match_anchor_position(target_pose, target_anchor_seqpos, fuzz_pose, fuzz_anchor)

def apply_match(target_pose, fuzz_pose, match):
    '''Apply a match to a target_pose and fuzz_pose.'''
    fuzz_res = match.fuzz_ball_matched_residue 
    target_seqpos = match.target_matched_residue
    
    rotamer_set = rosetta.core.pack.rotamer_set.bb_independent_rotamers( fuzz_pose.residue(fuzz_res).type(), True )
    set_inverse_rotamer(fuzz_pose, fuzz_res, rotamer_set[match.target_matched_rotamer])
    
    mutate_residues(target_pose, [target_seqpos], fuzz_pose.residue(fuzz_res).name3())

    replace_intra_residue_torsions(target_pose, target_seqpos, rotamer_set[match.target_matched_rotamer])

def anchor_is_good(target_pose, fuzz_pose, target_anchor, fuzz_anchor, ligand_residue):
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

    return True

def sc_heavy_atom_rmsd(residue1, residue2):
    '''Calculate the side chain function group RMSD between two residues.'''
    def side_chain_heavy_atom_xyzs(residue):
        side_chain_func_group_map = {
                'GLY':[],
                'ALA':['CB'],
                'PRO':['CB', 'CG', 'CD'],
                'VAL':['CB', 'CG1', 'CG2'],
                'LEU':['CB', 'CG', 'CD1', 'CD2'],
                'ILE':['CB', 'CG1', 'CG2', 'CD1'],
                'MET':['CG', 'SD', 'CE'],

                'PHE':['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
                'TYR':['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
                'TRP':['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],

                'SER':['CB', 'OG'],
                'THR':['CB', 'OG1', 'CG2'],
                'CYS':['CB', 'SG'],

                'ASP':['CB', 'CG', 'OD1', 'OD2'],
                'GLU':['CG', 'CD', 'OE1', 'OE2'],
                'ASN':['CB', 'CG', 'OD1', 'ND2'],
                'GLN':['CG', 'CD', 'OE1', 'NE2'],
                
                'LYS':['CD', 'CE', 'NZ'],
                'ARG':['NE', 'CZ', 'NH1', 'NH2'],
                'HIS':['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
                }
        
        #sc_heavy_atoms = list(range(residue.first_sidechain_atom(), residue.nheavyatoms() + 1))
        #return [residue.xyz(i) for i in sc_heavy_atoms]
        
        return [residue.xyz(a) for a in side_chain_func_group_map[residue.name3()]]

    def RMSD(xyzs1, xyzs2):
        assert(len(xyzs1) == len(xyzs2))
        diffs = [xyzs1[i] - xyzs2[i] for i in range(len(xyzs1))]
        return np.sqrt(sum(d.length_squared() for d in diffs) / len(diffs))

    xyzs1 = side_chain_heavy_atom_xyzs(residue1)
    xyzs2 = side_chain_heavy_atom_xyzs(residue2)
    return RMSD(xyzs1, xyzs2)

def find_matched_rotamers_for_residue_pairs(target_pose, target_seqpos, fuzz_pose, motif_res, rmsd_cutoff=2):
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
        
    for i in range(1, len(rotamer_set) + 1):
        replace_intra_residue_torsions(target_pose, target_seqpos, rotamer_set[i])
        rmsd = sc_heavy_atom_rmsd(target_pose.residue(target_seqpos), fuzz_pose.residue(motif_res))
        if rmsd < rmsd_cutoff:
            matched_rotamers.append((i, rmsd))
            #target_pose.dump_pdb('debug/target_matched_{0}_{1}_{2}.pdb'.format(target_seqpos, motif_res, i))###DEBUG
            #fuzz_pose.dump_pdb('debug/fuzz_matched_{0}_{1}_{2}.pdb'.format(target_seqpos, motif_res, i))#DEBUG

    if len(matched_rotamers) == 0: return None

    return min(matched_rotamers, key=lambda x : x[1])


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
                    motif_res)
            if matched_rotamer:
                matches.append(Match())
                
                matches[-1].fuzz_ball_anchor_residue = fuzz_anchor
                matches[-1].target_matched_residue = target_seqpos
                matches[-1].fuzz_ball_matched_residue = motif_res
                matches[-1].target_matched_rotamer = matched_rotamer[0]
                matches[-1].match_rmsd = matched_rotamer[1]

    return matches

def find_matched_rotamers_for_fuzz_ball(target_pose, target_matching_seqposes,
        fuzz_pose, ligand_residue, motif_residues):
    '''Find all rotamers that could match the fuzz ball to the target pose.
    Args:
        target_pose: The target pose.
        target_matching_seqposes: The residues on the target pose that are available for
            motif matching.
        fuzz_pose: The fuzz ball pose.
        ligand_residue: The id of the ligand residue.
        motif_residues: A list of motif residues.
        
    Return:
        A list of lists for matches of each anchor.
    '''
   
    # Set up the fold tree that roots on the ligand residue

    set_up_fold_tree_for_root_residue(fuzz_pose, ligand_residue)

    # Find the matches for each of the motif residues anchored
    
    matches = []

    for fuzz_anchor in motif_residues:
    #for fuzz_anchor in [10]: ###DEBUG

        rotamer_set = rosetta.core.pack.rotamer_set.bb_independent_rotamers( fuzz_pose.residue(fuzz_anchor).type(), True )

        for i in range(1, len(rotamer_set) + 1):
            
            set_inverse_rotamer(fuzz_pose, fuzz_anchor, rotamer_set[i])
            
            for target_anchor_seqpos in target_matching_seqposes:
                match_anchor_position(target_pose, target_anchor_seqpos, fuzz_pose, fuzz_anchor)
          
                if not anchor_is_good(target_pose, fuzz_pose, target_anchor_seqpos, fuzz_anchor, ligand_residue):
                    continue

                matches_for_anchor = find_matched_rotamers_for_anchored_fuzz_ball(target_pose, target_matching_seqposes,
                    fuzz_pose, ligand_residue, motif_residues, fuzz_anchor)
                for m in matches_for_anchor:
                    m.target_anchor_residue = target_anchor_seqpos 
                    m.fuzz_ball_anchor_rotamer = i
               
                matches += matches_for_anchor
        
                print len(matches_for_anchor)
                if len(matches_for_anchor) > 8: 
                    picked_matches = pick_non_clashing_lowest_rmsd_matches(target_pose, fuzz_pose, matches_for_anchor, ligand_residue)
                    print fuzz_anchor, i, target_anchor_seqpos, len(picked_matches)
              
                    if len(picked_matches) > 3:
                        output_path = os.path.join('debug', '{0}_{1}_{2}_{3}'.format(len(picked_matches) + 1, fuzz_anchor, i, target_anchor_seqpos))

                        if not os.path.exists(output_path):
                            os.mkdir(output_path)

                        dump_matches_for_an_anchor(target_pose, fuzz_pose, ligand_residue, picked_matches,
                                os.path.join(output_path, 'target_pose.pdb'), os.path.join(output_path, 'matched_fuzz_pose.pdb'))
                        exit()###DEBUG

                #fuzz_pose.dump_pdb('debug/test_fuzz_{0}_{1}.pdb'.format(fuzz_anchor, i)) ###DEBUG
                #target_pose.dump_pdb('debug/target.pdb') ###DEBUG

    
    print matches ###DEBUG
    return matches

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
    matched_fuzz_pose.append_residue_by_jump(fuzz_pose.residue(matches[0].fuzz_ball_anchor_residue).clone(), 1)
    for m in matches:
        matched_fuzz_pose.append_residue_by_jump(fuzz_pose.residue(m.fuzz_ball_matched_residue).clone(), 1)

    target_pose.dump_pdb(target_output_file)
    matched_fuzz_pose.dump_pdb(matched_fuzz_output_file)


if __name__ == '__main__':
    #pyrosetta.init(options='-extra_res_fa inputs/REN_no_charge_from_mol2.params')
    pyrosetta.init(options='-extra_res_fa inputs/REN_no_charge_from_mol2.params -mute all')
    
    fuzz_pose = rosetta.core.pose.Pose()
    #rosetta.core.import_pose.pose_from_file(fuzz_pose, 'inputs/binding_site_from_james_renumbered.pdb')
    #rosetta.core.import_pose.pose_from_file(fuzz_pose, 'inputs/binding_site_from_james_all.pdb')
    rosetta.core.import_pose.pose_from_file(fuzz_pose, 'inputs/binding_site_from_james_all_cleaned.pdb')
    fuzz_pose = clean_fuzz_pose(fuzz_pose, 1)
    fuzz_pose = filter_motif_residues(fuzz_pose, 1)

    target_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(target_pose, 'inputs/3tdn_barrel.pdb')
    matchable_positions_pdb = [130, 80, 171, 101, 48, 23, 5, 7, 9, 201, 144, 169, 126, 128, 103, 225, 224, 222, 78, 55, 50, 52]
    matchable_positions = [target_pose.pdb_info().pdb2pose('B', i) for i in matchable_positions_pdb]

    #fuzz_pose.dump_pdb('debug/test_fuzz.pdb') ###DEBUG
    #print fuzz_pose.size() ###DEBUG
    find_matched_rotamers_for_fuzz_ball(target_pose, matchable_positions, fuzz_pose, 1, list(range(2, fuzz_pose.size() + 1)))
