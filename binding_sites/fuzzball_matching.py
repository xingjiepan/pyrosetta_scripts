#!/usr/bin/env python2.7

import numpy as np

import pyrosetta
from pyrosetta import rosetta


class Match:
    '''A record of a match.'''
    def __init__(self):
        target_anchor_residue = 0
        fuzz_ball_anchor_residue = 0
        fuzz_ball_anchor_rotamer = 0
        
        target_matched_residue = 0
        fuzz_ball_matched_residue = 0
        fuzz_ball_matched_rotamer = 0

        match_rmsd = 999

    def __repr__(self):
        return str(self.__dict__)

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

def anchor_is_good(target_pose, fuzz_pose, target_anchor, fuzz_anchor, ligand_residue):
    '''Return true if an anchor is good.'''
    for res in range(1, target_pose.size()):
        if res == target_anchor: continue

        if residue_heavy_atom_clashes(fuzz_pose.residue(ligand_residue), target_pose.residue(res), res2_sc=False):
            return False

        if -1 <= target_anchor - res <=1: continue
        
        if residue_heavy_atom_clashes(fuzz_pose.residue(fuzz_anchor), target_pose.residue(res), res2_sc=False):
            return False

    return True

def find_matched_rotamers_for_residue_pairs(target_pose, target_seqpos, fuzz_pose, motif_res, rmsd_cutoff=2):
    '''Find all rotamers that matches the residue on the target pose to the motif residue.
    
    Return:
        The rotamer id that has lowest RMSD and the RMSD.
    '''
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

    def sc_heavy_atom_rmsd(residue1, residue2):
        xyzs1 = side_chain_heavy_atom_xyzs(residue1)
        xyzs2 = side_chain_heavy_atom_xyzs(residue2)
        return RMSD(xyzs1, xyzs2)

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
        fuzz_pose, ligand_residue, motif_residues, fuzz_anchor, fuzz_ball_anchor_rotamer_id):
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
                matches[-1].fuzz_ball_anchor_rotamer = fuzz_ball_anchor_rotamer_id
                matches[-1].target_matched_residue = target_seqpos
                matches[-1].fuzz_ball_matched_residue = motif_res
                matches[-1].fuzz_ball_matched_rotamer = matched_rotamer[0]
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

    # Find the matches for each of the motif residues achored
    
    matches = []

    for fuzz_anchor in motif_residues:

        rotamer_set = rosetta.core.pack.rotamer_set.bb_independent_rotamers( fuzz_pose.residue(fuzz_anchor).type(), True )

        for i in range(1, len(rotamer_set) + 1):
            
            set_inverse_rotamer(fuzz_pose, fuzz_anchor, rotamer_set[i])
            
            for target_anchor_seqpos in target_matching_seqposes:
                match_anchor_position(target_pose, target_anchor_seqpos, fuzz_pose, fuzz_anchor)
           
                if not anchor_is_good(target_pose, fuzz_pose, target_anchor_seqpos, fuzz_anchor, ligand_residue):
                    continue

                matches_for_achor.append(find_matched_rotamers_for_anchored_fuzz_ball(target_pose, target_matching_seqposes,
                    fuzz_pose, ligand_residue, motif_residues, fuzz_anchor, i))
               
                matches += matches_for_achor
               
                print matches_for_achor
                #fuzz_pose.dump_pdb('debug/test_fuzz_{0}_{1}.pdb'.format(fuzz_anchor, i)) ###DEBUG
                #target_pose.dump_pdb('debug/target.pdb') ###DEBUG
                #exit()

    
    print matches ###DEBUG
    return matches


if __name__ == '__main__':
    pyrosetta.init(options='-extra_res_fa inputs/REN_no_charge_from_mol2.params')
    
    fuzz_pose = rosetta.core.pose.Pose()
    #rosetta.core.import_pose.pose_from_file(fuzz_pose, 'inputs/binding_site_from_james_renumbered.pdb')
    rosetta.core.import_pose.pose_from_file(fuzz_pose, 'inputs/binding_site_from_james_all.pdb')

    #####OBSOLETE
    #target_pose = rosetta.core.pose.Pose()
    #rosetta.core.import_pose.pose_from_file(target_pose, 'inputs/2lvb_no_terms.pdb')
    #find_matched_rotamers_for_fuzz_ball(target_pose, 6, list(range(37, 50)) + list(range(63, 77)),
    #        fuzz_pose, 1, list(range(2, fuzz_pose.size() + 1)))
    ######OBSOLETE

    target_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(target_pose, 'inputs/3tdn_barrel.pdb')
    matchable_positions_pdb = [130, 80, 171, 101, 48, 23, 5, 7, 9, 201, 144, 169, 126, 128, 103, 225, 224, 222, 78, 55, 50, 52]
    matchable_positions = [target_pose.pdb_info().pdb2pose('B', i) for i in matchable_positions_pdb]

    #fuzz_pose.dump_pdb('debug/test_fuzz.pdb') ###DEBUG
    find_matched_rotamers_for_fuzz_ball(target_pose, matchable_positions, fuzz_pose, 1, list(range(2, fuzz_pose.size() + 1)))
