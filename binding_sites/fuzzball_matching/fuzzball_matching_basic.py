#!/usr/bin/env python2.7
'''Basic classes and  functions for fuzzball matching.
'''

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

        self.hbond_match = False

    def __repr__(self):
        return str(self.__dict__)

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

def residue_residue_total_energy(pose, res1, res2):
    '''Get the total interaction energy between two residues.'''
    e_edge = pose.energies().energy_graph().find_energy_edge(res1, res2)

    if e_edge:
        return e_edge.dot(pose.energies().weights())
    else:
        return 0

def residue_residue_hbond_energy(pose, res1, res2):
    '''Get the H-bond energy between two residues.'''
    e_edge = pose.energies().energy_graph().find_energy_edge(res1, res2)
    e_map = e_edge.fill_energy_map()

    return e_map[rosetta.core.scoring.hbond_sc]

def mutate_residues(pose, residues, res_name, keep_g_p=True):
    '''Mutate residues in a pose to a given type.'''
    mutater = rosetta.protocols.simple_moves.MutateResidue()
    mutater.set_res_name(res_name)

    for res in residues:
        if keep_g_p and (pose.residue(res).name3() in ['GLY', 'PRO']):
            continue
        mutater.set_target(res)
        mutater.apply(pose)

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

