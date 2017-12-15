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

def insert_alas(pose, position, length, insert_after=True, reset_fold_tree=True):
    '''Insert a poly-ALA peptide before or after a given position.,
    Set the fold tree to have a cutpoint before or after inserted residues.
    '''
    # Set the fold tree with a single cutpoint

    if reset_fold_tree:
        ft = rosetta.core.kinematics.FoldTree()
        
        if 1 < position < pose.size():
            cutpoint = position if insert_after else position - 1
            ft.add_edge(1, pose.size(), 1)
            ft.add_edge(1, cutpoint, -1)
            ft.add_edge(pose.size(), cutpoint + 1, -1)
        else:
            ft.add_edge(1, pose.size())
        
        pose.fold_tree(ft)

    # Append the residues

    residue_type_set = pose.residue_type_set_for_pose()
    new_rsd = rosetta.core.conformation.ResidueFactory.create_residue( residue_type_set.name_map("ALA") )
   
    for i in range(length):
        if insert_after:
            pose.conformation().safely_append_polymer_residue_after_seqpos(new_rsd, position + i, True)
            pose.set_omega(position + i, 180)
        else:
            pose.conformation().safely_prepend_polymer_residue_before_seqpos(new_rsd, position, True)
            pose.set_omega(position, 180)

    if insert_after and position + length < pose.size():
        rosetta.core.conformation.idealize_position(position + length, pose.conformation())
        rosetta.core.conformation.idealize_position(position + length + 1, pose.conformation())
    elif not insert_after and position > 1:
        rosetta.core.conformation.idealize_position(position - 1, pose.conformation())
        rosetta.core.conformation.idealize_position(position, pose.conformation())

def replace_intra_residue_torsions(pose, seqpos, ref_residue):
    '''Replace the intra residue torsions at seqpos
    by the torsions from the reference residue.
    '''
    ref_chis = ref_residue.chi()
    
    for i in range(1, pose.residue(seqpos).nchi() + 1):
        pose.set_chi(i, seqpos, ref_chis[i])

def insert_flanking_residues(pose, motif_residues):
    '''Insert flanking residues for the motif residues.
    For each residue, 3 residues are added to its front and
    back respectively.
    Return the motif residues after addition.
    '''
    motif_residues_sorted = sorted(motif_residues)

    for i in range(len(motif_residues_sorted)):
        insert_alas(pose, motif_residues_sorted[i], 3, insert_after=True, reset_fold_tree=False)
        insert_alas(pose, motif_residues_sorted[i], 3, insert_after=False, reset_fold_tree=False)

        motif_residues_sorted[i] += 3

        for j in range(i + 1, len(motif_residues_sorted)):
            motif_residues_sorted[j] += 6

    return motif_residues_sorted

if __name__ == '__main__':
    pyrosetta.init(options='-extra_res_fa inputs/LG1.params')

    pose = rosetta.core.pose.Pose()
    #rosetta.core.import_pose.pose_from_file(pose, 'inputs/ke07_active_site_no_substrate.pdb')
    rosetta.core.import_pose.pose_from_file(pose, 'inputs/ke07_active_site.pdb')
    remove_terminal_variants(pose)
    #rosetta.core.pose.correctly_add_cutpoint_variants(pose)

    insert_flanking_residues(pose, (1, 2, 3))

    pose.dump_pdb('test.pdb')
    exit()
    
    ft = rosetta.core.kinematics.FoldTree()
    ft.add_edge(4, 1, 1)
    ft.add_edge(4, 2, 2)
    ft.add_edge(4, 3, 3)
    ft.set_jump_atoms(1, 'C1', 'CZ2')
    ft.set_jump_atoms(2, 'C1', 'OE1')
    ft.set_jump_atoms(3, 'C1', 'CE1')

    pose.fold_tree(ft)

    rotamers = rosetta.core.pack.rotamer_set.bb_independent_rotamers( pose.residue(3).type(), True ) 

    print len(rotamers)

    for i in range(1, len(rotamers) + 1):
        replace_intra_residue_torsions(pose, 3, rotamers[i])

        pose.dump_pdb('test_{0}.pdb'.format(i))
