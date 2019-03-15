#!/usr/bin/env python3
'''Insert residues to a structure. Run as:
    ./insert_residues.py input.pdb output.pdb -s sequence
'''

from optparse import OptionParser

import pyrosetta
from pyrosetta import rosetta

def seq_to_aas(sequence):
    '''Convert a sequence of 1 letter code to a list of 3 letter codes'''
    aa_name_map = {'A':'ALA', 'P':'PRO', 'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET',
               'F':'PHE', 'Y':'TYR', 'W':'TRP', 'S':'SER', 'T':'THR', 'C':'CYS',
               'K':'LYS', 'R':'ARG', 'H':'HIS', 'D':'ASP', 'E':'GLU', 'N':'ASN',
               'Q':'GLN', 'G':'GLY'}

    return [aa_name_map[c] for c in sequence]


def insert_aas(pose, position, sequence, insert_after=True, reset_fold_tree=True, fold_tree_root=1):
    '''Insert residues before or after a given position.,
    Set the fold tree to have a cutpoint before or after inserted residues.
    '''
    aas = seq_to_aas(sequence)
    length = len(aas) 

    assert(1 <= position <= pose.size())

    # Set the fold tree with a single cutpoint

    def sub_fold_tree_add_edges_no_jump(ft, root, start, stop):
        '''Add edges to a sub-fold-tree that does not have
        and jumps.'''
        if start < root:
            ft.add_edge(root, start, -1)
        if stop > root:
            ft.add_edge(root, stop, -1)

    if reset_fold_tree:
        cutpoint = position if insert_after else position - 1
        ft = rosetta.core.kinematics.FoldTree()
        
        if fold_tree_root <= cutpoint and cutpoint < pose.size():
            sub_root = pose.size()
            ft.add_edge(fold_tree_root, sub_root, 1)
            sub_fold_tree_add_edges_no_jump(ft, sub_root, cutpoint + 1, pose.size())
            sub_fold_tree_add_edges_no_jump(ft, fold_tree_root, 1, cutpoint)
        
        elif fold_tree_root > cutpoint and cutpoint > 0:
            sub_root = 1
            ft.add_edge(fold_tree_root, sub_root, 1)
            sub_fold_tree_add_edges_no_jump(ft, sub_root, 1, cutpoint)
            sub_fold_tree_add_edges_no_jump(ft, fold_tree_root, cutpoint + 1, pose.size())

        else:
            sub_fold_tree_add_edges_no_jump(ft, fold_tree_root,  1, pose.size())
        
        pose.fold_tree(ft)

    # Append the residues

    residue_type_set = pose.residue_type_set_for_pose()
   
    for i in range(length):
        if insert_after:
            new_rsd = rosetta.core.conformation.ResidueFactory.create_residue( residue_type_set.name_map(aas[i]) )
            pose.conformation().safely_append_polymer_residue_after_seqpos(new_rsd, position + i, True)
            pose.set_omega(position + i, 180)
        else:
            new_rsd = rosetta.core.conformation.ResidueFactory.create_residue( residue_type_set.name_map(aas[len(aas) - i - 1]) )
            pose.conformation().safely_prepend_polymer_residue_before_seqpos(new_rsd, position, True)
            pose.set_omega(position, 180)

    if insert_after:
        rosetta.core.conformation.idealize_position(position + length, pose.conformation())
        
        if position + length + 1 <= pose.size():
            rosetta.core.conformation.idealize_position(position + length + 1, pose.conformation())
    else:
        if position - 1 > 0:
            rosetta.core.conformation.idealize_position(position - 1, pose.conformation())
        rosetta.core.conformation.idealize_position(position, pose.conformation())



if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-n", "--n_term", dest="n_term", action="store_true", default=False,
            help="Insert from the n terminus")
    parser.add_option("-s", "--sequence", dest="sequence", action="store", default=1,
            help="Insert sequence")

    (options, args) = parser.parse_args()

    ipdb = args[0]
    opdb = args[1]

    pyrosetta.init()
    
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, ipdb)

    seq = options.sequence

    if options.n_term: 
        insert_aas(pose, 1, seq, insert_after=False) 
    else:
        insert_aas(pose, pose.size(), seq) 


    pose.dump_pdb(opdb)
