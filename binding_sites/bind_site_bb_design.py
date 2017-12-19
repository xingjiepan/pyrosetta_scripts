#!/usr/bin/env python2.7

import numpy as np

import pyrosetta
from pyrosetta import rosetta


def normalize_np_v(v):
    '''Return a normalized numpy vector.'''
    return v / np.linalg.norm(v)

def xyzV_to_np_array(xyz):
    return np.array([xyz.x, xyz.y, xyz.z])

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

def insert_flanking_residues(pose, motif_residues, ligand_residue):
    '''Insert flanking residues for the motif residues.
    For each residue, 3 residues are added to its front and
    back respectively.
    Return the motif residues after addition and the updated
    ligand residue id.
    '''
    motif_residues_sorted = sorted(motif_residues)

    for i in range(len(motif_residues_sorted)):
        if ligand_residue > motif_residues_sorted[i]:
            ligand_residue += 6
        
        # The order is important

        insert_alas(pose, motif_residues_sorted[i], 3, insert_after=True, reset_fold_tree=False)
        insert_alas(pose, motif_residues_sorted[i], 3, insert_after=False, reset_fold_tree=False)

        motif_residues_sorted[i] += 3
        
        for seqpos in range(motif_residues_sorted[i] - 3, motif_residues_sorted[i] + 3):
            pose.conformation().declare_chemical_bond(seqpos, "C", seqpos + 1, "N")

        for j in range(i + 1, len(motif_residues_sorted)):
            motif_residues_sorted[j] += 6

    return motif_residues_sorted, ligand_residue

def set_up_fold_tree(pose, motif_residues, ligand_residue, ligand_anchor_atom, motif_anchor_atoms=None):
    '''Set the fold tree for the system.
    Also add the cutpoint variants to the residues.
    '''
    ft = rosetta.core.kinematics.FoldTree()
  
    if motif_anchor_atoms is None:
        motif_anchor_atoms = []
        for res in motif_residues:
            #motif_anchor_atoms.append(pose.residue(res).atom_name(pose.residue(res).nheavyatoms())) ###This setting fails sometimes for unkonwn reason
            motif_anchor_atoms.append(pose.residue(res).atom_name(pose.residue(res).nheavyatoms() - 2))

    for i, res in enumerate(motif_residues):
        ft.add_edge(ligand_residue, res, i + 1)
        ft.set_jump_atoms(i + 1, ligand_residue, ligand_anchor_atom,
                motif_residues[i], motif_anchor_atoms[i], True)
        ft.add_edge(res, res - 3, -1)
        ft.add_edge(res, res + 3, -1)

    pose.fold_tree(ft)

    for res in motif_residues:
         rosetta.core.pose.remove_variant_type_from_pose_residue(pose,
                rosetta.core.chemical.CUTPOINT_UPPER, res - 3)
         rosetta.core.pose.remove_variant_type_from_pose_residue(pose,
                rosetta.core.chemical.CUTPOINT_LOWER, res + 3)

    ### DEBUG
    #atom_pointer = rosetta.core.id.AtomID_Map_std_shared_ptr_core_kinematics_tree_Atom_t() 
    #rosetta.core.conformation.build_tree( ft, pose.conformation().const_residues(), atom_pointer );
   
    #aid = rosetta.core.id.AtomID(8, 33)
    #atom_pointer[aid].show()

    #for i in [ligand_residue] + motif_residues:
    #    for j in range(1, pose.residue(i).nheavyatoms() + 1):
    #        aid = rosetta.core.id.AtomID(j, i)
    #        print aid
    #        for k in range(atom_pointer[aid].n_children()):
    #            print '\t', atom_pointer[aid].child(k).id()
    
    #exit()
    #### DEBUG

def apply_ideal_helix(pose, start, stop):
    '''Make the residues from start to stop ideal helix.'''
    for i in range(start, stop + 1):
       pose.set_phi(i, -57)
       pose.set_psi(i, -47)
       pose.set_omega(i, 180)

def apply_ideal_strand(pose, start, stop):
    '''Make the residues from start to stop ideal strand.'''
    for i in range(start, stop + 1):
       pose.set_phi(i, -120)
       pose.set_psi(i, 117)
       pose.set_omega(i, 180)

def motif_res_helix_direction(pose, motif_residue):
    '''Return the helix direction of a motif residue.'''
    cos = [xyzV_to_np_array(pose.residue(i).xyz('O') - pose.residue(i).xyz('C'))
            for i in range(motif_residue - 2, motif_residue + 3)]

    return normalize_np_v(sum(cos))

def motif_res_strand_direction(pose, motif_residue):
    '''Return the strand direction of a motif residue.'''
    return normalize_np_v(xyzV_to_np_array(
        pose.residue(motif_residue + 1).xyz('CA') - pose.residue(motif_residue - 1).xyz('CA')))

def check_clashes(pose):
    '''Return true if there are no clashes in the pose.'''
    scale_factor = 0.36

    for i in range(1, pose.size() + 1):
        res_i = pose.residue(i)
        if res_i.is_virtual_residue(): continue
        
        for j in range(i + 1, pose.size() + 1):
            if i == j or i == j + 1 or i == j - 1:
                continue

            res_j = pose.residue(j)
            if res_j.is_virtual_residue(): continue

            for ai in range(1, res_i.nheavyatoms() + 1):
                for aj in range(1, res_j.nheavyatoms() + 1):

                    vi = res_i.xyz(ai)
                    vj = res_j.xyz(aj)
                    ri = res_i.atom_type(ai).lj_radius()
                    rj = res_j.atom_type(aj).lj_radius()

                    if (vi - vj).length_squared() < scale_factor * ((ri + rj) ** 2):
                        return False

    return True

def set_sses(pose, motif_residues, motif_ss_types):
    '''Set secondary structures for motif residues.'''
    for i, res in enumerate(motif_residues):
        if 'H' == motif_ss_types[i]:
            apply_ideal_helix(pose, res - 3, res + 3)
        else:
            apply_ideal_strand(pose, res - 3, res + 3)

def find_parallel_sses(pose, motif_residues, motif_ss_types):
    '''Find Secondary structures that support a given binding site,
    The secondary structures should be approximately parallel to each
    other.
    '''
    # Set the secondary structures 

    set_sses(pose, motif_residues, motif_ss_types)

    # Do combinitorial search

    def next_rotamer_combination_id(combination_id, rotamer_set_sizes):
        '''Change the combination id to the next rotamer combination id.
        Return False if it is the last combination.
        '''
        for i in range(len(combination_id)):
            combination_id[i] = (combination_id[i] + 1) % rotamer_set_sizes[i]
            if combination_id[i] != 0:
                return True

        return False

    rotamer_sets = [rosetta.core.pack.rotamer_set.bb_independent_rotamers( pose.residue(res).type(), True )
                    for res in motif_residues]
    rotamer_set_sizes = [len(rs) for rs in rotamer_sets]
    combination_id = [0] * len(rotamer_set_sizes)

    while True:
        
        # Apply a new set of conformers

        for i, res in enumerate(motif_residues): 
            rot_id = combination_id[i] + 1
            replace_intra_residue_torsions(pose, res, rotamer_sets[i][rot_id])

        # Get the SSE directions
    
        directions = []
        for i, res in enumerate(motif_residues):
            if 'H' == motif_ss_types[i]:
                directions.append(motif_res_helix_direction(pose, res))
            else:
                directions.append(motif_res_strand_direction(pose, res))
        
        # Check if the SSE directions are well aligned 
        
        sse_aligned = True

        for i in range(len(directions)):
            for j in range(i + 1, len(directions)):
                if np.absolute(np.dot(directions[i], directions[j])) < 0.6:
                    sse_aligned = False

        if sse_aligned:
            if check_clashes(pose):
                pose.dump_pdb('debug/test.{0}.{1}.pdb'.format('_'.join(motif_ss_types), '_'.join(str(i + 1) for i in combination_id)))


        if not next_rotamer_combination_id(combination_id, rotamer_set_sizes):
            break

def dump_all_rotamers(pose, residue, ligand_residue):
    '''Dump all rotamers for a residue.'''
    rotamer_set = rosetta.core.pack.rotamer_set.bb_independent_rotamers( pose.residue(residue).type(), True )
    
    seqposes = rosetta.utility.vector1_unsigned_long()
    seqposes.append(ligand_residue)
    for i in range(residue - 3, residue + 4):
        seqposes.append(i)
    
    for i in range(len(rotamer_set)):
        replace_intra_residue_torsions(pose, residue, rotamer_set[i + 1])

        local_pose = rosetta.core.pose.Pose() 
        rosetta.core.pose.pdbslice(local_pose, pose, seqposes)
        local_pose.pdb_info().obsolete(True)
        local_pose.dump_pdb('debug/test.{0}.{1}.pdb'.format(residue, i + 1))


if __name__ == '__main__':
    #pyrosetta.init(options='-extra_res_fa inputs/LG1.params')
    pyrosetta.init(options='-extra_res_fa inputs/REN_no_charge_from_mol2.params')

    pose = rosetta.core.pose.Pose()
    #rosetta.core.import_pose.pose_from_file(pose, 'inputs/ke07_active_site.pdb')
    rosetta.core.import_pose.pose_from_file(pose, 'inputs/binding_site_from_james_renumbered.pdb')
    #rosetta.core.import_pose.pose_from_file(pose, 'inputs/binding_site_from_james_renumbered_small.pdb')
    remove_terminal_variants(pose)

    #motif_residues, ligand_residue = insert_flanking_residues(pose, (1, 2, 3), 4)
    #set_up_fold_tree(pose, motif_residues, ligand_residue, 'C1')
    #find_parallel_sses(pose, motif_residues, ['H', 'E', 'H'])

    motif_residues, ligand_residue = insert_flanking_residues(pose, list(range(2, pose.size() + 1)), 1)
    set_up_fold_tree(pose, motif_residues, ligand_residue, 'C1')
    #find_parallel_sses(pose, motif_residues, ['E', 'H', 'H'])
    set_sses(pose, motif_residues, ['H'] * (pose.size() - 1))

    print pose.fold_tree()

    for res in motif_residues:
        dump_all_rotamers(pose, res, ligand_residue)

    #pose.dump_pdb('test.pdb')
    

