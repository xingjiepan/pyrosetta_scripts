#!/usr/bin/env python2.7

import pyrosetta
from pyrosetta import rosetta

from fuzzball_matching_basic import *


def get_bb_compatible_rotamers_for_pose(pose_original, positions, cutoff_repulsion_score=10):
    '''Get the set of rotamers that are compatible with
    the backbone of a given pose.
    Return:
        A dictionary with the format {position: {'AA':[rotamer_ids]}}
    '''
    def repulsion_score(pose):
        return pose.energies().total_energies()[rosetta.core.scoring.fa_rep] 
    
    AAs = ['VAL', 'LEU', 'ILE', 'MET',
           'PHE', 'TYR', 'TRP',
           'SER', 'THR', 'CYS',
           'ARG', 'LYS', 'HIS',
           'ASP', 'GLU', 'ASN', 'GLN']

    pose = pose_original.clone()
    mutate_residues(pose, range(1, pose.size() + 1), 'ALA')
   
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)
    ref_score = repulsion_score(pose)

    bb_compatible_rotamers = {}

    for seqpos in positions:
        bb_compatible_rotamers[seqpos] = {'ALA':[], 'GLY':[], 'PRO':[]}

        for aa in AAs:
            bb_compatible_rotamers[seqpos][aa] = []
            
            mutate_residues(pose, [seqpos], aa, keep_g_p=False)
            rotamer_set = rosetta.core.pack.rotamer_set.bb_independent_rotamers( pose.residue(seqpos).type(), True )

            for i in range(1, len(rotamer_set) +1):
                replace_intra_residue_torsions(pose, seqpos, rotamer_set[i])
                sfxn(pose)
                score_diff = repulsion_score(pose) - ref_score
                #print seqpos, aa, i, score_diff

                if score_diff < cutoff_repulsion_score:
                    bb_compatible_rotamers[seqpos][aa].append(i)

            mutate_residues(pose, [seqpos], 'ALA', keep_g_p=False)

    return bb_compatible_rotamers

def clean_fuzz_pose(fuzz_pose, ligand_residue):
    '''The fuzz pose might have a lot of residues that are
    nearly identical, remove the redundancy.
    After cleaning, the ligand will be the first residue.
    '''
    cleaned_pose = rosetta.core.pose.Pose(fuzz_pose, ligand_residue, ligand_residue)

    for i in range(1, fuzz_pose.size() + 1):    
        if i == ligand_residue: continue

        # Discard all ALA, GLY and PRO

        if fuzz_pose.residue(i).name3() in ['ALA', 'GLY', 'PRO']: continue

        # Check redundancy

        redundant = False
        for j in range(2, cleaned_pose.size() + 1):
            if fuzz_pose.residue(i).name3() == cleaned_pose.residue(j).name3():
                if sc_heavy_atom_rmsd(fuzz_pose.residue(i), cleaned_pose.residue(j)) < 0.2:
                    redundant = True

        if not redundant:
            cleaned_pose.append_residue_by_jump(fuzz_pose.residue(i).clone(), 1)

    return cleaned_pose

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

def load_cleaned_filtered_fuzz_pose(pdb_file, ligand_id):
    '''Load a fuzz pose from a pdb file. Clean it and
    filter out bad motifs.
    Return the pose.
    '''
    fuzz_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(fuzz_pose, pdb_file)
    
    fuzz_pose = clean_fuzz_pose(fuzz_pose, ligand_id)
    fuzz_pose = filter_motif_residues(fuzz_pose, ligand_id)

    return fuzz_pose

def find_interface_seqposes_noGP(pose, chain1, chain2, cutoff_distance=15):
    '''Return the sequence positions on the interface
    between two chains.
    The selected residues are within a cutoff CA distance and
    the CA-CB vectors are pointing to the other chain.
    GLYs and PROs are ignored
    '''
    chain1_residues = [i for i in range(1, pose.size() + 1) if pose.pdb_info().chain(i) == chain1]
    chain2_residues = [i for i in range(1, pose.size() + 1) if pose.pdb_info().chain(i) == chain2]

    interface_residues = set()

    for res1 in chain1_residues:
        if pose.residue(res1).name3() in ['GLY', 'PRO']: continue 
        ca1 = pose.residue(res1).xyz('CA')
        cb1 = pose.residue(res1).xyz('CB')

        for res2 in chain2_residues:
            if pose.residue(res2).name3() in ['GLY', 'PRO']: continue 
            ca2 = pose.residue(res2).xyz('CA')
            cb2 = pose.residue(res2).xyz('CB')

            if cb1.distance(cb2) > cutoff_distance:
                continue

            if (cb1 - ca1).dot(ca2 - ca1) > 0:
                interface_residues.add(res1)
            
            if (cb2 - ca2).dot(ca1 - ca2) > 0:
                interface_residues.add(res2)

    return list(interface_residues)

def print_pymol_selection_for_residues(pose, residues):
    '''Print the pymol selection command for the residues.'''
    res_commands = ['(c. {0} and res {1})'.format(pose.pdb_info().chain(i), pose.pdb_info().number(i))
                    for i in residues]

    print 'sele ' + ' or '.join(res_commands)

if __name__ == '__main__':
    pyrosetta.init(options='-extra_res_fa test/inputs/REN_no_charge_from_mol2.params')

    #fuzz_pose = load_cleaned_filtered_fuzz_pose('test/inputs/fuzz_balls/REN_0001-single_pose.pdb', 1)
    #print "Number of motif residues =", fuzz_pose.size()
    #fuzz_pose.dump_pdb('test/outputs/cleaned_filtered_fuzz_pose.pdb')

    #target_pose = rosetta.core.pose.Pose()
    #rosetta.core.import_pose.pose_from_file(target_pose, 'test/inputs/3tdn_barrel.pdb')
    #matchable_positions_pdb = [130, 80, 171, 101, 48, 23, 5, 7, 9, 201, 144, 169, 126, 128, 103, 225, 224, 222, 78, 55, 50, 52]
    #matchable_positions = [target_pose.pdb_info().pdb2pose('B', i) for i in matchable_positions_pdb]
    #bb_compatible_rotamers = get_bb_compatible_rotamers_for_pose(target_pose, matchable_positions)
    #print bb_compatible_rotamers

    target_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(target_pose, 'test/inputs/target_pdbs/1svx.pdb')

    matchable_positions = find_interface_seqposes_noGP(target_pose, 'A', 'B')
    print_pymol_selection_for_residues(target_pose, matchable_positions)
    print(len(matchable_positions))

