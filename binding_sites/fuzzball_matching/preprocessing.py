#!/usr/bin/env python2.7

import pyrosetta
from pyrosetta import rosetta


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

