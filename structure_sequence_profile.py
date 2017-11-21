#!/usr/bin/env python2.7
'''Functions to find the sequence profile of a
structure.
'''
import re

import numpy as np

import pyrosetta
from pyrosetta import rosetta


AAs = [aa for aa in 'ACDEFGHIKLMNPQRSTVWY']

def get_chunk_abego(chunk, abego_manager):
    '''Get the ABEGO sequence for a vall chunk.'''
    abego_list = []

    for i in range(1, chunk.size() + 1):
        abego_list.append(abego_manager.index2symbol(abego_manager.torsion2index_level1(
            chunk.at(i).phi(), chunk.at(i).psi(), chunk.at(i).omega())))

    return ''.join(abego_list)

def get_pose_abege(pose):
    '''Get the ABEGO sequence for a pose'''
    abego_manager = rosetta.core.sequence.ABEGOManager()
    abego_list = []
    
    for i in range(1, pose.size() + 1):
        abego_list.append(abego_manager.index2symbol(abego_manager.torsion2index_level1(
            pose.phi(i), pose.psi(i), pose.omega(i))))

    # Set the abego at the terminus to be its nearst neighbors ABEGO

    abego_list[0] = abego_list[1]
    abego_list[pose.size() - 1] = abego_list[pose.size() - 2]

    return ''.join(abego_list)

def find_sequences_for_abego_from_chunk(chunk, abego_pattern, abego_manager):
    '''Find sequences corresponding to a string of abego_pattern from a chunk of VALL.'''
    sequences = []

    chunk_sequence = chunk.get_sequence()
    chunk_abegeo = get_chunk_abego(chunk, abego_manager)

    match_positions = [m.start() for m in re.finditer(abego_pattern, chunk_abegeo)]
    
    for p in match_positions:
        sequences.append(chunk_sequence[p : p + len(abego_pattern)])

    return sequences

def find_sequences_for_abego_from_vall(vall_provider, abego_pattern):
    '''Find sequences corresponding to a string of abego_pattern from a VALL.'''
    sequences = []
    
    abego_manager = rosetta.core.sequence.ABEGOManager()

    chunk_ids = np.array(range(1, vall_provider.size() + 1))
    np.random.shuffle(chunk_ids)
    for i in chunk_ids:
        chunk = vall_provider.at(i)
        sequences += find_sequences_for_abego_from_chunk(chunk, abego_pattern, abego_manager)
        if len(sequences) > 1000: # Return if enough sequences are found
            return sequences

    return sequences

def get_pose_sequence_profile(pose, vall_path):
    '''Get the pose sequence profile from the vall database'''
    pose_abego = get_pose_abege(pose)
   
    vall_provider = rosetta.protocols.frag_picker.VallProvider()
    vall_provider.vallChunksFromLibrary(vall_path)

    # Get sequences for each position

    frag_length = 9
    sequences = []

    for i in range(pose.size() + 1 - frag_length):
        sequences.append(find_sequences_for_abego_from_vall(
            vall_provider, pose_abego[i : i + frag_length]))

    # Count the number of AA at each position

    profile = []
    for i in range(pose.size()):
        p = {}
        for aa in AAs:
            p[aa] = 0.1
        profile.append(p)

    for i, ss in enumerate(sequences):
        for s in ss:
            for j in range(frag_length):
                profile[i + j][s[j]] += 1

    # Calculate the AA frequency at each position

    for p in profile:
        s = sum(p[k] for k in p.keys())
        for k in p.keys():
            p[k] = p[k] / s
        
    return profile

def save_profile_to_pssm_file(profile, output_file):
    '''Save a sequence profile to a text file.'''
    with open(output_file, 'w') as f:
        f.write(' '.join(AAs) + '\n')
        for i, d in enumerate(profile):
            best_aa = 'A'
            best_prob = 0
            for k in AAs:
                if d[k] > best_prob:
                    best_aa = k
                    best_prob = d[k]

            frequencies = ['{0:.4f}'.format(d[k]) for k in AAs] 
            f.write('{0} {1} '.format(i + 1, best_aa) + '\t'.join(frequencies) + '\n')

def get_rosetta_profile_from_pssm(pssm_file):
    '''Load a rosetta profile from PSSM'''
    profile = rosetta.core.sequence.SequenceProfile()
    profile.read_from_checkpoint(pssm_file)

def profile_to_energy_profile(profile, kT=1):
    '''Convert the profile of frequencies to the profile
    of energies.
    '''
    energy_profile = []
    
    for p in profile:
        ep = {}

        for aa in AAs:
            ep[aa] = - kT * np.log(p[aa])

        energy_profile.append(ep)
    
    return energy_profile

def get_profile_constraint_list(profile):
    '''Get a list of residue type constraint for the profile.'''
    aa_name_map = {'A':'ALA', 'P':'PRO', 'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET',
                   'F':'PHE', 'Y':'TYR', 'W':'TRP', 'S':'SER', 'T':'THR', 'C':'CYS',
                   'K':'LYS', 'R':'ARG', 'H':'HIS', 'D':'ASP', 'E':'GLU', 'N':'ASN',
                   'Q':'GLN', 'G':'GLY'}
    
    energy_profile = profile_to_energy_profile(profile)

    constraints = []

    for i, ep in enumerate(energy_profile):
        for aa in AAs:
            constraints.append(rosetta.core.scoring.constraints.ResidueTypeConstraint(
                i + 1, aa, aa_name_map[aa], -ep[aa] ))

    return constraints

if __name__ == '__main__':
    pyrosetta.init()
    
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, 'inputs/input.pdb')

    #vall_path = '/home/xingjie/Softwares/Rosetta/githubRepo/tools/fragment_tools/vall.jul19.2011.gz'
    vall_path = '/home/xingjie/Softwares/Rosetta/githubRepo/main/database/sampling/small.vall.gz'

    profile = get_pose_sequence_profile(pose, vall_path)
    save_profile_to_pssm_file(profile, 'pssm.txt')

    get_profile_constraint_list(profile)

    #get_rosetta_profile_from_pssm(pyrosetta.rosetta.utility.file.FileName('pssm.txt'))

    

    #with open('test.fasta', 'w') as f:
    #    for i in range(min(999, len(sequences))):
    #        f.write('> a\n')
    #        f.write(sequences[i] + '\n')

    #print sequences
