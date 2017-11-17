#!/usr/bin/env python2.7
'''Functions to find the sequence profile of a
structure.
'''
import re

import numpy as np

import pyrosetta
from pyrosetta import rosetta


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
    AAs = ['A', 'P', 'V', 'L', 'I', 'M', 'G', 'F', 'Y', 'W', 
           'S', 'T', 'C', 'H', 'R', 'K', 'D', 'E', 'N', 'Q']
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

if __name__ == '__main__':
    pyrosetta.init()
    
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, 'inputs/input.pdb')

    #vall_path = '/home/xingjie/Softwares/Rosetta/githubRepo/tools/fragment_tools/vall.jul19.2011.gz'
    vall_path = '/home/xingjie/Softwares/Rosetta/githubRepo/main/database/sampling/small.vall.gz'


    get_pose_sequence_profile(pose, vall_path)

    #with open('test.fasta', 'w') as f:
    #    for i in range(min(999, len(sequences))):
    #        f.write('> a\n')
    #        f.write(sequences[i] + '\n')

    #print sequences
