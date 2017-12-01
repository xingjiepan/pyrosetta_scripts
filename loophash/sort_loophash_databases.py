#!/usr/bin/env python2.7

import re
import os
import shutil

import pyrosetta
from pyrosetta import rosetta

import lh_basic


def get_bb_seg_abego(abego_manager, backbone_seg):
    '''Get the ABEGO sequence of a backbone segment.'''
    abego_list = []

    for i in range(backbone_seg.length()):
        abego_list.append(abego_manager.index2symbol(abego_manager.torsion2index_level1(
            backbone_seg.phi()[i], backbone_seg.psi()[i], backbone_seg.omega()[i])))

    return ''.join(abego_list)


def check_loop_abego(lh_library, loop_id, prefix_pattern, loop_pattern, suffix_pattern, abego_manager):
    '''Check if the given loop in a hashmap has the
    desired ABEGO pattern.
    '''
    total_length = len(prefix_pattern) + len(loop_pattern) + len(suffix_pattern)
   
    hashmap = lh_library.gethash(len(loop_pattern))

    # Get the bb segment

    cp = hashmap.get_peptide( loop_id )
    
    bb_data = rosetta.protocols.loophash.BBData()
    extra_data = rosetta.protocols.loophash.BBExtraData()
    lh_library.backbone_database().get_protein( cp.index, bb_data )
    lh_library.backbone_database().get_extra_data(bb_data.extra_key, extra_data)

    protein_size = len(extra_data.sequence)
    seq_offset = cp.offset / 3
    
    if seq_offset <= len(prefix_pattern) or seq_offset + len(loop_pattern) + len(suffix_pattern) >= protein_size:
        return False

    backbone_seg = rosetta.protocols.loophash.BackboneSegment()
    lh_library.backbone_database().get_backbone_segment( cp.index, cp.offset - 3 * len(prefix_pattern), total_length, backbone_seg )

    abego = get_bb_seg_abego(abego_manager, backbone_seg)
    desired_pattern = prefix_pattern + loop_pattern + suffix_pattern

    if re.match(desired_pattern, abego):
        #print abego
        return True
    return False

def save_loops_with_given_abego_pattern_into_a_new_db(lh_library, prefix_pattern, loop_pattern, suffix_pattern, output_file):
    '''Save the loops that have the given abego pattern to a new database.'''
    abego_manager = rosetta.core.sequence.ABEGOManager()
    hashmap = lh_library.gethash(len(loop_pattern))

    selected_ids = []
    for i in range(hashmap.n_loops()):
    #for i in range(1000): ###DEBUG
        if check_loop_abego(lh_library, i, prefix_pattern, loop_pattern, suffix_pattern, abego_manager):
            selected_ids.append(i)
   
    if len(selected_ids) > 0: 
        with open(output_file, 'w') as f:
            for i in selected_ids:
                leap_index = hashmap.get_peptide(i)
                f.write('{0} {1} {2}\n'.format(leap_index.index, leap_index.offset, leap_index.key))

def sort_loophash_db_into_a_new_db(lh_db_path, loop_sizes, pattern_getter, output_path):
    '''Save a loophash database into a new database in which all loops
    have a desired ABEGO pattern.
    Args:
        pattern_getter: a function that convert the length of the loop to a tuple
            of (prefix_pattern, loop_pattern, suffix_pattern)
    '''
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    shutil.copy(os.path.join(lh_db_path, 'backbone.db'), os.path.join(output_path, 'backbone.db'))

    lh_library = lh_basic.load_lh_library(loop_sizes, lh_db_path)

    for loop_size in loop_sizes:
        hashmap = lh_library.gethash(loop_size)
        prefix_pattern, loop_pattern, suffix_pattern = pattern_getter(loop_size)
        
        save_loops_with_given_abego_pattern_into_a_new_db(lh_library, prefix_pattern, 
                loop_pattern, suffix_pattern, os.path.join(output_path, 'loopdb.{0}.db'.format(loop_size)))

def helix_pattern_getter(loop_size):
    '''Loop pattern getter for helices.'''
    return '', '.' + 'A' * (loop_size - 1), ''

def strand_pattern_getter(loop_size):
    '''Loop pattern getter for strands.'''
    return '', '.' + 'B' * (loop_size - 1), ''

def linker_HH_pattern_getter(loop_size):
    '''Loop pattern getter for a linker that connects two helices.'''
    return 'AAA', 'A' + '.' * (loop_size - 1), 'AAAA'

def linker_EE_pattern_getter(loop_size):
    '''Loop pattern getter for a linker that connects two strands.'''
    return 'BB', 'B' + '.' * (loop_size - 1), 'BBB'

def linker_HE_pattern_getter(loop_size):
    '''Loop pattern getter for a linker that connects a helix and a strand.'''
    return 'AAA', 'A' + '.' * (loop_size - 1), 'BBB'

def linker_EH_pattern_getter(loop_size):
    '''Loop pattern getter for a linker that connects a strand and a helix.'''
    return 'BB', 'B' + '.' * (loop_size - 1), 'AAAA'

if __name__ == '__main__':
    pyrosetta.init()

    loop_sizes = list(range(3, 15))
    lh_db_path = "/home/xingjie/DataBases/loophash_db"

    sort_loophash_db_into_a_new_db(lh_db_path, loop_sizes, helix_pattern_getter, 'lh_db_helix')
    #sort_loophash_db_into_a_new_db(lh_db_path, loop_sizes, strand_pattern_getter, 'lh_db_strand')
    #sort_loophash_db_into_a_new_db(lh_db_path, loop_sizes, linkerHH_pattern_getter, 'lh_db_linker_HH')
    #sort_loophash_db_into_a_new_db(lh_db_path, loop_sizes, linkerEE_pattern_getter, 'lh_db_linker_EE')
    #sort_loophash_db_into_a_new_db(lh_db_path, loop_sizes, linkerHE_pattern_getter, 'lh_db_linker_HE')
    #sort_loophash_db_into_a_new_db(lh_db_path, loop_sizes, linkerEH_pattern_getter, 'lh_db_linker_EH')
