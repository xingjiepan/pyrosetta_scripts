#!/usr/bin/env python3
'''Print the layer definition of each residue.
Usage:
    ./print_protein_layers.py pdb_file [--use_sc_neighbors]
'''

from optparse import OptionParser

import pyrosetta
from pyrosetta import rosetta


def get_layers(pose, use_sc_neighbors=False):
    '''Get the layers of a pose'''
    layer_selector = rosetta.core.select.residue_selector.LayerSelector() 
    layer_selector.set_use_sc_neighbors(use_sc_neighbors)
  
    # Select core residues 

    core_residues = []
    layer_selector.set_layers(True, False, False)
    layers = layer_selector.apply(pose)
    for i in range(1, len(layers) + 1):
        if 0 != layers[i]:
            core_residues.append(i)

    # Select boundary residues 

    boundary_residues = []
    layer_selector.set_layers(False, True, False)
    layers = layer_selector.apply(pose)
    for i in range(1, len(layers) + 1):
        if 0 != layers[i]:
            boundary_residues.append(i)

    # Select surface residues 

    surface_residues = []
    layer_selector.set_layers(False, False, True)
    layers = layer_selector.apply(pose)
    for i in range(1, len(layers) + 1):
        if 0 != layers[i]:
            surface_residues.append(i)

    return core_residues, boundary_residues, surface_residues

def get_pymol_script_for_one_layer(layer_name, layer_residues):
    '''Get the pymol script for one layer.'''
    return 'sele ' + layer_name + '_layer, res ' + ' res '.join([str(i) for i in layer_residues])

def get_pymol_scripts_for_layers(pose, use_sc_neighbors=False):
    '''Return the pymol script for selecting the layers.'''
    core_residues, boundary_residues, surface_residues = get_layers(pose, use_sc_neighbors)

    return '\n'.join([get_pymol_script_for_one_layer('core', core_residues),
                      get_pymol_script_for_one_layer('boundary', boundary_residues),
                      get_pymol_script_for_one_layer('surface', surface_residues),
                      'color magenta, core_layer and n. c*',
                      'color green, boundary_layer and n. c*',
                      'color cyan, surface_layer and n. c*'])

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-s", "--use_sc_neighbors", action="store_true", dest="use_sc_neighbors")
  
    (options, args) = parser.parse_args()

    pyrosetta.init()

    pdb_file = args[0]
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
            
    print(get_pymol_scripts_for_layers(pose, use_sc_neighbors=options.use_sc_neighbors))

