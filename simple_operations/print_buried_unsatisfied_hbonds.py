#!/usr/bin/env python3
'''Print buried unsatisfied hydrogen bonds of a protein.
Usage:
    ./print_buried_unsatisfied_hbonds.py pdb_file
'''

import sys
import json

import pyrosetta
from pyrosetta import rosetta


if __name__ == '__main__':
    pdb_file = sys.argv[1]

    #pyrosetta.init(options='-ignore_unrecognized_res true -ignore_waters false')
    pyrosetta.init(options='-ignore_unrecognized_res true')
    
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    # Get the buried unsats

    bupc = rosetta.protocols.simple_pose_metric_calculators.BuriedUnsatisfiedPolarsCalculator(
            'default', 'default')

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)
  
    metric_value = rosetta.basic.MetricValueBase()

    buhs_for_each_res = json.loads(bupc.get('residue_bur_unsat_polars', pose))

    # Print residues with buried unsats

    for i in range(pose.size()):
        if buhs_for_each_res[i] > 0:
            
            print('{0}{1}: {2}'.format(pose.pdb_info().chain(i + 1), 
                pose.pdb_info().number(i + 1), buhs_for_each_res[i]))
