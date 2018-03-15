#!/usr/bin/env python2.7
'''
Match a set of fuzz balls to a target interface.
Note that the interface should between chain A and chain B
in the target pdb.

Usage:
    ./run_jobs.py -f=<FBD> -t=<TPD> -p=<PF> -l=<LI> -o=<OD> -n=<NT>

Options:
    --fuzz_ball_dir=<FBD>, -f=<FBD>
        The directory of fuzz ball structures

    --target_pdb_dir=<TPD>, -t=<TPD> 
        The directory of target structures.

    --params_file=<PF>, -p=<PF>
        The params file of the ligand.

    --ligand_id=<LI>, -l=<LI>
        The ligand id of the fuzz ball pose

    --output_dir=<OD>, -o=<OD>
        The output directory

    --num_tasks=<NT>, -n=<NT>
        Number of jobs for parallel run.
'''

import os

import docopt

import pyrosetta
from pyrosetta import rosetta

import preprocessing
import fuzzball_matching
import pick_matches


def match(fuzzball_pdb, target_pdb, ligand_id, output_path):
    '''Find matches for a pair of target and fuzz ball.'''

    fuzz_pose = preprocessing.load_cleaned_filtered_fuzz_pose(fuzzball_pdb, ligand_id)
    motif_residues = [i for i in range(1, fuzz_pose.size() + 1) if i != ligand_id]
    
    print 'The size of the fuzzball is', fuzz_pose.size()

    target_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(target_pose, target_pdb)
    matchable_positions = preprocessing.find_interface_seqposes_noGP(target_pose, 'A', 'B', cutoff_distance=10)###DEBUG
    #matchable_positions = preprocessing.find_interface_seqposes_noGP(target_pose, 'A', 'B')
    
    print 'The number of matchable positions is', len(matchable_positions)
    
    bb_compatible_rotamers = preprocessing.get_bb_compatible_rotamers_for_pose(target_pose, matchable_positions)
    
    print 'Start matching.'

    fuzzball_matching.find_matched_rotamers_for_fuzz_ball(target_pose, matchable_positions, fuzz_pose, ligand_id, motif_residues, bb_compatible_rotamers)


if __name__ == '__main__':
    arguments = docopt.docopt(__doc__)

    try:
        task_id = int(os.environ['SGE_TASK_ID']) - 1
    except KeyError:
        task_id = 0

    pyrosetta.init(options='-extra_res_fa {0} -mute all'.format(arguments['--params_file']))

    # Get the tuples of (fuzzball_pdb, target_pdb, output_path)

    fuzz_ball_dir = arguments['--fuzz_ball_dir']
    target_pdb_dir = arguments['--target_pdb_dir']

    jobs = []

    for ff in os.listdir(fuzz_ball_dir):
        if not ff.endswith('.pdb'): continue

        for tf in os.listdir(target_pdb_dir):
            if not tf.endswith('.pdb'): continue
         
            job_output_dir = os.path.join(arguments['--output_dir'], ff[:-4] + '_' + tf[:-4])

            try:
                os.mkdir(job_output_dir)
            except OSError:
                pass

            jobs.append((os.path.join(fuzz_ball_dir, ff),
                         os.path.join(target_pdb_dir, tf),
                         job_output_dir))

    # Run jobs that belongs to this task

    for i, job in enumerate(jobs):
        if i % int(arguments['--num_tasks']) == task_id:
            print 'Start job {0}/{1}.'.format(i, len(jobs))
            match(job[0], job[1], int(arguments['--ligand_id']), job[2])
            print 'Finish job {0}/{1}.'.format(i, len(jobs))
