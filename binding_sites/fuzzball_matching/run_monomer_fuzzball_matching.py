#!/usr/bin/env python2.7
#$ -S /netapp/home/xingjiepan/.local/bin/python2  #-- the shell for the job                                                                                                          
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -l mem_free=3G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=1G,scratch=1G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=48:00:00                #-- runtime limit 

'''
Match a set of fuzz balls to a target monomer.

Usage:
    ./run_jobs.py -f=<FBD> -t=<TPD> -p=<PF> -l=<LI> -o=<OD> -n=<NT> [--run_test]

Options:
    --fuzz_ball_dir=<FBD>, -f=<FBD>
        The directory of fuzz ball structures

    --target_pdb_dir=<TPD>, -t=<TPD> 
        The directory of target structures. The directory should
        also have the .pos files that specifies the matchable positions
        on the target structure.

    --params_file=<PF>, -p=<PF>
        The params file of the ligand.

    --ligand_id=<LI>, -l=<LI>
        The ligand id of the fuzz ball pose

    --output_dir=<OD>, -o=<OD>
        The output directory

    --num_tasks=<NT>, -n=<NT>
        Number of jobs for parallel run.

    --run_test
        Run as a test job.
'''

import os
import time

import docopt
import numpy as np

import pyrosetta
from pyrosetta import rosetta

import preprocessing
import fuzzball_matching
import pick_matches


def match(fuzzball_pdb, target_pdb, pos_file, ligand_id, output_path, min_match_size=3, max_num_outputs=10):
    '''Find matches for a pair of target and fuzz ball.'''
    # Load inputs

    fuzz_pose = preprocessing.load_cleaned_filtered_fuzz_pose(fuzzball_pdb, ligand_id)
    ligand_id = 1 #Ligand id will be 1 after cleaning
    motif_residues = [i for i in range(1, fuzz_pose.size() + 1) if i != ligand_id]
    
    print 'The size of the fuzzball is', fuzz_pose.size()

    target_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(target_pose, target_pdb)
    matchable_positions = preprocessing.load_matchable_positions_from_pos_file(target_pose, pos_file)

    print 'The number of matchable positions is', len(matchable_positions)
    
    bb_compatible_rotamers = preprocessing.get_bb_compatible_rotamers_for_pose(target_pose, matchable_positions)
   
    # Find matches

    print 'Start matching.'

    matches = fuzzball_matching.find_matched_rotamers_for_fuzz_ball(target_pose, matchable_positions, fuzz_pose, ligand_id, motif_residues, bb_compatible_rotamers)

    print 'Found {0} raw matches.'.format(len(matches)) 
    if len(matches) > 0: print 'The max of matched motifs is {0}'.format(max(len(m) for m in matches))

    # Pick matches

    print 'Start picking matches.'

    all_picked_matches = []

    for matches_for_anchor in matches:
       
        # Only construct binding sites from the raw matches with more than min_match_size * 1.5 matches

        if len(matches_for_anchor) >= min_match_size * 1.5: 
            picked_matches = pick_matches.pick_lowest_score_matches_greedy(target_pose, fuzz_pose, matches_for_anchor, ligand_id)
            print len(picked_matches[0]), picked_matches[1]
            
            if len(picked_matches[0]) >= min_match_size:
                all_picked_matches.append(picked_matches)
    
    all_picked_matches = sorted(all_picked_matches, key=lambda x : x[1])

    # Dump the picked matches

    for match_id, picked_matches in enumerate(all_picked_matches[:max_num_outputs]):

        match_output_path = os.path.join(output_path, '{0}_{1}_{2}'.format(match_id, len(picked_matches[0]), int(picked_matches[1])))

        if not os.path.exists(match_output_path):
            os.mkdir(match_output_path)
    
        matched_motif_residues = '_'.join([str(m.target_matched_residue) for m in picked_matches[0]])

        pick_matches.dump_matches_for_an_anchor(target_pose, fuzz_pose, ligand_id, picked_matches[0],
                os.path.join(match_output_path, 'target_pose_M{0}M.pdb'.format(matched_motif_residues)), 
                os.path.join(match_output_path, 'matched_fuzz_pose.pdb'))


if __name__ == '__main__':
    arguments = docopt.docopt(__doc__)

    try:
        task_id = int(os.environ['SGE_TASK_ID']) - 1
    except KeyError:
        task_id = 0

    pyrosetta.init(options='-extra_res_fa {0} -ignore_unrecognized_res true -mute all'.format(arguments['--params_file']))

    # Get the tuples of (fuzzball_pdb, target_pdb, output_path)

    fuzz_ball_dir = arguments['--fuzz_ball_dir']
    target_pdb_dir = arguments['--target_pdb_dir']

    jobs = []

    for ff in os.listdir(fuzz_ball_dir):
        if not ff.endswith('.pdb'): continue

        for tf in os.listdir(target_pdb_dir):
            if not tf.endswith('.pdb'): continue
        
            # Find the position file
            
            pf_found = False
            for pf in os.listdir(target_pdb_dir):
                if pf.endswith('.pos') and pf.startswith(tf[:4]):
                    pf_found = True
                    break

            if not pf_found: continue 

            # Create the output directory

            job_output_dir = os.path.join(arguments['--output_dir'], tf[:-4] + '_' + ff[:-4])

            try:
                os.mkdir(job_output_dir)
            except OSError:
                pass

            jobs.append((os.path.join(fuzz_ball_dir, ff),
                         os.path.join(target_pdb_dir, tf),
                         os.path.join(target_pdb_dir, pf),
                         job_output_dir))

    # Run jobs that belongs to this task

    for i, job in enumerate(jobs):
        if i % int(arguments['--num_tasks']) == task_id:
            print 'Start job {0}/{1}.'.format(i, len(jobs))
            print 'Job defined as', job
            start_time = time.time()
           
            if arguments['--run_test']:
                match(job[0], job[1], job[2], int(arguments['--ligand_id']), job[3], min_match_size=1)
            else:
                match(job[0], job[1], job[2], int(arguments['--ligand_id']), job[3])
            
            end_time = time.time()
            print 'Finish job {0}/{1} in {2} seconds.'.format(i + 1, len(jobs), int(end_time - start_time))
