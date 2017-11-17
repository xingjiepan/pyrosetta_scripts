#!/usr/bin/env python2.7
'''Make ramachandran plot of a structure.
    Usage:
        ./plot_rama input.pdb
'''
import sys

import numpy as np
import matplotlib.pyplot as plt

import pyrosetta
from pyrosetta import rosetta


def plot_rama_score_heat_map(residue_type):
    '''Plot the ramachandran scores of a residue type
    as a heat map.
    '''
    aa = rosetta.core.chemical.aa_from_name(residue_type)
    rama = rosetta.core.scoring.ScoringManager.get_instance().get_Ramachandran()

    scores = []

    for i in range(-18, 18):
        scores.append([])
        for j in range(-18, 18):
            scores[-1].append(rama.eval_rama_score_residue(aa, 10 * i + 5, 10 * j + 5))

    fig, ax = plt.subplots()

    extent = [-180, 180, -180, 180]
    cax = ax.imshow(np.transpose(scores), origin='lower', interpolation='gaussian', extent=extent)
    cbar = fig.colorbar(cax)

    plt.xlabel('phi')
    plt.ylabel('psi')


if __name__ == '__main__':
    pyrosetta.init()

    input_pdb = sys.argv[1]

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, input_pdb)

    phis = [pose.phi(i) for i in range(2, pose.size())]
    psis = [pose.psi(i) for i in range(2, pose.size())]

    
    plot_rama_score_heat_map('VAL')
    plt.scatter(phis, psis)
    plt.axis([-180, 180, -180, 180])
    plt.show()
