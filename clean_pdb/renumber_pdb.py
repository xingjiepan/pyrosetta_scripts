#!/usr/bin/env python3
'''Renumber a pdb file. Run as:
    ./renumber_pdb.py input.pdb output.pdb -p params_file -s start_num
'''

from optparse import OptionParser

import pyrosetta
from pyrosetta import rosetta


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-p", "--params_file", dest="params_file", action="store", default=None,
            help="Parameter file.")
    parser.add_option("-s", "--start_num", dest="start_num", action="store", default=1,
            help="Number of the first residue.")

    (options, args) = parser.parse_args()

    ipdb = args[0]
    opdb = args[1]
    start_num = int(options.start_num)

    if options.params_file:
        pyrosetta.init(options='-extra_res_fa {0} -ignore_unrecognized_res true'.format(options.params_file))
    else:
        pyrosetta.init(options='-ignore_unrecognized_res true')

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, ipdb)

    for i in range(1, pose.size() + 1):
        pose.pdb_info().set_resinfo(i, 'A', i + start_num - 1)

        print(i, pose.pdb_info().pose2pdb(i))

    pose.dump_pdb(opdb)
