#!/bin/bash

../run_interface_fuzzball_matching.py -f inputs/fuzz_balls_small/ -t inputs/target_pdbs_interface/ -p inputs/REN_no_charge_from_mol2.params -l 1 -n 1 -o outputs/ --run_test
