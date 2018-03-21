#!/bin/bash

#../run_monomer_fuzzball_matching.py -f inputs/fuzz_balls_small/ -t inputs/target_pdbs_monomer/ -p inputs/REN_no_charge_from_mol2.params -l 1 -n 1 -o outputs/ --run_test

../run_monomer_fuzzball_matching.py -f inputs/fuzz_balls/ -t inputs/target_pdbs_monomer/ -p inputs/REN_no_charge_from_mol2.params -l 1 -n 1 -o outputs/ --run_test
