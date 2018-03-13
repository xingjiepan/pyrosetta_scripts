load target_pose.pdb
load matched_fuzz_pose.pdb

sele ligand, matched_fuzz_pose and res 1
sele fuzz_scs, matched_fuzz_pose and (not (n. ca+c+o+n)) and (not h.) and (not ligand)
sele matches,target_pose and (not resn ALA) and (not h.)

hide
show cartoon, target_pose
show sticks, ligand
show lines, fuzz_scs
show lines, matches
