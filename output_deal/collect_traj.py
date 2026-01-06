from ase.io import read, write

# Define the samples to collect into one file
# Here we collect the optimization steps for LiFePO4 structures when 100%, 75% and 0% Lithiated as well as the MD trajectory of these
filepath = r'/7.workflow/train_data/'

traj_to_collect = [
                   'pro_43_6_7_Mg10.traj',
    'pro_43_6_9_Mg1-cp2k.traj',
    'pro_43_6_9_Mg1-v2c.traj',
    'pro_43_6_17_Mg3-cp2k.traj',
    'pro_43_6_17_Mg3-v2c.traj',
    'pro_43_6_21_Mg4.traj',
    'pro_43_6_25_Mg6.traj',
    'pro_43_6_39_Mg8.traj',
    'pro_43_6_44_Mg2-cp2k.traj',
    'pro_43_6_44_Mg2-v2c.traj',

                   ]


# Read the structures and write them to a new file
traj = []
for traj_files in traj_to_collect:
    traj_file = filepath + traj_files
    atoms = read(traj_file,':')
    for atom in atoms:
        traj.append(atom)
write(filepath+'init_dataset.traj',traj)