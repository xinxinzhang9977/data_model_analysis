from ase.io import read, write

# Define the samples to collect into one file
# Here we collect the optimization steps for LiFePO4 structures when 100%, 75% and 0% Lithiated as well as the MD trajectory of these
filepath = r'E:\Project\57.MgTFSI2_DME_interface\7.workflow\train_data/'

traj_to_collect = [
                   # 'pro_43_6_21_p1.traj',
                   # 'pro_43_6_37_p1.traj',
                   # 'pro_43_6_37_p2.traj',
                   # 'pro_43_6_41_p1.traj',
                   # 'pro_43_6_41_p2.traj',
                   # 'pro_43_6_44.traj',
                   #'pro_43_6_48_p1.traj',
                   #'pro_43_6_48_p2.traj',
                   #'pro_43_6_7_p1.traj',
                   #'pro_43_6_9.traj',
                   #'interface_cluster_Mg1.traj',
                   #'interface_cluster_Mg2_p1.traj',
                   #'interface_cluster_Mg2_p2.traj',
                   #'interface_solvent_p1.traj',
                   #'interface_solvent_p2.traj',
                   'pro_43_6_21_Mg4.traj'
                   ]


# Read the structures and write them to a new file
traj = []
for traj_files in traj_to_collect:
    traj_file = filepath + traj_files
    atoms = read(traj_file,':')
    for atom in atoms:
        traj.append(atom)
write(filepath+'init_dataset.traj',traj)