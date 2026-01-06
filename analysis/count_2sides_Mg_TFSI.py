import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

filepath = r'/1.800_solvents/7.MD_init2_nvt_5V/'
topo_file = 'NVEa3.data'
trj_file = 'NVEe_1_5.lammpsdump'
sys = mda.Universe(filepath+topo_file,filepath+trj_file)
zmax = max(sys.atoms.positions[:,2])
zmin = min(sys.atoms.positions[:,2])
zmid = (zmax+zmin)/2
mg_resid_select = sys.select_atoms('resid 801 to 820')
tfsi_resid_select = sys.select_atoms('resid 821 to 860')
for ts in sys.trajectory[::10]:
    print(ts.frame)
    n_mg_bot = 0
    n_mg_top = 0
    n_tfsi_bot = 0
    n_tfsi_top = 0
    for atom in mg_resid_select:
        if atom.position[2] < zmid:
            n_mg_bot += 1
        if atom.position[2] >= zmid:
            n_mg_top += 1
    for atom in tfsi_resid_select:
        if atom.position[2] < zmid:
            n_tfsi_bot += 1
        if atom.position[2] >= zmid:
            n_tfsi_top += 1
    n_tfsi_bot = n_tfsi_bot/15
    n_tfsi_top = n_tfsi_top/15
    print(f'n_mg_bot: {n_mg_bot},n_mg_top:{n_mg_top}, n_tfsi_bot:{n_tfsi_bot}, n_tfsi_top:{n_tfsi_top}')

