import MDAnalysis as mda
import numpy as np
import MDAnalysis.transformations as mdt
import matplotlib.pyplot as plt

filepath = r'/1.800_solvents/6.MD_e_0V/'
topo_filename = 'NVEe.data'
trj_filenames = []

atoms_center = [801,818,805]
for i in range(1, 6):
    trj_filename = f'NVEe_1_{i}.lammpsdump'
    trj_filenames.append(trj_filename)
for j in range(2):
    dist_total = []
    for i,trj_filename in enumerate(trj_filenames):
        print(trj_filename)
        sys = mda.Universe(filepath+topo_filename, filepath+trj_filename)
        sys_atom = sys.atoms
        transform = mdt.wrap(sys_atom)
        sys.trajectory.add_transformations(transform)
        atom1_sel_str = f'resid {atoms_center[j]}'
        atom2_sel_str = f'resid {atoms_center[j+1]}'
        atom1_sel = sys.select_atoms(atom1_sel_str, updating=True)
        atom2_sel = sys.select_atoms(atom2_sel_str, updating=True)
        for ts in sys.trajectory[::10]:
            pos1 = atom1_sel.positions
            pos2 = atom2_sel.positions
            dis = np.linalg.norm(pos1-pos2)
            dist_total.append([ts.frame+i*1000,dis])
    dist_total = np.array(dist_total)
    np.savetxt(filepath+f'atom_distance_{atom1_sel_str}_{atom2_sel_str}.dat', dist_total, delimiter='\t')
    plt.plot(dist_total[:,0],dist_total[:,1])
plt.show()