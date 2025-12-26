import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def find_name(resid):
    if resid == 0:
        return "slab"
    if 0<resid<=800:
        return 'DME'
    if 801<=resid<=820:
        return 'Mg'
    if 821<=resid<=860:
        return 'TFSI'

filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\2.MD_e_2V/'
topo_filename = 'NVEa3.data'

dict = {'TFSI':1,"DME":2,"Mg":3,'slab':4}
for ion_now in range(801,821):
    selection = fr'resid {ion_now}' # Mg
    print(selection)
    trj_filenames = []
    solvation_total = []
    for i in range(1,6):
        trj_filename = f'NVEe_1_{i}.lammpsdump'
        trj_filenames.append(trj_filename)
    for i, trj_filename in enumerate(trj_filenames):
        sys = mda.Universe(filepath+topo_filename,filepath+trj_filename)
        print(trj_filename)
        center_atom = sys.select_atoms(selection)
        for ts in sys.trajectory[::10]:
            solvation_sel_str = 'around 2.8 group center_atom'
            solvation_sel = sys.select_atoms(solvation_sel_str,center_atom=center_atom,updating=True)
            solvation_now = [ts.frame+i*1000,0,0,0,0]
            for solvation_atom in solvation_sel.atoms:
                solvation_now[dict[find_name(solvation_atom.resid)]] += 1
            solvation_total.append(solvation_now)
    solvation_total = np.array(solvation_total)
    headline = 'timestep\t' + '\t'.join(dict.keys())
    np.savetxt(filepath+f'solvation_total_{selection}.dat',solvation_total,header=headline,delimiter='\t')
plt.plot(solvation_total[:,0],solvation_total[:,1],label='TFSI')
plt.plot(solvation_total[:,0],solvation_total[:,2],label='DME')
plt.plot(solvation_total[:,0],solvation_total[:,3],label='Mg')
plt.plot(solvation_total[:,0],solvation_total[:,4],label='Slab')
plt.legend()
plt.show()





