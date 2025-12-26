import MDAnalysis as mda
import numpy as np
import MDAnalysis.analysis as ma
import matplotlib.pyplot as plt


def get_z_pos(topo, trj, select, step):
    print(trj)
    sys = mda.Universe(topo, trj)
    sel = sys.select_atoms(select)
    print(sel.atoms.indices)
    n_sel = len(sel.atoms)
    z_pos = []
    for i, ts in enumerate(sys.trajectory[::step]):
        # sys.trajectory[i]
        data_now = sel.positions[:, 2]
        z_pos.append(data_now)
    return n_sel, z_pos


filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\3.MD_e_8V/'
topo_filename = 'NVEa3.data'

# 设置全局字体大小
plt.rcParams.update({'font.size': 18})

for ion_now in range(801, 802):
    trj_filenames = []
    z_total = []
    selection = f'index 15093 15260 10602 8087'  # N-TFSI 'type 39' Mg ions 'resid 801 to 820'
    print(selection)
    for i in range(1, 2):
        trj_filename = f'NVEe_1_{i}.lammpsdump'
        trj_filenames.append(trj_filename)
    for trj_filename in trj_filenames:
        n_sel, z_now = get_z_pos(filepath + topo_filename, filepath + trj_filename, selection, 1)
        n_z = len(z_now)
        for i in range(n_z):
            z_total.append(z_now[i])
    indexs = np.arange(len(z_total)) * 2000
    z_total = np.array(z_total)
    np.savetxt(filepath + f'z_total_{selection}.txt', z_total, delimiter='\t')

fig, axes = plt.subplots(1, 1, figsize=(12, 4))

# 加粗边框
for spine in axes.spines.values():
    spine.set_linewidth(2)  # 加粗边框线宽

# 刻度向内
axes.tick_params(direction='in', width=2, length=6)  # 刻度向内，加粗刻度线

# 绘图
for i in range(n_sel):
    axes.plot(indexs / 1000, z_total[:, i])

# 设置标题
axes.set_title(selection, fontsize=18)

# 设置图例，去掉边框
legend = axes.legend(['H-DME1', 'H-DME2', 'F-TFSI1', 'F-TFSI2'],
                     frameon=False,  # 去掉图例边框
                     fontsize=18,
                     bbox_to_anchor=(1.05, 1),
                     loc='upper left',
                     borderaxespad=0)  # 确保图例字体大小一致

# 设置坐标轴标签字体大小
axes.set_xlabel(axes.get_xlabel(), fontsize=18)
axes.set_ylabel(axes.get_ylabel(), fontsize=18)

# 设置刻度标签字体大小
axes.tick_params(axis='both', which='major', labelsize=18)

# 紧凑布局
plt.tight_layout()
plt.show()