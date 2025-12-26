from MDAnalysis.analysis.dihedrals import Dihedral
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np

# 确保角度在[-180, 180]度范围内
def wrap_to_180(angles):
    """将角度包装到[-180, 180]范围"""
    angles = (angles + 180) % 360 - 180
    return angles

filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\6.MD_e_0V/'
topo_filename = 'NVEe.data'
trj_filenames = []
angles = []
for i in range(1, 6):
    trj_filename = f'NVEe_1_{i}.lammpsdump'
    trj_filenames.append(trj_filename)
for i, trj_filename in enumerate(trj_filenames):
    print(trj_filename)
    # 加载轨迹
    u = mda.Universe(filepath+topo_filename, filepath+trj_filename)
    # 例如：'resid 806 and name C'
    atom_groups = [u.select_atoms('resid 806 808 812 816')] # 请替换为实际的原子选择]
    # 创建Dihedral分析对象
    dihedral_analysis = Dihedral(atom_groups).run()
    # 获取结果
    angles_now = dihedral_analysis.results.angles  # 弧度
    # 转换为角度
    angles_deg = np.degrees(angles_now)
    angles_deg_wrapped = wrap_to_180(angles_deg)
    angles.append(angles_deg_wrapped)
angles = np.array(angles)
# 绘制结果
plt.figure(figsize=(12,10))
for i in range(len(angles)):
    plt.plot(np.arange(0,len(u.trajectory))+i*1000, angles[i], label=f'二面角 {i+1}')
plt.xlabel('帧数')
plt.ylabel('二面角 (°)')
plt.ylim(-180, 180)  # 设置y轴范围
plt.axhline(y=180, color='r', linestyle='--', alpha=0.3)
plt.axhline(y=-180, color='r', linestyle='--', alpha=0.3)
plt.legend()
plt.grid(True)
plt.show()