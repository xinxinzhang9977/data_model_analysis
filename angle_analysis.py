import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt


# 计算角度（使用自定义函数）
def calculate_angle_from_positions(positions):
    """
    计算三个原子的角度
    输入: 3x3的numpy数组，包含三个原子的坐标
    输出: 角度（度）
    """
    v1 = positions[0] - positions[1]
    v2 = positions[2] - positions[1]

    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    return np.degrees(np.arccos(cos_angle))


filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\3.MD_e_8V/'
topo_filename = 'NVEa3.data'
trj_filenames = []
angles = []
for i in range(1, 6):
    trj_filename = f'NVEe_1_{i}.lammpsdump'
    trj_filenames.append(trj_filename)
for i, trj_filename in enumerate(trj_filenames):
# 加载轨迹
    u = mda.Universe(filepath+topo_filename, filepath+trj_filename)
    print(trj_filename)
    # 选择原子组
    atom_group1 = u.select_atoms('resid 806 808 812')
    atom_group2 = u.select_atoms('resid 808 812 816')
    # 分析轨迹
    for ts in u.trajectory[::10]:
        positions1 = atom_group1.positions
        positions2 = atom_group2.positions
        angle1 = calculate_angle_from_positions(positions1)
        angle2 = calculate_angle_from_positions(positions2)
        angles.append([ts.frame+i*1000,angle1,angle2])
angles = np.array(angles)
np.savetxt(filepath+'angle_analysis.dat', angles, delimiter='\t')
plt.plot(angles[:,0], angles[:,1],label='resid 806 808 812')
plt.plot(angles[:,0], angles[:,2],label='resid 808 812 816')
plt.legend()
plt.show()
print(f"角度统计:")
print(f"范围: {np.min(angles):.2f}° - {np.max(angles):.2f}°")
print(f"平均: {np.mean(angles):.2f}° ± {np.std(angles):.2f}°")