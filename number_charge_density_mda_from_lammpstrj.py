import MDAnalysis as mda
import MDAnalysis.analysis.lineardensity as mald
import numpy as np
import matplotlib.pyplot as plt
import time

filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\8.MD_init2_0V/'
trj_filename = 'NVEe.lammpsdump'
data_filename = 'NVEa3.data'
sys = mda.Universe(filepath + data_filename, filepath + trj_filename)  # 修正：添加轨迹文件
print(sys.dimensions)

slab_selection_str = 'resid 0'
slab_selection = sys.select_atoms(slab_selection_str)
solvent_select_str = 'resid 1 to 800'
solvent_selection = sys.select_atoms(solvent_select_str)
ions_select_str = 'resid 801 to 860'
ions_selection = sys.select_atoms(ions_select_str)
mg_selection_str = 'resid 801 to 820'
mg_selection = sys.select_atoms(mg_selection_str)
anion_selection_str = 'resid 821 to 860'
anion_selection = sys.select_atoms(anion_selection_str)

# 存储所有帧的结果
mg_density_all = []
anion_density_all = []
solvent_density_all = []

# 获取轨迹总帧数
n_frames = len(sys.trajectory)
print(f"总帧数: {n_frames}")

# 设置要分析的帧范围（前1000帧）
analysis_frames = min(1000, n_frames)

# 每隔10帧计算一次线性密度
for i in range(0, analysis_frames, 100):
    sys.trajectory[i]  # 跳到第i帧

    # 计算Mg的线性密度
    mg_density = mald.LinearDensity(mg_selection, binsize=0.2)
    mg_density.run()
    mg_density_all.append(mg_density.results.z.pos)

    # 计算阴离子的线性密度
    anion_density = mald.LinearDensity(anion_selection, binsize=0.2)
    anion_density.run()
    anion_density_all.append(anion_density.results.z.pos)

    solvent_density = mald.LinearDensity(solvent_selection, binsize=0.2)
    solvent_density.run()
    solvent_density_all.append(solvent_density.results.z.pos)
    print(f"已处理第 {i} 帧")

# 转换为numpy数组以便计算平均值
mg_density_all = np.array(mg_density_all)
anion_density_all = np.array(anion_density_all)
solvent_density_all = np.array(solvent_density_all)

# 计算平均值
mg_density_avg = np.mean(mg_density_all, axis=0)
anion_density_avg = np.mean(anion_density_all, axis=0)
solvent_density_avg = np.mean(solvent_density_all, axis=0)

print(f"共处理了 {len(mg_density_all)} 帧数据")

# 生成z轴坐标（假设所有帧的bins都是一样的）
z_pos = np.arange(len(mg_density_avg)) * 0.2

# 绘制平均密度
plt.figure(figsize=(10, 6))
plt.plot(z_pos, mg_density_avg, label='Mg', linewidth=2)
plt.plot(z_pos, anion_density_avg / 15, label='An', linewidth=2)
plt.plot(z_pos, solvent_density_avg / 16, label='Sol', linewidth=2)
plt.xlabel('Z position (Å)')
plt.ylabel('Density')
plt.legend()
plt.title(f'Average Linear Density ({len(mg_density_all)} frames, every 100 frames)')
plt.grid(True, alpha=0.3)
plt.show()

# 可选：保存平均数据

np.savetxt(filepath+'mg_density_average.txt', np.column_stack((z_pos, mg_density_avg)))
np.savetxt(filepath+'anion_density_average.txt', np.column_stack((z_pos, anion_density_avg)))
output_name = filepath + trj_filename + '_density_average.txt'
np.savetxt(output_name,
           np.column_stack((z_pos, mg_density_avg, anion_density_avg)),
           header='z_pos(A) mg_density_avg anion_density_avg',
           fmt='%.6f')