import MDAnalysis as mda
import MDAnalysis.analysis.lineardensity as mald
import numpy as np
import matplotlib.pyplot as plt
import time
import os
'''
用MDanalysis统计离子数和溶剂的密度分布
跑的比较慢，跑一个隔100帧的分轨迹大约需要40分钟
'''
t0 = time.perf_counter()  # 计时开始（更高精度）

filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\7.MD_init2_nvt_5V/'
trj_filename = 'NVEe_1_1.lammpsdump'
data_filename = 'NVEa3.data'

sys = mda.Universe(os.path.join(filepath, data_filename),
                   os.path.join(filepath, trj_filename))
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

analysis_frames = min(1000, n_frames)

step = 100  # 每隔多少帧取一次
binsize = 0.2

loop_t0 = time.perf_counter()
n_processed = 0

for i in range(0, analysis_frames, step):
    frame_t0 = time.perf_counter()

    sys.trajectory[i]

    mg_density = mald.LinearDensity(mg_selection, binsize=binsize)
    mg_density.run()
    mg_density_all.append(mg_density.results.z.pos)

    anion_density = mald.LinearDensity(anion_selection, binsize=binsize)
    anion_density.run()
    anion_density_all.append(anion_density.results.z.pos)

    solvent_density = mald.LinearDensity(solvent_selection, binsize=binsize)
    solvent_density.run()
    solvent_density_all.append(solvent_density.results.z.pos)

    n_processed += 1
    frame_t1 = time.perf_counter()
    print(f"已处理第 {i} 帧 | 本次耗时: {frame_t1 - frame_t0:.3f} s")

loop_t1 = time.perf_counter()
print(f"密度计算循环结束：共处理 {n_processed} 帧 | 循环耗时: {loop_t1 - loop_t0:.2f} s")

# 转换为numpy数组
mg_density_all = np.array(mg_density_all)
anion_density_all = np.array(anion_density_all)
solvent_density_all = np.array(solvent_density_all)

# 平均
mg_density_avg = np.mean(mg_density_all, axis=0)
anion_density_avg = np.mean(anion_density_all, axis=0)
solvent_density_avg = np.mean(solvent_density_all, axis=0)

print(f"共处理了 {len(mg_density_all)} 帧数据")

# z轴坐标（和binsize一致）
z_pos = np.arange(len(mg_density_avg)) * binsize

# 绘图
plt.figure(figsize=(10, 6))
plt.plot(z_pos, mg_density_avg, label='Mg', linewidth=2)
plt.plot(z_pos, anion_density_avg / 15, label='An', linewidth=2)
plt.plot(z_pos, solvent_density_avg / 16, label='Sol', linewidth=2)
plt.xlabel('Z position (Å)')
plt.ylabel('Density')
plt.legend()
plt.title(f'Average Linear Density ({len(mg_density_all)} frames, every {step} frames)')
plt.grid(True, alpha=0.3)

# ===== 保存图像（新增） =====
fig_name = os.path.join(filepath, f"{trj_filename}_avg_linear_density.png")
plt.savefig(fig_name, dpi=300, bbox_inches='tight')
print(f"图像已保存: {fig_name}")

#plt.show()

# 保存数据
np.savetxt(os.path.join(filepath, 'mg_density_average.txt'),
           np.column_stack((z_pos, mg_density_avg)))
np.savetxt(os.path.join(filepath, 'anion_density_average.txt'),
           np.column_stack((z_pos, anion_density_avg)))

output_name = os.path.join(filepath, f"{trj_filename}_density_average.txt")
np.savetxt(output_name,
           np.column_stack((z_pos, mg_density_avg, anion_density_avg)),
           header='z_pos(A) mg_density_avg anion_density_avg',
           fmt='%.6f')

# ===== 总运行时间（新增） =====
t1 = time.perf_counter()
print(f"总运行时间: {t1 - t0:.2f} s")
