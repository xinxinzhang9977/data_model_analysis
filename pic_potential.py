import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# 设置全局字体大小
plt.rcParams.update({'font.size': 14})

filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\7.MD_init2_nvt_5V/'
filename1 = 'NVEe_1_1_element_labeled.lammpstrj_potential_0_100_10.dat'
filename2 = 'NVEe_1_5_element_labeled.lammpstrj_potential_900_1000_10.dat'
data1 = np.loadtxt(filepath+filename1,skiprows=1)
data2 = np.loadtxt(filepath+filename2,skiprows=1)

# 创建图形和坐标轴
fig, ax = plt.subplots(figsize=(12, 4))

# 绘制数据
ax.plot(data1[:,0], data1[:,3], '-', label='Potential init', linewidth=2)
ax.plot(data2[:,0], data2[:,3], '-', label='Potential final', linewidth=2)

# 设置刻度朝内
ax.tick_params(axis='both', direction='in', length=8, width=2, labelsize=14,
               top=True, right=True, which='major')

# 加粗边框
for spine in ax.spines.values():
    spine.set_linewidth(2)

# 设置图例，加大字号
ax.legend(fontsize=14, frameon=False, framealpha=1, edgecolor='black')

# 设置坐标轴标签（如果有的话）
ax.set_xlabel('Z(A)', fontsize=16, fontweight='bold')
ax.set_ylabel('Potential(V)', fontsize=16, fontweight='bold')

# 设置标题（如果需要）
ax.set_title('Potential Comparison', fontsize=18, fontweight='bold')

# 调整布局
plt.tight_layout()

# 显示图形
plt.show()