import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import os

# 设置全局字体大小
plt.rcParams.update({'font.size': 14})
# 读取数据
filepaths = [
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\8.MD_init2_0V/',
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\10.MD_init2_1V/',
    r"E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\11.MD_init2_2V/",
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\9.MD_init2_3V/',
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\12.MD_init2_4V/',
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\7.MD_init2_nvt_5V/'
]

filename1 = 'NVEe_1_5_element_labeled.lammpstrj_potential_900_1000_10.dat'
#filename2 = 'NVEe_element_labeled_1_5.lammpstrj_potential_900_1000_10.dat'
fig, ax = plt.subplots(figsize=(12, 4))
for filepath in filepaths:
    data_now = np.loadtxt(filepath + filename1, skiprows=1)

    # 创建图形和坐标轴

    folder_name = os.path.basename(os.path.normpath(filepath.rstrip('/\\')))
    label_now = folder_name.split('_')[-1]
    # 绘制数据
    ax.plot(data_now[:,0], data_now[:,3], '-', label=label_now, linewidth=2)


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
#ax.set_title('Potential Comparison', fontsize=18, fontweight='bold')

# 调整布局
plt.tight_layout()

# 显示图形
plt.show()