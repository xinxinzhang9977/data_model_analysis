import numpy as np
import matplotlib.pyplot as plt

# 请填写您的文件路径和文件名
filepath = r'/1.800_solvents/6.MD_e_0V/'  # 例如: r'C:\data\'
filename = 'NVEe_element_labeled_1_5.lammpstrj_broaden_density_900_1000_10.dat'   # 例如: 'data.txt'

# 读取第一行作为标签
with open(filepath + filename, 'r', encoding='utf-8') as f:
    first_line = f.readline().strip()
    labels = first_line.split('\t')

# 读取数据（跳过第一行）
data = np.loadtxt(filepath + filename, skiprows=1,delimiter='\t')

# 获取数据列数（第一列通常是x轴，所以y轴数据列数为n_data-1）
n_cols = data.shape[1]

# 创建图形和坐标轴
fig, axes = plt.subplots(1, 1, figsize=(12,4))

# 绘制每条曲线（从第2列开始，因为第1列是x轴）
for i in range(1, n_cols-1):
    if i-1 < len(labels):  # 确保标签索引不越界
        axes.plot(data[:, 0], data[:, i], label=labels[i])
    else:
        axes.plot(data[:, 0], data[:, i], label=f'Series {i}')

# 设置图形样式
# 1. 加粗边框
for spine in axes.spines.values():
    spine.set_linewidth(2)

# 2. 增大字体
plt.rcParams.update({'font.size': 18})
axes.tick_params(axis='both', which='major', labelsize=18)
axes.tick_params(axis='both', which='minor', labelsize=18)

# 3. 刻度向内
axes.tick_params(direction='in', which='both', length=6, width=1.5)

# 设置坐标轴标签
if len(labels) > 0:
    axes.set_xlabel(labels[0], fontsize=18, fontweight='bold')
axes.set_ylabel('Y Value', fontsize=18, fontweight='bold')

# 添加图例并设置样式
legend = axes.legend(['Mg','N-TFSI'],loc='best', fontsize=18, frameon=False)
legend.get_frame().set_linewidth(1.5)
legend.get_frame().set_edgecolor('black')

# 添加网格
axes.grid(True, linestyle='--', alpha=0.7)

# 自动调整布局
plt.tight_layout()

# 显示图形
plt.show()