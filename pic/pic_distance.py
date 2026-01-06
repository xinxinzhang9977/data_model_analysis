import os
import numpy as np
import matplotlib.pyplot as plt

# 设置全局字体和样式参数
plt.rcParams.update({
    'font.size': 12,  # 字体大小
    'font.weight': 'bold',  # 字体加粗
    'axes.linewidth': 2,  # 坐标轴线宽（边框加粗）
    'axes.labelsize': 14,  # 坐标轴标签字体大小
    'axes.labelweight': 'bold',  # 坐标轴标签加粗
    'axes.titlesize': 16,  # 标题字体大小
    'axes.titleweight': 'bold',  # 标题加粗
    'xtick.direction': 'in',  # x轴刻度朝内
    'ytick.direction': 'in',  # y轴刻度朝内
    'xtick.major.width': 2,  # x轴主刻度线宽
    'ytick.major.width': 2,  # y轴主刻度线宽
    'xtick.major.size': 8,  # x轴主刻度长度
    'ytick.major.size': 8,  # y轴主刻度长度
    'xtick.minor.width': 1.5,  # x轴次刻度线宽
    'ytick.minor.width': 1.5,  # y轴次刻度线宽
    'xtick.minor.size': 4,  # x轴次刻度长度
    'ytick.minor.size': 4,  # y轴次刻度长度
    'legend.fontsize': 12,  # 图例字体大小
    'legend.frameon': True,  # 显示图例边框
    'legend.framealpha': 0.8,  # 图例边框透明度
    'legend.edgecolor': 'black',  # 图例边框颜色
})

filepath = r'/1.800_solvents/6.MD_e_0V/'
filenames = ['atom_distance_resid 806_resid 808.dat',
             'atom_distance_resid 808_resid 812.dat',
             'atom_distance_resid 812_resid 816.dat', ]

# 初始化 data_total 为 None
data_total = None

for filename in filenames:
    data_now = np.loadtxt(filepath + filename, delimiter='\t')

    # 修改：检查 data_total 是否为 None
    if data_total is None:
        data_total = data_now
    else:
        data_total = np.hstack((data_total, data_now))

fig, ax = plt.subplots(figsize=(10, 6), dpi=100)

# 定义颜色和线型
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # 专业配色
line_styles = ['-', '-', '-']
line_widths = [2, 2, 2]

for i in range(0, len(filenames)):
    # 简化图例标签（可选）
    label = filenames[i].replace('atom_distance_resid ', 'Res ').replace('_resid ', '-Res ').replace('.dat', '')

    # 使用面向对象的方式绘图
    ax.plot(data_total[:, i * 2],
            data_total[:, i * 2 + 1],
            label=label,
            color=colors[i % len(colors)],
            linestyle=line_styles[i % len(line_styles)],
            linewidth=line_widths[i % len(line_widths)])

# 设置坐标轴标签（请根据实际数据修改）
ax.set_xlabel('Time (ps)', fontsize=14, fontweight='bold')
ax.set_ylabel('Distance (Å)', fontsize=14, fontweight='bold')
ax.set_title('Interatomic Distance Analysis', fontsize=16, fontweight='bold', pad=15)

# 设置刻度参数
ax.tick_params(axis='both', which='major',
               length=8, width=2,
               direction='in',
               labelsize=12)
ax.tick_params(axis='both', which='minor',
               length=4, width=1.5,
               direction='in')

# 添加网格（可选）
ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

# 添加图例
legend = ax.legend(loc='best', frameon=True,
                   fancybox=True, shadow=True,
                   edgecolor='black', facecolor='white')
legend.get_frame().set_linewidth(1.5)

# 调整布局
plt.tight_layout()

# 保存高清图片（可选）
plt.savefig('distance_analysis.png', dpi=300, bbox_inches='tight', facecolor='white')

plt.show()