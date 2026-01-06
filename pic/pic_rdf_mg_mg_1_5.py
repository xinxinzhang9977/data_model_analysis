import os
import numpy as np
import matplotlib.pyplot as plt

# 设置中文字体（如果需要显示中文）
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 读取数据
data_paths = [
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\6.MD_e_0V/',
    r"E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\2.MD_e_2V/",
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\4.MD_e_4V/',
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\5.MD_e_6V/',
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\3.MD_e_8V/'
]

data_files = ['NVEe_1_1.lammpsdump_average_rdf_results.dat', 'NVEe_1_5.lammpsdump_average_rdf_results.dat']

# 创建图窗和子图
fig, axes = plt.subplots(len(data_paths), 1, figsize=(8, 3 * len(data_paths)), dpi=100)

# 如果只有一个子图，将axes转换为列表以便迭代
if len(data_paths) == 1:
    axes = [axes]

# 设置全局字体大小
plt.rcParams.update({
    'font.size': 13,
    'axes.titlesize': 14,
    'axes.labelsize': 13,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12
})

# 定义线条样式
line_styles = ['-', '-']
colors = ['#1f77b4', '#ff7f0e']  # 两种不同的颜色

# 统一设置y轴范围为[0, 25)
ylim_range = [0, 25]

# 读取所有数据文件并绘制
for i, data_path in enumerate(data_paths):
    ax = axes[i]

    # 设置子图标题
    path_parts = data_path.rstrip('/').split('/')
    title_text = path_parts[-1] if path_parts else data_path
    ax.set_title(title_text, fontsize=14, fontweight='bold', pad=10)

    # 绘制每个数据文件
    for j, data_file in enumerate(data_files):
        file_path = os.path.join(data_path, data_file)

        try:
            # 读取数据文件
            data = np.loadtxt(file_path)

            # 检查数据维度
            if data.shape[1] >= 4:
                # 绘制第一列和第四列的数据
                ax.plot(data[:, 0], data[:, 3],
                        linewidth=2.5,  # 加粗图线
                        linestyle=line_styles[j % len(line_styles)],  # 不同线条样式
                        color=colors[j % len(colors)],  # 不同颜色
                        label=data_file)
            else:
                print(f"警告: {file_path} 列数不足4列，实际列数: {data.shape[1]}")

        except Exception as e:
            print(f"读取文件 {file_path} 时出错: {e}")
            continue

    # 设置刻度向内并加粗
    ax.tick_params(axis='both', which='both', direction='in',
                   top=True, right=True,
                   length=6, width=1.5,  # 加粗刻度线
                   labelsize=12)

    # 加粗子图边框
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)

    # 统一设置y轴范围
    ax.set_ylim(ylim_range)

    # 设置轴标签
    if i == len(data_paths) - 1:  # 最后一个子图显示x轴标签
        ax.set_xlabel('Distance (Å)', fontsize=13, fontweight='bold', labelpad=8)
    ax.set_ylabel('g(r)', fontsize=13, fontweight='bold', labelpad=8)

    # 添加图例，放在左上角
    legend = ax.legend(frameon=True, framealpha=0.9, loc='upper left')
    if legend:
        legend.get_frame().set_linewidth(1.5)

    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--')

    # 设置x轴范围（根据数据自动调整）
    ax.autoscale(enable=True, axis='x', tight=True)

    # 设置刻度密度
    if i == len(data_paths) - 1:
        ax.xaxis.set_major_locator(plt.MaxNLocator(6))

# 调整布局
plt.tight_layout()

# 显示图形
plt.show()

# 可选：保存图形
# fig.savefig('RDF_comparison_all_voltages.png', dpi=300, bbox_inches='tight')
# print("图形已保存为 'RDF_comparison_all_voltages.png'")