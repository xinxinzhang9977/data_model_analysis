import os
import numpy as np
import matplotlib.pyplot as plt

# 设置字体为微软雅黑
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'Arial']  # 使用微软雅黑
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 读取数据
data_paths = [
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\6.MD_e_0V/',
    r"E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\2.MD_e_2V/",
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\4.MD_e_4V/',
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\5.MD_e_6V/',
    r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\3.MD_e_8V/'
]

data_file = 'NVEe_1_1.lammpsdump_average_rdf_results.dat'
rdf_data_list = []  # 存储每个文件的(r, g(r))数据

# 读取所有数据文件
for data_path in data_paths:
    file_path = os.path.join(data_path, data_file)

    try:
        # 读取数据文件，假设数据文件有四列，我们取第一列和第四列
        data = np.loadtxt(file_path)
        if data.shape[1] >= 4:
            rdf_data = np.column_stack((data[:, 0], data[:, 3]))  # 第一列和第四列
            rdf_data_list.append(rdf_data)
            print(f"成功读取: {file_path}, 数据形状: {rdf_data.shape}")
        else:
            print(f"警告: {file_path} 列数不足4列，实际列数: {data.shape[1]}")
    except Exception as e:
        print(f"读取文件 {file_path} 时出错: {e}")
        continue

# 检查是否读取到数据
if not rdf_data_list:
    print("错误: 未读取到任何数据文件!")
    exit()

# 将所有数据垂直堆叠
rdf_mg_mg_total = np.vstack(rdf_data_list)

# 设置全局字体大小 - 增大字号
plt.rcParams.update({
    'font.size': 18,           # 增大全局字体大小
    'axes.titlesize': 18,      # 增大标题字体大小
    'axes.labelsize': 18,      # 增大坐标轴标签字体大小
    'xtick.labelsize': 18,     # 增大X轴刻度标签字体大小
    'ytick.labelsize': 18,     # 增大Y轴刻度标签字体大小
    'legend.fontsize': 18,     # 增大图例字体大小
    'axes.linewidth': 2.5,     # 加粗坐标轴线宽（边框）
})

# 创建图形
plt.figure(figsize=(8,6), dpi=100)

# 绘制每个数据集的曲线
for i, (data_path, rdf_data) in enumerate(zip(data_paths, rdf_data_list)):
    # 从路径中提取标签（最后一个文件夹名，按_分割的最后一个值）
    path_parts = os.path.normpath(data_path).split(os.sep)
    folder_name = path_parts[-1] if path_parts else data_path
    label_parts = folder_name.split('_')
    label_now = label_parts[-1] if label_parts else folder_name

    # 绘制曲线，加粗线宽
    plt.plot(rdf_data[:, 0], rdf_data[:, 1]+8*i,
             linewidth=2.5,  # 加粗图线
             label=f'{label_now}' if label_now.replace('V', '').isdigit() else label_now)

# 设置刻度 - 加粗刻度线并增大刻度标签字体
plt.tick_params(
    axis='both',
    which='both',
    direction='in',
    top=True,
    right=True,
    width=2,          # 加粗刻度线宽度
    length=6,         # 增大刻度线长度
    labelsize=18,     # 增大刻度标签字体大小
    pad=8            # 增加刻度标签与刻度线的间距
)

# 获取当前坐标轴并加粗边框
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2.5)  # 加粗坐标轴边框

# 添加标题和轴标签 - 使用更大更粗的字体
plt.title('Mg-Mg RDF at Different Voltages', fontsize=18, fontweight='bold', pad=20)
plt.xlabel('Distance (Å)', fontsize=18, fontweight='bold')
plt.ylabel('g(r)', fontsize=18, fontweight='bold')

# 添加图例 - 使用更大的字体，去除边框
plt.legend(
    frameon=False,  # 移除图例边框
    fontsize=14,
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    borderaxespad=0
)

# 设置网格
plt.grid(True, alpha=0.3, linestyle='--', linewidth=1)

# 调整布局
plt.tight_layout()

# 显示图形
plt.show()

# 可选：保存图形
# plt.savefig('Mg-Mg_RDF_comparison.png', dpi=300, bbox_inches='tight')
# print("图形已保存为 'Mg-Mg_RDF_comparison.png'")