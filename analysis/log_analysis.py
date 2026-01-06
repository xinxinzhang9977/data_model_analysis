import lammps_logfile
import matplotlib.pyplot as plt
import numpy as np

# 读取数据
data_path = r"E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\7.MD_init2_nvt/"
log = lammps_logfile.File(data_path+'log.ele')
run_num = 0
x = log.get("Step", run_num)
y = log.get("c_qbot", run_num)

# 从第10个数据点开始
x_data = x[1:]
y_data = y[1:]

print(f"平均电荷值: {y_data.mean()}")

# 保存数据到txt文件
data_save_path = "lammps_charge_data.txt"
# 将数据组合成两列：Step 和 c_qbot
data_to_save = np.column_stack((x_data, y_data))
# 保存数据，添加表头
np.savetxt(data_path+data_save_path, data_to_save,
           fmt='%.6f',  # 控制数据精度
           delimiter='\t',  # 制表符分隔
           header='Step\tc_qbot',  # 表头
           comments='')  # 去掉注释符号

print(f"数据已保存到: {data_save_path}")
print(f"数据形状: {data_to_save.shape}")
print(f"数据范围: Step从 {x_data.min()} 到 {x_data.max()}, c_qbot从 {y_data.min():.6f} 到 {y_data.max():.6f}")

# 设置全局绘图参数
plt.rcParams.update({
    'font.size': 14,
    'font.weight': 'bold',
    'axes.linewidth': 2,
    'lines.linewidth': 2,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.major.size': 8,
    'ytick.major.size': 8,
})

# 创建图形和坐标轴
fig, ax = plt.subplots(figsize=(10, 6))

# 绘制数据
ax.plot(x_data, y_data, color='blue', alpha=0.8)

# 设置标题和轴标签
ax.set_title('LAMMPS Simulation: Bottom Charge vs Time Step',
             fontsize=16, fontweight='bold', pad=20)
ax.set_xlabel('Time Step', fontsize=14, fontweight='bold')
ax.set_ylabel('Bottom Charge (c_qbot)', fontsize=14, fontweight='bold')

# 设置刻度参数
ax.tick_params(axis='both', which='major', labelsize=12, length=6, width=1.5)

# 添加网格
ax.grid(True, alpha=0.3, linestyle='--')

# 添加平均值的水平参考线
mean_val = y_data.mean()
ax.axhline(y=mean_val, color='red', linestyle='--', alpha=0.8,
           label=f'Average: {mean_val:.4f}')

# 添加图例
ax.legend(fontsize=12, frameon=True, fancybox=True, shadow=True)

# 调整布局
plt.tight_layout()

# 显示图形
plt.show()

# 可选：保存图片
# plt.savefig('lammps_charge_plot.png', dpi=300, bbox_inches='tight')