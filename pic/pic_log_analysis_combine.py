import lammps_logfile
import matplotlib.pyplot as plt
import numpy as np

# 读取数据
data_paths = [r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\8.MD_init2_0V/',
              r"E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\10.MD_init2_1V/",
              r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\11.MD_init2_2V/',
              r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\9.MD_init2_3V/',
              r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\12.MD_init2_4V/',
              r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\7.MD_init2_nvt_5V/'


              ]
charge_total = []
head = 'index\t'
for data_path in data_paths:
    path_last = data_path.split('\\')[-1]
    head = head+path_last+'\t'
    log = lammps_logfile.File(data_path+'log.ele')
    print('Dealing with '+path_last)
    run_num = 0
    x = log.get("Step", run_num)
    y = log.get("c_qbot", run_num)

    # 从第10个数据点开始
    x_data = x[1:]
    y_data = y[1:]

    print(f"平均电荷值: {y_data.mean()}")
    if charge_total==[]:
        charge_total = np.vstack((x_data, y_data))
    else:
        charge_total = np.vstack((charge_total, y_data))
charge_total = charge_total.T


# 保存数据到txt文件
data_save_path = r"E:\Project\57.MgTFSI2_DME_interface\1.800_solvents/"+"charge_data_total.txt"

# 保存数据，添加表头
np.savetxt(data_save_path, charge_total,
           fmt='%.6f',  # 控制数据精度
           delimiter='\t',  # 制表符分隔
           header='Step\t'+head,  # 表头
           comments='')  # 去掉注释符号

print(f"数据已保存到: {data_save_path}")
print(f"数据形状: {charge_total.shape}")
print(f"数据范围: Step从 {x_data.min()} 到 {x_data.max()}, c_qbot从 {y_data.min():.6f} 到 {y_data.max():.6f}")

# 设置全局绘图参数
plt.rcParams.update({
    'font.size': 18,
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
fig, ax = plt.subplots(figsize=(8,6))

for i in range(1,charge_total.shape[1]):
    # 绘制数据
    ax.plot(charge_total[:,0]/1000, charge_total[:,i], alpha=0.8)

# 设置标题和轴标签
ax.set_title('LAMMPS Simulation: Bottom Charge vs Time Step',
             fontsize=16, fontweight='bold', pad=20)
ax.set_xlabel('Time(ps)', fontsize=14, fontweight='bold')
ax.set_ylabel('Bottom Slab Charge(q)', fontsize=18, fontweight='bold')

# 设置刻度参数
ax.tick_params(axis='both', which='major', labelsize=16, length=6, width=1.5)

# 添加网格
ax.grid(True, alpha=0.3, linestyle='--')

# 添加平均值的水平参考线
# mean_val = y_data.mean()
# ax.axhline(y=mean_val, color='red', linestyle='--', alpha=0.8,
#            label=f'Average: {mean_val:.4f}')

# 添加图例
ax.legend(['0V','1V','2V','3V','4V','5V'],fontsize=18, frameon=False,bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

# 调整布局
plt.tight_layout()

# 显示图形
plt.show()

# 可选：保存图片
# plt.savefig('lammps_charge_plot.png', dpi=300, bbox_inches='tight')