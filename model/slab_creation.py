from ase import Atoms
from ase.build import bulk
import numpy as np
'''
创建Mg slab的pdb文件
'''
# 定义晶格参数
a1 = 3.20940
a2 = 5.55884
a3 = 5.21050

# 定义基元原子位置（分数坐标）
basis_positions = [
    [0.25, 0.3333, 0.25],
    [0.25, 0.6667, 0.75],
    [0.75, 0.8333, 0.25],
    [0.75, 0.1667, 0.75]
]

# 定义超胞尺寸
vXSIZE = 10*a1  # 15 a1
vYSIZE = 6*a2   # 9 a2

# 计算复制次数
n_x = int(round(vXSIZE / a1))
n_y = int(round(vYSIZE / a2))
n_z_bot = int(round(10 / a3))  # -12到-2，厚度10Å
n_z_top = int(round(10 / a3))  # 82到92，厚度10Å

print(f"复制次数: x={n_x}, y={n_y}, z_bot={n_z_bot}, z_top={n_z_top}")

# 创建原胞
# 使用正交晶系，设置晶格向量
cell = np.array([[a1, 0, 0],
                 [0, a2, 0],
                 [0, 0, a3]])

# 将分数坐标转换为笛卡尔坐标
positions = []
for basis in basis_positions:
    cartesian_pos = np.dot(basis, cell)
    positions.append(cartesian_pos)

# 创建原胞原子对象
primitive_cell = Atoms('Mg4', positions=positions, cell=cell, pbc=[True, True, True])

# 创建底部平板 (-12到-2 Å)
bottom_slab = primitive_cell * (n_x, n_y, n_z_bot)
# 平移到底部位置
bottom_slab.translate([0, 0, -12])

# 创建顶部平板 (82到92 Å)
top_slab = primitive_cell * (n_x, n_y, n_z_top)
# 平移到顶部位置
top_slab.translate([0, 0, 182])

# 合并两个平板
final_system = bottom_slab + top_slab

# 设置最终系统的盒子尺寸
final_system.set_cell([[vXSIZE, 0, 0],
                       [0, vYSIZE, 0],
                       [0, 0, 200]])  # z方向从-15到105，总高度120Å

print(f"系统信息:")
print(f"原子总数: {len(final_system)}")
print(f"底部平板原子数: {len(bottom_slab)}")
#print(f"顶部平板原子数: {len(top_slab)}")
print(f"系统尺寸: {final_system.get_cell().lengths()} Å")

# 可选：保存为文件
filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\13.MD_init3_long_5V/'
final_system.write(filepath+'mg_slab_model.xyz')  # 保存为XYZ格式
final_system.write(filepath+'mg_slab_model.pdb')  # 保存为PDB格式

# 可选：可视化（需要安装ase-gui）
from ase.visualize import view
view(final_system)