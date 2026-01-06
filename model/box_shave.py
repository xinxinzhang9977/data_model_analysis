import MDAnalysis as mda
import numpy as np
import os

filepath = r'/4.pro_43_6/'

# 遍历所有pdb文件
for filename in os.listdir(filepath):
    if filename.endswith('.pdb'):
        # 加载体系
        sys = mda.Universe(os.path.join(filepath, filename))

        # 计算新维度和中心
        new_dim = sys.dimensions[0:3] / 2
        new_center = new_dim / 2

        # 平移原子
        cluster_center = sys.select_atoms('all').center_of_geometry()
        sys.atoms.translate(new_center - cluster_center)

        # 更新维度并保存
        sys.dimensions = [*new_dim, 90, 90, 90]
        output_name = filename.replace('.pdb', '_shaved.pdb')
        sys.atoms.write(os.path.join(filepath, output_name))