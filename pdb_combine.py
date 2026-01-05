from ase.io import read, write
from ase import Atoms

'''
把液相和电极的pdb结构合到一起去
电极slab的pdb文件来自slab_creation.py
液相的pdb文件来自packmol
'''
filepath = r'E:\Project\57.MgTFSI2_DME_interface\7.workflow\intial_structure/'
filename1 = 'mg_slab_model.pdb'
filename2 = 'MgTFSI2_DME_30.pdb'
# 读取两个PDB文件
atoms1 = read(filepath+filename1)
atoms2 = read(filepath+filename2)

# 平移第二个结构 (0, 0, 12)
atoms2.translate([0, 0, 0])

# 合并两个结构
combined = atoms1 + atoms2

# 可选：设置晶胞参数（根据实际情况调整）
# combined.set_cell([a, b, c, alpha, beta, gamma])
# combined.set_pbc(True)  # 如果需要周期性边界条件

# 将合并后的结构输出为PDB文件
write(filepath+'combined.pdb', combined)