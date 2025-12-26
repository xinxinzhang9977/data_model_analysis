from ase.io import read, write

# 读取vasprun.xml文件（所有帧）
filepath = r'E:\Project\57.MgTFSI2_DME_interface\4.pro_43_6\2.pro_43_6_9_Mg1\1.vasp/'


# 使用二进制模式读取文件
with open(filepath+'vasprun.xml', 'rb') as f:
    trajectory = read(f, index=':', format='vasp-xml')

# 写入xyz文件
write(filepath+'output.xyz', trajectory, format='xyz')