from ase.io import read, write

# 读取vasprun.xml文件（所有帧）
filepath = r'/4.pro_43_6/7.pro_43_6_17_Mg3/'


# 使用二进制模式读取文件
with open(filepath+'vasprun.xml', 'rb') as f:
    trajectory = read(f, index=':', format='vasp-xml')
trajectory_10 = trajectory[::10]
# 写入xyz文件
write(filepath+'pro_43_6_17_Mg3-pos-1.xyz', trajectory_10, format='xyz')