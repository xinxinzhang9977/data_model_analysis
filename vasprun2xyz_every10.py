from ase.io import read, write
import os

# 读取vasprun.xml文件（所有帧）
filepath = r'E:\Project\57.MgTFSI2_DME_interface\4.pro_43_6\7.pro_43_6_17_Mg3/'

# 确保路径正确
xml_file = os.path.join(filepath, 'vasprun.xml')

# 方法1：使用二进制模式读取
with open(xml_file, 'rb') as f:
    # 读取所有帧
    trajectory = read(f, index=':', format='vasp-xml')

# 方法2：或者指定编码为UTF-8
# with open(xml_file, 'r', encoding='utf-8') as f:
#     trajectory = read(f, index=':', format='vasp-xml')

# 设置文件名前缀
filename_pre = 'pro_43_6_17_Mg3'

# 确保输出目录存在
if not os.path.exists(filepath):
    os.makedirs(filepath)

# 循环每一帧，每隔10帧保存一帧
for i in range(0, len(trajectory), 10):
    # 获取当前帧的原子结构
    frame = trajectory[i]

    # 生成文件名
    filename = f"{filename_pre}-{i}.xyz"

    # 完整文件路径
    full_path = os.path.join(filepath, filename)

    # 将当前帧写入xyz文件
    write(full_path, frame, format='xyz')

    print(f"已保存: {full_path}")

print(f"总共保存了 {len(range(0, len(trajectory), 10))} 帧")