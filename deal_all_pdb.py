##调用脚本
import os
import subprocess

# 要转换的PDB文件列表

def get_pdb_files(directory="file"):
    """获取指定目录下所有PDB文件"""
    pdb_files = []

    # 检查目录是否存在
    if not os.path.exists(directory):
        print(f"错误：目录 '{directory}' 不存在")
        return pdb_files

    # 获取目录下所有文件
    all_files = os.listdir(directory)

    # 筛选出.pdb文件（不区分大小写）
    pdb_files = [f for f in all_files if f.lower().endswith('.pdb')]

    # 添加完整路径
    pdb_files_with_path = [os.path.join(directory, f) for f in pdb_files]

    return pdb_files_with_path


# 使用示例
pdb_list = get_pdb_files(r"E:\Project\57.MgTFSI2_DME_interface\6.interface\2.interface_small_cluster\1.Mg1/")
print(f"找到 {len(pdb_list)} 个PDB文件:")
for file in pdb_list:
    print(f"  {file}")

# 或者自动获取所有PDB文件
# pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]

for pdb_file in pdb_list:
    if os.path.exists(pdb_file):
        out_name = os.path.basename(pdb_file)
        output_file = f"POSCAR_{out_name.replace('.pdb', '')}"
        print(f"Converting {pdb_file} to {output_file}...")

        # 调用转换脚本
        subprocess.run(["python", "pdb2poscar.py", pdb_file, pdb_list+output_file])