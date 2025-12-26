from ase.io import read, write
from ase.visualize import view
import os
import glob

filepath = r'E:\Project\57.MgTFSI2_DME_interface\7.workflow\train_data/'

# 查找所有符合 pattern 的 vasprun.xml 文件
pattern = os.path.join(filepath, '*_vasprun.xml')
xml_files = glob.glob(pattern)

print(f"找到 {len(xml_files)} 个 vasprun.xml 文件")

for xml_file in xml_files:
    try:
        # 提取文件名（不包含路径和扩展名）
        base_name = os.path.basename(xml_file)
        filename = base_name.replace('_vasprun.xml', '')

        print(f"正在处理: {base_name} -> {filename}.traj")

        # 读取vasprun.xml文件
        with open(xml_file, 'r', encoding='utf-8') as f:
            atoms = read(f, index=':')  # 从文件对象读取

        # 生成输出文件路径
        output_file = os.path.join(filepath, f'{filename}.traj')

        # 保存为xyz格式
        write(output_file, atoms)
        print(f"成功转换: {output_file}")

    except Exception as e:
        print(f"处理文件 {xml_file} 时出错: {e}")

print("所有文件处理完成！")
atoms_from_traj = read(filepath+'interface_cluster_Mg1.traj')
print("晶胞应力:", atoms_from_traj.get_stress())
print("总能:", atoms_from_traj.get_potential_energy())
print(atoms_from_traj.get_forces())
