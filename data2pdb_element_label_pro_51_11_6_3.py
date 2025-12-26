import pandas as pd
from collections import Counter

for i in range(1, 2):
    pdb_name = f"DME10.pdb"
    filepath = r'E:\Project\57.MgTFSI2_DME_interface\6.interface\1.interface_solvent/'
    filename = filepath + pdb_name
    n_slab = 72
    n_dme = 10
    n_mg = 0
    n_tfsi = 0
    name_slab = ['Mg']
    name_dme = 'c h h h o c c h h h h o c h h h'
    name_dme = name_dme.upper().split()
    na_dme = len(name_dme)

    name_mg = ['Mg']
    name_li = ['Li']
    name_tfsi = 'n o o o o c c f f f f f f s s'
    name_tfsi = name_tfsi.upper().split()
    na_tfsi = len(name_tfsi)

    # 读取文件
    with open(filename, 'r') as file:
        lines = file.readlines()

    # 处理文件内容
    i_line = 2  # 从第三行开始（索引为2），重命名变量避免与循环变量冲突
    while i_line < len(lines):
        line = lines[i_line].strip()

        # 跳过空行和非数据行
        if not line or not line.startswith('ATOM') and not line.startswith('HETATM'):
            i_line += 1
            continue

        # 分割行数据，PDB文件使用固定宽度格式
        # 标准PDB格式：原子名称在13-16列，残基序号在23-26列
        if len(line) >= 26:  # 确保行足够长
            try:
                # 提取原子名称和残基序号
                atom_name = line[12:16].strip()
                resid_str = line[22:26].strip()
                resid = int(resid_str) if resid_str else 0
            except (ValueError, IndexError):
                i_line += 1
                continue

            # 根据resid范围修改原子名称
            if resid == 0:  # slab
                new_atom_name = name_slab[0].rjust(3) if len(name_slab[0]) < 4 else name_slab[0]
                # 修改最后一个值（元素符号）
                if len(line) >= 78:  # 确保行有足够的长度
                    # 最后一个值通常在76-78列（元素符号）
                    new_line = line[:12] + new_atom_name.ljust(4) + line[16:76] + name_slab[0].ljust(2) + line[78:]
                else:
                    new_line = line[:12] + new_atom_name.ljust(4) + line[16:]
                lines[i_line] = new_line + '\n'
                i_line += 1

            elif resid <= n_dme:  # DME分子
                for j in range(na_dme):
                    if i_line + j < len(lines):
                        current_line = lines[i_line + j].strip()
                        if current_line and (current_line.startswith('ATOM') or current_line.startswith('HETATM')):
                            # 修改原子名称（第13-16列）
                            # 保持PDB格式，原子名称通常右对齐
                            new_atom_name = name_dme[j].rjust(3) if len(name_dme[j]) < 4 else name_dme[j]
                            # 修改最后一个值（元素符号）
                            if len(current_line) >= 78:
                                new_line = current_line[:12] + new_atom_name.ljust(4) + current_line[16:76] + name_dme[
                                    j].ljust(2) + current_line[78:]
                            else:
                                new_line = current_line[:12] + new_atom_name.ljust(4) + current_line[16:]
                            lines[i_line + j] = new_line + '\n'
                i_line += na_dme

            elif resid <= n_dme + n_mg:  # mg原子
                # 修改原子名称
                new_atom_name = name_mg[0].rjust(3) if len(name_mg[0]) < 4 else name_mg[0]
                # 修改最后一个值（元素符号）
                if len(line) >= 78:
                    new_line = line[:12] + new_atom_name.ljust(4) + line[16:76] + name_mg[0].ljust(2) + line[78:]
                else:
                    new_line = line[:12] + new_atom_name.ljust(4) + line[16:]
                lines[i_line] = new_line + '\n'
                i_line += 1

            elif resid <= n_dme + n_mg + n_tfsi:  # TFSI分子
                for j in range(na_tfsi):
                    if i_line + j < len(lines):
                        current_line = lines[i_line + j].strip()
                        if current_line and (current_line.startswith('ATOM') or current_line.startswith('HETATM')):
                            # 修改原子名称
                            new_atom_name = name_tfsi[j].rjust(3) if len(name_tfsi[j]) < 4 else name_tfsi[j]
                            # 修改最后一个值（元素符号）
                            if len(current_line) >= 78:
                                new_line = current_line[:12] + new_atom_name.ljust(4) + current_line[16:76] + name_tfsi[
                                    j].ljust(2) + current_line[78:]
                            else:
                                new_line = current_line[:12] + new_atom_name.ljust(4) + current_line[16:]
                            lines[i_line + j] = new_line + '\n'
                i_line += na_tfsi

            else:
                i_line += 1
        else:
            i_line += 1

    # 保存修改后的文件
    output_filename = filename.replace('.pdb', '_modified.pdb')
    with open(output_filename, 'w') as file:
        file.writelines(lines)

    print(f"文件已保存为: {output_filename}")
    print(f"处理了 {len(lines)} 行数据")