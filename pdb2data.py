#!/usr/bin/env python3
"""
将 PDB 文件转换为 LAMMPS data 文件的脚本 - 在 masses 注释中显示元素类型
需要安装 ASE: pip install ase
"""

import argparse
from ase.io import read, write
from ase import Atoms
import os
import numpy as np
from ase.data import atomic_masses, atomic_numbers


def pdb_to_data_fixed(pdb_file, data_file, atom_style='full'):
    """
    将 PDB 文件转换为 LAMMPS data 文件 - 在 masses 注释中显示元素类型

    参数:
        pdb_file: 输入的 PDB 文件路径
        data_file: 输出的 LAMMPS data 文件路径
        atom_style: LAMMPS 原子样式 (full, atomic, charge, molecular 等)
    """

    # 检查输入文件是否存在
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB 文件不存在: {pdb_file}")

    # 使用 ASE 读取 PDB 文件
    print(f"正在读取 PDB 文件: {pdb_file}")

    try:
        atoms = read(pdb_file)
    except Exception as e:
        print(f"读取 PDB 文件时出错: {e}")
        # 尝试使用不同的方法读取
        atoms = read(pdb_file, format='pdb')

    # 输出原子系统信息
    print(f"系统包含 {len(atoms)} 个原子")
    print(f"晶胞向量: {atoms.cell}")

    # 手动创建类型映射
    symbols = atoms.get_chemical_symbols()
    unique_symbols = sorted(set(symbols))
    element_to_type = {elem: i + 1 for i, elem in enumerate(unique_symbols)}
    type_to_element = {v: k for k, v in element_to_type.items()}  # 反向映射

    print(f"检测到元素: {list(element_to_type.keys())}")
    print(f"元素到类型的映射: {element_to_type}")

    # 创建类型到质量的映射
    type_to_mass = {}
    for elem, type_id in element_to_type.items():
        atomic_number = atomic_numbers[elem]
        type_to_mass[type_id] = atomic_masses[atomic_number]

    print(f"类型到质量的映射: {type_to_mass}")

    # 创建一个新的 Atoms 对象，确保所有属性正确设置
    new_atoms = Atoms(
        symbols=symbols,
        positions=atoms.positions,
        cell=atoms.cell,
        pbc=atoms.pbc
    )

    # 手动设置原子类型
    types = [element_to_type[sym] for sym in symbols]
    new_atoms.set_array('types', np.array(types, dtype=int))

    # 设置电荷（如果没有）
    if not atoms.has('charges'):
        charges = np.zeros(len(atoms))
        #new_atoms.set_charges(charges)
    else:
        new_atoms.set_charges(atoms.get_charges())

    # 设置质量（基于元素）
    masses = []
    for sym in symbols:
        atomic_number = atomic_numbers[sym]
        masses.append(atomic_masses[atomic_number])
    new_atoms.set_masses(masses)

    # 写入 LAMMPS data 文件 - 使用更明确的参数
    print(f"正在写入 LAMMPS data 文件: {data_file}")

    try:
        # 方法1: 直接写入
        write(data_file, new_atoms, format='lammps-data', atom_style=atom_style)
    except Exception as e:
        print(f"方法1失败: {e}")
        try:
            # 方法2: 使用不同的参数
            write(data_file, new_atoms, format='lammps-data',
                  atom_style=atom_style, units='real', specorder=unique_symbols)
        except Exception as e2:
            print(f"方法2失败: {e2}")
            # 方法3: 使用最简化的方式
            write(data_file, new_atoms, format='lammps-data',
                  atom_style='atomic')  # 使用最简单的atomic样式

    # 如果 ASE 没有正确写入 masses，我们可以手动添加
    if not check_masses_in_data(data_file):
        print("检测到 data 文件缺少 masses 部分，手动添加...")
        add_masses_to_data(data_file, type_to_mass, type_to_element)
    else:
        # 如果 ASE 已经写了 masses，但注释是 type1, type2 等，我们需要替换它们
        replace_masses_comments(data_file, type_to_element)

    print("转换完成!")


def check_masses_in_data(data_file):
    """检查 data 文件是否包含 masses 部分"""
    if not os.path.exists(data_file):
        return False

    with open(data_file, 'r') as f:
        content = f.read()

    # 检查是否包含 "Masses" 部分
    return "Masses" in content


def add_masses_to_data(data_file, type_to_mass, type_to_element):
    """手动将 masses 信息添加到 data 文件，注释中显示元素类型"""
    # 读取原始文件内容
    with open(data_file, 'r') as f:
        lines = f.readlines()

    # 找到 Atoms 部分的位置
    atoms_index = -1
    for i, line in enumerate(lines):
        if "Atoms" in line and "#" not in line:
            atoms_index = i
            break

    if atoms_index == -1:
        print("警告: 无法找到 Atoms 部分，无法添加 masses")
        return

    # 创建新的文件内容
    new_lines = []

    # 添加 masses 部分
    new_lines.append("\nMasses\n\n")
    for type_id, mass in sorted(type_to_mass.items()):
        element = type_to_element.get(type_id, "Unknown")
        new_lines.append(f"{type_id} {mass:.6f}  # {element}\n")

    # 插入 masses 部分到 Atoms 部分之前
    updated_lines = lines[:atoms_index] + new_lines + lines[atoms_index:]

    # 写入更新后的文件
    with open(data_file, 'w') as f:
        f.writelines(updated_lines)

    print("已手动添加 masses 部分到 data 文件，注释中显示元素类型")


def replace_masses_comments(data_file, type_to_element):
    """替换 masses 部分的注释为元素类型"""
    # 读取原始文件内容
    with open(data_file, 'r') as f:
        lines = f.readlines()

    # 找到 Masses 部分的位置
    masses_start = -1
    masses_end = -1
    in_masses_section = False

    for i, line in enumerate(lines):
        if "Masses" in line and "#" not in line:
            masses_start = i
            in_masses_section = True
            continue

        if in_masses_section:
            # 检查是否到达下一个部分或文件末尾
            if line.strip() == "" or (i + 1 < len(lines) and
                                      (lines[i + 1].startswith("\n") or
                                       lines[i + 1].startswith("Atoms") or
                                       lines[i + 1].startswith("Velocities") or
                                       lines[i + 1].startswith("Bonds") or
                                       lines[i + 1].startswith("Angles") or
                                       lines[i + 1].startswith("Dihedrals"))):
                masses_end = i
                break

    if masses_start == -1:
        print("警告: 无法找到 Masses 部分")
        return

    # 替换 Masses 部分的注释
    for i in range(masses_start + 1, masses_end + 1):
        line = lines[i]
        if line.strip() and not line.startswith("#") and not line.startswith("\n"):
            # 提取类型ID
            parts = line.split()
            if len(parts) >= 2:
                try:
                    type_id = int(parts[0])
                    element = type_to_element.get(type_id, "Unknown")
                    # 保留质量值，替换注释
                    mass_value = parts[1]
                    lines[i] = f"{type_id} {mass_value}  # {element}\n"
                except ValueError:
                    # 如果无法解析类型ID，跳过这一行
                    pass

    # 写入更新后的文件
    with open(data_file, 'w') as f:
        f.writelines(lines)

    print("已更新 masses 部分的注释为元素类型")


def pdb_to_data_simple(pdb_file, data_file):
    """
    简化的转换方法，避免复杂的属性设置，但确保包含 masses 和正确的元素注释
    """
    print(f"使用简化方法转换: {pdb_file}")

    # 读取PDB文件
    atoms = read(pdb_file)

    # 创建新的原子系统，只保留基本信息
    new_atoms = Atoms(
        symbols=atoms.get_chemical_symbols(),
        positions=atoms.positions,
        cell=atoms.cell,
        pbc=atoms.pbc
    )

    # 手动添加 masses
    symbols = new_atoms.get_chemical_symbols()
    unique_symbols = sorted(set(symbols))
    element_to_type = {elem: i + 1 for i, elem in enumerate(unique_symbols)}
    type_to_element = {v: k for k, v in element_to_type.items()}  # 反向映射

    # 创建类型到质量的映射
    type_to_mass = {}
    for elem, type_id in element_to_type.items():
        atomic_number = atomic_numbers[elem]
        type_to_mass[type_id] = atomic_masses[atomic_number]

    # 设置原子类型
    types = [element_to_type[sym] for sym in symbols]
    new_atoms.set_array('types', np.array(types, dtype=int))

    # 设置质量
    masses = []
    for sym in symbols:
        atomic_number = atomic_numbers[sym]
        masses.append(atomic_masses[atomic_number])
    new_atoms.set_masses(masses)

    # 写入LAMMPS data文件
    write(data_file, new_atoms, format='lammps-data')

    # 检查并确保 masses 部分存在且注释正确
    if not check_masses_in_data(data_file):
        print("检测到 data 文件缺少 masses 部分，手动添加...")
        add_masses_to_data(data_file, type_to_mass, type_to_element)
    else:
        # 如果 ASE 已经写了 masses，但注释是 type1, type2 等，我们需要替换它们
        replace_masses_comments(data_file, type_to_element)

    print(f"简化转换完成: {data_file}")


def main():
    parser = argparse.ArgumentParser(description='将 PDB 文件转换为 LAMMPS data 文件')
    parser.add_argument('pdb_file', help='输入的 PDB 文件')
    parser.add_argument('-o', '--output', default='data.lammps',
                        help='输出的 LAMMPS data 文件 (默认: data.lammps)')
    parser.add_argument('--atom_style', default='full',
                        choices=['atomic', 'charge', 'molecular', 'full'],
                        help='LAMMPS 原子样式 (默认: full)')
    parser.add_argument('--simple', action='store_true',
                        help='使用简化转换方法')

    args = parser.parse_args()

    try:
        if args.simple:
            pdb_to_data_simple(args.pdb_file, args.output)
        else:
            pdb_to_data_fixed(args.pdb_file, args.output, args.atom_style)
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()