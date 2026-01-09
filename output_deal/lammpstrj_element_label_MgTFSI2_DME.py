import pandas as pd
from collections import Counter
import os

filepath = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\12.MD_init2_4V/'
filename = filepath + 'NVEe_pl_1_10.lammpstrj'
base_name, ext = os.path.splitext(filename)
n_slab = 2160 #960 #2160
n_dme = 800
n_mg = 20
n_tfsi = 40

name_slab = ['Mg']
name_dme = 'c h h h o c c h h h h o c h h h'
name_dme = name_dme.upper().split()
na_dme = len(name_dme)

name_mg = ['Mg']

name_tfsi = 'n o o o o c c f f f f f f s s'
name_tfsi = name_tfsi.upper().split()
na_tfsi = len(name_tfsi)

n_total = n_slab + n_dme * na_dme + n_tfsi * na_tfsi + n_mg
print(f"预期总原子数: {n_total}")


# 定义函数：处理一行数据并返回处理后的数据列表
def process_line(line, element_name):
    """
    输入一行字符串，将第三个量改成相应的元素名
    返回处理后的数据列表
    """
    parts = line.strip().split()
    if len(parts) >= 3:
        # 替换第三个元素为指定的元素名
        parts[2] = element_name
        return parts
    else:
        # 如果行格式不正确，返回原始行
        return line.strip().split()


# 主程序
frame_count = 0  # 用于统计处理的帧数

with open(filename, 'r') as f:
    output_file_path = base_name + '_element_labeled' + ext
    output_file = open(output_file_path, 'w')

    while True:
        line = f.readline()
        if not line:
            print('处理完成！')
            break

        # 写入当前行到输出文件
        output_file.write(line)

        # 检查是否到达原子数据部分
        if line.strip() == 'ITEM: ATOMS id mol element q x y z':
            frame_count += 1
            print(f"\n正在处理第 {frame_count} 帧...")

            # 创建空的DataFrame来存储原子数据
            atom_data = []

            # 读取slab原子
            for i in range(n_slab):
                line = f.readline()
                processed_data = process_line(line, name_slab[0])
                atom_data.append(processed_data)

            # 读取DME原子
            for i in range(n_dme):
                for j in range(na_dme):
                    line = f.readline()
                    processed_data = process_line(line, name_dme[j])
                    atom_data.append(processed_data)

            # 读取Mg原子
            for i in range(n_mg):
                line = f.readline()
                processed_data = process_line(line, name_mg[0])
                atom_data.append(processed_data)

            # 读取TFSI原子
            for i in range(n_tfsi):
                for j in range(na_tfsi):
                    line = f.readline()
                    processed_data = process_line(line, name_tfsi[j])
                    atom_data.append(processed_data)

            # 创建DataFrame
            df = pd.DataFrame(atom_data)

            # 统计各元素数量
            if len(df.columns) >= 3:
                element_counts = df[2].value_counts().sort_index()
                total_atoms = len(df)

                print(f"第 {frame_count} 帧元素统计:")
                print("-" * 30)
                for element, count in element_counts.items():
                    percentage = (count / total_atoms) * 100
                    print(f"{element}: {count} 个 ({percentage:.2f}%)")
                print(f"总原子数: {total_atoms}")
                print("-" * 30)

                # 验证原子数量是否符合预期
                if total_atoms != n_total:
                    print(f"警告: 实际原子数 ({total_atoms}) 与预期原子数 ({n_total}) 不一致!")

            # 按第三列（元素名）排序
            if len(df.columns) >= 3:
                df_sorted = df.sort_values(by=2)  # 按第三列排序
            else:
                df_sorted = df

            # 重新排序原子ID（从1开始连续编号），但保留原始分子ID
            # 创建新的原子ID列
            new_atom_ids = range(1, len(df_sorted) + 1)

            # 将排序后的DataFrame写入输出文件
            for index, (_, row) in enumerate(df_sorted.iterrows()):
                # 创建新的行数据：新原子ID + 原始分子ID + 其余数据
                new_row_data = [str(new_atom_ids[index])] + list(row[1:])
                line_str = ' '.join(new_row_data) + '\n'
                output_file.write(line_str)

    output_file.close()

    # 最终统计
    print(f"\n{'=' * 50}")
    print(f"处理完成总结:")
    print(f"共处理了 {frame_count} 帧数据")
    print(f"输出文件已保存至: {output_file_path}")
    print(f"{'=' * 50}")