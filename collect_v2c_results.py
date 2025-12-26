import os

filepath = r'E:\Project\57.MgTFSI2_DME_interface\4.pro_43_6\2.pro_43_6_9_Mg1\2.vasp2cp2k\xyz_files/'
filename_pre = 'pro_43_6_9_Mg1-'

# 创建并写入energy文件的第一行
energy_filename = filename_pre + '1.ener'
energy_file = filepath + energy_filename
with open(energy_file, 'w') as f:
    f.write(
        '#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]\n')

# 创建并写入stress文件的第一行
stress_filename = filename_pre + '1.stress'
stress_file = filepath + stress_filename
with open(stress_file, 'w') as f:
    f.write(
        '#   Step   Time [fs]            xx [bar]            xy [bar]            xz [bar]            yx [bar]            yy [bar]            yz [bar]            zx [bar]            zy [bar]            zz [bar]\n')

force_filename = filename_pre + 'frc-1.xyz'
force_file = filepath + force_filename
# 如果force_file存在则删除，然后创建新文件
if os.path.exists(force_file):
    os.remove(force_file)

# 创建新文件
open(force_file, 'w').close()
print(f"已准备好文件: {force_file}")

# 处理多个输出文件
for i in range(0, 6):
    frame_now = i * 10
    outfile_name = filename_pre + f'{frame_now}.out'
    outfile = filepath + outfile_name

    # 检查文件是否存在
    if not os.path.exists(outfile):
        print(f"文件 {outfile} 不存在，跳过")
        continue

    elements = []
    energy_now = None
    force_now = []
    stress_now = []

    try:
        with open(outfile, 'r') as f:
            lines = f.readlines()

        line_index = 0
        while line_index < len(lines):
            line = lines[line_index].strip()

            # 处理Hirshfeld Charges部分
            if 'Hirshfeld Charges' in line:
                # 跳过两行
                line_index += 3
                while line_index < len(lines) and lines[line_index].strip() != '':
                    line_now = lines[line_index].strip().split()
                    if len(line_now) > 1:
                        elements.append(line_now[1])
                    line_index += 1

            # 处理能量
            elif line.startswith('ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]'):
                parts = line.split()
                if parts:
                    energy_now = parts[-1]

            # 处理力
            elif 'FORCES|   Atom     x               y               z               |f|' in line:
                line_index += 1
                while (line_index < len(lines) and
                       not lines[line_index].strip().startswith('FORCES| Sum')):
                    line_now = lines[line_index].strip().split()
                    if len(line_now) >= 5:
                        try:
                            f_x = float(line_now[2])
                            f_y = float(line_now[3])
                            f_z = float(line_now[4])
                            force_now.append([f_x, f_y, f_z])
                        except ValueError:
                            pass
                    line_index += 1

            # 处理应力张量
            elif 'STRESS| Analytical stress tensor [bar]' in line:
                line_index += 1
                for _ in range(4):
                    if line_index < len(lines):
                        line_now = lines[line_index].strip().split()
                        if len(line_now) >= 5:
                            try:
                                # 提取xx, xy, xz分量
                                stress_row = [float(line_now[2]), float(line_now[3]), float(line_now[4])]
                                stress_now.append(stress_row)
                            except ValueError:
                                pass
                        line_index += 1

            line_index += 1

        # 输出当前帧的结果
        print(f"\n处理文件: {outfile_name}")
        print(f"元素列表: {elements}")
        print(f"能量: {energy_now}")
        print(f"力的数量: {len(force_now)}")
        print(f"应力张量: {stress_now}")

        # 这里可以添加代码将数据写入文件
        # 例如，将力写入xyz格式文件
        if force_now:
            with open(force_file, 'a') as f:
                f.write(f"{len(force_now)}\n")
                f.write(f"Frame {frame_now}\n")
                for j, force in enumerate(force_now):
                    if j < len(elements):
                        element = elements[j]
                    else:
                        element = "X"  # 默认元素符号
                    f.write(f"{element} {force[0]:12.6f} {force[1]:12.6f} {force[2]:12.6f}\n")

        # 将能量写入energy文件
        if energy_now is not None:
            with open(energy_file, 'a') as f:
                f.write(
                    f"{frame_now:12d} {0:15.6f} {0:15.6f} {0:15.6f} {float(energy_now):15.6f} {0:15.6f} {0:15.6f}\n")

        # 将应力写入stress文件
        if stress_now and len(stress_now) == 3:
            with open(stress_file, 'a') as f:
                f.write(f"{frame_now:8d} {0:12.6f}")
                for i in range(3):
                    for j in range(3):
                        if i < len(stress_now) and j < len(stress_now[i]):
                            f.write(f" {stress_now[i][j]:20.6f}")
                        else:
                            f.write(" 0.000000")
                f.write("\n")

    except Exception as e:
        print(f"处理文件 {outfile_name} 时出错: {e}")