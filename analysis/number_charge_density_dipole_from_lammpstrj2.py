import MDAnalysis as mda
from MDAnalysis.analysis import lineardensity
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import os

# 设置全局绘图样式
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fontsize'] = 10


def line2float(line):
    """
    把读入的原子坐标行转化成float数据类型
    :param line:
    :return:
    """
    line_now = line.split()
    list_now = []

    for x in line_now:
        try:
            # 尝试将字符串转换为float
            list_now.append(float(x))
        except ValueError:
            # 如果转换失败，保留原始字符串
            list_now.append(x)

    return list_now


def gaussian(x, x0, hw):
    """
    计算高斯函数的值
    x: 自变量数组
    x0: 高斯函数中心
    hw: 高斯函数半宽
    """
    return np.exp(-0.5 * ((x - x0) / hw) ** 2)


class DensityFromLammpstrj():
    def __init__(self, path, trj_file, start=0, end=-1, step=1):
        """

        :param file:
        :param start: 从start帧开始
        :param end: 到end帧结束，-1为最后一帧

        """
        self.path = path
        self.trj_file = self.path + trj_file
        self.start_index = start * LINE_NUMBER_PER_FRAME + 1  # 要算的第一帧的开始行数
        # start end是帧数，从0开始
        if end > 0:
            self.end_index = end * LINE_NUMBER_PER_FRAME + 1  # 要算的最后一帧的开始行数
        else:
            self.end_index = -1
        self.line_number = 1  # 行数从1开始（和txt统一）
        self.step_index = LINE_NUMBER_PER_FRAME * step  # 要算的帧间隔的行数

    def get_density(self, data, bin=0.1):
        z_size = 120 #200 #120
        z_pos = np.arange(-20, z_size, bin)
        num_of_pos = len(z_pos)
        num_of_atom = len(data)
        charge_dens = np.zeros(num_of_pos)
        ion_dict = {'Mg': 0, 'N': 1, 'H': 2}
        ion_dens = np.zeros((num_of_pos, 3))
        data = data.sort_values(by='z')
        index = 0
        for i in range(1, num_of_pos):
            while z_pos[i - 1] < data.iloc[index]['z'] <= z_pos[i]:
                charge_dens[i] += data.iloc[index]['q']
                if data.iloc[index]['element'] in ion_dict.keys() and data.iloc[index]['mol']!=0:
                    ion_dens[i][ion_dict[data.iloc[index]['element']]] += 1
                if index < num_of_atom - 1:
                    index += 1
                else:
                    break
        return z_pos, charge_dens, ion_dens

    def get_dipole(self, data, bin=0.1, polar=False):
        z_size = 120 #200 #120
        e = 1.6e-19
        z_pos = np.arange(-20, z_size, bin)
        num_of_pos = len(z_pos)
        num_of_atom = len(data)
        dipole_dens = np.zeros(num_of_pos)
        z_center = data['z'].sum() / num_of_atom
        data = data.sort_values(by='z')
        index = 0

        if polar:
            pass
        if not polar:
            for i in range(1, num_of_pos):
                while z_pos[i - 1] < data.iloc[index]['z'] <= z_pos[i]:
                    if data.iloc[index]['element'] == 'WL':
                        # 计算
                        dipole_dens[i] += data.iloc[index]['q'] * e * (data.iloc[index]['z'] - z_center) * 1e-10
                    else:
                        pass
                    if index < num_of_atom - 1:
                        index += 1
                    else:
                        break
        return z_pos, dipole_dens

    def get_potential(self):
        x_size = 48.141 #32.094 #48.141
        y_size = 50.03 #33.353 #50.03
        e = 1.6e-19
        epsilon = 8.85e-12
        # 读取电荷密度文件，跳过第一行的列名
        charge_density_ave = np.loadtxt(self.path + CHARGE_NAME, skiprows=1)
        bin = charge_density_ave[1][0] - charge_density_ave[0][0]
        efield = np.zeros_like(charge_density_ave[:, 0])
        potential = np.zeros_like(charge_density_ave[:, 0])
        for i in range(1, len(charge_density_ave)):
            efield[i] = efield[i - 1] + charge_density_ave[i][1] * e / (
                    x_size * y_size * bin *epsilon * 1e-30) * bin * 1e-10

            # potential[i] = potential[i - 1] - (efield[i - 1] * bin + charge_density_ave[i][2] * bin) # charge_density第三列是偶极矩密度
            potential[i] = potential[i - 1] - efield[i - 1] * bin * 1e-10  # charge_density第三列是偶极矩密度
        combined = np.column_stack((charge_density_ave[:, 0:2], efield, potential))

        # 保存电势文件，添加列名
        header = "z_position(A)\tcharge_density\telectric_field(V/m)\telectric_potential(V)"
        np.savetxt(self.path + POTENTIAL_NAME, combined, delimiter='\t', header=header, comments='')

        return charge_density_ave[:, 0], potential

    def broaden_density(self, hw=1):
        # 读取电荷密度文件，跳过第一行的列名
        number_density_ave = np.loadtxt(self.path + CHARGE_NAME, skiprows=1)
        number_density_broaden = np.zeros_like(number_density_ave[:, 3:6])
        num_of_pos = len(number_density_ave)
        for i in range(num_of_pos):
            for j in range(3, 6):
                if number_density_ave[i][j] > 0:
                    line_now = gaussian(number_density_ave[:, 0], number_density_ave[i, 0], hw)
                    number_density_broaden[:, j - 3] += line_now
        combined = np.column_stack((number_density_ave[:, 0], number_density_broaden))

        # 保存展宽密度文件，添加列名
        header = "z_position(A)\tMg_broaden_density\tN_broaden_density\tCl_broaden_density"
        np.savetxt(self.path + BROADEN_NAME, combined, delimiter='\t', header=header, comments='')

        return number_density_ave[:, 0], number_density_broaden

    def get_charge(self):
        columns = ['id', 'mol', 'element', 'q', 'x', 'y', 'z']
        trj_file = open(self.trj_file, 'r')
        end_flag = False
        number_of_frame = 0
        line_now = trj_file.readline()
        charge_density_total = []
        dipole_density_total = []
        ion_density_total = []
        while self.start_index < self.end_index or (
                self.end_index == -1 and not end_flag):  # 如果还没读到最后一帧的起始位置或者没有规定最后位置且没读完
            while self.line_number < self.start_index and not end_flag:
                line_now = trj_file.readline()
                if line_now == '':
                    end_flag = True
                    break
                self.line_number += 1
            if not end_flag:
                number_of_frame += 1  # 记录总共统计了多少帧
                print('number of frame is ' + str(number_of_frame))
                print('start_index is ' + str(self.start_index))
                this_frame_line_number = 1
                this_frame_data_frame = pd.DataFrame(columns=columns)
                while this_frame_line_number <= LINE_NUMBER_PER_FRAME:  # 处理这一帧数据
                    #  从第10行开始是原子信息
                    if this_frame_line_number >= 10:
                        list_now = line2float(line_now)
                        this_line_series = pd.Series(list_now, index=columns)
                        this_frame_data_frame = this_frame_data_frame.append(this_line_series, ignore_index=True)
                    else:  # 跳过前面的信息
                        pass
                    line_now = trj_file.readline()
                    this_frame_line_number += 1
                    self.line_number += 1
                z_pos, charge_density_this_frame, ion_density_this_frame = self.get_density(this_frame_data_frame)
                z_pos, dipole_density_this_frame = self.get_dipole(this_frame_data_frame)
                # plt.plot(z_pos, dipole_density_this_frame)
                # plt.title("frame" + str(number_of_frame))
                # plt.show()
                if len(charge_density_total) == 0:
                    charge_density_total = charge_density_this_frame
                    dipole_density_total = dipole_density_this_frame
                    ion_density_total = ion_density_this_frame
                else:
                    charge_density_total += charge_density_this_frame
                    dipole_density_total += dipole_density_this_frame
                    ion_density_total += ion_density_this_frame
                self.start_index += self.step_index  # 下一帧开始的行数
        charge_density_ave = charge_density_total / number_of_frame
        dipole_density_ave = dipole_density_total / number_of_frame
        ion_density_ave = ion_density_total / number_of_frame
        combined = np.column_stack((z_pos, charge_density_ave, dipole_density_ave, ion_density_ave))

        # 保存电荷密度文件，添加列名
        header = "z_position(A)\tcharge_density\tdipole_density\tMg_density\tN_density\tCl_density"
        np.savetxt(self.path + CHARGE_NAME, combined, delimiter='\t', header=header, comments='')

        trj_file.close()
        return z_pos, ion_density_ave


def create_subplots(z_pos, ion_density_ave, z_potential, potential, z_broaden, number_density_broaden, sys_name):
    """
    创建三个子图的函数
    """
    # 创建图形和子图
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
    fig.suptitle(f'Analysis Results for {sys_name}', fontsize=16, fontweight='bold')

    # 第一个子图：离子密度
    ax1.plot(z_pos, ion_density_ave[:, 0], label='Mg$^{2+}$', linewidth=2, color='blue')
    ax1.plot(z_pos, ion_density_ave[:, 1], label='TFSI$^-$', linewidth=2, color='red')
    #ax1.plot(z_pos, ion_density_ave[:, 2], label='H-DME', linewidth=2, color='orange')
    #ax1.set_xlim([8, 70])
    ax1.set_ylim([0, 5])
    ax1.set_xlabel('Z coordinate (Å)', fontweight='bold')
    ax1.set_ylabel('Ion Density', fontweight='bold')
    ax1.set_title('(a) Ion Density Distribution', fontweight='bold')
    ax1.legend(loc='best', framealpha=0.9)
    ax1.grid(True, alpha=0.3)

    # 设置边框和刻度
    for spine in ax1.spines.values():
        spine.set_linewidth(2)

    # 第二个子图：电势
    ax2.plot(z_potential, potential, label='Electric Potential', linewidth=2, color='green')
    ax2.set_xlabel('Z coordinate (Å)', fontweight='bold')
    ax2.set_ylabel('Potential (V)', fontweight='bold')
    ax2.set_title('(b) Electric Potential Profile', fontweight='bold')
    ax2.legend(loc='best', framealpha=0.9)
    ax2.grid(True, alpha=0.3)

    # 设置边框和刻度
    for spine in ax2.spines.values():
        spine.set_linewidth(2)

    # 第三个子图：展宽密度
    ax3.plot(z_broaden, number_density_broaden[:, 0], label='Mg$^{2+}$', linewidth=2, color='blue')
    ax3.plot(z_broaden, number_density_broaden[:, 1], label='N (TFSI)', linewidth=2, color='red')
    #ax3.plot(z_broaden, number_density_broaden[:, 2], label='H-DME', linewidth=2, color='orange')
    #ax3.set_xlim([8, 70])
    ax3.set_xlabel('Z coordinate (Å)', fontweight='bold')
    ax3.set_ylabel('Broadened Density', fontweight='bold')
    ax3.set_title('(c) Broadened Ion Density', fontweight='bold')
    ax3.legend(loc='best', framealpha=0.9)
    ax3.grid(True, alpha=0.3)

    # 设置边框和刻度
    for spine in ax3.spines.values():
        spine.set_linewidth(2)

    # 调整布局
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)

    # 保存图片
    plt.savefig(f'{sys_name}_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()


def load_existing_data(f_path):
    """加载已存在的数据文件"""
    try:
        # 读取电荷密度文件
        charge_data = np.loadtxt(f_path + CHARGE_NAME, skiprows=1)
        z_pos = charge_data[:, 0]
        ion_density_ave = charge_data[:, 3:6]  # Mg, N, Cl密度

        # 读取电势文件
        potential_data = np.loadtxt(f_path + POTENTIAL_NAME, skiprows=1)
        z_potential = potential_data[:, 0]
        potential = potential_data[:, 3]  # 电势在第四列

        # 读取展宽密度文件
        broaden_data = np.loadtxt(f_path + BROADEN_NAME, skiprows=1)
        z_broaden = broaden_data[:, 0]
        number_density_broaden = broaden_data[:, 1:4]  # 展宽后的Mg, N, Cl密度

        return z_pos, ion_density_ave, z_potential, potential, z_broaden, number_density_broaden

    except Exception as e:
        print(f"读取数据文件时出错: {e}")
        return None


if __name__ == '__main__':
    strat_time = time.time()
    title_dict = {
                  'mgtfsi2_dme_800_5V': r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\7.MD_init2_nvt_5V/',
                  'mgtfsi2_dme_800_2V':r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\11.MD_init2_2V/',
                  'mgtfsi2_dme_800_4V': r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\12.MD_init2_4V/',
                  'mgtfsi2_dme_800_3V': r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\9.MD_init2_3V/',
                  'mgtfsi2_dme_800_1V': r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\10.MD_init2_1V/',
                  'mgtfsi2_dme_800_0V': r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\8.MD_init2_0V/',
                  'mgtfsi2_dme_800_5V_scale_0.5':r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\14.MD_init2_charge_scale_0.5_5V/',
                  'mgtfsi2_dme_800_5V_long': r"E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\13.MD_init3_long_5V/",
                  'mgtfsi2_dme_800_5V_scale_0.9':r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\16.MD_init2_charge_scale_0.9_5V/',
                  'mgtfsi2_dme_800_5V_scale_0.7': r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\15.MD_init2_charge_scale_0.7_5V/',

    }

    sys_name = 'mgtfsi2_dme_800_5V'
    LINE_NUMBER_PER_FRAME = 15589 # 14389 #15589  # p每帧的行数
    LAMMPSTRJ_NAME = 'NVEe_pl_1_10_element_labeled.lammpstrj'

    for start in [900]:
        end = start + 100
        step = 10
        CHARGE_NAME = LAMMPSTRJ_NAME + f'_total_density_{start}_{end}_{step}.dat'
        POTENTIAL_NAME = LAMMPSTRJ_NAME + f'_potential_{start}_{end}_{step}.dat'
        BROADEN_NAME = LAMMPSTRJ_NAME + f'_broaden_density_{start}_{end}_{step}.dat'
        F_DENSITY = False  # 设置为False来测试读取现有数据

        f_path = title_dict[sys_name]

        # 检查文件是否存在
        charge_file_exists = os.path.exists(f_path + CHARGE_NAME)
        potential_file_exists = os.path.exists(f_path + POTENTIAL_NAME)
        broaden_file_exists = os.path.exists(f_path + BROADEN_NAME)

        if F_DENSITY:
            # 计算新数据
            density = DensityFromLammpstrj(f_path, LAMMPSTRJ_NAME, start, end, step)
            z_pos, ion_density_ave = density.get_charge()
            z_potential, potential = density.get_potential()
            z_broaden, number_density_broaden = density.broaden_density(hw=0.5)
        else:
            # 检查文件是否存在
            if charge_file_exists and potential_file_exists and broaden_file_exists:
                print("数据文件存在，读取现有数据...")
                result = load_existing_data(f_path)
                if result is not None:
                    z_pos, ion_density_ave, z_potential, potential, z_broaden, number_density_broaden = result
                else:
                    print("读取数据失败，程序退出")
                    exit()
            else:
                print("数据文件不存在，请检查文件路径或设置F_DENSITY=True重新计算")
                missing_files = []
                if not charge_file_exists:
                    missing_files.append(CHARGE_NAME)
                if not potential_file_exists:
                    missing_files.append(POTENTIAL_NAME)
                if not broaden_file_exists:
                    missing_files.append(BROADEN_NAME)
                print(f"缺失的文件: {missing_files}")
                exit()

        # 创建子图
        create_subplots(z_pos, ion_density_ave, z_potential, potential, z_broaden, number_density_broaden, sys_name)

    end_time = time.time()
    print('Running time is ' + str((end_time - strat_time) / 60) + 'minutes')