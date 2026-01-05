from ase import Atoms
import ase.io
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator


def parse_cp2k_ener_file(ener_file):
    """解析CP2K的.ener文件"""

    with open(ener_file, 'r') as f:
        lines = f.readlines()

    # 找到数据开始的行
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip() and not line.startswith('#'):
            data_start = i
            break

    # 解析列标题
    headers = lines[data_start-1].strip().split()
    headers = ['Step Nr.','Time[fs]','Kin.[a.u.]','Temp[K]','Pot.[a.u.]','Cons Qty[a.u.]','UsedTime[s]']
    # 提取数据
    data = {header: [] for header in headers}

    for line in lines[data_start::10]:
        if line.strip() and not line.startswith('#'):
            values = line.split()
            for header, value in zip(headers, values):
                data[header].append(float(value))

    return data

def parse_cp2k_stress_file(stress_file):
    """解析CP2K的.ener文件"""

    with open(stress_file, 'r') as f:
        lines = f.readlines()

    # 找到数据开始的行
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip() and not line.startswith('#'):
            data_start = i
            break

    # 解析列标题
    headers = lines[data_start-1].strip().split()
    headers = ['Step','Time [fs]','xx [bar]','xy [bar]','xz [bar]','yx [bar]','yy [bar]','yz [bar]','zx [bar]','zy [bar]','zz [bar]']
    # 提取数据
    data = {header: [] for header in headers}

    for line in lines[data_start:]:
        if line.strip() and not line.startswith('#'):
            values = line.split()
            for header, value in zip(headers, values):
                data[header].append(float(value))

    return data

def create_aimd_trajectory(pos_file, frc_file=None, ener_file=None,stress_file=None,
                           cell=None, pbc=True):
    """
    从CP2K输出创建轨迹

    参数:
    pos_pattern: 坐标文件模式，如 'pos-*.xyz'
    frc_pattern: 受力文件模式，如 'frc-*.xyz'
    ener_file: 能量文件路径
    cell: 晶胞向量（如果有）
    pbc: 周期性边界条件
    """

    # 获取所有坐标文件
    atoms = ase.io.read(pos_file, format='xyz', index=':')
    ase.io.write(filepath+'cp2k_aimd_xyz.traj', atoms,format='traj')
    forces = ase.io.read(frc_file, format='xyz', index=':')
    n_frames = len(atoms)


    # 读取能量数据
    energy_data = None
    if ener_file:
        energy_data = parse_cp2k_ener_file(ener_file)

    stress_data = None
    if stress_file:
        stress_data = parse_cp2k_stress_file(stress_file)

    trajectory = []

    for j in range(n_frames):
        # 设置晶胞
        if cell is not None:
            atoms[j].set_cell(cell)
            atoms[j].set_pbc(pbc)
        # 添加受力
        # 添加能量
        # 添加位力
        stress_now = [stress_data['xx [bar]'][j],stress_data['yy [bar]'][j],stress_data['zz [bar]'][j],
                      stress_data['yz [bar]'][j],stress_data['xz [bar]'][j],stress_data['xy [bar]'][j]]
        calc = SinglePointCalculator(atoms=atoms[j],
                                     forces=forces[j].get_positions(),
                                     energy=energy_data['Cons Qty[a.u.]'][j],
                                     stress=stress_now)
        atoms[j].calc = calc
    trajectory.append(atoms)
    return atoms


# 使用示例
if __name__ == "__main__":
    filepath = r'E:\Project\57.MgTFSI2_DME_interface\4.pro_43_6\2.pro_43_6_9_Mg1\2.cp2k/'
    program_name = 'pro_43_6_9_Mg1'
    pos_filename = program_name+'-pos-1.xyz'
    frc_filename = program_name+'-frc-1.xyz'
    ener_filename = program_name+'-1.ener'
    stress_filename = program_name+'-1.stress'

    # 示例2：使用模式匹配
    traj = create_aimd_trajectory(
        pos_file=filepath+pos_filename,
        frc_file=filepath+frc_filename,
        ener_file=filepath+ener_filename,
        stress_file=filepath+stress_filename,
        cell=[[40.06, 0, 0], [0, 40.06, 0], [0, 0, 40.06]],  # 晶胞参数
        pbc=True
    )

    # 保存轨迹
    ase.io.write(filepath+program_name+'.traj', traj,format='traj')

    # 读取和查看轨迹
    read_traj = ase.io.read(filepath+program_name+'.traj', index=':')
    print(f"轨迹包含 {len(read_traj)} 帧")

    # 提取特定信息
    energies = [atoms.get_potential_energy() for atoms in read_traj]
    stresses = [atoms.get_stress() for atoms in read_traj]
    energies = np.array(energies)
    stresses = np.array(stresses)
    print(f"能量范围: {energies.min():.3f} 到 {energies.max():.3f} eV")
    print(f'应力范围: {stresses.min():.3f} 到 {stresses.max():.3f} bar')
