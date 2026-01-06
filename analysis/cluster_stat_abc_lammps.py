# coding=utf-8
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import GRO, XTC
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import shutil
from subprocess import Popen, PIPE, STDOUT

'''
统计团簇的长宽高
- 联动bash脚本，运行Multiwfn计算团簇长宽高

stat
统计团簇组分并计算长宽高
输出cluster_abc.dat

pic
可以单独运行
读取cluster_abc.dat中的数据
出团簇组分分布图，和cluster_stat一样
其他出图方式在abc_stat.py中
可调整限值，出现次数过少的团簇不显示，一般0.1合适
输出数据文件cluster_freq.dat (包括限值以下的数据）

get_abc
计算团簇长宽高
不单独运行

copy_multiwfn
将父文件夹中的Multiwfn程序，get_abc.txt, get_abc.bat拷到数据文件夹中
每个实验第一次运行时启用


'''


class Cluster_stat():

    def __init__(self, file_path, center, step=100):
        self.filepath = file_path
        self.out_file = ABC_NAME
        self.sys = mda.Universe(file_path + DATA_NAME, file_path + XYZ_NAME)
        sys_tmp = mda.Universe(file_path + DATA_NAME)
        self.sys.dimensions = sys_tmp.dimensions

        self.traj = self.sys.trajectory[::step]
        self.frames = len(self.traj)
        self.center_sel = center
        self.center = self.sys.select_atoms(center)
        self.q_id = list()  # 生成空队列
        self.q_name = list()
        self.head = 0
        self.tail = 0
        self.exist = list()
        self.n_Free = 0
        self.n_CIP = 0
        self.n_AGG = 0
        self.cluster_stat = dict()

    def copy_multiwfn(self):
        files = ['Multiwfn.exe', 'get_abc.bat', 'get_abc.txt']
        for file in files:
            father = os.path.dirname(os.path.dirname(self.filepath))
            shutil.copy(father + '/' + file, self.filepath + file)
            # os.system('copy {:s} {:s}'.format(os.path.dirname(self.filepath)+file, self.filepath))

    def get_resname(self, resid, sysname):
        dict = {'normal': [800, 820, 840, 900],
                'less': [800, 820, 830, 880],
                'more': [800, 820, 860, 940],
                'concen': [400, 420, 440, 500],
                'dilute': [1200, 1220, 1240, 1300]}
        if resid <= dict[sysname][0]:
            return 'THF'
        elif resid <= dict[sysname][1]:
            return "Mg"
        elif resid <= dict[sysname][2]:
            return 'Li'
        elif resid <= dict[sysname][3]:
            return 'Cl'

    def in_queue(self, id, name):
        self.q_id.append(id)
        self.q_name.append(name)
        self.exist.append(id)
        self.tail += 1

    def free_cip_agg(self):
        cluster_file = open(self.filepath + self.out_file, 'r')
        line = cluster_file.readline()
        while line:
            cluster_data = dict()
            line_n = line.split()
            l_n = len(line_n)
            for i in range(1, l_n, 2):
                cluster_data[line_n[i]] = int(line_n[i + 1])
            '''
            if cluster_data['M3C'] == 1 and cluster_data['TFSI'] > 0:
                self.n_CIP += 1
            elif cluster_data['M3C'] > 1:
                self.n_AGG += 1
            else:
                self.n_Free += 1
            '''
            '''
            if 0 < cluster_data['TFSI'] <= 2:
                if cluster_data['M3C'] == 0:
                    self.n_Free += cluster_data['TFSI']
                elif cluster_data['M3C'] > 0:
                    self.n_CIP += cluster_data['TFSI']
            elif cluster_data['TFSI'] > 2:
                self.n_AGG += cluster_data['TFSI']
            '''
            if cluster_data['M3C'] == 1:
                if cluster_data['TFSI'] == 0:
                    self.n_Free += cluster_data['M3C']
                elif cluster_data['TFSI'] > 0:
                    self.n_CIP += cluster_data['M3C']
            elif cluster_data['M3C'] > 1:
                self.n_AGG += cluster_data['M3C']
            line = cluster_file.readline()
        total = [self.n_Free, self.n_CIP, self.n_AGG]
        sum_cluster = sum(total) / 100
        print('---------Free CIP AGG Report---------')
        print('{:12s} {:12s} {:12s} {:12s}'.format(' ', 'Free', 'CIP', 'AGG'))
        print('{:12s} {:12d} {:12d} {:12d}'.format('n', self.n_Free, self.n_CIP, self.n_AGG))
        print('{:12s} {:12.2f} {:12.2f} {:12.2f}'.
              format('%', self.n_Free / sum_cluster, self.n_CIP / sum_cluster, self.n_AGG / sum_cluster))
        cluster_file.close()

    def get_abc(self):
        dim = self.sys.dimensions
        box_center = dim[0:3] / 2
        resids = ''
        for id in self.q_id:
            resids = resids + str(id) + ' '
        cluster = self.sys.select_atoms(f'resid {resids}')
        cluster_center = cluster.center_of_geometry(pbc=True)
        print('origin center of geometry at', cluster_center)
        cluster.write(self.filepath + 'before.pdb')
        while sum(abs(box_center - cluster_center)) > 0.3:
            print(sum(abs(box_center - cluster_center)))
            self.sys.atoms.translate(box_center - cluster_center)
            cluster.wrap(center='cog')
            # coord = cluster.pack_into_box(inplace=True)
            # print(coord)
            cluster_center = cluster.center_of_geometry()
            print('new center of geometry as ', cluster_center)
        cluster.write(self.filepath + 'centered.pdb')
        time.sleep(1)
        os.chdir(self.filepath)
        p = Popen("cmd.exe /c" + self.filepath + 'get_abc.bat')
        time.sleep(1)
        abc_file = open(self.filepath + 'get_abc.out', 'r')
        abc = ['0', '0', '0']
        abc_data = abc_file.readline()
        while abc_data:
            if 'Length of the three sides' in abc_data:
                abc = abc_data.split()[5:8]
                break
            abc_data = abc_file.readline()
        abc_file.close()
        return abc

    def stat(self):
        out_file = open(self.filepath + self.out_file, 'w')

        for ts in self.traj:
            print(ts.frame)
            self.exist = []
            for atom in self.center.atoms:
                if atom.resid not in self.exist:
                    self.q_id = []
                    self.q_name = []
                    self.head = 0
                    self.tail = 0
                    self.in_queue(atom.resid, self.get_resname(atom.resid, sys_name))
                    while self.head < self.tail:
                        id_now = self.q_id[self.head]
                        name_now = self.q_name[self.head]
                        self.head += 1
                        center_now = self.sys.select_atoms(f'resid {id_now}')
                        selection = f'{self.center_sel} and around 2.7 group center_now'  # 4.0 for Na, 2.7 for Li
                        new_item = self.sys.select_atoms(selection, center_now=center_now, updating=True, periodic=False)
                        for atom_new in new_item:
                            id_new = atom_new.resid
                            name_new = self.get_resname(atom_new.resid, sys_name)
                            if id_new not in self.exist:
                                self.in_queue(id_new, name_new)

                    print(self.q_id)
                    abc = self.get_abc()
                    print(abc)
                    cluster_data = dict()
                    resname_data = self.center_sel.split()[1:]
                    dict_ion = {5: 'Mg', 6: 'Li', 7: 'Cl'}
                    for name in resname_data:
                        cluster_data[dict_ion[int(name)]] = self.q_name.count(dict_ion[int(name)])
                    out_str = ''
                    for name in cluster_data.keys():
                        out_str = out_str + ' ' + name + ' ' + str(cluster_data[name]) + ' '
                    if out_str in self.cluster_stat.keys():
                        self.cluster_stat[out_str] += 1
                    else:
                        self.cluster_stat[out_str] = 1
                    out_data = '{:6d} {:30s} {:10s} {:10s} {:10s} \n'.format(ts.frame * 5, out_str, *abc)
                    out_file.write(out_data)
        out_file.close()

    def pic(self):
        cluster_file = open(self.filepath + self.out_file, 'r')
        freq_file = open(self.filepath + FREQ_NAME, 'w')
        cluster_dict = dict()
        line = cluster_file.readline()
        while line:
            combine = ' '.join(line.split()[1:-3])

            if combine in cluster_dict.keys():
                cluster_dict[combine] += 1
            else:
                cluster_dict[combine] = 1
            line = cluster_file.readline()
        cluster_dict = sorted(cluster_dict.items(), key=lambda x: x[1])
        label = []
        num = []
        for cluster in cluster_dict:
            num_n = int(cluster[1]) / self.frames
            freq_file.write('{:30s} {:8.3f}\n'.format(cluster[0], num_n))
            if num_n >= 0.1:
                label.append(cluster[0])
                num.append(num_n)
        freq_file.close()
        plt.figure(figsize=(10, 10))
        plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.05)
        b = plt.barh(range(len(label)), num, 0.1, align='edge', tick_label=label)
        for rect in b:
            w = rect.get_width()
            plt.text(w + 0.02, rect.get_y() + rect.get_height() / 2, '{:.2f}'.format(w), ha='left', va='center')
        plt.title(self.filepath)
        plt.savefig(self.filepath + "cluster_stat.png")
        print(self.filepath + 'finished')
        # plt.show()
        cluster_file.close()

    def main(self):
        self.copy_multiwfn()
        self.stat()
        # self.free_cip_agg()
        self.pic()


if __name__ == "__main__":
    ABC_NAME = 'cluster_abc_nopbc.dat'
    DATA_NAME = 'NVTp.data'
    XYZ_NAME = 'NVTp.xyz'
    FREQ_NAME = 'cluster_freq_nopbc.dat'

    start = time.time()
    # filepaths = ['D:/Project/cluster_sampling/MgCl2_LiCl_size_limit/1/',
    # 'D:/Project/cluster_sampling/MgCl2_LiCl_size_limit/2/',
    # 'D:/Project/cluster_sampling/MgCl2_LiCl_size_limit/3/'
    # ]
    title_dict = {'less': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\21/',
                  'normal': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\11/',
                  'concen': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\Concen/',
                  'dilute': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\Dilute/',
                  'more': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\12/',
                  'MNCC-bulk': r'D:\Project\54.MgCl2_NaCl\1.bulk/',
                  'MNCC-efield': r'D:\Project\54.MgCl2_NaCl\2.slab_e_up/',
                  'wall': r'D:\Project\51.cluster_interface\3.wall_vs_slab\1.wall_E\1.MLCC_11/'
                  }
    for sys_name in ['normal', 'less', 'more', 'concen', 'dilute']:
        f_path = title_dict[sys_name]
        step = 100  # 隔100帧统计一次（隔10帧太慢了）
        for i in range(1, 6):
            path = f_path + str(i) + '/'
            print(path)
            cluster = Cluster_stat(path, 'type 5 6 7', step)
            cluster.main()
        end = time.time()
        print('运行时间{:.2f}min'.format((end - start) / 60))
