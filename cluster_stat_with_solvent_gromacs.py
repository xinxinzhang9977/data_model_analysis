import MDAnalysis as mda
from MDAnalysis.tests.datafiles import GRO, XTC
import numpy as np
import matplotlib.pyplot as plt
import time
import os


class Cluster_stat():

    def __init__(self, file_path, center):
        self.filepath = file_path
        self.sys = mda.Universe(file_path+'NVT_MD.gro', file_path+'NVT_MD.xtc')
        self.frames = len(self.sys.trajectory[::100])
        self.center_sel = center
        self.center = self.sys.select_atoms(center)
        self.q_id = list() # 生成空队列
        self.q_name = list()
        self.head = 0
        self.tail = 0
        self.exist = list()
        self.n_Free = 0
        self.n_CIP = 0
        self.n_AGG = 0
        self.cluster_stat = dict()

    def in_queue(self, id, name):
        self.q_id.append(id)
        self.q_name.append(name)
        self.exist.append(id)
        self.tail += 1

    def stat(self):
        out_file = open(self.filepath+'cluster_py.dat', 'w')
        out_counter = 0
        for ts in self.sys.trajectory[::100]:
            print(ts.frame)
            self.exist = []
            for atom in self.center.atoms:
                if atom.resid not in self.exist:
                    self.q_id = []
                    self.q_name = []
                    self.head = 0
                    self.tail = 0
                    self.in_queue(atom.resid, atom.resname)
                    while self.head < self.tail:
                        id_now = self.q_id[self.head]
                        center_now = self.sys.select_atoms(f'resid {id_now}')
                        name_now = self.q_name[self.head]
                        self.head += 1
                        if center_now.resnames[0] in ['MG', 'TFSI']:
                            if center_now.resnames[0] == 'MG':  # 是Mg,扩展寻找，但Mg周围的所有东西都入队
                                selection = f'around 2.8 group center_now'  # 所有东西都入队
                            elif center_now.resnames[0] == 'TFSI':
                                selection = f'{self.center_sel} and around 2.8 group center_now'  # 是TFSI， 扩展寻找，只有周围的Mg or TFSI才入队
                            new_item = self.sys.select_atoms(selection, center_now=center_now, updating=True)
                            for atom_new in new_item:
                                id_new = atom_new.resid
                                name_new = atom_new.resname
                                if id_new not in self.exist:
                                    self.in_queue(id_new, name_new)

                    print(self.q_id)
                    out_selction_str = 'resid '+' '.join(map(str,self.q_id))
                    out_selction = self.sys.select_atoms(out_selction_str)
                    box_center = self.sys.dimensions[0:3]/2
                    out_sel_center = out_selction.center_of_geometry(pbc=True)
                    while sum(abs(box_center - out_sel_center)) > 0.3:
                        self.sys.atoms.translate(box_center-out_sel_center)
                        out_selction.wrap(center='cog')
                        out_sel_center = out_selction.center_of_geometry()

                    cluster_data = {'MG': 0, 'TFSI': 0, 'DME': 0, 'MEL': 0}
                    resname_data = self.center_sel.split()[1:]
                    for name in cluster_data.keys():
                        cluster_data[name] = self.q_name.count(name)
                    out_str = ''
                    for name in cluster_data.keys():
                        out_str = out_str + ' ' + name + ' ' + str(cluster_data[name]) + ' '
                    if cluster_data['MG'] == 1 and cluster_data['TFSI'] == 2 and cluster_data['DME'] == 2:
                        a = [str(x) for x in self.q_id]
                        print(' '.join(a))
                    a = [str(x) for x in self.q_id]
                    #print(' '.join(a))
                    #print('\n')
                    valence = cluster_data['MG'] * 2 - cluster_data['TFSI']
                    if valence == -2:
                        # 构建目标文件夹路径
                        target_dir = os.path.join(self.filepath, f'cluster_charge_{valence}')

                        # 检查文件夹是否存在，如果不存在则创建
                        if not os.path.exists(target_dir):
                            os.makedirs(target_dir)
                            print(f"创建文件夹: {target_dir}")

                        # 写入输出文件
                        out_counter += 1
                        out_sel_filename = os.path.join(target_dir, f'pro_43_6_charge_{valence}_{out_counter}.pdb')
                        out_selction.write(out_sel_filename)
                    if out_str in self.cluster_stat.keys():
                        self.cluster_stat[out_str] += 1
                    else:
                        self.cluster_stat[out_str] = 1
                    out_str = out_str + ' ' + str(valence)
                    out_file.write(str(ts.frame) + out_str+'\n')
        out_file.close()

    def pic(self):
        cluster_file = open(self.filepath+'cluster_py.dat', 'r')
        cluster_dict = dict()
        line = cluster_file.readline()
        total_cluster = 0
        while line:
            total_cluster += 1
            combine = line.split(' ', 1)[1]
            if combine in cluster_dict.keys():
                cluster_dict[combine] += 1
            else:
                cluster_dict[combine] = 1
            line = cluster_file.readline()
        cluster_dict = sorted(cluster_dict.items(), key=lambda x: x[1])
        label = []
        num = []
        for cluster in cluster_dict:
            num_n = int(cluster[1]) / total_cluster * 100

            if num_n >= 1:
                label.append(cluster[0])
                num.append(num_n)
        plt.figure(figsize=(10, 10))
        plt.subplots_adjust(left=0.22, right=0.95, top=0.95, bottom=0.05)
        b = plt.barh(range(len(label)), num, 0.1, align='edge', tick_label=label)
        for rect in b:
            w = rect.get_width()
            plt.text(w, rect.get_y()+rect.get_height()/2, '{:.2f}%'.format(w), ha='left', va='center')
        plt.title(self.filepath)
        plt.show()
        cluster_file.close()

    def main(self):
        self.stat()
        self.pic()




if __name__ == "__main__":
    start = time.time()
    filepaths = [
                 r'E:\Project\57.MgTFSI2_DME_interface\4.pro_43_6/',
                 #r'D:\Project\43.Melm\7.MgTFSI2_DME_Melm_v11_scaled/',
                 #r'D:\Project\43.Melm\5.MgTFSI2_DME_Melm_v21_scaled/'
                 #r'D:\Project\43.Melm\1.MgHMDS2_THF_scaled/',
                 #r'D:\Project\43.Melm\9.MgHMDS2_THF_Melm_v11_scaled/'
                 ]
    for filepath in filepaths:
        cluster = Cluster_stat(filepath, 'resname MG TFSI')
        cluster.main()
    end = time.time()
    print('运行时间{:.2f}min'.format((end - start)/60))