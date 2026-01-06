import MDAnalysis as mda
from Bio.Data.IUPACData import atom_weights
from MDAnalysis.tests.datafiles import GRO, XTC
import numpy as np
import matplotlib.pyplot as plt
import time


class Cluster_stat():

    def __init__(self, file_path, center):
        self.filepath = file_path
        self.sys = mda.Universe(file_path+'NVTp.data', file_path+'NVTp_1_1.xyz')
        self.sys.dimensions = np.array([42.63228,42.63228,42.63228,90,90,90])
        n_atoms = len(self.sys.atoms)
        for i in range(n_atoms):
            self.sys.atoms[i].resname = self.find_name(self.sys.atoms[i].resid)
        self.frames = len(self.sys.trajectory[::1000])
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
    def find_name(self, id):
        if id <=400:
            name = 'DME'
        elif id <=440:
            name = 'LI'
        else:
            name = 'TFSI'
        return name

    def free_cip_agg(self):
        cluster_file = open(self.filepath + 'cluster_py.dat', 'r')
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
            if cluster_data['MG'] == 1:
                if cluster_data['TFSI'] == 0:
                    self.n_Free += cluster_data['MG']
                elif cluster_data['TFSI'] > 0:
                    self.n_CIP += cluster_data['MG']
            elif cluster_data['MG'] > 1:
                self.n_AGG += cluster_data['MG']
            line = cluster_file.readline()
        total = [self.n_Free, self.n_CIP, self.n_AGG]
        sum_cluster = sum(total)/100
        print('---------Free CIP AGG Report---------')
        print('{:12s} {:12s} {:12s} {:12s}'.format(' ', 'Free', 'CIP', 'AGG'))
        print('{:12s} {:12d} {:12d} {:12d}'.format('n', self.n_Free, self.n_CIP, self.n_AGG))
        print('{:12s} {:12.2f} {:12.2f} {:12.2f}'.
              format('%', self.n_Free/sum_cluster, self.n_CIP/sum_cluster, self.n_AGG/sum_cluster))
        cluster_file.close()

    def stat(self):
        out_file = open(self.filepath+'cluster_py.dat', 'w')
        out_counter = 0
        for ts in self.sys.trajectory[::1000]:
            print(ts.frame)
            self.exist = []
            for atom in self.center.atoms:
                if atom.resid not in self.exist:
                    self.q_id = []
                    self.q_name = []
                    self.head = 0
                    self.tail = 0
                    self.in_queue(atom.resid, self.find_name(atom.resid))
                    while self.head < self.tail:
                        id_now = self.q_id[self.head]
                        center_now = self.sys.select_atoms(f'resid {id_now}')
                        name_now = self.q_name[self.head]
                        self.head += 1
                        if 401 < center_now.resids[0]<=480:
                            if 401<center_now.resids[0]<= 440:  # 是Mg,扩展寻找，但Mg周围的所有东西都入队
                                selection = f'around 2.8 group center_now'  # 所有东西都入队
                            elif 441<center_now.resids[0] <=480:
                                selection = f'{self.center_sel} and around 2.8 group center_now'  # 是TFSI， 扩展寻找，只有周围的Mg or TFSI才入队
                            new_item = self.sys.select_atoms(selection, center_now=center_now, updating=True)
                            for atom_new in new_item:
                                id_new = atom_new.resid
                                name_new = self.find_name(atom_new.resid)
                                if id_new not in self.exist:
                                    self.in_queue(id_new, name_new)

                    print(self.q_id)
                    out_selction_str = 'resid '+' '.join(map(str,self.q_id))
                    out_selction = self.sys.select_atoms(out_selction_str)
                    box_center = np.array([42.63228,42.63228,42.63228]) / 2
                    out_sel_center = out_selction.center_of_geometry()
                    while sum(abs(box_center - out_sel_center)) > 0.3:
                        self.sys.atoms.translate(box_center-out_sel_center)
                        out_selction.wrap(center='cog')
                        out_sel_center = out_selction.center_of_geometry()

                    cluster_data = {'LI': 0, 'TFSI': 0, 'DME': 0, 'MEL': 0}
                    resname_data = self.center_sel.split()[1:]
                    for name in cluster_data.keys():
                        cluster_data[name] = self.q_name.count(name)
                    out_str = ''
                    for name in cluster_data.keys():
                        out_str = out_str + ' ' + name + ' ' + str(cluster_data[name]) + ' '
                    if cluster_data['LI'] == 1 and cluster_data['TFSI'] == 2 and cluster_data['DME'] == 2:
                        a = [str(x) for x in self.q_id]
                        print(' '.join(a))
                    a = [str(x) for x in self.q_id]
                    #print(' '.join(a))
                    #print('\n')
                    valence = cluster_data['LI'] - cluster_data['TFSI']
                    if valence == 0:
                        out_counter += 1
                        out_sel_filename = self.filepath + f'pro_51_11_6_3_{out_counter}.pdb'
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
        #self.free_cip_agg()
        self.pic()




if __name__ == "__main__":
    start = time.time()
    filepaths = [
                 r'E:\Project\57.MgTFSI2_DME_interface\5.pro_51_11_6_3/',
                 #r'D:\Project\43.Melm\7.MgTFSI2_DME_Melm_v11_scaled/',
                 #r'D:\Project\43.Melm\5.MgTFSI2_DME_Melm_v21_scaled/'
                 #r'D:\Project\43.Melm\1.MgHMDS2_THF_scaled/',
                 #r'D:\Project\43.Melm\9.MgHMDS2_THF_Melm_v11_scaled/'
                 ]
    for filepath in filepaths:
        cluster = Cluster_stat(filepath, 'resid 401 to 480')
        cluster.main()
    end = time.time()
    print('运行时间{:.2f}min'.format((end - start)/60))