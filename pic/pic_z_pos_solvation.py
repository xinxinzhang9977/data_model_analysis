import os
import numpy as np
import matplotlib.pyplot as plt
'''
绘制上下两个子图，上图是镁离子的z坐标变化，下图是镁离子配体配位数量的变化
保存图片pic_z_pos_solvation_{selection}.jpg
'''
for i in range(801,821):
    selection = f'resid {i}'
    print(selection)
    filepath = r'/1.800_solvents/6.MD_e_0V/'
    z_pos_filename = f'z_total_{selection}.txt'
    solvation_filename = f'solvation_total_{selection}.dat'
    z_pos_data = np.loadtxt(filepath+z_pos_filename, delimiter='\t')
    solvation_data = np.loadtxt(filepath+solvation_filename, delimiter='\t')
    n_z = len(z_pos_data)
    timestep_z = np.arange(n_z)
    fig,axes = plt.subplots(2,1)

    axes[0].plot(timestep_z,z_pos_data)
    axes[0].set_xlabel('time')
    axes[0].set_ylabel('z')

    axes[1].plot(solvation_data[:,0],solvation_data[:,1],label='TFSI')
    axes[1].plot(solvation_data[:,0],solvation_data[:,2],label='DME')
    axes[1].plot(solvation_data[:,0],solvation_data[:,3],label='Mg')
    axes[1].set_xlabel('time')
    axes[1].set_ylabel('coordination number')
    axes[1].legend()

    plt.title(f'{selection}')
    #plt.show()
    plt.savefig(filepath+f'pic_z_pos_solvation_{selection}.jpg')

