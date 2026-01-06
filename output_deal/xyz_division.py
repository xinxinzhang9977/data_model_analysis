import os


def split_file(input_file, output_prefix, output_suffix="xyz", steps=100):
    input_counter = 0
    output_counter = 0
    with open(input_file, 'r') as input_file:
        output_file_path = f"{output_prefix}_{output_counter + 1}.{output_suffix}"
        output_file = open(output_file_path, 'w')
        print(output_file_path + " writing!")
        first_line = input_file.readline()
        input_file.seek(0)
        for line in input_file:
            if line == first_line:
                input_counter += 1
                if input_counter > steps:
                    output_file.close()
                    output_counter += 1
                    input_counter = 0
                    output_file_path = f"{output_prefix}_{output_counter + 1}.{output_suffix}"
                    output_file = open(output_file_path, 'w')
                    print(output_file_path + " writing!")
            output_file.write(line)

        output_file.close()


if __name__ == "__main__":
    title_dict = {'less': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\21/',
                  'normal': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\11/',
                  'concen': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\Concen/',
                  'dilute': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\Dilute/',
                  'more': r'D:\Project\51.cluster_interface\11.lammps5\1.parameter\12/',
                  'MNCC-bulk': r'D:\Project\54.MgCl2_NaCl\1.bulk/',
                  'MNCC-efield': r'D:\Project\54.MgCl2_NaCl\2.slab_e_up/',
                  'slab-np1': r'D:\Project\51.cluster_interface\12.lammps-interface\8.nodrude_electrode\6.large_init1\1.0V/',
                  'slab-np2':r'D:\Project\51.cluster_interface\12.lammps-interface\8.nodrude_electrode\7.large_init2\1.0V/',
                  'slab-p1': r'D:\Project\51.cluster_interface\12.lammps-interface\9.drude_electrode\2.init1_0V/',
                  'slab-p2': r'D:\Project\51.cluster_interface\12.lammps-interface\9.drude_electrode\3.init2_0V/',
                  'slab-p2_5': r'D:\Project\51.cluster_interface\12.lammps-interface\9.drude_electrode\4.init2_5V_vs_PZC/',
                  'slab-p2_2':r'D:\Project\51.cluster_interface\12.lammps-interface\9.drude_electrode\5.init2_2V_vs_PZC/',
                  'slab-np1_2': r'D:\Project\51.cluster_interface\12.lammps-interface\8.nodrude_electrode\6.large_init1\3.2V_vs_PZC/',
                  'slab-np1_5': r'D:\Project\51.cluster_interface\12.lammps-interface\8.nodrude_electrode\6.large_init1\2.5V_vs_PZC/',
                  'litfsi':r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\1.LiTFSI_DME/',
                  'litfsi40':r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\3.LiTFSI_DME_40/',
                  'litfsi60': r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\10.LiTFSI_DME_60/',
                  'litfsi80': r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\4.LiTFSI_DME_80/',
                  'litfsi100': r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\5.LiTFSI_DME_100/',
                  'liclo': r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\6.LiClO4_THF/',
                  'liclo40': r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\7.LiClO4_THF_40/',
                  'liclo60': r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\11.LiClO4_THF_60/',
                  'liclo80': r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\8.LiClO4_THF_80/',
                  'liclo100': r'D:\Project\51.cluster_interface\11.lammps5\6.model_system\9.LiClO4_THF_100/',
                  'mgtfsi2_dme_800':r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents/'

                  }
    sys_name = 'mgtfsi2_dme_800'
    format = 'lammpstrj'
    input_prefix = "NVEe"
    input_filename = input_prefix + "." + format
    steps = 1000
    f_path = title_dict[sys_name]
    f_path = r'/1.800_solvents/12.MD_init2_4V/'
    paths = []
    for i in [1]:
        path = f_path
        input_file = path + input_filename  # 替换为你的输入文件路径
        if os.path.exists(input_file):
            output_file_prefix = path + input_prefix + "_" + str(i)  # 替换为你想要的输出文件前缀
            split_file(input_file, output_file_prefix, format, steps)
        else:
            print(input_file + ' not found')
            continue
