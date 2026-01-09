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
    title_dict = {
                  'mgtfsi2_dme_800':r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents/'

                  }
    sys_name = 'mgtfsi2_dme_800'
    format = 'lammpstrj'
    input_prefix = "NVEe_pl"
    input_filename = input_prefix + "." + format
    steps = 1000
    f_path = title_dict[sys_name]
    f_path = r'E:\Project\57.MgTFSI2_DME_interface\1.800_solvents\7.MD_init2_nvt_5V/'
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
