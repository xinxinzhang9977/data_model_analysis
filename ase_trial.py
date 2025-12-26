import ase.io
import numpy as np
import re

filepath = r'E:\Project\57.MgTFSI2_DME_interface\4.pro_43_6\cp2k/'
program_name = 'pro_43_6_21_Mg4'
pos_filename = program_name+'-pos-1.xyz'
frc_filename = program_name+'-frc-1.xyz'
ener_filename = program_name+'-1.ener'
stress_filename = program_name+'-1.stress'

pos_file = filepath+pos_filename
atoms = ase.io.read(pos_file, format='xyz',index=':')
print(atoms)