import sys
sys.path.insert(0, "/Users/rodrigo/google_drive/code/")

import argparse
import numpy as np
from astropy.table import Table

parser = argparse.ArgumentParser(
                    prog='process_dsim_output',
                    description='This script makes the output of dsim readable by Aladin.')

parser.add_argument('io_dir', help='Path of input/output directory.')
# parser.add_argument('--infilename', default='output.txt', required=False, 
#                     help='Name of the file containing the dsim output. Defaults to output.txt')

args = parser.parse_args()
print(args)
io_dir = args.io_dir


columns = ['OBJNAME', 'RA', 'DEC', 'EQX', 'MAG', 'band', 'PCODE', 'LIST', 'SEL?', 'PA', 'L1', 'L2']
dtype = []

for column in columns:
    if column=='SEL?':
        dtype.append('i4')
    else:
        dtype.append('<U12')

def read_target_file(filename):
    with open(filename, 'r') as targets:
        lines = targets.readlines()

    align = []
    targets = []
    for line in lines:
        line_data = line.strip().split()
        if len(line_data)==9:
            align.append(line_data)
        elif len(line_data)==10:
            targets.append(line_data)
    return align, targets


infilename = "{}output.txt".format(io_dir)
out_selected = "{}selected.csv".format(io_dir)
out_not_selected = "{}not_selected.csv".format(io_dir)
out_align = "{}align.csv".format(io_dir)

align, targets = read_target_file(infilename)

df_targets = Table(np.array(targets), names=columns[:-2], dtype=dtype[:-2])
df_align = Table(np.array(align[:-1]), names=columns[:-3], dtype=dtype[:-3])

df_targets[df_targets['SEL?']==1].write(out_selected, format='ascii.csv', overwrite=True)
df_targets[df_targets['SEL?']==0].write(out_not_selected, format='ascii.csv', overwrite=True)
df_align.write(out_align, format='ascii.csv', overwrite=True)

df_targets[df_targets['SEL?']==1]