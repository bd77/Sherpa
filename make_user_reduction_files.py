'''
Created on 23 Aug 2017

@author: degraba
'''

output_folder = 'D:/SHERPA/FUA112/reduction_input_files/perSNAPandPrecursor/'
precursor_list = ['NOx', 'NMVOC', 'NH3', 'PPM', 'SOx']
snap_list = range(1, 11)
header = 'POLL'
for snap in snap_list:
    header = header + '\tMS' + str(snap)
header += '\n'    

# write the 50 files with a reduction for a specific SNAP - precursor combination
for snap_file in range(1, 11):
    for precursor_file in precursor_list:
        fname = 'user_reduction_SNAP' + str(snap_file) + '_' + precursor_file + '.txt'
        f = open(output_folder + fname, 'w')
        f.write(header)
        for precursor in precursor_list:
            f.write(precursor)
            for snap in snap_list:
                if snap_file == snap and precursor_file == precursor:
                    f.write('\t' + str(50))
                else:
                    f.write('\t' + str(0))
            f.write('\n')
        f.close()

# write the 5 files with a reduction for a specific precursor
output_folder = 'D:/SHERPA/FUA112/reduction_input_files/perPercursor/'
for precursor_file in precursor_list:
    fname = 'user_reduction_SNAPall_' + precursor_file + '.txt'
    f = open(output_folder + fname, 'w')
    f.write(header)
    for precursor in precursor_list:
        f.write(precursor)
        for snap in snap_list:
            if precursor_file == precursor:
                f.write('\t' + str(50))
            else:
                f.write('\t' + str(0))
        f.write('\n')
    f.close()

# write the 10 files with a reduction for a specific SNAP sector
output_folder = 'D:/SHERPA/FUA112/reduction_input_files/perSNAP/'
for snap_file in range(1, 11):
    fname = 'user_reduction_SNAP' + str(snap_file) + '_all.txt'
    f = open(output_folder + fname, 'w')
    f.write(header)
    for precursor in precursor_list:
        f.write(precursor)
        for snap in snap_list:
            if snap_file == snap:
                f.write('\t' + str(50))
            else:
                f.write('\t' + str(0))
        f.write('\n')
    f.close()

# write the file with a reduction for all SNAP sectors and precursors
output_folder = 'D:/SHERPA/FUA112/reduction_input_files/all/'
fname = 'user_reduction_SNAPall_all.txt'
f = open(output_folder + fname, 'w')
f.write(header)
for precursor in precursor_list:
    f.write(precursor)
    for snap in snap_list:
        f.write('\t' + str(50))
    f.write('\n')
f.close()

if __name__ == '__main__':
    pass