#!/usr/bin/python
import sys
if len(sys.argv) < 5:
    print('\nNot enough arguments entered.\n')
    print('Usage: dir of species level motifs in transfac format, seq file, number of matrix permutations, bg file, sp, dir of cichlid level motifs in transfac format\n')
    sys.exit()

print('\nImporting...')
import os
import math
import subprocess
import time

print ('Working with spcies motif dir:\t' + sys.argv[1])
print ('Working with cichlid motif dir:\t' + sys.argv[2])

path = '/tgac/workarea/Research-Groups/RG-cichlids/CSPSSMs/'
print ('Working in dir:\t' + path)

# try:
#     os.mkdir(path)
# except FileExistsError:
#     print('Dir ' + path + 'already exists, writing there... ')
name_list = []

for fil in os.listdir(sys.argv[1]):
    print(sys.argv[1] + fil)
    name_list.append([sys.argv[1], fil])

print('\n')

for fil in os.listdir(sys.argv[6]):
    print(sys.argv[6] + fil)
    name_list.append([sys.argv[6], fil])

# sys.exit()
# motif_dict  = {}
# with open(sys.argv[1]) as motif_file:
#     lines = motif_file.read()
#     for motif in lines.split('//'):
#
#         motif.split('\n').append('//')
#         # print (motif)
#         try:
#             motif_dict[motif.split('\n')[1].split()[1]] = motif
#         except IndexError:
#             try:
#                 motif_dict[motif.split('\n')[0].split()[1]] = motif
#             except IndexError:
#                 pass
#
#
# for name, motif in motif_dict.items():
#     name_list.append('temp_' + name + '.transfac')
#     out_file = open(path + 'temp_' + name + '.transfac', 'w')
#     out_file.write(motif.strip())
#     out_file.write('\n//\n')
#     out_file.close()

print('\nSubmitting ' + str(len(name_list)) + ' jobs to the cluster...\n')

for name in name_list:

    if name[1].split('.')[0] == sys.argv[5]:
        out_file_name = path + 'sub_scripts/' + name[1] + '.sh'
        out_file      = open(out_file_name, 'w')

        out_file.write('#!/bin/bash\n#SBATCH -N 1\n#SBATCH -p ei-long\n#SBATCH -c 1\n#SBATCH -n 1\n#SBATCH --mem 500\n#SBATCH -t 10-00:00\n#SBATCH -o ' + path + 'logs/' + name[1] + '.STDOUT\n#SBATCH -e ' + path + 'logs/' + name[1] + '.STDERR\n#SBATCH -J ' + name[1] + '\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=will.j.nash@gmail.com\n\nexport PATH=/tgac/workarea/group-eg/project_cichlids/src/rsat/python-scripts:/tgac/workarea/group-eg/project_cichlids/src/rsat/perl-scripts:/tgac/workarea/group-eg/project_cichlids/src/rsat/bin:/tgac/software/testing/bin/core/../..//make/4.1/x86_64/bin:/tgac/software/testing/binutils/2.25.1/x86_64/bin/:/tgac/software/testing/bin/core/../..//weblogo/2.8.2/x86_64/bin:/tgac/software/testing/bin/core/../..//gcc/5.2.0/x86_64/bin:/tgac/software/testing/libraries/bin/core/../..//openblas/0.2.15/x86_64/bin:/tgac/software/testing/libraries/bin/core/../..//numpy/1.12.0/x86_64/bin:/tgac/software/production/bin/core/../..//hpccore/5/x86_64/bin:/tgac/software/testing/bin/core/../..//python/3.5.1/x86_64/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/tgac/software/production/bin:/tgac/software/production/libraries/bin:/tgac/software/testing/bin:/tgac/software/testing/libraries/bin:/opt/cis/bin:/opt/cis/help:/usr/users/EI_ga011/nashw/.local/bin:/usr/users/EI_ga011/nashw/bin:/tgac/software/production/bin:/tgac/software/production/libraries/bin:/tgac/software/testing/bin:/tgac/software/testing/libraries/bin:/opt/cis/bin:/opt/cis/help:/usr/users/EI_ga011/nashw/.local/bin:/usr/users/EI_ga011/nashw/bin:/usr/users/EI_ga011/nashw/.local/bin:/usr/users/EI_ga011/nashw/bin\n\nsource perl-5.22.1\n\n\nsrun perl /tgac/workarea/group-eg/project_cichlids/src/rsat/perl-scripts/matrix-quality -v 0 -ms ' + ''.join(name) + ' -matrix_format transfac -pseudo 1 -seq positive_set ' + sys.argv[2] + ' -seq_format fasta -perm positive_set ' + sys.argv[3] + ' -bgfile ' + sys.argv[4] + ' -bg_format oligo-analysis  -o ' + path + 'run1_' + name[1] + '_\n')

        # subprocess.run(['sbatch ' + out_file_name], shell = True, check = True)
        # subprocess.run(['rm -f ' + out_file_name], shell = True, check = True)

        # out_file.close()

    else:
        out_file_name = path + 'sub_scripts/' + name[1] + '_' + sys.argv[5] + '.sh'
        out_file      = open(out_file_name, 'w')

        out_file.write('#!/bin/bash\n#SBATCH -N 1\n#SBATCH -p ei-long\n#SBATCH -c 1\n#SBATCH -n 1\n#SBATCH --mem 500\n#SBATCH -t 10-00:00\n#SBATCH -o ' + path + 'logs/' + name[1] + '_' + sys.argv[5] + '.STDOUT\n#SBATCH -e ' + path + 'logs/' + name[1] + '_' + sys.argv[5] + '.STDERR\n#SBATCH -J ' + name[1] + '_' + sys.argv[5] + '\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=will.j.nash@gmail.com\n\nexport PATH=/tgac/workarea/group-eg/project_cichlids/src/rsat/python-scripts:/tgac/workarea/group-eg/project_cichlids/src/rsat/perl-scripts:/tgac/workarea/group-eg/project_cichlids/src/rsat/bin:/tgac/software/testing/bin/core/../..//make/4.1/x86_64/bin:/tgac/software/testing/binutils/2.25.1/x86_64/bin/:/tgac/software/testing/bin/core/../..//weblogo/2.8.2/x86_64/bin:/tgac/software/testing/bin/core/../..//gcc/5.2.0/x86_64/bin:/tgac/software/testing/libraries/bin/core/../..//openblas/0.2.15/x86_64/bin:/tgac/software/testing/libraries/bin/core/../..//numpy/1.12.0/x86_64/bin:/tgac/software/production/bin/core/../..//hpccore/5/x86_64/bin:/tgac/software/testing/bin/core/../..//python/3.5.1/x86_64/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/tgac/software/production/bin:/tgac/software/production/libraries/bin:/tgac/software/testing/bin:/tgac/software/testing/libraries/bin:/opt/cis/bin:/opt/cis/help:/usr/users/EI_ga011/nashw/.local/bin:/usr/users/EI_ga011/nashw/bin:/tgac/software/production/bin:/tgac/software/production/libraries/bin:/tgac/software/testing/bin:/tgac/software/testing/libraries/bin:/opt/cis/bin:/opt/cis/help:/usr/users/EI_ga011/nashw/.local/bin:/usr/users/EI_ga011/nashw/bin:/usr/users/EI_ga011/nashw/.local/bin:/usr/users/EI_ga011/nashw/bin\n\nsource perl-5.22.1\n\n\nsrun perl /tgac/workarea/group-eg/project_cichlids/src/rsat/perl-scripts/matrix-quality -v 0 -ms ' + ''.join(name) + ' -matrix_format transfac -pseudo 1 -seq positive_set ' + sys.argv[2] + ' -seq_format fasta -perm positive_set ' + sys.argv[3] + ' -bgfile ' + sys.argv[4] + ' -bg_format oligo-analysis  -o ' + path + 'run1_' + name[1] + '_' + sys.argv[5] + '_\n')

    out_file.close()

    subprocess.run(['sbatch ' + out_file_name], shell = True, check = True)
    subprocess.run(['rm -f ' + out_file_name], shell = True, check = True)

print('\nDone.\n\n')
