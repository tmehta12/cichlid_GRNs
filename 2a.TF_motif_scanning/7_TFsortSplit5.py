#!/usr/bin/python
#=====================================================================================================================================================
'''
To be used following the TFBSextrapolation_v1(or 2).py script in the TFBS exprapolation from model to non-model genomes workflow
This script will sort the output of the previous script by TF and then print a FASTA file per TF pritned with each of the
TFBS sequences as a new entry. These TF files can then be passed to the info-gibbs tool from the RSAT suite to be converted into matrices.

'''

import sys

if len(sys.argv) < 3:
    print('\nNot enough arguments entered.\nUsage: TFBSextrapolation output, orthology file , TFBSextrapolation input\n')
    sys.exit()

print('===================================================================================\n')
print('TFsortSplit5.py\n')
print("\nImporting...\n")
import os
import csv
import time
# import math
# import statisitcs
import subprocess

from math import ceil
from statistics import median

print('Done.\n')

print('Working with extroaploated TFBS data from:\t' + str(sys.argv[1]) + '\n')

TFBS_dict = {}
print('Sorting extroaploated TFBS...')
with open(sys.argv[1]) as TFBS_extrap:
    for TFBS_entry in csv.reader(TFBS_extrap, delimiter = ';'):
        if TFBS_entry[0] != 'tf':                                                           # get rid of header line
            if '-' not in TFBS_entry[2]:                                                    # get rid of any gaps that have snuck through
                if TFBS_entry[0] in TFBS_dict:
                    TFBS_dict[TFBS_entry[0]].append((TFBS_entry[1], TFBS_entry[2], 'Hs'))
                else:
                    TFBS_dict[TFBS_entry[0]] = [(TFBS_entry[1], TFBS_entry[2], 'Hs')]


print('\n\n\nA total of ' + str(len(TFBS_dict)) + ' TFs extroaploated to the novel species.\n')

ortho_dict = {}
with open(sys.argv[2]) as ortho_infoHs:
    for ortho in csv.reader(ortho_infoHs, delimiter = ';'):
        if ortho[1] not in ortho_dict:
            ortho_dict[ortho[1]] = [ortho[0]]
        else:
            ortho_dict[ortho[1]].append(ortho[0])



input_TFBS_dict = {}
print('Sorting input TFBS...')
with open(sys.argv[3]) as TFBS_extrap:
    for TFBS_entry in csv.reader(TFBS_extrap, delimiter = '\t'):
        if TFBS_entry[0] != 'tf':
            # print(ortho_dict[TFBS_entry[0]])
            if TFBS_entry[0] in input_TFBS_dict:
                input_TFBS_dict[TFBS_entry[0]].append(len(TFBS_entry[2]))
            else:
                input_TFBS_dict[TFBS_entry[0]] = [len(TFBS_entry[2])]

print('Done.\n')

print('===================================================================================\n')

print('\nCreating output directory:\t' + '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_fasta')
try:
    os.mkdir('/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_fasta')
except(FileExistsError):
    print('!!!\tFASTA dir already exists, writing there...')

print('\nCreating output directory:\t' + '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_ig')
try:
    os.mkdir('/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_ig')
except(FileExistsError):
    print('!!!\tig dir already exists, writing there...')

print('\nCreating output directory:\t' + '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_transfac')
try:
    os.mkdir('/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_transfac')
except(FileExistsError):
    print('!!!\ttransfac dir already exists, writing there...')


print('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n')
#
fail_ig_count    = 0
fail_cm_count    = 0
complete_count   = 0
TFs_too_few_tfbs = []
count            = 0
for TF, TFBS_set in TFBS_dict.items():
    # print(TF + '\t' + str(ortho_dict[TF]))
    count = count + 1
    if len(TFBS_set) > 3:

        print('\n\n' + TF)

        parp = []
        for papa in ortho_dict[TF]:
            try:
                # print(str(median(input_TFBS_dict[papa])) + '\t' + str(len(input_TFBS_dict[papa])))
                parp = parp + input_TFBS_dict[papa]
            except KeyError:
                pass
        print('\tMedian motif length in input data is:\t' + str(ceil(median(parp))) + '\n')
        motif_length = ceil(median(parp))

        out_file = open('/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_fasta/' + TF + '.fasta', 'w')
        prat = []
        for TFBS in TFBS_set:
            # print('>' + TFBS[0] + '_' + TFBS[2])
            out_file.write('>' + TFBS[0] + '_' + TFBS[2] + '\n')
            # print(TFBS[1].upper())
            out_file.write(TFBS[1].upper() + '\n')
            prat.append(len(TFBS[1]))
        out_file.close()
        print('\tMedian of extrap sites is:\t\t' + str(ceil(median(prat))) + '\n')

        with open('./temp_sub_' + str(count) + '.sh', 'w') as outfile:
            outfile.write('#!/bin/bash\n#SBATCH -N 1\n#SBATCH -p ei-medium\n#SBATCH -c 1\n#SBATCH -n 1\n#SBATCH --mem 10000\n#SBATCH -t 5-00:00\n#SBATCH -o Scr/logs/' + TF + '_sub.STDOUT\n#SBATCH -e Scr/logs/' + TF + '_sub.STDERR\n#SBATCH -J ' + TF + '\n' + '#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=will.j.nash@gmail.com\n\nsource python-3.5.1\nsource hpccore-5\nsource biopython-1.66\nsource perl-5.22.1\n\nsrun Cic/src/rsat/contrib/info-gibbs/info-gibbs   -w ' + str(motif_length) + ' -e 1 -r 5 -n 1000 -i ' +  '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_fasta/' + TF + '.fasta > '  +  '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_ig/' + TF + '.ig\n\nsrun perl Cic/src/rsat/perl-scripts/convert-matrix -from infogibbs -to transfac -i ' +  '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_ig/' + TF + '.ig -pseudo 1 -multiply 1 -decimals 1 -perm 0 -bg_pseudo 0.01 -return counts,consensus,parameters -to transfac -o '  +  '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_transfac/' + TF + '.transfac\n\n')

        subprocess.run(['sbatch ./temp_sub_' + str(count) + '.sh'], shell = True, check = True)

        subprocess.run(['rm -f ./temp_sub_' + str(count) + '.sh'], shell = True, check = True)

        print('\n==================================================================================~\n')
        #
        # print('\tRunning RSAT info-gibbs-python on\t' + '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_fasta/' + TF + '.aln.fasta using a motif length of\t' + str(motif_length) + '\tcalculated from the median length of input sites' )
        # try:
        #     subprocess.run([ 'Cic/src/rsat/contrib/info-gibbs/info-gibbs   -w ' + str(motif_length) + ' -e 1 -r 5 -n 1000 -i ' +  '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_fasta/' + TF + '.fasta > '  +  '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_ig/' + TF + '.ig' ], shell = True, check = True)
        #
        #     print('\n==================================================================================~\n')
        #
        # except(subprocess.CalledProcessError):
        #     print('\t!!!\tNOT POSSIBLE TO RUN info-gibbs-python script on\t' + '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_fasta/' + TF + '.fasta\n\n\n')
        #     fail_ig_count = fail_ig_count + 1
        #     continue
        #
        # print('\tRunning RSAT convert-matrix script on ' + '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_ig/' + TF + '.ig')
        #
        # try:
        #     subprocess.run([ 'perl Cic/src/rsat/perl-scripts/convert-matrix -from infogibbs -to transfac -i ' +  '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_ig/' + TF + '.ig -pseudo 1 -multiply 1 -decimals 1 -perm 0 -bg_pseudo 0.01 -return counts,consensus,parameters -to transfac -o '  +  '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_transfac/' + TF + '.transfac' ], shell = True, check = True)
        #
        #     print('\n==================================================================================~\n')
        #     complete_count = complete_count + 1
        #     print('Done.\n\n\n\n\n')
        #
        # except(subprocess.CalledProcessError):
        #     print('\t!!!\tNOT POSSIBLE TO RUN convert-matrix script on\t' + '/'.join(sys.argv[1].split('/')[0:len(sys.argv[1].split('/'))-1]) + '/TF_fasta/' + TF + '.aln.fasta\n\n\n')
        #     fail_cm_count = fail_cm_count + 1
        #     continue

    else:
        print('\t\tNot enough sites to calculate matrix...\n\t\tPassing to ortho group stage.\n\n\t\t' + TF)
        TFs_too_few_tfbs.append(TF)
        continue
#
print('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Complete.')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
# print(str(fail_mafft_count) + '\tprocesses fail at mafft')
# print(str(fail_cons_count)  + '\tprocesses fail at cons')
# print(str(fail_ig_count)    + '\tprocesses fail at ig')
# print(str(fail_cm_count)    + '\tprocesses fail at cm')
# print(str(complete_count)   + '\tprocesses complete successfully')
# print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')

with open(sys.argv[1].split('.')[0] + '.fails.out', 'w') as out_log:
    for entry in TFs_too_few_tfbs:
        out_log.write(entry + '\n')


print('Fails log written.\n\nAll outputs closed.')

sys.exit()
