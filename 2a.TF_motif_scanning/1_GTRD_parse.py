#!/usr/bin/python

#=======================================================================================================
# Exptrapolation extension written for paddy
# aim: to construct the tf:tg:tfbs output file from the GTRD database that has been cleaned using bedtools
#intersect to remove overlapping peaks and then bedtools closest to rmeove all non adjacent genes wthin a
#10kb window of peaks, should tka plus and minus strand files and merge them, and then add sequence.
#

###Runs fastest if depolyed over a job array
#use split -l # -a # -d /in /out to chunk the GTRD data then
# SBATCH --array=0:n in the batch script with %A(process) %a(job id) and $SLURM_ARRAY_TASK_ID in the call


#=======================================================================================================

import sys
if len(sys.argv) < 3:
    print('\nNot enough arguments entered...\n')
    print('Usage: plus strand reduced GTRD file, minus strand reduced GTRD file, genome FASTA, tfClass2ensembl.txt\n')
    sys.exit()

print('===================================================================================\n')
print('GTRD_parse.py\n')
print('requires biopython\n')
print("Importing...\n")
import os
import csv
import time
# import collections

from tqdm import tqdm
from Bio.Seq import Seq
from collections import OrderedDict

def readFastaGenerator(seqin):
    for line in seqin:
        yield line.strip()[1:], next(seqin).strip()

print('Writing output to:\t' + str(sys.argv[1].split('.')[0]) + '.parse1.TgTfTfbs.in\n')

out_file1 = open(sys.argv[1].split('.')[0] + '.parse1.TgTfTfbs.in', 'w')

print('===================================================================================\n\n')

print('Working with sequences from from:\t' + str(sys.argv[3] + '\n'))

print('Constructing sequence dictionary...')
seq_dict = {}


with open(sys.argv[3]) as seq_file:
    seq_gen = readFastaGenerator(seq_file)
    for _id, _seq in seq_gen:
        # print(_id)
        seq_dict[_id.split(' ')[0].replace('>','')] = _seq

if len(seq_dict) > 0:
    print('Done. ' + 'Loaded ' + str(len(seq_dict.keys())) + ' chromosomes/scaffolds.\n')
else:
    print('\n\n\t--- Break due to no sequences loaded ---')

# time.sleep(2)
print('===================================================================================\n\n')

print('Loading tfClass2ensembl...')
tfC2ens_dict = {}
with open(sys.argv[4]) as tfClass2ensembl:
    for line in csv.reader(tfClass2ensembl, delimiter = '\t'):
        if line[2] == 'Homo sapiens':
            # print('\t'.join(line))
            if line[0] not in tfC2ens_dict:
                tfC2ens_dict[line[0]] = [line[1]]
            else:
                # print('\nDupliacte tf identity in list\n')
                tfC2ens_dict[line[0]].append(line[1])

# for ip, bip in tfC2ens_dict.items():
#     print(ip)
#     for chip in bip:
#         print('\t' + chip)


print('===================================================================================\n\n')
# sys.exit()


print('Working with GTRD TFBS on the + strand from:\t' + str(sys.argv[1] + '\n'))
print('Working with GTRD TFBS on the - strand from:\t' + str(sys.argv[2] + '\n'))
print('Constructing GTRD dictionary...\n')

tfbs_dict = {}
with open(sys.argv[1]) as tfbs_info:
    for tfbs in tqdm(csv.reader(tfbs_info, delimiter = '\t')):

        for entry in tfC2ens_dict[tfbs[4]]:
            if entry not in tfbs_dict:
                tfbs_dict[entry] = {tfbs[20].split(';')[0].split('=')[1]: [(tfbs[0], tfbs[1], tfbs[2], tfbs[18])]}

            else:
                if tfbs[14].split(';')[0].split('=')[0] not in tfbs_dict[entry]:
                    tfbs_dict[entry][tfbs[20].split(';')[0].split('=')[1]] = [(tfbs[0], tfbs[1], tfbs[2], tfbs[18])]

                else:
                    tfbs_dict[entry][tfbs[20].split(';')[0].split('=')[1]].append((tfbs[0], tfbs[1], tfbs[2], tfbs[18]))

with open(sys.argv[2]) as tfbs_info:
    for tfbs in tqdm(csv.reader(tfbs_info, delimiter = '\t')):
        for entry in tfC2ens_dict[tfbs[4]]:
            if entry not in tfbs_dict:
                tfbs_dict[entry] = {tfbs[20].split(';')[0].split('=')[1]: [(tfbs[0], tfbs[1], tfbs[2], tfbs[18])]}

            else:
                if tfbs[14].split(';')[0].split('=')[0] not in tfbs_dict[entry]:
                    tfbs_dict[entry][tfbs[20].split(';')[0].split('=')[1]] = [(tfbs[0], tfbs[1], tfbs[2], tfbs[18])]

                else:
                    tfbs_dict[entry][tfbs[20].split(';')[0].split('=')[1]].append((tfbs[0], tfbs[1], tfbs[2], tfbs[18]))

for number, (tf, targets) in enumerate(tqdm(tfbs_dict.items())):
    # print(tf)
    # time.sleep(1)
    for tg, peaks in targets.items():
        # print('\t' + tg)
        for tfbs in peaks:
            if tfbs[3] == '+':
                out_file1.write(tf + '\t' + tg + '_(+)' + '\t' + str(seq_dict[tfbs[0].replace('chr', '')][int(tfbs[1])-1:int(tfbs[2])]) + '\n')
            elif tfbs[3] == '-':
                out_file1.write(tf + '\t' + tg + '_(-)' + '\t' + str(Seq(seq_dict[tfbs[0].replace('chr', '')][int(tfbs[1])-1:int(tfbs[2])]).reverse_complement()) + '\n')
            else:
                print('\n\n\n --- Error in dict construction -- line l07 --- \n\n\n')
                sys.exit()

    with open(sys.argv[1].split('.')[0] + '_' + str(number), 'w') as individual_tf_outfile:
        for tg, peaks in targets.items():
            # print('\t' + tg)
            for tfbs in peaks:
                if tfbs[3] == '+':
                    individual_tf_outfile.write(tf + '\t' + tg + '_(+)' + '\t' + str(seq_dict[tfbs[0].replace('chr', '')][int(tfbs[1])-1:int(tfbs[2])]) + '\n')
                elif tfbs[3] == '-':
                    individual_tf_outfile.write(tf + '\t' + tg + '_(-)' + '\t' + str(Seq(seq_dict[tfbs[0].replace('chr', '')][int(tfbs[1])-1:int(tfbs[2])]).reverse_complement()) + '\n')
                else:
                    print('\n\n\n --- Error in dict construction -- line l20 --- \n\n\n')
                    sys.exit()




out_file1.close()
print('\nOutput file closed.\n')
sys.exit()
