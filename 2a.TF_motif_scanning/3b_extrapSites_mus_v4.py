#!/usr/bin/python

#=======================================================================================================
#Aim: Substitutes binding sites extracted from ChipSeq databases (GTRD, ORegAnno, JASPAR, HTRI, TRUUST) with
#the core motif as described in the JASPAR2018 db. If this cannot be found, this will replace the binding
#site seq with the top kmer from the Uniprobe db. If this is also not possible the consensus motif will
#be created from the total sites using a method

###Runs fastest if depolyed over a job array
#use split -l # -a # -d /in /out to chunk the GTRD data then
# SBATCH --array=0:n in the batch script with %A(process) %a(job id) and $SLURM_ARRAY_TASK_ID in the call


#=======================================================================================================

import sys
if len(sys.argv) < 6:
    print('\nNot enough arguments entered...\n')
    print('Usage: JASPAR hg38/mm10 infile, JASPAR meta, JASPAR .sites dir, db infile (non-redundant TF:TG:TFBS), Uniprobe top Kmer infile, PARSED hocomoco data\n')
    sys.exit()

print('===================================================================================\n')
print('extrapSites_hum_v1.py\n')
# print('requires biopython\n')
print("Importing...\n")
import os
import csv
import time
import subprocess

from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# from tqdm import tqdm
# from collections import OrderedDict

def readFastaGenerator(seqin):
    for line in seqin:
        yield line.strip()[1:], next(seqin).strip()

# def hitJudge(nMM, lenCut, scoreCut):
#     print('1mm -- score threshold == ' + str((len(motifin)*2) - scoreCut))
#     print('\tmin length == ' + str((len(motifin)- nMM)))
#
#     if a[4] - a[3] >= len(motifin) - lenCut:
#         if a[2] >= (len(motifin)*2) - scoreCut:
#             print(format_alignment(*a))
#             # time.sleep(1)
#             if motifin not in motifout:
#                 return motifin
#                 # line[3].append(motifin)
#         else:
#             print('low score @:\t' + str(a[2]) + '\n')
#     return


def motifMatcher(seqin, motifIn, motifout):
    motifin = str(motifIn).strip()
    if len(motifin) >= 5:
        for a in pairwise2.align.localms(str(seqin), str(motifin), 2, -2, -10, -1):
            if 5 <= len(motifin) < 6:

                # motifOut = hitJudge(1, 0, 4)
                #
                # return motifOut

                print('1mm -- score threshold == ' + str((len(motifin)*2)-8))
                print('\tmin length == ' + str((len(motifin)-2)))

                if a[4] - a[3] >= len(motifin):
                    if a[2] >= (len(motifin)*2) - 8:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                elif a[4] - a[3] >= len(motifin) - 2:
                    if a[2] >= (len(motifin)*2):
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')
                else:
                    print('REJECTED: alignment score: ' + str(a[2]) + '\talignment length: ' + str(a[4]-a[3]) + '\n')

            elif 6 <= len(motifin) < 12:
                print('2mm -- score threshold == ' + str((len(motifin)*2)-12))
                print('\tmin length == ' + str((len(motifin)-3)))


                if a[4] - a[3] >= len(motifin):
                    if a[2] >= (len(motifin)*2) - 12:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                elif a[4] - a[3] >= len(motifin) - 2:
                    if a[2] >= (len(motifin)*2) - 8:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                elif a[4] - a[3] >= len(motifin) - 3:
                    if a[2] >= (len(motifin)*2):
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)

                    else:
                        print('low score @:\t' + str(a[2]) + '\n')
                else:
                    print('REJECTED: alignment score: ' + str(a[2]) + '\talignment length: ' + str(a[4]-a[3]) + '\n')

            elif 12 <= len(motifin) < 18:
                print('3mm --  score threshold == ' + str((len(motifin)*2)-16))
                print('\tmin length     == ' + str((len(motifin)-4)) + '\n')

                # print(format_alignment(*a))

                # print(a)

                if a[4] - a[3] >= len(motifin):
                    if a[2] >= (len(motifin)*2)-16:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                elif a[4] - a[3] >= len(motifin) - 2:
                    if a[2] >= (len(motifin)*2)-12:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                elif a[4] - a[3] >= len(motifin) - 3:
                    if a[2] >= (len(motifin)*2)-8:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                elif a[4] - a[3] >= len(motifin) - 4:
                    if a[2] >= (len(motifin)*2):
                        # print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')
                else:
                    print('REJECTED: alignment score: ' + str(a[2]) + '\talignment length: ' + str(a[4]-a[3]) + '\n')

            elif len(motifin) >= 18:
                print('\t4mm -- score threshold == ' + str((len(motifin)*2)-20))
                print('\tmin length == ' + str((len(motifin)-5)))

                if a[4] - a[3] >= len(motifin):
                    if a[2] >= (len(motifin)*2)-20:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                if a[4] - a[3] >= len(motifin) - 2:
                    if a[2] >= (len(motifin)*2)-16:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                elif a[4] - a[3] >= len(motifin) - 3:
                    if a[2] >= (len(motifin)*2)-12:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                elif a[4] - a[3] >= len(motifin) - 4:
                    if a[2] >= (len(motifin)*2)-8:
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')

                elif a[4] - a[3] >= len(motifin) - 5:
                    if a[2] >= (len(motifin)*2):
                        print(format_alignment(*a))
                        # time.sleep(1)
                        if motifin not in motifout:
                            return motifin
                            # line[3].append(motifin)
                    else:
                        print('low score @:\t' + str(a[2]) + '\n')
                else:
                    print('REJECTED: alignment score: ' + str(a[2]) + '\talignment length: ' + str(a[4]-a[3]) + '\n')
    return

print('===================================================================================\n')

print('Making Ens2JASPAR dictionary...')
ens2jas_dict = {}

with open(sys.argv[1]) as ens2jas_in:
    for motif in csv.reader(ens2jas_in, delimiter = '\t'):
        stored_ver = 1
        # print(motif)
        with open(sys.argv[2]) as jas_meta:
            for entry in csv.reader(jas_meta, delimiter = ','):
                if entry[2] == 'Mus musculus':
                    if entry[0].split('.')[0] == motif[0]:
                        # print('\t' + entry[0])
                        for pap in os.listdir(sys.argv[3]):
                            if entry[0] in pap:
                                # print('\t\t' + pap)
                                if int(entry[0].split('.')[1]) > stored_ver:
                                    stored_ver = int(entry[0].split('.')[1])
                        # print('\t' + str(entry[0].split('.')[0]) + '\t' +  str(stored_ver))

        for mot in motif[3].split(';'):
            # print(mot)
            if mot not in ens2jas_dict:
                ens2jas_dict[mot] = [[(motif[2], motif[0] + '.' + str(stored_ver))]]
            else:
                ens2jas_dict[mot][0].append((motif[2], motif[0] + '.' + str(stored_ver)))


print('Parsing Uniprobe top Kmers and adding...')
# uniprobe_dict = {}
with open(sys.argv[5]) as uniprobTopKmers:
    for line in csv.reader(uniprobTopKmers, delimiter = '\t'):
        # print(line)
        for entry in line[2].split(';'):
            if entry in ens2jas_dict:
                ens2jas_dict[entry].append(line[3])

print('Done.\n')

print('===================================================================================\n')

print('Parsing HOCOMOCO data...\n')
hocoParse = {}
with open(sys.argv[6]) as hocomoco:
    for line in csv.reader(hocomoco, delimiter = '\t'):
        hocoParse[line[0]] = line[1]

print('File locations for ' + str(len(hocoParse)) + ' HOCOMOCO dinucleatide word alignments parsed.\n')
# sys.exit()
print('===================================================================================\n\n')

print('Parsing TF:TGTFBS in file:\t' + str(sys.argv[4] + '\n'))

counter           = 1
count_list        = []
no_site_name_list = []
with open(sys.argv[4]) as db_tftgtfbs_file, open(sys.argv[4].split('.')[0] + '.sites3.out', 'w') as out_file1:
    for line in csv.reader(db_tftgtfbs_file, delimiter = '\t'):
        line.append([])
        print('Target:\t' + str(line))
        print(line)
        # time.sleep(5)

        if line[0] in hocoParse:
            print('\nHOCOMOCO data for this TF! Loading...\n')
            print('Cic/data/motifs/extrapolation/hocomoco/mouse_words/' + hocoParse[line[0]] + '.words open...\n')
            # time.sleep(5)
            with open('Cic/data/motifs/extrapolation/hocomoco/mouse_words/' + hocoParse[line[0]] + '.words') as hoco_in:
                for word in hoco_in:
                    print('\tTFBS: ' + line[2])
                    print('\tSite: ' + str(word.strip()))

                    print('\tSite len: ' + str(len(word.strip())) + '\n')

                    mot = motifMatcher(line[2], word, line[3])

                    if mot is not None:
                        print('\n\n' + str(mot) + '\n\n')
                        count_list.append(mot)
                        line[3].append(mot)

            for ele in line[3]:
                print(ele)
                if len(ele) == len(max(line[3], key=len)):
                    # print(line[2])
                    print(line[0]+'\t'+'\t'.join(line[1:2]) + '\t' + ele + '\n')
                    # time.sleep(10)
                    out_file1.write(line[0]+'\t'+'\t'.join(line[1:2]) + '\t' + ele + '\n')
                    print('\a')
        else:
            print('\nNo HOCOMOCO data.\n')

            # print('\n' + str(len(count_list)) + '\n')

            if line[0] in ens2jas_dict:
                if len(ens2jas_dict[line[0]][0]) == 1:
                    print('\t' + str(ens2jas_dict[line[0]]))

                # sys.exit()
                    try:
                        with open(sys.argv[3] + ens2jas_dict[line[0]][0][0][1] + '.sites_singleLine') as sites_file:
                            for site_id, site_seq in readFastaGenerator(sites_file):
                                # print('\t\t' + site_id)
                                # print('\t\t' + site_seq)
                                # time.sleep(.1)

                                seq_to_match = ''.join([base for base in site_seq if base.isupper()])

                                print('\tTFBS: ' + line[2])
                                print('\tSite: ' + str(seq_to_match))

                                print('\tSite len: ' + str(len(seq_to_match)))

                                mot = motifMatcher(line[2], seq_to_match, line[3])

                                if mot is not None:
                                    print('\n\n' + str(mot) + '\n\n')
                                    count_list.append(mot)
                                    line[3].append(mot)


                                # if len(seq_to_match) >= 5:
                                #     for a in pairwise2.align.localms(str(line[2]), str(seq_to_match), 2, -2, -10, -1):
                                #         if 5 <= len(seq_to_match) < 6:
                                #             print('1mm -- score threshold == ' + str((len(seq_to_match)*2)-4))
                                #
                                #             if a[2] >= (len(seq_to_match)*2)-4:
                                #                 print(format_alignment(*a))
                                #                 # time.sleep(1)
                                #                 if seq_to_match not in line[3]:
                                #                     line[3].append(seq_to_match)
                                #             else:
                                #                 print('low score @:\t' + str(a[2]))
                                #
                                #         elif 6 <= len(seq_to_match) < 12:
                                #             print('2mm -- score threshold == ' + str((len(seq_to_match)*2)-8))
                                #
                                #             if a[2] >= (len(seq_to_match)*2)-8:
                                #                 print(format_alignment(*a))
                                #                 # time.sleep(1)
                                #                 if seq_to_match not in line[3]:
                                #                     line[3].append(seq_to_match)
                                #             else:
                                #                 print('low score @:\t' + str(a[2]))
                                #
                                #         elif 12 <= len(seq_to_match) < 18:
                                #             print('3mm -- score threshold == ' + str((len(seq_to_match)*2)-12))
                                #
                                #             if a[2] >= (len(seq_to_match)*2)-12:
                                #                 print(format_alignment(*a))
                                #                 # time.sleep(1)
                                #                 if seq_to_match not in line[3]:
                                #                     line[3].append(seq_to_match)
                                #             else:
                                #                 print('low score @:\t' + str(a[2]))
                                #
                                #         elif len(seq_to_match) >= 18:
                                #             print('4mm -- score threshold == ' + str((len(seq_to_match)*2)-16))
                                #
                                #             if a[2] >= (len(seq_to_match)*2)-16:
                                #                 print(format_alignment(*a))
                                #                 # time.sleep(1)
                                #                 if seq_to_match not in line[3]:
                                #                     line[3].append(seq_to_match)
                                #
                                #             else:
                                #                 print('low score @:\t' + str(a[2]))

                    except FileNotFoundError:
                        print('\tNo JASPAPAR2018 sites file available for ' + str(ens2jas_dict[line[0]][0][0][1]) + '...')
                        if len(ens2jas_dict[line[0]]) == 2:
                            print('\t\tRecovered Uniprobe top kmers for ' + str(ens2jas_dict[line[0]][0][0][1]) + '...\n')
                            time.sleep(5)
                            print('\t\tScanning site using top kmers:')
                            for topKmer in ens2jas_dict[line[0]][1].split(','):
                                time.sleep(1)

                                mot = motifMatcher(line[2], topKmer, line[3])

                                if mot is not None:
                                    print('\n\n' + str(mot) + '\n\n')
                                    count_list.append(mot)
                                    line[3].append(mot)

                                        # if len(topKmer) >= 5:
                                    # for a in pairwise2.align.localms(str(line[2]), str(topKmer), 2, -2, -10, -1):
                                    #     if 5 <= len(topKmer) <= 6:
                                    #         print('\t\t1mm -- score threshold == ' + str((len(topKmer)*2)-2))
                                    #
                                    #         if a[2] >= (len(topKmer)*2)-2:
                                    #             print(format_alignment(*a))
                                    #             # time.sleep(1)
                                    #             if topKmer not in line[3]:
                                    #                 line[3].append(topKmer)
                                    #         else:
                                    #             print('low score @:\t' + str(a[2]))
                                    #
                                    #     elif 7 <= len(topKmer) <= 12:
                                    #         print('\t\t2mm -- score threshold == ' + str((len(topKmer)*2)-4))
                                    #
                                    #         if a[2] >= (len(topKmer)*2)-4:
                                    #             print(format_alignment(*a))
                                    #             # time.sleep(1)
                                    #             if topKmer not in line[3]:
                                    #                 line[3].append(topKmer)
                                    #         else:
                                    #             print('low score @:\t' + str(a[2]))
                                    #
                                    #     elif 13 <= len(topKmer) <= 18:
                                    #         print('\t\t3mm -- score threshold == ' + str((len(topKmer)*2)-6))
                                    #
                                    #         if a[2] >= (len(topKmer)*2)-6:
                                    #             print(format_alignment(*a))
                                    #             # time.sleep(1)
                                    #             if topKmer not in line[3]:
                                    #                 line[3].append(topKmer)
                                    #         else:
                                    #             print('low score @:\t' + str(a[2]))
                                    #
                                    #     elif len(topKmer) >= 19:
                                    #         print('\t\t4mm -- score threshold == ' + str((len(topKmer)*2)-8))
                                    #
                                    #         if a[2] >= (len(topKmer)*2)-8:
                                    #             print(format_alignment(*a))
                                    #             # time.sleep(1)
                                    #             if topKmer not in line[3]:
                                    #                 line[3].append(topKmer)
                                    #
                                    #         else:
                                    #             print('\t\t\tlow score @:\t' + str(a[2]))

                        else:
                            print('\t\tNo Uniprobe top kmers available either...')
                            continue


                    if len(line[3]) == 0:
                        print('\n\n\tNo sites found in exisiting JASPAR2018 sites file...') #GOES TO UNIPROBE
                        print('\tNo matching topKmers found in Uniprobe db ...\n\n') #GOES TO UNIPROBE

                        continue

                    print(line)
                    print('\n\n\n')
                    print(line[0:3])

                    # sys.exit()


                    for ele in line[3]:
                        print(ele)
                        if len(ele) == len(max(line[3], key=len)):
                            # print(line[2])
                            print('\t'.join(line[0:2]) + '\t' + ele + '\n')
                            # time.sleep(10)
                            out_file1.write('\t'.join(line[0:2]) + '\t' + ele + '\n')
                        # else:
                        #     print(line[2])
                        #     print(' -----  nub discarded')


                else:

                    tf_suffix = ''
                    for motif_variant in ens2jas_dict[line[0]][0]:

                        print('\t' + str(motif_variant))
                        if motif_variant[0].endswith('(var.2)'):
                            tf_suffix = '_v2'
                            print('\t\tMotif suffix updated.\n')

                        try:
                            with open(sys.argv[3] + motif_variant[1] + '.sites_singleLine') as sites_file:
                                for site_id, site_seq in readFastaGenerator(sites_file):

                                    seq_to_match = ''.join([base for base in site_seq if base.isupper()])

                                    print('\tTFBS: '+line[2])
                                    print('\tSite: '+str(seq_to_match))

                                    print('\tSite len: '+str(len(seq_to_match)))

                                    mot = motifMatcher(line[2], seq_to_match, line[3])

                                    if mot is not None:
                                        print('\n\n' + str(mot) + '\n\n')
                                        count_list.append(mot)
                                        line[3].append(mot)

                                    # if len(seq_to_match) >= 5:
                                    #     for a in pairwise2.align.localms(str(line[2]), str(seq_to_match), 2, -2, -10, -1):
                                    #         if 5 <= len(seq_to_match) < 6:
                                    #             print('1mm -- score threshold == ' + str((len(seq_to_match)*2)-2))
                                    #
                                    #             if a[2] >= (len(seq_to_match)*2)-2:
                                    #                 print(format_alignment(*a))
                                    #                 # time.sleep(1)
                                    #                 if seq_to_match not in line[3]:
                                    #                     line[3].append(seq_to_match)
                                    #             else:
                                    #                 print('low score @:\t' + str(a[2]))
                                    #
                                    #         elif 6 <= len(seq_to_match) < 12:
                                    #             print('2mm -- score threshold == ' + str((len(seq_to_match)*2)-4))
                                    #
                                    #             if a[2] >= (len(seq_to_match)*2)-4:
                                    #                 print(format_alignment(*a))
                                    #                 # time.sleep(1)
                                    #                 if seq_to_match not in line[3]:
                                    #                     line[3].append(seq_to_match)
                                    #             else:
                                    #                 print('low score @:\t' + str(a[2]))
                                    #
                                    #         elif 12 <= len(seq_to_match) < 18:
                                    #             print('3mm -- score threshold == ' + str((len(seq_to_match)*2)-6))
                                    #
                                    #             if a[2] >= (len(seq_to_match)*2)-6:
                                    #                 print(format_alignment(*a))
                                    #                 # time.sleep(1)
                                    #                 if seq_to_match not in line[3]:
                                    #                     line[3].append(seq_to_match)
                                    #             else:
                                    #                 print('low score @:\t' + str(a[2]))
                                    #
                                    #         elif len(seq_to_match) >= 18:
                                    #             print('4mm -- score threshold == ' + str((len(seq_to_match)*2)-8))
                                    #
                                    #             if a[2] >= (len(seq_to_match)*2)-8:
                                    #                 print(format_alignment(*a))
                                    #                 # time.sleep(1)
                                    #                 if seq_to_match not in line[3]:
                                    #                     line[3].append(seq_to_match)
                                    #
                                    #             else:
                                    #                 print('low score @:\t' + str(a[2]))

                        except FileNotFoundError:
                            print('\tNo JASPAPAR2018 sites file available for ' + str(motif_variant[1]) + '...')
                            if len(motif_variant) == 2:
                                print('\t\tRecovered Uniprobe top kmers for ' + str(motif_variant[1]) + '...\n')
                                time.sleep(5)
                                print('\t\tScanning site using top kmers:')
                                for topKmer in ens2jas_dict[line[0]][1].split(','):

                                    mot = motifMatcher(line[2], topKmer, line[3])

                                    if mot is not None:
                                        print('\n\n' + str(mot) + '\n\n')
                                        count_list.append(mot)
                                        line[3].append(mot)

                                    # if len(topKmer) >= 5:
                                    #     for a in pairwise2.align.localms(str(line[2]), str(topKmer), 2, -2, -10, -1):
                                    #         if 5 <= len(topKmer) <= 6:
                                    #             print('\t\t1mm -- score threshold == ' + str((len(topKmer)*2)-2))
                                    #
                                    #             if a[2] >= (len(topKmer)*2)-2:
                                    #                 print(format_alignment(*a))
                                    #                 # time.sleep(1)
                                    #                 if topKmer not in line[3]:
                                    #                     line[3].append(topKmer)
                                    #             else:
                                    #                 print('low score @:\t' + str(a[2]))
                                    #
                                    #         elif 7 <= len(topKmer) <= 12:
                                    #             print('\t\t2mm -- score threshold == ' + str((len(topKmer)*2)-4))
                                    #
                                    #             if a[2] >= (len(topKmer)*2)-4:
                                    #                 print(format_alignment(*a))
                                    #                 # time.sleep(1)
                                    #                 if topKmer not in line[3]:
                                    #                     line[3].append(topKmer)
                                    #             else:
                                    #                 print('low score @:\t' + str(a[2]))
                                    #
                                    #         elif 13 <= len(topKmer) <= 18:
                                    #             print('\t\t3mm -- score threshold == ' + str((len(topKmer)*2)-6))
                                    #
                                    #             if a[2] >= (len(topKmer)*2)-6:
                                    #                 print(format_alignment(*a))
                                    #                 # time.sleep(1)
                                    #                 if topKmer not in line[3]:
                                    #                     line[3].append(topKmer)
                                    #             else:
                                    #                 print('low score @:\t' + str(a[2]))
                                    #
                                    #         elif len(topKmer) >= 19:
                                    #             print('\t\t4mm -- score threshold == ' + str((len(topKmer)*2)-8))
                                    #
                                    #             if a[2] >= (len(topKmer)*2)-8:
                                    #                 print(format_alignment(*a))
                                    #                 # time.sleep(1)
                                    #                 if topKmer not in line[3]:
                                    #                     line[3].append(topKmer)
                                    #
                                    #             else:
                                    #                 print('\t\t\tlow score @:\t' + str(a[2]))

                            else:
                                print('\t\tNo Uniprobe top kmers available either...')
                                continue

                        if len(line[3]) == 0:
                            print('\n\n\tNo sites found in exisiting JASPAR2018 sites file...') #GOES TO UNIPROBE
                            print('\tNo matching topKmers found in Uniprobe db ...\n\n') #GOES TO UNIPROBE

                            continue
                            # time.sleep(30)

                        print(line)
                        print('\n\n\n')
                        print(line[0:3])


                        for ele in line[3]:
                            print(ele)
                            if len(ele) == len(max(line[3], key=len)):
                                # print(line[2])
                                print(line[0]+tf_suffix+'\t'+'\t'.join(line[1:2]) + '\t' + ele + '\n')
                                # time.sleep(10)
                                out_file1.write(line[0]+tf_suffix+'\t'+'\t'.join(line[1:2]) + '\t' + ele + '\n')
                            # else:
                            #     print(line[2])
                            #     print(' -----  nub discarded')


                        # time.sleep(30)

            else:
                print('no sites...') #GOES TO UNIPROBE
                # print(str(ens2jas_dict[line[0]]))
                try:
                    if len(ens2jas_dict[line[0]]) == 2:
                        print('OMG UNIPROBE TO THE RESCUE')
                        sys.exit()
                except KeyError:
                    print('...not even in the ref dict')
                    print('\nAdding tfbs to fasta\n\n\n')

                    # try:
                    #     os.mkdir(sys.argv[3] + 'tf_m2_' + line[0])
                    # except FileExistsError:
                    #     print('Dir ' + sys.argv[3] + 'tf_m2_' + line[0] + ' already exists, writing there... ')
                    #
                    #
                    # if 'tf_m2_' + line[0] + '.fasta' in os.listdir(sys.argv[3] + 'tf_m2_' + line[0]):
                    #     with open(sys.argv[3] + 'tf_m2_' + line[0] + '/tf_m2_' + line[0] +  '.fasta', 'a') as tf_out_fasta:
                    #         tf_out_fasta.write('>' + str(counter).strip() + '_' + line[1].strip() + '\n')
                    #         tf_out_fasta.write(line[2].strip() + '\n')
                    # else:
                    #     with open(sys.argv[3] + 'tf_m2_' + line[0] + '/tf_m2_' + line[0] +  '.fasta', 'w') as tf_out_fasta:
                    #         tf_out_fasta.write('>' + str(counter).strip() + '_' + line[1].strip() + '\n')
                    #         tf_out_fasta.write(line[2].strip() + '\n')
                    #
                    # counter = counter + 1
                    if line[0] not in no_site_name_list:
                        no_site_name_list.append(line[0])


print('\n' + str(len(count_list)) + '\n')

print('\n\n\nEnd of HOCOMOCO/JASPAR/UNIPROBE sets reached.\n\n\n')
if len(no_site_name_list) > 0:
    for ent in no_site_name_list:
        print(str(ent) + ' shoud have a tf_m2_...\n\n\n')

sys.exit()


for lol in no_site_name_list:

    # print('Aligning using mafft...')
    # subprocess.run(['mafft --localpair --maxiterate 1000 --quiet ' + sys.argv[3] + 'tf_m2_' + lol + '/tf_m2_' + lol + '.fasta > ' + sys.argv[3] + 'tf_m2_' + lol + '/tf_m2_' + lol + '_aln.fasta'], shell = True, check = True)
    #
    # print('Done.\n')

    print('Scanning for motifs using MEME...\n')
    subprocess.run(['meme ' + sys.argv[3] + 'tf_m2_' + lol + '/tf_m2_' + lol + '.fasta -dna -oc ' + sys.argv[3] + 'tf_m2_' + lol + '/ -minw 5 -maxw 30 -maxsize 100000000 -nostatus'], shell = True, check = True)
    print('Done.\n')


    filer = open(sys.argv[3] + 'tf_m2_' + lol + '/meme.txt')
    liner = filer.readlines()

    switch     = 0
    switch2    = 0
    store_dict = {}
    for f in liner:
        if 'Motif 1 sites sorted by position p-value' in f:
            switch = switch + 1
        if switch == 1:
            if f[0].isdigit():
                store_dict[f.split()[0].split('_')[1]] = f.split()[1:]
            elif f.startswith('-'):
                switch2 = switch2 + 1

        if switch2 == 3:
            break

    core_motif_list = []
    for stick, stack in store_dict.items():
        print(stick + '\t' + str(stack))
        if len(stack) == 5:
            if stack[3].strip('-') not in core_motif_list:
                core_motif_list.append(stack[3].strip('-'))
        elif len(stack) == 3:
                if stack[2].strip('-') not in core_motif_list:
                    core_motif_list.append(stack[2].strip('-'))
        else:
            print('\n\n\n\tBreak due to no non-conventional motif --- line 464\n\n\n')
            sys.exit()

    longest_core = []
    print('\n')
    for ele in core_motif_list:
        if len(ele) == len(max(core_motif_list, key=len)):
            longest_core.append(ele)
            print('\t\t\t\t\t\t   ' + ele)

    print('\n\nScript has found ' + str(len(longest_core)) + ' core motifs for ' + str(lol) + '\n\n\n')

    with open(sys.argv[3] + 'tf_m2_' + lol + '/tf_m2_' + lol +  '.memeSites', 'w') as tf_out_fasta:
        for kek, pek in store_dict.items():

            if len(pek[3].strip('-')) == len(longest_core[0]):
                tf_out_fasta.write('>' + str(kek) + '\n')
                tf_out_fasta.write(pek[3].strip('-') + '\n')


with open(sys.argv[4]) as db_tftgtfbs_file, open(sys.argv[4].split('.')[0] + '.sites3.out', 'a') as out_file2:
    for line in csv.reader(db_tftgtfbs_file, delimiter = '\t'):
        line.append([])
        print('Target:\t' + str(line))
        print('Opening:\t' + sys.argv[3] + 'tf_m2_' + line[0] + '/tf_m2_' + line[0] + '.memeSites')
        try:
            with open(sys.argv[3] + 'tf_m2_' + line[0] + '/tf_m2_' + line[0] + '.memeSites') as sites_file:
                for site_id, site_seq in readFastaGenerator(sites_file):
                    print('\t\t' + site_id)
                    print('\t\t' + site_seq)
                    # time.sleep(.1)

                    seq_to_match = ''.join([base for base in site_seq if base.isupper()])

                    print('\tTFBS: '     + line[2])
                    print('\tSite: '     + str(seq_to_match))

                    print('\tSite len: ' + str(len(seq_to_match)))


                    mot = motifMatcher(line[2], seq_to_match, line[3])

                    if mot is not None:
                        print('\n\n' + str(mot) + '\n\n')
                        count_list.append(mot)
                        line[3].append(mot)


                    print('\n')
                    print(line)
                    print('\n')
                    print(line[0:3])
                    print('\n')



                if len(line[3]) != 0:
                    for ele in line[3]:
                        print(ele)
                        if len(ele) == len(max(line[3], key=len)):
                            # print(line[2])
                            print(line[0]+'\t'+'\t'.join(line[1:2]) + '\t' + ele + '\n')
                            # time.sleep(10)
                            out_file2.write(line[0]+'\t'+'\t'.join(line[1:2]) + '\t' + ele + '\n')

        except FileNotFoundError:
            print('\tmemeSites not generated for ' + line[0] + '. May have JASPAR or Uniprobe sites.')
            continue
            # sys.exit()


out_file1.close()
out_file2.close()
print('\n\nAll outputs closed.\n\n')

sys.exit()
