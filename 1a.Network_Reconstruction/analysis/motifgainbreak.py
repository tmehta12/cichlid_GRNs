#!/usr/bin/env python

########################################################################################################################################################################################################
## 1.	Categorize TFBS mutations as breaking/gaining

# 1. replace all SNPs overlapping motifs in the promoter with the ALT genotype, then
# 2. DONE in /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/SNP_overlap.sh: re-scan for motifs (this tests whether the SNP breaks the site as flanking is unchanged) and compare the p-value for the original and mutated hit, then
# 3. DONE in /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/motifgainbreak_pt2.py: Categorise according to motif gain (from ancestor REF SNP), loss (from ancestor from REF SNP)


# strategy
# first file are your fasta files with >gene IDs followed by promoter sequence
# second file is a tabulated index file with all SNP details
# then:
	# a. read fasta and index file: match ID
	# b. read along fasta seq and stop when get to index col1 position in promoter
	# c. read col2 and check that position nt in promoter seq matches REF
	# d. read col3 and change that position nt in promoter seq to ALT
		# c and d. if fasta_pos=index_col2, fasta_pos=index_col3
########################################################################################################################################################################################################

# # Files to use
# /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta
# /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta

# /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap
# /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/*_*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt


## Test files:
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed

# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed

# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/mz11_pn1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/pn1_mz11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt


########################################################################################################################################################################################################

#### 0. Load in the necessary modules
import numpy as np # this is for any statistical calculations that we may do
import csv
import sys # pass system arguments
import time # this is used for time delays
import re # this does string substitution
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from more_itertools import unique_everseen # this is used for getting unique lines
# from tqdm import tqdm # this is to time each mapping

# defined a function to read Fasta file and store as a dictionary (this is far better than SeqIO)
def readFastaGenerator(seqin):
    for line in seqin:
        yield line.strip()[1:], next(seqin).strip()

# an argument for whether the number of command input files is too short
if len(sys.argv) < 4:
    print('\nNot enough arguments entered.\nUsage: 1. Promoter_FASTA  2. SNP_VCF_file 3. Promoter_BED 4. Outfile e.g. MzProm_pn1SNP.fasta\n')
    sys.exit()

print('----------\nOriginal 5kb promoter FASTA file is:\t' + sys.argv[1])
print('----------\nSNP VCF file is:\t' + sys.argv[2])
print('----------\nOriginal 5kb promoter BED file is:\t' + sys.argv[3])

#### 1. Open each file as sysargs and create dictionaries of each

# A. Dictionary of Promoter FASTA file using function above where geneID is key and sequence is value
with open(sys.argv[1]) as prom_fasta:
    prom_dict = {}
    pf = readFastaGenerator(prom_fasta)
    for _id, _seq in pf:
        prom_dict[_id] = _seq
print('----------\nPromoter FASTA loaded')

# B. Dictionary of the SNPs where each geneID is the key and every SNP line is the value
with open(sys.argv[2], 'r') as f:
    SNPdict = {}
    for line in csv.reader(f, delimiter = '\t'):
        if line[16] in SNPdict:
            SNPdict[line[16]].append(line)
        else:
            SNPdict[line[16]] = [line]
    # print(SNPdict)
print('----------\nSNPs loaded')

# with open('test.txt', 'r') as f:
#     SNPdict = {}
#     for line in csv.reader(f, delimiter = '\t'):
#         if line[16] in SNPdict:
#             SNPdict[line[16]].append(line)
#         else:
#             SNPdict[line[16]] = [line]


# load in the promoter bed so as to work out the SNP position in the promoter sequence
with open(sys.argv[3], 'r') as f:
    # Prom_BED = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list
    PromBED = {}
    for line in csv.reader(f, delimiter = '\t'):
        PromBED[line[4]] = line
    # print(PromBED)
    # PromBED = {gene_id:line.strip('\n').split('\t')}
print('----------\nPromoter BED file loaded')

## Test files
# mz_prom = open("Mz_Prom_SAMPLE.fasta")
# mz_prom_dict = SeqIO.to_dict(SeqIO.parse(mz_prom, "fasta"))
# with open('mz11_pn1_subst_SAMPLE.txt', 'r') as f:
#     mz11_pn1_SNPs = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list
# # load in the promoter bed so as to work out the SNP position in the promoter sequence
# with open('Mz_PromBED.txt', 'r') as f:
#     Mz_Prom_BED = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list

print('\n1. All files successfully loaded and turned into dictionaries\n----------')

print('\n2. Looping over dictionaries and now running the matching strategy\n----------')

#### 2. Carry out the matching strategy as so - this runs ~5000 SNP lines using 2GB memory in <5 seconds

with open(sys.argv[4], 'w') as outfile:
    for m, seq in prom_dict.items(): # iterate through the key (geneID) and value (seq) of the prom dict
        if m in SNPdict: # if the geneID is in the SNPdict
            for snp in SNPdict[m]: # and for each instance/line (snp line) of the geneID in the prom dict
                outfile.write('>' + m + "_SNPid:" + snp[0] + '\n') # create a new file that first appends a header line e.g. >ab.gene.s1.1_SNPid:1
                outfile.write(seq[0:(int(snp[2]) - int(PromBED[m][1]))] + snp[6] + seq[(int(snp[3]) - int(PromBED[m][2])):len(seq)] + '\n') # then in the line line, output the geneID promoter seqence up to the point of the SNP, then the SNP nucleotide, then the rest of the sequence (this is clearer in the extended examples below)

## to time with tqdm, use below
# with open(sys.argv[4], 'w') as outfile:
#     for m, seq in tqdm(prom_dict.items()):
#         if m in SNPdict:
#             for snp in tqdm(SNPdict[m]):
#                 outfile.write('>' + m + "_SNPid:" + snp[0] + '\n')
#                 outfile.write(seq[0:(int(snp[2]) - int(PromBED[m][1]))] + snp[6] + seq[(int(snp[3]) - int(PromBED[m][2])):len(seq)] + '\n')

# print(SNPs_prom)
print('\n3. Matching strategy DONE and written to file\n----------')

print('\n4. File successfully created and outputted to:\t' + sys.argv[4])

# Faster version of below - takes approx. 4 mins for each SNP row using 4GB memory

# prom = open(sys.argv[1])
# prom_dict = SeqIO.to_dict(SeqIO.parse(prom, "fasta"))
# Using list comprehension, open and read each file line by line, adding each line to a list and each element of the first list to another list (list of list)
# with open(sys.argv[2], 'r') as f:
    # SNPs = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list
# with open(sys.argv[3], 'r') as f:
    # Prom_BED = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list
#
# SNPs_prom=[]
# for m in prom_dict: # loop through each entry in the promoter fasta file dictionary
#     for b in Prom_BED: # loop through the promoter BED file
#         for row in SNPs: # loop through each row of the SNPs file
#             if row[16] == b[4] and prom_dict[m].id == row[16]: # a. perform function if all gene IDs in FASTA, BED and SNP files match
#                 SNPs_prom.append('>' + prom_dict[m].id + "_SNPid:" + row[0])
#                 SNPs_prom.append(prom_dict[m].seq[0:(int(row[2]) - int(b[1]))] + row[6] + prom_dict[m].seq[(int(row[3]) - int(b[2])):len(prom_dict[m].seq)])



# Below is a slow version and takes approx. 8 mins for each SNP row using 4GB memory
# a. match ID in fasta and index file
# b. read along fasta seq and stop when get to index col1 position in promoter
	# you can use the SNP pos to select a pos in the fasta list e.g. y[x[1]] > this needs to be 1-based to select position.
# c. read col2 and check that position nt in promoter seq matches REF
# d. read col3 and change that position nt in promoter seq to ALT
		# c and d. if fasta_pos=index_col2, fasta_pos=index_col3

# SNPs_prom=[]
# for m in prom_dict: # loop through each entry in the promoter fasta file dictionary
#     sequence_key1 = (prom_dict[m].id) # a. assign sequence ID in fasta to sequence_key1
#     # print(sequence_key1)
#     sequence = (prom_dict[m].seq) # b. assign sequences in fasta to sequence key
#     # print(sequence)
#     # time.sleep(10)
#     # sequence = sequence.tomutable()
#     for b in Prom_BED: # loop through the promoter BED file
#         sequence_key2 = b[4] # c. assign promoter gene ID in BED file to sequence_key2
#         prom_start = b[1] # d. assign promoter start pos to key
#         prom_end = b[2] # e. assign promoter end pos to key
#         for row in SNPs: # loop through each row of the SNPs file
#             sequence_key = row[16] # f. assign gene ID of overlapping SNP to sequence_key
#             # print(sequence_key)
#             if sequence_key == sequence_key2 and sequence_key1 == sequence_key: # g. perform function if all gene IDs in FASTA, BED and SNP files match
#                 # motif_start = row[12]
#                 SNP_start = row[2]
#                 SNP_end = row[3]
#                 SNP_pos_start = int(SNP_start) - int(prom_start)
#                 SNP_pos_end = int(SNP_end) - int(prom_start)
#                 SNP_alt = row[6]
#                 SNP_ID = row[0]
#                 # print(row, "\n------- SNP pos =", SNP_pos_start,"-",SNP_pos_end, "\n------- SNP nucleotide is:", str(prom_dict[sequence_key1][SNP_pos_start:SNP_pos_end].seq),"\n------- Gene promoter is:", sequence_key1,"\nSummary: Position",SNP_pos_start,"-",SNP_pos_end,"is the nucleotide","-",str(prom_dict[sequence_key1][SNP_pos_start:SNP_pos_end].seq),"-","in",sequence_key1,"gene promoter to be replaced with the alternative SNP nucleotide -",SNP_alt, "\n------- Gene promoter fasta is:","\n",sequence,"\n")
#                 # time.sleep(10)
#                 # print(row, "\n------- SNP pos =", SNP_pos_start,"-",SNP_pos_end, "\n------- SNP nucleotide is:", str(prom_dict[sequence_key1][SNP_pos_start:SNP_pos_end].seq),"\n------- Gene promoter is:", sequence_key1,"\nSummary: Position",SNP_pos_start,"-",SNP_pos_end,"is the nucleotide","-",str(prom_dict[sequence_key1][SNP_pos_start:SNP_pos_end].seq),"-","in",sequence_key1,"gene promoter to be replaced with the alternative SNP nucleotide -",SNP_alt,"\n")
#                 # time.sleep(10)
#                 # print("original sequence is:\n" + prom_dict[m].format("fasta"))
#                 # print("amended sequence is:\n" + '>' + prom_dict[m].id)
#                 # print(sequence[0:SNP_pos_start] + "--" + SNP_alt + "--" + sequence[SNP_pos_end:len(sequence)])
#                 # time.sleep(5)
#                 SNPs_prom.append('>' + prom_dict[m].id + "_SNPid:" + SNP_ID)
#                 SNPs_prom.append(sequence[0:SNP_pos_start] + SNP_alt + sequence[SNP_pos_end:len(sequence)])

## For the test files
# mz11_pn1_SNPs_prom=[]
# for m in mz_prom_dict:
#     sequence_key1 = (mz_prom_dict[m].id) # a. match ID in fasta and index file
#     # print(sequence_key1)
#     sequence = (mz_prom_dict[m].seq) # pull out just the sequences to amend
#     # print(sequence)
#     # time.sleep(10)
#     # sequence = sequence.tomutable()
#     for b in Mz_Prom_BED:
#         sequence_key2 = b[4] # a. match ID in fasta and index file
#         prom_start = b[1]
#         prom_end = b[2]
#         for row in mz11_pn1_SNPs:
#             sequence_key = row[16] # a. match ID in fasta and index file
#             # print(sequence_key)
#             if sequence_key == sequence_key2 and sequence_key1 == sequence_key:
#                 # motif_start = row[12]
#                 SNP_start = row[2]
#                 SNP_end = row[3]
#                 SNP_pos_start = int(SNP_start) - int(prom_start)
#                 SNP_pos_end = int(SNP_end) - int(prom_start)
#                 SNP_alt = row[6]
#                 SNP_ID = row[0]
#                 # print(row, "\n------- SNP pos =", SNP_pos_start,"-",SNP_pos_end, "\n------- SNP nucleotide is:", str(mz_prom_dict[sequence_key1][SNP_pos_start:SNP_pos_end].seq),"\n------- Gene promoter is:", sequence_key1,"\nSummary: Position",SNP_pos_start,"-",SNP_pos_end,"is the nucleotide","-",str(mz_prom_dict[sequence_key1][SNP_pos_start:SNP_pos_end].seq),"-","in",sequence_key1,"gene promoter to be replaced with the alternative SNP nucleotide -",SNP_alt, "\n------- Gene promoter fasta is:","\n",sequence,"\n")
#                 # time.sleep(10)
#                 # print(row, "\n------- SNP pos =", SNP_pos_start,"-",SNP_pos_end, "\n------- SNP nucleotide is:", str(mz_prom_dict[sequence_key1][SNP_pos_start:SNP_pos_end].seq),"\n------- Gene promoter is:", sequence_key1,"\nSummary: Position",SNP_pos_start,"-",SNP_pos_end,"is the nucleotide","-",str(mz_prom_dict[sequence_key1][SNP_pos_start:SNP_pos_end].seq),"-","in",sequence_key1,"gene promoter to be replaced with the alternative SNP nucleotide -",SNP_alt,"\n")
#                 # time.sleep(10)
#                 # print("original sequence is:\n" + mz_prom_dict[m].format("fasta"))
#                 # print("amended sequence is:\n" + '>' + mz_prom_dict[m].id)
#                 # print(sequence[0:SNP_pos_start] + "--" + SNP_alt + "--" + sequence[SNP_pos_end:len(sequence)])
#                 # time.sleep(5)
#                 mz11_pn1_SNPs_prom.append('>' + mz_prom_dict[m].id + "_SNPid:" + SNP_ID)
#                 mz11_pn1_SNPs_prom.append(sequence[0:SNP_pos_start] + SNP_alt + sequence[SNP_pos_end:len(sequence)])

# print('\n3. Matching strategy DONE, now writing to file\n----------')

# 3. Write the nested list out to a tab-separated file, with trailing newlines

## Test output:
# with open("test.fasta", 'w') as file: # open a file for writing
#     file.writelines(''.join(i) + '\n' for i in mz11_pn1_SNPs_prom) # join each element as tab-separated and trailing newlines for each line in the list
# file.close() # close the file

# with open(sys.argv[4], 'w') as file: # open a file for writing
#     file.writelines(''.join(i) + '\n' for i in SNPs_prom) # join each element as tab-separated and trailing newlines for each line in the list
# file.close() # close the file

# print('\n4. File successfully created and outputted to:\t' + sys.argv[4])
