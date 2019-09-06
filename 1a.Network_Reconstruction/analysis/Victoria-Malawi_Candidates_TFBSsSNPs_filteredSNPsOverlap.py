#!/usr/bin/env python

########################################################################################################################################################################################################
# This script does an intersect of the Lake Malawi and Lake Victoria VCF candidate overlaps (filtered by Will for overlap with promoter alignments) to find sites that segregate in both lakes

# How to run the script, by running on the terminal:
# python Victoria-Malawi_Candidates_TFBSsSNPs_filteredSNPsOverlap.py VCFFileA VCFFileB VCFFileA_speciesID outfile
# the script will output lines where a SNP is found to be divergent and segregating in the two VCFs compared, outputting lines of the second file (VCFFileB above)

# the script can run when FileA and FileB are interchanged (and needs to be run both ways) e.g python Victoria-Malawi_Candidates_TFBSsSNPs_filteredSNPsOverlap.py VCFFileB VCFFileA outfile
# this will output matching lines of VCFFileA
########################################################################################################################################################################################################

# cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/

# Columns in each file are:
#
# 0.	Pairwise_SNP_Ref_chr - Except, NOTE: this is the Pairwise_SNP_Alt_chr in [species]_pn1 cases e.g. nb1_pn1 - this is the chr for pn1
# 1.	Pairwise_SNP_Ref_start (0-based) - Except, NOTE: this is the Pairwise_SNP_Alt_start in [species]_pn1 cases e.g. nb1_pn1 - this is the start for pn1
# 2.	Pairwise_SNP_Ref_end (0-based) - Except, NOTE: this is the Pairwise_SNP_Alt_end in [species]_pn1 cases e.g. nb1_pn1 - this is the end for pn1 > this is the 1-based position of the SNP
# 3.	Pairwise_SNP_comparison_species
# 4.	Pairwise_SNP_quality - default:100
# 5.	Pairwise_SNP_strand
# 6.	Pairwise_SNP_Ref_allele
# 7.	Pairwise_SNP_Alt_allele
# 8.	Pairwise_SNP_Alt_chr - Except, NOTE: this is the Pairwise_SNP_Ref_chr in [species]_pn1 cases e.g. nb1_pn1 - this is the chr for nb1
# 9.	Pairwise_SNP_Alt_position (1-based) - Except, NOTE: this is the Pairwise_SNP_Ref_position in [species]_pn1 cases e.g. nb1_pn1 - this is the pos for nb1
# 10.	WG_motif_chr
# 11.	WG_motif_start
# 12.	WG_motif_end
# 13.	cichlid_geneID
# 14.	motif_genesymbolDr
# 15.	motif_genesymbolGa
# 16.	motif_genesymbolSp
# 17.	motif_strand
# 18.	TFmotif_ID
# 19.	TFmotif_cichlidID
# 20.	TFmotif_gene_symbol
# 21.	fimo_score
# 22.	fimo_pval
# 23.	fimo_qval
# 24.	motif_seq
# 25.	conf_level
# 26.	conf_score
# 27.	MalawiVCF_SNP_Ref_chr
# 28.	MalawiVCF_SNP_Ref_position
# 29.	MalawiVCF_SNP_Ref_ID
# 30.	MalawiVCF_SNP_Ref_allele
# 31.	MalawiVCF_SNP_Alt_allele
# 32.	MalawiVCF_SNP_Quality
# 33.	MalawiVCF_SNP_filter
# 34.	MalawiVCF_SNP_info
# 35.	MalawiVCF_SNP_format
# 36 to End of columns. genotype in Lake Malawi species

########################### create a base python version for the matching


#### 0. Load in the necessary modules
import numpy as np # this is for any statistical calculations that we may do
import sys # pass system arguments
import time # this is used for time delays
from more_itertools import unique_everseen # this is used for getting unique lines

if len(sys.argv) < 4:
    print('\nNot enough arguments entered.\nUsage: VCFFileA VCFFileB VCFFileA_speciesID outfile\n')
    sys.exit()

# 1. Open each file as sysargs and create a list of each row, and then each element as a list (list of a list)

# Using list comprehension, open and read each file line by line, adding each line to a list and each element of the first list to another list (list of list)
with open(sys.argv[1], 'r') as f:
    rowlistA = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list
# with open('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2', 'r') as f:
#     rowlistA = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list

with open(sys.argv[2], 'r') as f:
    rowlistB = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list
# with open('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed', 'r') as f:
#     rowlistB = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list

# # test dummy files created to run a test first
# with open('A-test', 'r') as f:
#     rowlistA = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list
#
# with open('B-test', 'r') as f:
#     rowlistB = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list

print('----------\nVCF FileA is:\t' + sys.argv[1])
print('VCF FileB is:\t' + sys.argv[2])
print('\n1. All files successfully loaded and turned into rows of lists and elements of lists (list of lists)\n----------')


# 2. Split the pairwise species column (col3) and store as two elements at the start of list of list - do for both files
# this is required for the matching in both files
rowlistAb=[] # create an empty list to append to on each loop
for line in rowlistA: # read each line
    rowlistA_sp = line[3] # take only the third element (species pairwise)
    rowlistAa = rowlistA_sp.split('_') # split the species pairwise
    # print(rowlistAa)
    rowlistAb.append(rowlistAa + line) # append the new species pairwise as separate elements to each line, maintaining the list of list structure by using '+'

rowlistBb=[] # create an empty list to append to on each loop
for line in rowlistB: # read each line
    rowlistB_sp = line[3] # take only the third element (species pairwise)
    rowlistBa = rowlistB_sp.split('_') # split the species pairwise
    # print(rowlistBa)
    rowlistBb.append(rowlistBa + line) # append the new species pairwise as separate elements to each line, maintaining the list of list structure by using '+'

print('\n2. Pairwise species column split and stored as two elements at start of list of list\n----------')

# This is the format of the list now - rowlistAb[0] accesses first line, [2] second and so on
# rowlistAb[0][0] == Ref_sp
# rowlistAb[0][1] == Alt_sp
# rowlistAb[0][2] == Ref_chr
# rowlistAb[0][3] == Ref_start
# rowlistAb[0][4] == Ref_end
# rowlistAb[0][5] == Pairwise_comparison
# rowlistAb[0][6] == Score
# rowlistAb[0][7] == SNP_strand
# rowlistAb[0][8] == Ref_allele
# rowlistAb[0][9] == Alt_allele
# rowlistAb[0][10] == Alt_chr
# rowlistAb[0][11] == Alt_position

# 3. loop over both new list of lists and do the matching, then append to another list - this will output the lines in FileA that overlap with those sites of FileB
### A KEY POINT HERE:
# 1. As an example for FileA if FileA is Mz-Malawi there are cases of e.g. mz11_on11 AND on11_mz11.
# 2. Where mz11 is the ref [0], the chr is [2], start is [3] and end is [4]
# 3. Where mz11 is the alt [1], the chr is also [2], start is [3] and end is [4] NOT chr [10] and pos [11]
# the comaprisons will be:
# pn_mz vs mz_ab
# pn_mz vs mz_nb
# pn_mz vs mz_on
# mz_pn vs pn_ab
# mz_pn vs pn_nb
# mz_pn vs pn_on

# at command line input, input the species (mz11, pn1, ab1, nb1, on11) of FileA as the third command line argument e.g. if FileA is Mz, input mz11 as third command line argument
matchingA_vs_B=[]
for x in rowlistAb:
    for y in rowlistBb:
        if x[0]==sys.argv[3] and y[1]==sys.argv[3] and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]:
            matchingA_vs_B.append(x)
        elif x[1]==sys.argv[3] and y[1]==sys.argv[3] and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]:
            matchingA_vs_B.append(x)
        elif x[0]==sys.argv[3] and y[0]==sys.argv[3] and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]:
            matchingA_vs_B.append(x)
        elif x[1]==sys.argv[3] and y[0]==sys.argv[3] and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]:
            matchingA_vs_B.append(x)
matchingA_vs_B_unique=list(unique_everseen(matchingA_vs_B)) # only keep unique lines

# # this is the test
# matchingA_vs_B=[]
# for x in rowlistAb:
#     for y in rowlistBb:
#         if x[0]=="mz11" and y[1]=="mz11" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]:
#             matchingA_vs_B.append(x)
#         elif x[1]=="mz11" and y[1]=="mz11" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]:
#             matchingA_vs_B.append(x)
#         elif x[0]=="mz11" and y[0]=="mz11" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]:
#             matchingA_vs_B.append(x)
#         elif x[1]=="mz11" and y[0]=="mz11" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]:
#             matchingA_vs_B.append(x)
# matchingA_vs_B_unique=list(unique_everseen(matchingA_vs_B)) # only keep unique lines
#
# x[0]=="mz11" and y[1]=="mz11" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]
# x[1]=="mz11" and y[1]=="mz11" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]
# x[0]=="mz11" and y[0]=="mz11" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]
# x[1]=="mz11" and y[0]=="mz11" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]
#
# x[0]=="pn1" and y[1]=="pn1" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]
# x[1]=="pn1" and y[1]=="pn1" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]
# x[0]=="pn1" and y[0]=="pn1" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]
# x[1]=="pn1" and y[0]=="pn1" and x[2]==y[10] and y[11]>=x[3] and y[11]<=x[4]

print('\n3. Looped over both list of lists and matching of FileA and FileB DONE\n----------')

# 4. Write the nested list out to a tab-separated file, with trailing newlines
with open(sys.argv[4], 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in matchingA_vs_B_unique) # join each element as tab-separated and trailing newlines for each line in the list
file.close() # close the file

print('\n4. Matched lines of FileA successfully outputted to:\t' + sys.argv[4])
print('\n5. The species that was assessed was:\t' + sys.argv[3])

####################################################################################################################################

########### Can also run the intersect using pandas but would recommend against this as certain functions are acting oddly when matching
# the script works, only thing that needs amending is the format of the output file

# #### 0. Load in the necessary modules
# import pandas as pd # pandas is the module for manipulating dfs
# import numpy as np # this is for any statistical calculations that we may do
# import sys # pass system arguments
# import sys
#
# if len(sys.argv) < 3:
#     print('\nNot enough arguments entered.\nUsage: VCFFileA VCFFileB outfile\n')
#     sys.exit()
#
# # 1. Read, by providing as sys args, the mz/pn file
# fileA = pd.read_csv(sys.argv[1], sep='\t', header = None) # e.g. mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2
# # fileA = pd.read_csv('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2', sep='\t', header = None)
# # fileA = pd.read_csv('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed', sep='\t', header = None)
#
# # 2. Read, by providing as sys args, the pn/mz file
# fileB = pd.read_csv(sys.argv[2], sep='\t', header = None) # e.g. pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed
# # fileB = pd.read_csv('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed', sep='\t', header = None)
# # fileB = pd.read_csv('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2', sep='\t', header = None)
#
# print('----------\nVCF FileA is:\t' + sys.argv[1])
# print('VCF FileB is:\t' + sys.argv[2])
# print('\n1. All files successfully loaded as pandas\n----------')
#
#
# # 3. Convert all inputs to a dataframe for splitting and joining columns
# fileADF = pd.DataFrame(data=fileA)
# fileBDF = pd.DataFrame(data=fileB)
#
# print('\n2. All pandas successfully converted to dataframes for splitting and joining columns\n----------')
#
# # 4. Define a panda function to do the matching:
# # For this, do all the following intersects of coordinates, then sort -k1,1 -k2,2n to identify presence in all species comparisons
# # Compare like so:
# # FileA = PnVCF
# # FileB = MzVCF
# # 1. Pairwise_SNP_comparison_species: pn_mz vs mz_ab (.split('_')[2] == .split('_')[1]) & 2. Pairwise_SNP_Alt_chr (Mz in FileA) == Pairwise_SNP_Ref_chr (Mz in FileB) & 3. Pairwise_SNP_Alt_position (FileA) >= Pairwise_SNP_Ref_start & <= Pairwise_SNP_Ref_end (FileB)
# # 1. Pairwise_SNP_comparison_species: pn_mz vs mz_nb (.split('_')[2] == .split('_')[1]) & 2. Pairwise_SNP_Alt_chr (Mz in FileA) == Pairwise_SNP_Ref_chr (Mz in FileB) & 3. Pairwise_SNP_Alt_position (FileA) >= Pairwise_SNP_Ref_start & <= Pairwise_SNP_Ref_end (FileB)
# # 1. Pairwise_SNP_comparison_species: pn_mz vs mz_on (.split('_')[2] == .split('_')[1]) & 2. Pairwise_SNP_Alt_chr (Mz in FileA) == Pairwise_SNP_Ref_chr (Mz in FileB) & 3. Pairwise_SNP_Alt_position (FileA) >= Pairwise_SNP_Ref_start & <= Pairwise_SNP_Ref_end (FileB)
# # 1. Pairwise_SNP_comparison_species: mz_pn vs pn_ab (.split('_')[2] == .split('_')[1]) & 2. Pairwise_SNP_Alt_chr (Pn in FileB) == Pairwise_SNP_Ref_chr (Pn in FileA) & 3. Pairwise_SNP_Alt_position (FileB) >= Pairwise_SNP_Ref_start & <= Pairwise_SNP_Ref_end (FileA)
# # 1. Pairwise_SNP_comparison_species: mz_pn vs pn_nb (.split('_')[2] == .split('_')[1]) & 2. Pairwise_SNP_Alt_chr (Pn in FileB) == Pairwise_SNP_Ref_chr (Pn in FileA) & 3. Pairwise_SNP_Alt_position (FileB) >= Pairwise_SNP_Ref_start & <= Pairwise_SNP_Ref_end (FileA)
# # 1. Pairwise_SNP_comparison_species: mz_pn vs pn_on (.split('_')[2] == .split('_')[1]) & 2. Pairwise_SNP_Alt_chr (Pn in FileB) == Pairwise_SNP_Ref_chr (Pn in FileA) & 3. Pairwise_SNP_Alt_position (FileB) >= Pairwise_SNP_Ref_start & <= Pairwise_SNP_Ref_end (FileA)
#
# # 5. split the Pairwise_SNP_comparison_species column of each dataframe
# fileADF2 = fileADF.iloc[:,3].str.split('_', expand=True) # this splits the column
# fileBDF2 = fileBDF.iloc[:,3].str.split('_', expand=True) # this splits the column
#
# print('\n3. All split the Pairwise_SNP_comparison_species column of each dataframe successfully DONE\n----------')
#
# # 6. join this on to the rest of the dataframe
# fileADF3 = fileADF2.join(fileADF, lsuffix='_sp', rsuffix='_rest')
# fileBDF3 = fileBDF2.join(fileBDF, lsuffix='_sp', rsuffix='_rest')
#
# print('\n4. All joining of splits above to the rest of the dataframe successfully DONE\n----------')
#
# # # do a test for the matching - created some dummy files where there will be some and no matches
# # mztest = pd.read_csv("mz.test", sep='\t') # read in the file
# # pntest = pd.read_csv("pn.test", sep='\t') # read in the file
# #
# # # this is just a test - prepared in such a way that mz vs pn and pn vs mz comparison can be done by swapping the input files of vcf_row and other_vcf
# # vcf_row = mztest
# # otherVCF = pntest
# # print(vcf_row)
# # print(otherVCF)
# # matching_rows = otherVCF[
# #        (vcf_row['1_sp'] == otherVCF['0_sp']) & # match the comparison species in fileA to the ref species in fileB
# #        (vcf_row['8'] == otherVCF['0_rest']) & # match the Pairwise_SNP_Alt_chr in fileA to Pairwise_SNP_Ref_chr in fileB
# #        (vcf_row['9'] >= otherVCF['1_rest']) & # ensure that the Pairwise_SNP_Alt_position in fileA is >= Pairwise_SNP_Ref_start in fileB, and
# #        (vcf_row['9'] <= otherVCF['2']) # <= Pairwise_SNP_Ref_end in fileB
# # ]
# # print(matching_rows)
#
# print('\n5. Now creating the function for doing the matching......\n----------')
#
# # 7. Carry out the matching
#
# # this is just a test
# # testA = pd.read_csv('A-test', sep='\t', header=None)
# # testB = pd.read_csv('B-test', sep='\t', header=None)
# # testADF = pd.DataFrame(data=testA)
# # testBDF = pd.DataFrame(data=testB)
# # testADF2 = testADF.iloc[:,3].str.split('_', expand=True) # this splits the column
# # testBDF2 = testBDF.iloc[:,3].str.split('_', expand=True) # this splits the column
# # testADF3 = testADF2.join(testADF, lsuffix='_sp', rsuffix='_rest')
# # testBDF3 = testBDF2.join(testBDF, lsuffix='_sp', rsuffix='_rest')
# #
# # # this works but need to apply as a function to run on DFs with different number of rows..
# # vcf_row = testADF3.iloc[0,:] # this only applies on row0
# # print(vcf_row)
# # matching_rows = testBDF3[
# # (vcf_row['1_sp'] == testBDF3['0_sp']) & # match the comparison species in fileA to the ref species in fileB
# # (vcf_row[8] == testBDF3['0_rest']) & # match the Pairwise_SNP_Alt_chr in fileA to Pairwise_SNP_Ref_chr in fileB
# # (vcf_row[9] >= testBDF3['1_rest']) & # ensure that the Pairwise_SNP_Alt_position in fileA is >= Pairwise_SNP_Ref_start in fileB, and
# # (vcf_row[9] <= testBDF3[2]) # <= Pairwise_SNP_Ref_end in fileB
# # ]
# # print(matching_rows)
# #
# # def find_otherVCFtest(vcfrow):
# #     # for a single region in the liftover, find matching rows in the VCF file
# #     matching_rows2 = testBDF3[
# #     (vcfrow['1_sp'] == testBDF3['0_sp']) & # match the comparison species in fileA to the ref species in fileB
# #     (vcfrow[8] == testBDF3['0_rest']) & # match the Pairwise_SNP_Alt_chr in fileA to Pairwise_SNP_Ref_chr in fileB
# #     (vcfrow[9] >= testBDF3['1_rest']) & # ensure that the Pairwise_SNP_Alt_position in fileA is >= Pairwise_SNP_Ref_start in fileB, and
# #     (vcfrow[9] <= testBDF3[2]) # <= Pairwise_SNP_Ref_end in fileB
# #     ]
# #     # return the rows that adhere to the above conditionals
# #     # print(matching_rows2)
# #     # print('')
# #     # print(matching_rows2.shape[0])
# #     if matching_rows2.shape[0] > 0: # as pandas apply applies it to the whole dataframe, in instance where a row does not match, it returns an empty dataframe and so, we want to select where the row (shape) is not empty >0
# #         return matching_rows2.iloc[0,:] # this will return all the rows, including those that do not match as empty frames
# #
# # # this will apply the function above
# # intersectVCF = testADF3.apply(
# #     find_otherVCFtest,
# #     axis=1) # this will return all the rows, including those that do not match as empty frames
# #
# # intersectVCF2 = intersectVCF.dropna() # this will remove the empty rows
#
#
# # this is the function for the larger frame
# def find_otherVCF(vcfrow):
#     # for a single region in the liftover, find matching rows in the VCF file
#     matching_rows = fileBDF3[
#     (vcfrow['1_sp'] == fileBDF3['0_sp']) & # match the comparison species in fileA to the ref species in fileB
#     (vcfrow[8] == fileBDF3['0_rest']) & # match the Pairwise_SNP_Alt_chr in fileA to Pairwise_SNP_Ref_chr in fileB
#     (vcfrow[9] >= fileBDF3['1_rest']) & # ensure that the Pairwise_SNP_Alt_position in fileA is >= Pairwise_SNP_Ref_start in fileB, and
#     (vcfrow[9] <= fileBDF3[2]) # <= Pairwise_SNP_Ref_end in fileB
#     ]
#     # return the rows that adhere to the above conditionals
#     # print(matching_rows)
#     # print('')
#     # print(matching_rows.shape[0])
#     if matching_rows.shape[0] > 0: # as pandas apply applies it to the whole dataframe, in instance where a row does not match, it returns an empty dataframe and so, we want to select where the row (shape) is not empty >0
#         return matching_rows.iloc[0,:] # this will return all the rows, including those that do not match as empty frames
#
# print('\n6. Function for matching being applied......\n----------')
#
# # This applies the function but again, did not work ?!?!
# intersectVCF = fileADF3.apply(
#     find_otherVCF,
#     axis=1) # this will return all the rows, including those that do not match as empty frames
# #axis is applying the function to either 0 = cols, or here 1 = rows
# intersectVCF2 = intersectVCF.dropna() # this will remove the empty rows
#
# print('\n7. Matching DONE > file will be outputted as:\t' + sys.argv[3])
#
# # write out to file
# intersectVCF2.to_csv(sys.argv[3], encoding='utf-8', sep='\t', index=False)
