#!/usr/bin/env python

####################### create a base python version for amending species genotypes to actual alleles in the filtered SNP bed files

#### 0. Load in the necessary modules
import numpy as np # this is for any statistical calculations that we may do
import sys # pass system arguments
import time # this is used for time delays
from more_itertools import unique_everseen # this is used for getting unique lines
import re # regex matches


#### 1. Open each file as sysargs and create a list of each row, and then each element as a list (list of a list)

# Using list comprehension, open and read each file line by line, adding each line to a list and each element of the first list to another list (list of list)
with open('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2', 'r') as f:
    mz = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list

with open('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed', 'r') as f:
    pn = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list

#### 2. loop over list of lists and amend elements [36] to the end (species genotypes) using ref [30] and alt [31] to convert the genotypes to alleles
# how to interpret each species column:
# / - non-phased
# | - phased
# Example - if Ref = C and Alt = A
# 0/0 - the sample is homozygous reference = C/C
# 0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles = C/A
# 1/0 - the sample is heterozygous, carrying 1 copy of each of the ALT and REF alleles = A/C
# 1/1 - the sample is homozygous alternate = A/A
# create the different genotypes
homorefphase='0|0'
homorefnophase='0/0'
heterophase='0|1'
heteronophase='0/1'
heterophase2='1|0'
heteronophase2='1/0'
homoaltphase='1|1'
homoaltnophase='1/1'

startcol = 36
mznumcols = len(mz[1]) # get the total number of columns in each file
pnnumcols = len(pn[1]) # get the total number of columns in each file

# mz is phased so uses '|' but put an or statement of '|' or '/'
# create a counter runner then loop over that and assess and restructure each element to the actual allelic genotype by adding each to new list of lists
mzcounter=list(range(startcol, mznumcols))
mzlistcounter=list(range(0,len(mzcounter)))
mzspeciesgenotype=[]
for (x,y) in zip(mzcounter,mzlistcounter): # iterate over the two lists at the same time so not repeated
    mzspeciesgenotype.append([]) # create lists of lists that can be appended to for each column
    # print(x,y)
    for z in mz:
        # print(x,z[x])
        # time.sleep(1)
        if z[x].startswith(homorefphase) or z[x].startswith(homorefnophase):
            # print(z[36] + '\t' + "I'M A HOMO REF" + '\t' + z[30] + "|" + z[30])
            mzspeciesgenotype[y].append(z[30] + "|" + z[30])
        elif z[x].startswith(heterophase) or z[x].startswith(heteronophase):
            # print(z[36] + '\t' + "I'M A HETERO" + '\t' + z[30] + "|" + z[31])
            mzspeciesgenotype[y].append(z[30] + "|" + z[31])
        elif z[x].startswith(heterophase2) or z[x].startswith(heteronophase2):
            # print(z[36] + '\t' + "I'M A HETERO2" + '\t' + z[31] + "|" + z[30])
            mzspeciesgenotype[y].append(z[31] + "|" + z[30])
        elif z[x].startswith(homoaltphase) or z[x].startswith(homoaltnophase):
            # print(z[36] + '\t' + "I'M A HOMO ALT" + '\t' + z[31] + "|" + z[31])
            mzspeciesgenotype[y].append(z[31] + "|" + z[31])
        else:
            mzspeciesgenotype[y].append("NA")
# join the allelic genotypes onto the end of the original list table
# to each element [0-134] of mzspeciesgenotype, join each line [0-5712] of the original list e.g. mz[0]
mzspeciesgenotype2=[list(x) for x in zip(*mzspeciesgenotype)] # transpose the list so that it is 5712 lines and 135 columns along (from the 135 lines and 5712 columns along)
mzcounter2=list(range(0, len(mz))) #0-5711
mznew = [a + b for a, b in zip(mz,mzspeciesgenotype2)] # join the two lists

# then, amend the reference and alternative alleles so that they are represented as homozygous
for x in mznew: # reference alleles
    if x[30]=='A':
        x[30]='A|A'
    elif x[30]=='T':
        x[30]='T|T'
    elif x[30]=='C':
        x[30]='C|C'
    elif x[30]=='G':
        x[30]='G|G'

for x in mznew: # alternative alleles
    if x[31]=='A':
        x[31]='A|A'
    elif x[31]=='T':
        x[31]='T|T'
    elif x[31]=='C':
        x[31]='C|C'
    elif x[31]=='G':
        x[31]='G|G'

# PN is non-phased so uses '/' but put an or statement of '|' or '/'
# create a counter runner then loop over that and assess and restructure each element to the actual allelic genotype by adding each to new list of lists
pncounter=list(range(startcol, pnnumcols))
pnlistcounter=list(range(0,len(pncounter)))
pnspeciesgenotype=[]
for (x,y) in zip(pncounter,pnlistcounter): # iterate over the two lists at the same time so not repeated
    pnspeciesgenotype.append([]) # create lists of lists that can be appened to for each column
    # print(x,y)
    for z in pn:
        # print(x,z[x])
        # time.sleep(1)
        if z[x].startswith(homorefphase) or z[x].startswith(homorefnophase):
            # print(z[36] + '\t' + "I'M A HOMO REF" + '\t' + z[30] + "|" + z[30])
            pnspeciesgenotype[y].append(z[30] + "/" + z[30])
        elif z[x].startswith(heterophase) or z[x].startswith(heteronophase):
            # print(z[36] + '\t' + "I'M A HETERO" + '\t' + z[30] + "|" + z[31])
            pnspeciesgenotype[y].append(z[30] + "/" + z[31])
        elif z[x].startswith(heterophase2) or z[x].startswith(heteronophase2):
            # print(z[36] + '\t' + "I'M A HETERO2" + '\t' + z[31] + "|" + z[30])
            pnspeciesgenotype[y].append(z[31] + "|" + z[30])
        elif z[x].startswith(homoaltphase) or z[x].startswith(homoaltnophase):
            # print(z[36] + '\t' + "I'M A HOMO ALT" + '\t' + z[31] + "|" + z[31])
            pnspeciesgenotype[y].append(z[31] + "/" + z[31])
        else:
            pnspeciesgenotype[y].append("NA")
# join the allelic genotypes onto the end of the original list table
# to each element of pnspeciesgenotype, join each line of the original list e.g. pn[0]
pnspeciesgenotype2=[list(x) for x in zip(*pnspeciesgenotype)] # transpose the list
pncounter2=list(range(0, len(pn))) #0-17014
pnnew = [a + b for a, b in zip(pn,pnspeciesgenotype2)] # join the two lists

# then, amend the reference and alternative alleles so that they are represented as homozygous
for x in pnnew: # reference alleles
    if x[30]=='A':
        x[30]='A/A'
    elif x[30]=='T':
        x[30]='T/T'
    elif x[30]=='C':
        x[30]='C/C'
    elif x[30]=='G':
        x[30]='G/G'

for x in pnnew: # alternative alleles
    if x[31]=='A':
        x[31]='A/A'
    elif x[31]=='T':
        x[31]='T/T'
    elif x[31]=='C':
        x[31]='C/C'
    elif x[31]=='G':
        x[31]='G/G'

#### 3. Write out the new lists to tab-separated files, with trailing newlines
with open('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2', 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in mznew) # join each element as tab-separated and trailing newlines for each line in the list
file.close() # close the file

with open('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2', 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in pnnew) # join each element as tab-separated and trailing newlines for each line in the list
file.close() # close the file

####################################################################################################################################
