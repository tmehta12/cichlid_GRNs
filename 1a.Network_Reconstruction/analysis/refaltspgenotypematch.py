#!/usr/bin/env python

####################### create a base python version for amending species genotypes to actual alleles in the filtered SNP bed files

#### 0. Load in the necessary modules
import numpy as np # this is for any statistical calculations that we may do
import sys # pass system arguments
import time # this is used for time delays
from more_itertools import unique_everseen # this is used for getting unique lines
import re # regex matches

#### 1. Open each file and create a list of each row, and then each element as a list (list of a list)

# Using list comprehension, open and read each file line by line, adding each line to a list and each element of the first list to another list (list of list)
with open('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2', 'r') as f:
    mz_malawiGT = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list

with open('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2', 'r') as f:
    pn_victoriaGT = [line.strip('\n').split('\t') for line in f] # open the file as f, and for each line, strip by newline, split by tab and assign to list

#### 2. loop over list of lists, create two tables for each - one for the ref [30] and another for the alt [31] comparison where, in each lake species column, input 1 for same, 0 for different

# M. zebra - create a counter runner then loop over that and assess and restructure each element to the actual allelic genotype by adding each to new list of lists
startcol = 171
mz_malawiGTnumcols=len(mz_malawiGT[1])
mzcounter3=list(range(startcol, mz_malawiGTnumcols)) #cols 171-305
mz_malawiGTlistcounter=list(range(0,len(mzcounter3)))

# Mz Ref Match
mzrefmatch=[]
for (x,y) in zip(mzcounter3,mz_malawiGTlistcounter): # iterate over the two lists at the same time so not repeated
    mzrefmatch.append([]) # create lists of lists that can be appended to for each column
    # print(x)
    # print(y)
    for z in mz_malawiGT:
        if z[30]==z[x]:
            mzrefmatch[y].append("1")
        else:
            mzrefmatch[y].append("0")
# join the matches onto the end of the original list table
# to each element of mzrefmatch, join each line of the original list
mzrefmatch2=[list(x) for x in zip(*mzrefmatch)] # transpose the list so that it is 5712 lines and 135 columns along (from the 135 lines and 5712 columns along)
mzcounter4=list(range(0, len(mz_malawiGT))) #0-5711
mz_malawiGTnew1 = [a + b for a, b in zip(mz_malawiGT,mzrefmatch2)] # join the two lists

# Mz Alt Match
mzaltmatch=[]
for (x,y) in zip(mzcounter3,mz_malawiGTlistcounter): # iterate over the two lists at the same time so not repeated
    mzaltmatch.append([]) # create lists of lists that can be appended to for each column
    for z in mz_malawiGT:
        if z[31]==z[x]:
            mzaltmatch[y].append("1")
        else:
            mzaltmatch[y].append("0")
# join the matches onto the end of the original list table
# to each element of mzrefmatch, join each line of the original list
mzaltmatch2=[list(x) for x in zip(*mzaltmatch)] # transpose the list so that it is 5712 lines and 135 columns along (from the 135 lines and 5712 columns along)
mzcounter4=list(range(0, len(mz_malawiGT))) #0-5711
mz_malawiGTnew2 = [a + b for a, b in zip(mz_malawiGT,mzaltmatch2)] # join the two lists

# P. nyererei - create a counter runner then loop over that and assess and restructure each element to the actual allelic genotype by adding each to new list of lists
startcol = 139
pn_victoriaGTnumcols = len(pn_victoriaGT[1]) # get the total number of columns in each file
pncounter3=list(range(startcol, pn_victoriaGTnumcols)) #139-241
pn_victoriaGTlistcounter=list(range(0,len(pncounter3)))

# Pn Ref Match
pnrefmatch=[]
for (x,y) in zip(pncounter3,pn_victoriaGTlistcounter): # iterate over the two lists at the same time so not repeated
    pnrefmatch.append([]) # create lists of lists that can be appended to for each column
    for z in pn_victoriaGT:
        if z[30]==z[x]:
            pnrefmatch[y].append("1")
        else:
            pnrefmatch[y].append("0")
# join the matches onto the end of the original list table
# to each element of pnrefmatch, join each line of the original list
pnrefmatch2=[list(x) for x in zip(*pnrefmatch)] # transpose the list
pncounter4=list(range(0, len(pn_victoriaGT)))
pn_victoriaGTnew1 = [a + b for a, b in zip(pn_victoriaGT,pnrefmatch2)] # join the two lists

# Pn Alt Match
pnaltmatch=[]
for (x,y) in zip(pncounter3,pn_victoriaGTlistcounter): # iterate over the two lists at the same time so not repeated
    pnaltmatch.append([]) # create lists of lists that can be appended to for each column
    for z in pn_victoriaGT:
        if z[31]==z[x]:
            pnaltmatch[y].append("1")
        else:
            pnaltmatch[y].append("0")
# join the matches onto the end of the original list table
# to each element of pnrefmatch, join each line of the original list
pnaltmatch2=[list(x) for x in zip(*pnaltmatch)] # transpose the list
pncounter4=list(range(0, len(pn_victoriaGT)))
pn_victoriaGTnew2 = [a + b for a, b in zip(pn_victoriaGT,pnaltmatch2)] # join the two lists

#### 3. Write out the new lists to tab-separated files, with trailing newlines
with open('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2', 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in mz_malawiGTnew1) # join each element as tab-separated and trailing newlines for each line in the list
file.close() # close the file

with open('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2', 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in mz_malawiGTnew2) # join each element as tab-separated and trailing newlines for each line in the list
file.close() # close the file

with open('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2', 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in pn_victoriaGTnew1) # join each element as tab-separated and trailing newlines for each line in the list
file.close() # close the file

with open('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2', 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in pn_victoriaGTnew2) # join each element as tab-separated and trailing newlines for each line in the list
file.close() # close the file
