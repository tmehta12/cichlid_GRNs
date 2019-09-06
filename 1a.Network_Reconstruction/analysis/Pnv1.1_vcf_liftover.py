#!/usr/bin/env python

cd /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/

# create a python script to do the following
# if $6 in liftover == ($1 in vcf && $2 in vcf < (less than) $8 in liftover) then $2 in vcf = $2 in vcf - $9 in liftover


## Here, we use normal python coding to tackle the problem - don't use this, use the pandas method below

nano Pnv1.1_vcf_liftover-normal.py

#!/usr/bin/env python

## Here, we use normal python coding to tackle the problem

# liftover file and assign to a variable
liftover = open("/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/secondSet_RegionsOfInterest_forTarang_withNewPos.txt")

# a. first, add columns of the liftover file to a list
next(liftover) # this skips the first line (header)
lift = [] # add the different cols to a list
for line in liftover:
    scaffold = line.split('\t')[0]
    chr_newRef = line.split('\t')[5]
    start_newRef = line.split('\t')[6]
    end_newRef = line.split('\t')[7]
    diff_ref_newRef = line.split('\t')[8].rstrip('\n')
    lift.append((scaffold,chr_newRef, start_newRef, end_newRef, diff_ref_newRef))
#lift

# A. Find where the vcf columns headers are
lines=open('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf')
line_number = 0
for line in lines: # read the vcf, line by line, looping over it
    if not line.startswith('##'): # look for the point where the ## stops - this is the marker for finding the header that is the line below
#        print(line, line_number) # print the line and the line number that this occurs (linenumber will get amended)
        break # then use the 'break' function to stop the loop at the point this happens otherwise it will continue looping through the lines
#    line_number = line_number + 1 # print the linenumber + 1

# open the vcf file, skip the vcf rows as per above, plus the topmost header and assign to a variable
# also, save the colheader to a variable
with open('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf') as f:
    vcf = f.readlines()[6911:]
with open('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf') as f:
    vcfheader=f.readlines()[6910:6911]
# print(vcfheader)

# this the short version
with open('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf.short') as f:
    vcf = f.readlines()[6911:]
with open('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf.short') as f:
    vcfheader=f.readlines()[6910:6911]

# here, we do the actual match and output a file
list_of_lines = []
for line in vcf:
    vcf_chr = line.split('\t')[0]
    vcf_end = int(line.split('\t')[1])
    #vcf_rest = line.split('\t')[2:] # instead, create a list of lists each containing a single column?
    for region in lift:
        (scaffold,lift_chr_newRef, lift_start_newRef, lift_end_newRef, lift_diff_ref_newRef) = region
        if str(vcf_chr) == str(lift_chr_newRef) and int(vcf_end) >= int(lift_start_newRef) and int(vcf_end) <= int(lift_end_newRef):
            new_vcf_end = int(vcf_end) - int(lift_diff_ref_newRef)


            ### THIS PART STILL NEEDS SORTING OUT - WANT TO ADD THE REST OF THE VCF ONTO THE END, UNCHANGED - BUT WHEN ADD AS A LIST IT OUTPUTS IT AS A SINGLE STRING
            #list_of_lines.append((scaffold,new_vcf_end,line.split('\t')[2:]))
            list_of_lines.append('\t'.join(scaffold,new_vcf_end,'\t'.join(line.split('\t')[2:]))) # ORIGINAL LINE
            ######## TypeError: join() takes exactly one argument (3 given) ^^^^ error in above line - check

            print(scaffold,new_vcf_end,'\t'.join(line.split('\t')[2:]))
            #vcf_amended = str(scaffold) + str(new_vcf_end) + str(vcf_rest)

vcf_amended = open("tarangRegions.secondSet_Pnv1.1liftover.vcf", "w")
#vcf_amended = open("xx.vcf", "w")
vcf_amended.write(str('\n'.join(list_of_lines)))

# remember to close the file
vcf_amended.close()

# create a batch script

nano Pnv1.1_vcf_liftover-normal-sbatch.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 250000 # memory pool for all cores
#SBATCH -t 0-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

ml python/3.5

python3 Pnv1.1_vcf_liftover-normal.py

# run the above:
sbatch Pnv1.1_vcf_liftover-normal-sbatch.sh

#### Things to amend above still
# check the output as it appears the very first line is skipped
# the nt output, if correct, looking on the browser, should be CCGATGGCCTCGGCCAGAAC..



######### Here, we use the pandas module to tackle the problem

cd /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/
# head -100000 tarangRegions.secondSet.vcf > tarangRegions.secondSet.vcf.short # this is just to test first

### since the vcf file is particularly large, best to split it up into smaller chunks to run in parallel using bash
### then run as an array by inputting as sysargs in your python and batch script
### then rejoin the output files

mkdir splitvcf

nano splitvcf.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 100000 # memory pool for all cores
#SBATCH -t 0-10:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

split -d -l 100000 -a 3 tarangRegions.secondSet.vcf tarangRegions.secondSet.vcf_ # this splits every 100,000 lines and assigns three digit prefix to assigned filename

# then need to process the very first file so that you split based on the header rows (up to 6912)
awk 'NR < 6912 { print >> "splitvcf/tarangRegions.secondSet.vcf-header"; next } {print >> "splitvcf/tarangRegions.secondSet.vcf_START" }' splitvcf/tarangRegions.secondSet.vcf_0000
rm splitvcf/tarangRegions.secondSet.vcf_0000
mv splitvcf/tarangRegions.secondSet.vcf_START splitvcf/tarangRegions.secondSet.vcf_0000

# Save the single line colheader to a file so that it can be catted to each split file as the python panda script needs the identifier
awk 'FNR == 6911 {print;exit}' tarangRegions.secondSet.vcf > tarangRegions.secondSet.vcf.colhead
mkdir splitvcf/newfiles # make a new dir for the new files to cat the colheader
# cat to the top of each file
cd splitvcf/
for i in tarangRegions.secondSet.vcf_* ; do cat ../tarangRegions.secondSet.vcf.colhead $i > newfiles/$i ; done
rm tarangRegions.secondSet.vcf_*
mv newfiles/* .
rm -r newfiles

mkdir splitvcf/output # make a directory for output files of the python script

# run the above
sbatch splitvcf.sh #{DONE}


nano Pnv1.1_vcf_liftover-pandas.py

#!/usr/bin/env python

########################################################################################################################################################################################################
# The Pn vcf file was generated by mapping to an unpublished PacBio Pn Genome_sequences
# Since we want to use coordinates from the old Pn (v1.1) assembly, we need to change the file
# This script changes the Pn vcf file (that are candidate gene TFBSs in promoters overlapping a vcf file of SNPs in Lake Victoria sepcies) so that the coordinates match the old Pn (v1.1) assembly.

# How to run the script, by running on the terminal:
# liftover_file = ('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/secondSet_RegionsOfInterest_forTarang_withNewPos.txt')
# input_vcf = ('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf') # this has been split into lines of 100000 instead to run as an array
# python Pnv1.1_vcf_liftover-pandas.py liftover_file input_vcf outfile

########################################################################################################################################################################################################

#### 0. Load in the necessary modules
import pandas as pd # pandas is the module for manipulating dfs
import numpy as np # this is for any statistical calculations that we may do
import sys # pass system arguments
#import seaborn as sns # seaborn is a visualisation module
#import matplotlib.pyplot as plt # we use this for plotting
# %matplotlib inline # this is required for viewing plots, inline, in Jupyter only
# pd.__version__ # view pandas version


# 1. Read in the liftover file
# regions = pd.read_csv('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/secondSet_RegionsOfInterest_forTarang_withNewPos.txt', sep='\t')
# regions # view the format of the file
regions = pd.read_csv(sys.argv[1], sep='\t')

# 2. Read the lines of the vcf file and find out where the actual header starts
# lines = open('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf').readlines()
#
# # A. Find where the vcf columns headers are
# line_number = 0
# for line in lines: # read the vcf, line by line, looping over it
#     if not line.startswith('##'): # look for the point where the ## stops - this is the marker for finding the header that is the line below
#         print(line, line_number) # print the line and the line number that this occurs (linenumber will get amended)
#         break # then use the 'break' function to stop the loop at the point this happens otherwise it will continue lopping through the lines
#     line_number = line_number + 1 # print the linenumber + 1

# 3. Read in the vcf file at the point where the actual header starts
# vcfpanda = pd.read_csv('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf', skiprows=6910, sep='\s') # space is the separator
# vcfpanda = pd.read_csv('/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf.short', skiprows=6910, sep='\t') # this was to use as a test
vcfpanda = pd.read_csv(sys.argv[2], sep='\t')

# 4. Define a panda function to do the matching:
# if col6 (actually col5 in python as starts from 0) in liftover == (col1 in vcf && col2 in vcf (is >col7 and <col8 in liftover) then col2 in vcf = col2 in vcf - col9 in liftover and col1 in vcf = col1 in liftover

# this is just a test
# vcf_row = vcfpanda.iloc[0,:]
# print(vcf_row)
# matching_rows = regions[
#        (vcf_row['#CHROM'] == regions['region_newGenome']) & # match the chr in the VCF to the chr in the liftover file
#        (vcf_row['POS'] >= regions['start_newGenome']) & # match only lines where the pos in the VCF is > the start of the newRef in the liftover
#        (vcf_row['POS'] <= regions['end_newGenome']) # match only lines where the pos in the VCF is < the start of the newRef in the liftover
# ]
# print(matching_rows)

def find_regions(vcf_row):
    # for a single region in the liftover, find matching rows in the VCF file
    matching_rows = regions[
    (vcf_row['#CHROM'] == regions['region_newGenome']) & # match the chr in the VCF to the chr in the liftover file
    (vcf_row['POS'] >= regions['start_newGenome']) & # match only lines where the pos in the VCF is > the start of the newRef in the liftover
    (vcf_row['POS'] <= regions['end_newGenome']) # match only lines where the pos in the VCF is < the start of the newRef in the liftover
    ].iloc[0,:]
    row_vcf = vcf_row.copy()

    # update the positions in the VCF file as per above
    row_vcf['POS'] = row_vcf['POS'] - matching_rows['diff_ref_newRef']

    # update the chromosome column in the vcf file
    row_vcf['#CHROM'] = matching_rows['scaffold_punv1']
    return row_vcf

# This applies the function
newvcfpanda = vcfpanda.apply(
    find_regions,
    axis=1)
#axis is applying the function to either 0 = cols, or here 1 = rows

# write out to file
newvcfpanda.to_csv(sys.argv[3], encoding='utf-8', sep='\t', index=False)

#### Run the above using an sbatch array
# in total there are 3407 files that this needs to be ran on in splitvcf/
# this does not include splitvcf/tarangRegions.secondSet.vcf-header
# files are splitvcf/tarangRegions.secondSet.vcf_* # 0000-3406

nano Pnv1.1_vcf_liftover-pandas-sbatch_array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-3406
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 splitvcf/tarangRegions.secondSet.vcf_* > vcffiles # create a list of all split vcf files to use as input
mapfile -t vcffiles < vcffiles # assign as elements to $vcffiles variable

sed 's/tarang/output\/tarang/g' vcffiles | sed 's/$/.out/g' > vcfouts # create a list of all output files required
mapfile -t vcfouts < vcfouts # assign as elements to $vcfouts variable

# ml python/3.5
# won't be loading any python module and instead running your user installed version by just calling python - Python 3.6.5 |Anaconda, Inc.|

liftover_file=(/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/secondSet_RegionsOfInterest_forTarang_withNewPos.txt)
python Pnv1.1_vcf_liftover-pandas.py $liftover_file ${vcffiles[${SLURM_ARRAY_TASK_ID}]} ${vcfouts[${SLURM_ARRAY_TASK_ID}]}

# run the above
sbatch Pnv1.1_vcf_liftover-pandas-sbatch_array.sh #{DONE}

### several files did not run and hence we need to run those
ls -1 splitvcf/output/*.out | sort > vcfout_ran
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' vcfout_ran vcfouts | grep NA | cut -f1 > vcfouts2 # these are the outfiles that were not created
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' vcfout_ran vcfouts | grep NA | cut -f1 | sed 's/output\///g' | sed 's/.out//g' > vcffiles2 # these are the infiles that need re-running

nano Pnv1.1_vcf_liftover-pandas-sbatch_array_rerun.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-217
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t vcffiles2 < vcffiles2 # assign as elements to $vcffiles variable
mapfile -t vcfouts2 < vcfouts2 # assign as elements to $vcfouts variable

# ml python/3.5
# won't be loading any python module and instead running your user installed version by just calling python - Python 3.6.5 |Anaconda, Inc.|

liftover_file=(/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/secondSet_RegionsOfInterest_forTarang_withNewPos.txt)
python Pnv1.1_vcf_liftover-pandas.py $liftover_file ${vcffiles2[${SLURM_ARRAY_TASK_ID}]} ${vcfouts2[${SLURM_ARRAY_TASK_ID}]}

# run the above
sbatch Pnv1.1_vcf_liftover-pandas-sbatch_array_rerun.sh #{DONE}



# remove the very large files once you are sure of the outputs
rm tarangRegions.secondSet.vcf # DONE - a zipped version is backed up in group-vh
rm slurm.t* # DONE
rm splitvcf/tarangRegions.secondSet.vcf_* # DONE

# Then, move over to the script found below to work on the outfiles in /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/splitvcf/output/*.out
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/Identify_EMSAcand.sh




# ### Then, create a script to merge the split files back to one file - don't do this yet since the file will be very large to run a bedtools intersect! Keep working on the chunks!
#
# nano catfiles.sh
#
# #!/bin/bash -e
# #SBATCH -p tgac-medium # partition (queue)
# #SBATCH -N 1 # number of nodes
# #SBATCH -n 1 # number of tasks
# #SBATCH --mem 100000 # memory pool for all cores
# #SBATCH -t 0-23:59 # time (D-HH:MM)
# #SBATCH -o slurm.%N.%j.out # STDOUT
# #SBATCH -e slurm.%N.%j.err # STDERR
# #SBATCH --mail-type=ALL # notifications for job done & fail
# #SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#
# cd splitvcf/output
# for i in *.out ; do sed -i '1d' $i ; done # first, remove the colheader for all files, in place
#
# cat *.out | sort -k1,1 -k2,2n >> /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.Pnv1.1liftover.vcf
# cat /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf.colhead /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.Pnv1.1liftover.vcf >> /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.Pnv1.1liftover.vcf2
#
# # run the above
# sbatch catfiles.sh #DONE - but check output since the sorting was added after the script was running..(sorted file required for bedtools intersect)
