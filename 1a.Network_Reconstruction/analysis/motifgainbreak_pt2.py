#!/usr/bin/env python3

########################################################################################################################################################################################################
## Categorize TFBS mutations as breaking/gaining

# to be ran from: /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss

# 1. replaced all SNPs overlapping motifs in the promoter with the ALT genotype (in /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/motifgainbreak.py)
# 2. re-scanned for motifs (this tests whether the SNP breaks the site as flanking is unchanged) and compare the p-value for the original and mutated hit (in /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/SNP_overlap.sh)
# 3. HERE, in this script we do the bulk and categorise each SNP overlapping motif according to motif gain (from ancestor REF SNP) and loss (from ancestor from REF SNP)

# Strategy
# this will require:
	# 1. for each snp comparison e.g. mz vs nb; does the ancestral species (nb) have the same motif in its orthologous promoter and in the same location - this can be seen in the SNP overlapping motif files of the ancestral species (nb vs mz) as the SNP will be the same (REF vs ALT and ALT vs REF) in both species as this is based upon multiple alignment
        # /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.OGID
        # /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.OGID
        # /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.OGID
        # /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.OGID
        # /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.OGID
        # NOTE: All Num files have been deleted and hence to correlate to the SNP in the SNP file, using enumerate, count to the nth line of the SNPID in the fasta e.g. _SNPid:5 means the fifth SNP line
        # /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/{Mz,Pn,Ab,Nb,On}Prom_{mz,pn,ab,nb,on}SNP.5kbProm-CS_CW_JP.txt - mutant promoter filtered fimo scan files of CS, CW and JASPAR matrices
	# 2. if not - motif gain as a result of SNP in most extant level 1
	# 3. if so
		# A. compare the p-value for the motif hit in 1. original species motif hit e.g. mz; 2. ancestral species motif hit e.g. nb; and 3. SNP altered motif hit e.g. mz vs nb
		# B. if:
			# i. in the most extant species, p-val of 1 > 2 and 3 - motif gain as a result of SNP in most extant level 2
			# ii. motif present in ancestral species only - motif loss as a result of SNP in most extant level 1
			# iii. in the most extant species, p-val of 1 < 2 and 3 - motif loss as a result of SNP in most extant level 2

########################################################################################################################################################################################################

#### 0. Load in the necessary modules
import csv # read in files
import sys # pass system arguments
import time # this is used for time delays
from more_itertools import unique_everseen # this is used for getting unique lines
from tqdm import tqdm # this is to time each mapping


#### 1. arguments
# an argument for whether the number of command input files is too short
if len(sys.argv) < 7:
    print('\nNot enough arguments entered.\nUsage: 1. MainSp_SNP_VCF_file 2. AncSp_SNP_VCF_file 3. Mutant_promoter_fimo(mainsp_vs_ancspSNP) 4. Main species ID e.g. Mz 5. Ancestral species ID e.g. Pn 6. Outfile e.g. MzProm_pn1SNP_Motif_gain_loss.txt\n')
    sys.exit()

# print('----------\nOrthogrouping file is:\t' + sys.argv[1])
print('----------\nMost extant species SNP VCF file is:\t' + sys.argv[1] + ' (species is: ' + sys.argv[4] + ')')
print('----------\nAncestral species SNP file is:\t' + sys.argv[2] + ' (species is: ' + sys.argv[5] + ')')
print('----------\nMotif scan in mutant promoters (most extant vs ancestral) is:\t' + sys.argv[3])
print('----------\n>>>>>>>>>>')

#### 2. Load in the files as dictionaries
# # Load the orthogrouping as a dictionary - put each orthogroup as the key and then the OGID and cichlid gene IDS and gene symbols as the value
# with open("OGIDS.txt5", 'r') as f: #sys.argv[1]
#     OGIDdict = {}
#     for line in csv.reader(f, delimiter = '\t'):
#         if line[0] in OGIDdict:
#             OGIDdict[line[0]].append(line)
#         else:
#             OGIDdict[line[0]] = [line]
#     # print(OGIDSdict)
# print('----------\nOGIDS loaded as dictionary - orthogroup[key] and rest[values]')

# a. load in MainSp_SNP_VCF_file as a dictionary - put each target gene OGID as the key and then the rest of the information as values (that way, you can match the OGID first and not have to loop through everything)
with open(sys.argv[1], 'r') as f: #sys.argv[1] #mz11_pn1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.OGID.cut
    MsSNP = {}
    for line in csv.reader(f, delimiter = '\t'):
        if line[23] in MsSNP: # OGID column
            MsSNP[line[23]].append(line)
        else:
            MsSNP[line[23]] = [line]
    # print(MSSNP)
print('----------\nMain Species SNP VCF file loaded as dictionary - TG orthogroup[key] and rest[values]')

# c. load in AncSp_SNP_VCF_file as a dictionary - put each target gene OGID as the key and then the rest of the information as values (that way, you can match the OGID first and not have to loop through everything)
with open(sys.argv[2], 'r') as f: #sys.argv[2] #pn1_mz11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.OGID.cut
    AnSNP = {}
    for line in csv.reader(f, delimiter = '\t'):
        if line[23] in AnSNP: # OGID column
            AnSNP[line[23]].append(line)
        else:
            AnSNP[line[23]] = [line]
    # print(AnSNP)
print('----------\nAncestral Species SNP VCF file loaded as dictionary - TG orthogroup[key] and rest[values]')

# d. load in the processed Mutant_promoter_fimo(mainsp_vs_ancspSNP) as a dictionary - put each target gene OGID as the key and then the rest of the information as values (that way, you can match the OGID first and not have to loop through everything)
with open(sys.argv[3], 'r') as f: #sys.argv[3] #MzProm_pnSNP.5kbProm-CS_fimo.txt3.6.cut
    MtFIMO = {}
    MtFIMOheader = next(f) # skip the header and store in a separate variable (just in case it is required)
    for line in csv.reader(f, delimiter = '\t'):
        if line[5] in MtFIMO: # change to OGID column
            MtFIMO[line[5]].append(line)
        else:
            MtFIMO[line[5]] = [line]
    # print(MtFIMO)
print('----------\nMutant promoter fimo (mainsp_vs_ancspSNP) file loaded as dictionary - TG orthogroup[key] and rest[values]')
print('----------\n>>>>>>>>>>')

#### 3. Run the matching strategy and output in outfile

print('----------\nMatching strategy started and will be outputted to:\t' + sys.argv[6])

with open(sys.argv[6], 'w') as outfile:
    for k,v in MsSNP.items():
        if k in AnSNP: # shorten the loop by only looking in orthologous TGs
            for num1, hits1 in enumerate(MsSNP[k],1): # loop through each line and add a variable for line number to compare to mutant promoter in main species
                for num2, hits2 in enumerate(AnSNP[k],1): # loop through each line and add a variable for line number to compare to mutant promoter in ancestral species
                    if hits1[23]==hits2[23] and hits1[16]==hits2[16] and hits1[0]==hits2[9] and hits1[2]==hits2[10] and hits1[23] in MtFIMO: # match the TG, then the TF, then the SNP chr, then the SNP pos to identify whether the main species and ancestral species has the same motif in its orthologous promoter
                        # print(hits1,hits2,'-->\tSame motif and position in orthologous promoter: Compare p-value in mutant promoter, line number',num1,'(Main-Sp) and',num2,'(Anc-Sp)')
                        # time.sleep(5)
                        for MtFIMOhits in MtFIMO[hits1[23]]:
                            # print(num1,MtFIMOhits[0])
                            # time.sleep(5)
                            if str(num1)==str(MtFIMOhits[0]): # compare p-val for motif hit in 1. Main-sp e.g. mz; 2. Anc-sp e.g. pn; and 3. mutant promoter motif hit e.g. mz vs pn, and if
                                # print(num1,MtFIMOhits)
                                # time.sleep(5)
                                if hits1[26] > hits2[26] and hits1[26] > MtFIMOhits[14]: # i. in the most extant species, p-val of 1 > 2 and 3 - motif gain as a result of SNP in most extant level 2
                                    # print('-------Motif gain in main species (Mz) as a result of SNP from ancestral species (Pn) - Motif Gain Level 2' + '\nMAIN_SPECIES:\t' + '\t'.join(hits1) + '\t|||\tANCESTRAL_SPECIES:\t' + '\t'.join(hits2) + '\t|||\tMUTANT_5KB_PROMOTER_WITH_ANCESTRAL_ALLELIC_SNP:\t' + '\t'.join(MtFIMOhits) + '\n')
                                    # time.sleep(5)
                                    outfile.write('-------Motif gain in main species (' + sys.argv[4] + ') as a result of SNP from ancestral species(' + sys.argv[5] + ') - Motif Gain Level 2' + '\nMAIN_SPECIES:\t' + '\t'.join(hits1) + '\t|||\tANCESTRAL_SPECIES:\t' + '\t'.join(hits2) + '\t|||\tMUTANT_5KB_PROMOTER_WITH_ANCESTRAL_ALLELIC_SNP:\t' + '\t'.join(MtFIMOhits) + '\n')
                                elif hits1[26] < hits2[26] and hits1[26] < MtFIMOhits[14]: # ii. in the most extant species, p-val of 1 < 2 and 3 - motif loss as a result of SNP in most extant level 2
                                    # print('-------Motif loss in main species (Mz) as a result of SNP from ancestral species (Pn) - Motif Loss Level 2' + '\nMAIN_SPECIES:\t' + '\t'.join(hits1) + '\t|||\tANCESTRAL_SPECIES:\t' + '\t'.join(hits2) + '\t|||\tMUTANT_5KB_PROMOTER_WITH_ANCESTRAL_ALLELIC_SNP:\t' + '\t'.join(MtFIMOhits) + '\n')
                                    # time.sleep(5)
                                    outfile.write('-------Motif loss in main species (' + sys.argv[4] + ') as a result of SNP from ancestral species(' + sys.argv[5] + ') - Motif Loss Level 2' + '\nMAIN_SPECIES:\t' + '\t'.join(hits1) + '\t|||\tANCESTRAL_SPECIES:\t' + '\t'.join(hits2) + '\t|||\tMUTANT_5KB_PROMOTER_WITH_ANCESTRAL_ALLELIC_SNP:\t' + '\t'.join(MtFIMOhits) + '\n')
                    elif hits2[23]!=hits1[23] and hits2[16]!=hits1[16] and hits2[0]!=hits1[9] and hits2[2]!=hits1[10]: # iii. motif present in ancestral species only - motif loss as a result of SNP in most extant level 1
                        # print('-------Motif loss in main species (Mz) as a result of SNP from ancestral species (Pn) - Motif Loss Level 1' + '\nMAIN_SPECIES:\tNA' + '\t|||\tANCESTRAL_SPECIES:\t' + '\t'.join(hits2) + '\n')
                        # time.sleep(5)
                        outfile.write('-------Motif loss in main species (' + sys.argv[4] + ') as a result of SNP from ancestral species(' + sys.argv[5] + ') - Motif Loss Level 1' + '\nMAIN_SPECIES:\tNA' + '\t|||\tANCESTRAL_SPECIES:\t' + '\t'.join(hits2) + '\n')
                    elif hits1[23]!=hits2[23] and hits1[16]!=hits2[16] and hits1[0]!=hits2[9] and hits1[2]!=hits2[10]: # iv. motif present in main species only - motif gain as a result of SNP in most extant level 1
                        # print('-------Motif gain in main species (Mz) as a result of SNP from ancestral species (Pn) - Motif Gain Level 1' + '\nMAIN_SPECIES:\t' + '\t'.join(hits1) + '\t|||\tANCESTRAL_SPECIES:\tNA' + '\n')
                        # time.sleep(5)
                        outfile.write('-------Motif gain in main species (' + sys.argv[4] + ') as a result of SNP from ancestral species(' + sys.argv[5] + ') - Motif Gain Level 1' + '\nMAIN_SPECIES:\t' + '\t'.join(hits1) + '\t|||\tANCESTRAL_SPECIES:\tNA' + '\n')
