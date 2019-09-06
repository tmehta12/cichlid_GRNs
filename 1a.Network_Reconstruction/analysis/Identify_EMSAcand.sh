#!/bin/sh

##########################################################################################################################
#
# Identifying EMSA candidates in 5 cichlids (O. niloticus, P. nyererei, A. burtoni, M .zebra, N. brichardi)
#
#
# By Tarang Mehta, Earlham, UK
# Version 1.0 2018
#
##########################################################################################################################

# Structure
# 1. Focus on candidate genes
  # a. Overlap pairwise SNPs with candidate gene promoters TFBSs
  # b. Overlap the above with Lake species VCF to pull out segregating TFBSs in East African lake species
# 2. Focus on all genes
  # a. From regulon tables (matrix), pull out absent TFBSs
  # b. Overlap absent TFBSs with presence of SNP in SNP overlaps of TFBSs
  # c. Overlap all candidate genes with absent TFBSs with SNP overlaps


##########################################################################################################################
#
# 1.	Focus on candidate genes
#
##########################################################################################################################

########################################################################################
# ~~~~~~~ # a. Overlap pairwise SNPs with candidate gene promoters TFBSs
########################################################################################

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap
# DONE in SNP_overlap.sh - files are *-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.txt

########################################################################################
# ~~~~~~~ # b. Overlap the above with Lake species VCF to pull out segregating TFBSs in East African lake species
########################################################################################
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap

# Lake Malawi VCF here:
/tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf

# Lake Victoria VCF here (these are regions in the vcf overlapping pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.txt):
/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.Pnv1.1liftover.vcf

# to run a bedtools intersect, you need to convert the SNP overlap of candidates to a BED file
for h in ab mz nb on pn ; do
  for i in *-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.txt ; do
    awk -F':' '{print $1,$2}' OFS='\t' $i |
    awk '{print $2,$3,$4,$1,$8,$16,$6,$7,$11,$12,$13,$14,$15,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' OFS='\t' |
    sed "s/_substitutions_1217.bed_$h-motifenr_mergedmap2d_motifgenomeannot.txt//g" |
    grep -v 'si' |
    grep -v 'im' |
    grep -v '_danRer7' |
    grep -v '_gasAcu1' |
    grep -v '_oryLat2' |
    grep -v 'zgc' |
    sort -k1,1 -k2,2n > "$(basename "$i" .txt).bed"
  done
done


##### NOTE: AS OF 05/07/2018 YOU NEED TO CHECK WHETHER YOU SHOULD ADD THE STRAND FROM THE MOTIF ($16) ABOVE IN THE POSITION YOU DO. THIS IMPLIES THAT THE SNP IS THE SAME STRAND AS THE MOTIF, WHICH IS INCORRECT!!!!!!!
##### INSTEAD, WHAT YOU NEED TO DO, IS IF THE SNP HAS BEEN OUTPUTTED FROM A -VE STRAND ALIGNMENT, AND THE MOTIF IS CALLED IN THE +VE STRAND, THEN THE REF ALLELE AND THE ALT ALLELE NEED TO BE REVERSE COMPLEMENTED E.G. A>T, T>A ETC.
##### WHAT SHOULD BE DONE, SO AS TO NOT RE-RUN EVERYTHING, AT THE STAGE WHERE WILL IS FILTERING ALL THE TFBS-SNP OVERLAPS OF CANDIDATE GENE MOTIFS WITH WHETHER THEY ARE PRESENT IN THE PROMOTER ALIGNMENTS (BELOW), IS TO:
  # 1. READ IN THE ORIGINAL MAF ALIGNMENT FILES
  # 2. FIND OUT WHETHER THE SNP CAME FROM A -VE STRAND ALIGN (THEN ADD THE STRAND OF THE SNP TO THE FILE - REPLACE COL6)
  # 3. IF SO, AND THE CORRESPONDING OVERLAPPING MOTIF IS IN THE +VE STRAND, REVERSE COMPLEMENT THE REF AND ALT ALLELES
  # WILL has corrected for this in his filtering script



### in each script, what you can also do is take the reverse SNPs e.g nb1_mz11 but obviously the mz11 coordinates to do a VCF overlaps
### that is how the nr2c2 TFBS will be picked up in the VCF

# ###### THE BELOW NEEDS CHECKING AS THE NR2C2/RXRB SNP IN N.B SWS1 VS MZ IS NOT PICKED UP IN THE BEDTOOLS INTERSECT, this being checked with below
#
# # original
# grep 'nb.gene.s1.386' nb-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed | grep 'nb1_mz11' | head
# scaffold_1	19181358	19181359	nb1_mz11	100	-	C	T	scaffold_102	2340622	scaffold_1	19181348	19181362	nb.gene.s1.386	opn1sw1	opn1sw1	q2l6a1_oryla	+	motif_208	nb.gene.s44.59	nr2c2	22.4429	9.72e-09	0.00603	CAGGGTCAGAGGTCA	2c	0.115
#
# # changed to
# grep 'nb.gene.s1.386' nb-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed | head
# scaffold_102	2340622	2340623	nb1_mz11	100	-	C	T	scaffold_102	2340622	scaffold_1	19181348	19181362	nb.gene.s1.386	opn1sw1	opn1sw1	q2l6a1_oryla	+	motif_208	nb.gene.s44.59	nr2c2	22.4429	9.72e-09	0.00603	CAGGGTCAGAGGTCA	2c	0.115
#
# # first check where the SNP is in the vcf - will the bedtools intersect (below) actually pick this SNP up?!?!?
# bedtools intersect -a pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed
#
# ### the problem could be the coordinate, instead of doing:
# scaffold_102	2340622	2340623
# ### do you need to do:
# scaffold_102	2340622	2340622
#
# nano test.sh
#
# #!/bin/bash -e
# #SBATCH -p tgac-medium # partition (queue)
# #SBATCH -N 1 # number of nodes
# #SBATCH -n 1 # number of tasks
# #SBATCH --mem 48000 # memory pool for all cores
# #SBATCH -t 0-5:59 # time (D-HH:MM)
# #SBATCH -o slurm.%N.%j.out # STDOUT
# #SBATCH -e slurm.%N.%j.err # STDERR
# #SBATCH --mail-type=ALL # notifications for job done & fail
# #SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#
# head -n 3123 /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf > Malinsky_et_al_2017_filteredSNPs_head.vcf
# cat Malinsky_et_al_2017_filteredSNPs_head.vcf > scaf102_2340622_MalawifilteredSNPs.vcf
# cat Malinsky_et_al_2017_filteredSNPs_head.vcf > scaf102_MalawifilteredSNPs.vcf
# grep -wiF scaffold_102 /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf | grep '2340622' >> scaf102_2340622_MalawifilteredSNPs.vcf
# grep -wiF scaffold_102 /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf  >> scaf102_MalawifilteredSNPs.vcf
#
# # then do a bedtools intersect with the above output file to check if works
#
# ml bedtools/2.25.0
# ml zlib
#
# # do the original one first to check
# bedtools intersect -a nb-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b scaf102_2340622_MalawifilteredSNPs.vcf -wa -wb
# bedtools intersect -a nb-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b scaf102_MalawifilteredSNPs.vcf -wa -wb
#
# grep nr2c2 nb-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed > test1.bed
# # then change to scaffold_102	2340622	2340622
# bedtools intersect -a test1.bed -b scaf102_2340622_MalawifilteredSNPs.vcf -wa -wb
#
#
#
# #### main note: consider changing to $10-1,$10 in the bed files since when you put this in the browser, you get the correct nt e.g.
# scaffold_102	2340622	2340622 # insetad of
# scaffold_102	2340622	2340623
#
# ### first orignal file (mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed) had 6245 lines
# ###### THE BELOW NEEDS CHECKING AS THE NR2C2/RXRB SNP IN N.B SWS1 VS MZ IS NOT PICKED UP IN THE BEDTOOLS INTERSECT, this being checked with above

####### This is for M. zebra candidates - overlapping with the Lake Malawi VCF

nano mzTFBS-VCFoverlap.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

ml bedtools/2.25.0
ml zlib

bedtools intersect -a mz-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2
# then, take the reverse where mz11 is not the ref but the alt
# this is to create Mz as reference in other species SNPs
for i in ab nb on pn ; do
  grep '_mz11' $i-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed |
  grep -v 'si' |
  grep -v 'im' |
  grep -v '_danRer7' |
  grep -v '_gasAcu1' |
  grep -v '_oryLat2' |
  grep -v 'zgc' |
  sort -k1,1 -k2,2n |
  awk '{print $9,$10-1,$10,$4,$5,$6,$7,$8,$1,$2,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27}' OFS='\t' > $i-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed
done
bedtools intersect -a ab-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2
bedtools intersect -a nb-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2
bedtools intersect -a on-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2
bedtools intersect -a pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2

# run the above
qsub -q Test -l select=1:mem=500GB:ncpus=1 mzTFBS-VCFoverlap.sh # {DONE}

### OR, SLURM version
nano mzTFBS-VCFoverlap-bash.sh

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 250000 # memory pool for all cores
#SBATCH -t 0-23:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

ml bedtools/2.25.0
ml zlib

bedtools intersect -a mz-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2
# then, take the reverse where mz11 is not the ref but the alt
# this is to create Mz as reference in other species SNPs
for i in ab nb on pn ; do
  grep '_mz11' $i-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed |
  grep -v 'si' |
  grep -v 'im' |
  grep -v '_danRer7' |
  grep -v '_gasAcu1' |
  grep -v '_oryLat2' |
  grep -v 'zgc' |
  sort -k1,1 -k2,2n |
  awk '{print $9,$10-1,$10,$4,$5,$6,$7,$8,$1,$2,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27}' OFS='\t' > $i-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed
done
bedtools intersect -a ab-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2
bedtools intersect -a nb-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2
bedtools intersect -a on-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2
bedtools intersect -a pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.MZref.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2

# Columns are:
#
# 1.	Pairwise_SNP_Ref_chr - NOTE: this is the Pairwise_SNP_Alt_chr in [species]_mz11 cases e.g. nb1_mz11 - this is the chr for mz11
# 2.	Pairwise_SNP_ Ref_start (0-based) - NOTE: this is the Pairwise_SNP_Alt_start in [species]_mz11 cases e.g. nb1_mz11 - this is the start for mz11
# 3.	Pairwise_SNP_ Ref_end (0-based) - NOTE: this is the Pairwise_SNP_Alt_end in [species]_mz11 cases e.g. nb1_mz11 - this is the end for mz11 > this is the 1-based position of the SNP
# 4.	Pairwise_SNP_comparison_species
# 5.	Pairwise_SNP_quality
# 6.	Pairwise_SNP_strand
# 7.	Pairwise_SNP_Ref_allele
# 8.	Pairwise_SNP_Alt_allele
# 9.	Pairwise_SNP_Alt_chr - NOTE: this is the Pairwise_SNP_Ref_chr in [species]_mz11 cases e.g. nb1_mz11 - this is the chr for nb1
# 10.	Pairwise_SNP_Alt_position (1-based) - NOTE: this is the Pairwise_SNP_Ref_position in [species]_mz11 cases e.g. nb1_mz11 - this is the pos for nb1
# 11.	WG_motif_chr
# 12.	WG_motif_start
# 13.	WG_motif_end
# 14.	cichlid_geneID
# 15.	motif_genesymbolDr
# 16.	motif_genesymbolGa
# 17.	motif_genesymbolSp
# 18.	motif_strand
# 19.	TFmotif_ID
# 20.	TFmotif_cichlidID
# 21.	TFmotif_gene_symbol
# 22.	fimo_score
# 23.	fimo_pval
# 24.	fimo_qval
# 25.	motif_seq
# 26.	conf_level
# 27.	conf_score
# 28.	MalawiVCF_SNP_Ref_chr
# 29.	MalawiVCF_SNP_Ref_position
# 30.	MalawiVCF_SNP_Ref_ID
# 31.	MalawiVCF_SNP_Ref_allele
# 32.	MalawiVCF_SNP_Alt_allele
# 33.	MalawiVCF_SNP_Quality
# 34.	MalawiVCF_SNP_filter
# 35.	MalawiVCF_SNP_info
# 36.	MalawiVCF_SNP_format
# 37 to End of columns. genotype in Lake Malawi species

### Will is then filtering the above file for whether the SNP is actually present in the promoter alignments (since SNPs generated based on a multiple genome alignment)
cp mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed2 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap # Has 9067 lines
# The comparison done is whether the chr + coordinate in Will's file + the ref allele (correcting for reverse complement in cases of -ve strand) correspond to the chr (col1) + coordinates (col2,3) + allele (col7)

### the file can be found here - 5712 out of 9067 lines retained after filtering:
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 #5712/9067
/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 #5712/9067

### create two files from above:
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
# 1. Alternative allele same in pairwise species promoter comparison and Malawi VCF ($8==$32) # total is 1981
awk '$8==$32' mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718_sameALTallele.bed
    # total is 1981
# 2. Alternative allele different in pairwise species promoter comparison and Malawi VCF ($8!=$32) - this means that the SNP is totally segregated in Malawi species from Mz and other one of the four comparison cichlid representatives (Pn, Ab, Nb, On)
awk '$8!=$32' mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718_diffALTallele.bed
    # total is 3731


####### This is for P. nyererei candidates - overlapping with the Lake Victoria VCF

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap

# Lake Victoria VCF has been amended (for liftover) and split into chunks to run parallel jobs
# The Lake Victoria VCF (chunks) are regions in their (Ole's) large vcf overlapping pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.txt
# Liftover, chunked and amended files prepped in this script: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/Pnv1.1_vcf_liftover.py
# Files here (3407 in total): /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/splitvcf/output/*.out

### Run using a SLURM array

# Before, running prep some required files and folders
mkdir pnTFBSVCFoverlap_out # make an output directory for all the vcf chunks to be processed

# Take the reverse where pn1 is not the ref but the alt
# this is to create Pn as reference in other species SNPs
for i in ab nb on mz ; do
  grep '_pn1' $i-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed |
  grep -v 'si' |
  grep -v 'im' |
  grep -v '_danRer7' |
  grep -v '_gasAcu1' |
  grep -v '_oryLat2' |
  grep -v 'zgc' |
  sort -k1,1 -k2,2n |
  awk '{print $9,$10-1,$10,$4,$5,$6,$7,$8,$1,$2,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27}' OFS='\t' > $i-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.PNref.bed
done # DONE

nano catheader.sh

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

#### you need to write the vcf header to each input file so that it recognises as vcf
ls -1 /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/splitvcf/output/*.out > vcffilesorig # create a list of all split vcf files to use as input
mapfile -t vcffilesorig < vcffilesorig # assign as elements to $vcffilesorig variable

header=(/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/splitvcf/tarangRegions.secondSet.vcf-header)
cat $header ${vcffilesorig[${SLURM_ARRAY_TASK_ID}]} | awk '$4!="." && $5!="."' > ${vcffilesorig[${SLURM_ARRAY_TASK_ID}]}-2 ### Also, we filter the files so that there is actually a SNP (since it contains regions where there isn't since it is an overlap)

### if several files did not run then they will need to be re-ran - this will do the check and can be used as the file for the next array
ls -1 /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/splitvcf/output/*.out-2 | sed 's/out-2/out/g' > vcffiles_out_ran
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' vcffiles_out_ran vcffilesorig | grep NA | cut -f1 > vcffiles_out_ran2 # these are the infiles for which outfiles were not created

# run the above
sbatch catheader.sh # DONE


nano catheader_rerun.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-193
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t vcffiles_out_ran2 < vcffiles_out_ran2 # assign as elements to $vcffilesorig variable

header=(/tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/splitvcf/tarangRegions.secondSet.vcf-header)
cat $header ${vcffiles_out_ran2[${SLURM_ARRAY_TASK_ID}]} | awk '$4!="." && $5!="."' > ${vcffiles_out_ran2[${SLURM_ARRAY_TASK_ID}]}-2


nano pnTFBS-VCFoverlap-bash_array.sh

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

ml bedtools/2.25.0
ml zlib

#### then set up the array for bedtools intersect using the new generated files above
ls -1 /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/splitvcf/output/*.out-2 > vcffiles # create a list of all split vcf files to use as input
mapfile -t vcffiles < vcffiles # assign as elements to $vcffiles variable

sed 's/^.*output\///g' vcffiles | sed 's/^/pnTFBSVCFoverlap_out\//g' | sed 's/$/.candIntersect/g' > vcfouts # create a list of all output files required
mapfile -t vcfouts < vcfouts # assign as elements to $vcfouts variable

bedtools intersect -a pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed -b ${vcffiles[${SLURM_ARRAY_TASK_ID}]} -wa -wb > ${vcfouts[${SLURM_ARRAY_TASK_ID}]}
# Take the reverse where pn1 is not the ref but the alt
bedtools intersect -a ab-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.PNref.bed -b ${vcffiles[${SLURM_ARRAY_TASK_ID}]} -wa -wb >> ${vcfouts[${SLURM_ARRAY_TASK_ID}]}
bedtools intersect -a nb-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.PNref.bed -b ${vcffiles[${SLURM_ARRAY_TASK_ID}]} -wa -wb >> ${vcfouts[${SLURM_ARRAY_TASK_ID}]}
bedtools intersect -a on-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.PNref.bed -b ${vcffiles[${SLURM_ARRAY_TASK_ID}]} -wa -wb >> ${vcfouts[${SLURM_ARRAY_TASK_ID}]}
bedtools intersect -a mz-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.PNref.bed -b ${vcffiles[${SLURM_ARRAY_TASK_ID}]} -wa -wb >> ${vcfouts[${SLURM_ARRAY_TASK_ID}]}
# NOTE: you will get several ': ambiguous redirect' errors in cases where nothing can be added


### if several files did not run then they will need to be re-ran - this will do the check and can be used as the file for the next array
ls -1 pnTFBSVCFoverlap_out/ | sed 's/.candIntersect//g' | sed 's/^/\/tgac\/scratch\/mehtat\/Cichlid_GRNs\/functional_variant_analysis\/Lake_Victoria_Data_Ole2018\/2.fullTFBSNPs_overlap_May2018\/splitvcf\/output\//g' > vcffiles_overlap_ran
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' vcffiles_overlap_ran vcffiles | grep NA | cut -f1 > vcffiles_overlap_ran2 # these are the infiles for which outfiles were not created

nano pnTFBS-VCFoverlap-bash_array_rerun.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-133
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ml bedtools/2.25.0
ml zlib

### use the new rerun files as array IDs
mapfile -t vcffiles_overlap_ran2 < vcffiles_overlap_ran2 # assign as elements to $vcffiles variable

sed 's/^.*output\///g' vcffiles_overlap_ran2 | sed 's/^/pnTFBSVCFoverlap_out\//g' | sed 's/$/.candIntersect/g' > vcffiles_overlap_ran2outs # create a list of all output files required
mapfile -t vcffiles_overlap_ran2outs < vcffiles_overlap_ran2outs # assign as elements to $vcffiles_overlap_ran2outs variable

bedtools intersect -a pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed -b ${vcffiles_overlap_ran2[${SLURM_ARRAY_TASK_ID}]} -wa -wb > ${vcffiles_overlap_ran2outs[${SLURM_ARRAY_TASK_ID}]}
# Take the reverse where pn1 is not the ref but the alt
bedtools intersect -a ab-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.PNref.bed -b ${vcffiles_overlap_ran2[${SLURM_ARRAY_TASK_ID}]} -wa -wb >> ${vcffiles_overlap_ran2outs[${SLURM_ARRAY_TASK_ID}]}
bedtools intersect -a nb-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.PNref.bed -b ${vcffiles_overlap_ran2[${SLURM_ARRAY_TASK_ID}]} -wa -wb >> ${vcffiles_overlap_ran2outs[${SLURM_ARRAY_TASK_ID}]}
bedtools intersect -a on-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.PNref.bed -b ${vcffiles_overlap_ran2[${SLURM_ARRAY_TASK_ID}]} -wa -wb >> ${vcffiles_overlap_ran2outs[${SLURM_ARRAY_TASK_ID}]}
bedtools intersect -a mz-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.PNref.bed -b ${vcffiles_overlap_ran2[${SLURM_ARRAY_TASK_ID}]} -wa -wb >> ${vcffiles_overlap_ran2outs[${SLURM_ARRAY_TASK_ID}]}
# NOTE: you will get several ': ambiguous redirect' errors in cases where nothing can be added

# Columns are:
#
# 1.	Pairwise_SNP_Ref_chr - NOTE: this is the Pairwise_SNP_Alt_chr in [species]_pn1 cases e.g. nb1_pn1 - this is the chr for pn1
# 2.	Pairwise_SNP_ Ref_start (0-based) - NOTE: this is the Pairwise_SNP_Alt_start in [species]_pn1 cases e.g. nb1_pn1 - this is the start for pn1
# 3.	Pairwise_SNP_ Ref_end (0-based) - NOTE: this is the Pairwise_SNP_Alt_end in [species]_pn1 cases e.g. nb1_pn1 - this is the end for pn1 > this is the 1-based position of the SNP
# 4.	Pairwise_SNP_comparison_species
# 5.	Pairwise_SNP_quality
# 6.	Pairwise_SNP_strand
# 7.	Pairwise_SNP_Ref_allele
# 8.	Pairwise_SNP_Alt_allele
# 9.	Pairwise_SNP_Alt_chr - NOTE: this is the Pairwise_SNP_Ref_chr in [species]_pn1 cases e.g. nb1_pn1 - this is the chr for nb1
# 10.	Pairwise_SNP_Alt_position (1-based) - NOTE: this is the Pairwise_SNP_Ref_position in [species]_pn1 cases e.g. nb1_pn1 - this is the pos for nb1
# 11.	WG_motif_chr
# 12.	WG_motif_start
# 13.	WG_motif_end
# 14.	cichlid_geneID
# 15.	motif_genesymbolDr
# 16.	motif_genesymbolGa
# 17.	motif_genesymbolSp
# 18.	motif_strand
# 19.	TFmotif_ID
# 20.	TFmotif_cichlidID
# 21.	TFmotif_gene_symbol
# 22.	fimo_score
# 23.	fimo_pval
# 24.	fimo_qval
# 25.	motif_seq
# 26.	conf_level
# 27.	conf_score
# 28.	MalawiVCF_SNP_Ref_chr
# 29.	MalawiVCF_SNP_Ref_position
# 30.	MalawiVCF_SNP_Ref_ID
# 31.	MalawiVCF_SNP_Ref_allele
# 32.	MalawiVCF_SNP_Alt_allele
# 33.	MalawiVCF_SNP_Quality
# 34.	MalawiVCF_SNP_filter
# 35.	MalawiVCF_SNP_info
# 36.	MalawiVCF_SNP_format
# 37 to End of columns. genotype in Lake Malawi species

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap
# check that you have 3407 files
ls -1 pnTFBSVCFoverlap_out/ | wc -l # DONE - ok

### then remove the original outfiles, retaining the ones with a header and amended for actual presence of a SNP (*.out-2)
rm /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/splitvcf/output/*.out # DONE
rm slurm.t* # remove all slurm out files # DONE

##### Then, you need to merge to make one large file and sort based on col1 and col2 sort -k1,1 -k2,2n

nano merge_pnTFBSVCFoverlap.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 0-00:45 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

for i in pnTFBSVCFoverlap_out/*.out-2.candIntersect ; do
  cat $i | sort -k1,1 -k2,2n >> pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs.bed
done

# run the above
sbatch merge_pnTFBSVCFoverlap.sh


### Will is then filtering the above file for whether the SNP is actually present in the promoter alignments (since SNPs generated based on a multiple genome alignment)
cp pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs.bed /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap #23018
# The comparison done is whether the chr + coordinate in Will's file + the ref allele (correcting for reverse complement in cases of -ve strand) correspond to the chr (col1) + coordinates (col2,3) + allele (col7)

### the file can be found here - total of 17015 lines retained from 23018
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates
/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed
/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed

### create two files from above:
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
# 1. Alternative allele same in pairwise species promoter comparison and Malawi VCF ($8==$32) # total is 3806
awk '$8==$32' pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718_sameALTallele.bed
    # total is 3806
# 2. Alternative allele different in pairwise species promoter comparison and Malawi VCF ($8!=$32) - this means that the SNP is totally segregated in Malawi species from Mz and other one of the four comparison cichlid representatives (Pn, Ab, Nb, On)
awk '$8!=$32' pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718_diffALTallele.bed
    # total is 346



#### Another step, change the species genotypes to actual nucleotide alleles
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates
# created a python script to do this: addspgenotypes.py
python addspgenotypes.py # no sys args, has complete paths so run directly
# add column headers
printf 'Pairwise_SNP_Ref_chr[Altchr_in_mz11]\tPairwise_SNP_Ref_start[Altstart_in_mz11]\tPairwise_SNP_Ref_end[Altend_in_mz11]\tPairwise_SNP_comparison_species\tPairwise_SNP_quality\tPairwise_SNP_strand\tPairwise_SNP_Ref_allele\tPairwise_SNP_Alt_allele\tPairwise_SNP_Alt_chr[Refchr_in_mz11]\tPairwise_SNP_Alt_position[Refpos_in_mz11]\tWG_motif_chr\tWG_motif_start\tWG_motif_end\tcichlid_geneID\tmotif_genesymbolDr\tmotif_genesymbolGa\tmotif_genesymbolSp\tmotif_strand\tTFmotif_ID\tTFmotif_cichlidID\tTFmotif_gene_symbol\tfimo_score\tfimo_pval\tfimo_qval\tmotif_seq\tconf_level\tconf_score\tMalawiVCF_SNP_Ref_chr\tMalawiVCF_SNP_Ref_position\tMalawiVCF_SNP_Ref_ID\tMalawiVCF_SNP_Ref_allele\tMalawiVCF_SNP_Alt_allele\tMalawiVCF_SNP_Quality\tMalawiVCF_SNP_filter\tMalawiVCF_SNP_info\tMalawiVCF_SNP_format\t' > Mz_Malawi_colheader.1
printf 'Pairwise_SNP_Ref_chr[Altchr_in_pn1]\tPairwise_SNP_Ref_start[Altstart_in_pn1]\tPairwise_SNP_Ref_end[Altend_in_pn1]\tPairwise_SNP_comparison_species\tPairwise_SNP_quality\tPairwise_SNP_strand\tPairwise_SNP_Ref_allele\tPairwise_SNP_Alt_allele\tPairwise_SNP_Alt_chr[Refchr_in_pn1]\tPairwise_SNP_Alt_position[Refpos_in_pn1]\tWG_motif_chr\tWG_motif_start\tWG_motif_end\tcichlid_geneID\tmotif_genesymbolDr\tmotif_genesymbolGa\tmotif_genesymbolSp\tmotif_strand\tTFmotif_ID\tTFmotif_cichlidID\tTFmotif_gene_symbol\tfimo_score\tfimo_pval\tfimo_qval\tmotif_seq\tconf_level\tconf_score\tVictoriaVCF_SNP_Ref_chr\tVictoriaVCF_SNP_Ref_position\tVictoriaVCF_SNP_Ref_ID\tVictoriaVCF_SNP_Ref_allele\tVictoriaVCF_SNP_Alt_allele\tVictoriaVCF_SNP_Quality\tVictoriaVCF_SNP_filter\tVictoriaVCF_SNP_info\tVictoriaVCF_SNP_format\t' > Pn_Victoria_colheader.1
cat Mz_Malawi_colheader.1 Malinsky_et_al_2017_filteredSNPs.header.vcf | tr -d "\n" > Malinsky_et_al_2017_filteredSNPs.header.vcf1
cat Pn_Victoria_colheader.1 tarangRegions.secondSet.Victoria.colhead.vcf | tr -d "\n" > tarangRegions.secondSet.Victoria.colhead.vcf1
# created some extra colheaders just in notepad to add other headers of GT > *colhead.vcf2
cat Malinsky_et_al_2017_filteredSNPs.header.vcf2 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2 > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3 ; rm mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2 ; mv mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2
cat tarangRegions.secondSet.Victoria.colhead.vcf2 pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2 > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3 ; rm pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2 ; mv pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3 pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2

#### After, this, then added 1 or 0 for whether the Ref or Alt allele are the same as the other lake species
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates
# created a python script to do this: refaltspgenotypematch.py
python refaltspgenotypematch.py # no sys args, has complete paths so run directly
# final files:
# mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2
# mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2
# pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2
# pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2
# add column headers
# created some extra colheaders just in notepad to add other headers of GT and GT-PA (Presence/Absence)
sed -i.bak '1d' mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2
sed -i.bak '1d' mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2
sed -i.bak '1d' pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2
sed -i.bak '1d' pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2
cat Malinsky_et_al_2017_filteredSNPs.header.vcf3 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2 > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3 ; rm mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2 ; mv mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2
cat Malinsky_et_al_2017_filteredSNPs.header.vcf3 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2 > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3 ; rm mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2 ; mv mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2
cat tarangRegions.secondSet.Victoria.colhead.vcf3 pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2 > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3 ; rm pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2 ; mv pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3 pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed2
cat tarangRegions.secondSet.Victoria.colhead.vcf3 pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2 > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3 ; rm pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2 ; mv pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3 pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed2


#### Then, filter the above files out all of the above for the 'candidates' that are from your own work, in other those that are marked as 'SNP_related_' in:
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/9.Candidates/Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt
# the genes are these: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/9.Candidates/Candidates_IDs_SNPrelated.txt2
for i in *-Candidates_TFBSsSNPs_*_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed2 ; do
  grep -wiFvf /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/9.Candidates/Candidates_IDs_SNPrelated.txt2 $i > "$(basename "$i" .bed2).bed3" ;
done

for i in *-Candidates_TFBSsSNPs_*_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.*SpMatch.bed2 ; do
  grep -wiFvf /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/9.Candidates/Candidates_IDs_SNPrelated.txt2 $i > "$(basename "$i" .bed2).bed3" ;
done

mv mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed
for i in *-Candidates_TFBSsSNPs_*_filteredSNPs_SNPsInPromAln.0718.bed ; do
  grep -wiFvf /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/9.Candidates/Candidates_IDs_SNPrelated.txt2 $i > "$(basename "$i" .bed).bed2" ;
done

###### final files:
# mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2
# pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed2

# mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3
# pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3

# mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3
# mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3
# pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3
# pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3


#### Then, do an intersect of the Lake Malawi and Lake Victoria VCF candidate overlaps (filtered by Will for overlap with promoter alignments) to find sites that segregate in both lakes

# For this, do all the following intersects of coordinates, then sort -k1,1 -k2,2n to identify presence in all species comparisons
# pn_mz vs mz_ab
# pn_mz vs mz_nb
# pn_mz vs mz_on
# mz_pn vs pn_ab
# mz_pn vs pn_nb
# mz_pn vs pn_on

# created a python script to do this
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/Victoria-Malawi_Candidates_TFBSsSNPs_filteredSNPsOverlap.py
# /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/Victoria-Malawi_Candidates_TFBSsSNPs_filteredSNPsOverlap.py

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates
# run the script
python Victoria-Malawi_Candidates_TFBSsSNPs_filteredSNPsOverlap.py mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3 pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3 mz11 Mz-Pn_BothLakes_Candidates_TFBSsSNPs_Malawi-Victoria_filteredSNPs_SNPsInPromAln.0718.bed
python Victoria-Malawi_Candidates_TFBSsSNPs_filteredSNPsOverlap.py pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.bed3 pn1 Pn-Mz_BothLakes_Candidates_TFBSsSNPs_Victoria-Malawi_filteredSNPs_SNPsInPromAln.0718.bed
cut -f17-19,23 Mz-Pn_BothLakes_Candidates_TFBSsSNPs_Malawi-Victoria_filteredSNPs_SNPsInPromAln.0718.bed | less
cut -f17-19,23 Pn-Mz_BothLakes_Candidates_TFBSsSNPs_Victoria-Malawi_filteredSNPs_SNPsInPromAln.0718.bed | less
# there are just many foxp2 and sox2 directed TFBSs?
# there is also foxm1 directed to zmiz2

#### Analysis (using Pandas) of individual SNP-TFBS-Promoter-VCF filt files
# created a jupyter-notebook for this: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/SNP-TFBS-Promoter-VCF_PandasAnalysis


#### Looking at examples of segregating sites in species:
cut -f15-17,21,307-441 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3 > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonly.bed3
cut -f15-17,21,243-345 pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3 > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonly.bed3
head -1 mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonly.bed3 | sed $'s/$/\tsum/g' > mz.colhead
head -1 pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonly.bed3 | sed $'s/$/\tsum/g' > pn.colhead
awk '{for(i=1;i<=NF;i++){NUM=NUM?NUM+$i:$i};$(NF+1)=NUM;NUM=""} 1' OFS='\t' mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonly.bed3 | sort -k140,140 -rn > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3 # this sums along the row, sorts high (conserved SNP) to low (segregating SNPs)
awk '{for(i=1;i<=NF;i++){NUM=NUM?NUM+$i:$i};$(NF+1)=NUM;NUM=""} 1' OFS='\t' pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonly.bed3 | sort -k108,108 -rn > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3 # this sums along the row, sorts high (conserved SNP) to low (segregating SNPs)
cat mz.colhead mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3 > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3a; rm mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3 ; mv mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3a mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3
cat pn.colhead pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3 > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3a ; rm pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3 ; mv pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3a pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.PAonlySUM.bed3

#### Looking at candidate genes in individual SNP-TFBS-Promoter-VCF filt files
### In Mz-Malawi, look at sws1, rh2, rho, neurod1, dlx5
## sws1 - focus on mz vs nb variants
grep -wiF 'opn1sw1' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 | cut -f4,15-17,21 | sort -u
# nb1_mz11	opn1sw1	opn1sw1	q2l6a1_oryla	fosl1 - present in both nb and mz
# nb1_mz11	opn1sw1	opn1sw1	q2l6a1_oryla	nkx2-1 - present in both nb and mz
# nb1_mz11	opn1sw1	opn1sw1	q2l6a1_oryla	nr2c2 - absent in mz, present in nb < target this, repressor of nuclear receptor signaling (rxr)
# nb1_mz11	opn1sw1	opn1sw1	q2l6a1_oryla	runx2 - present in both nb and mz
# nb1_mz11	opn1sw1	opn1sw1	q2l6a1_oryla	rxrb - absent in mz, present in nb < target this, photoreceptor neuron protection
  # nr2c2/rxrb: A|A: 3 (2%); A|G: 3 (2%); G|G: 129 (96%) (Nb = G/G, Mz = A|A)
  # scaffold_102	2340621	2340622	nb1_mz11	100	+	C	T	scaffold_1	19181358	scaffold_1	19181348	19181362	nb.gene.s1.386	opn1sw1	opn1sw1	q2l6a1_oryla	+	motif_208	nb.gene.s44.59	nr2c2	22.4429	9.72e-09	0.00603	CAGGGTCAGAGGTCA	2c	0.115	scaffold_102	2340622	.	T	C
  # scaffold_102	2340621	2340622	nb1_mz11	100	+	C	T	scaffold_1	19181358	scaffold_1	19181348	19181361	nb.gene.s1.386	opn1sw1	opn1sw1	q2l6a1_oryla	+	motif_541	nb.gene.s79.100	rxrb	20.3519	7.97e-08	0.0205	AGGGTCAGAGGTCA	2c	0.115	scaffold_102	2340622	.	T	C
# select some higher confidence motifs for presence in Mz only
grep -wiF 'opn1sw1' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 | grep '1b' | cut -f4,15-17,21 | sort -u
# mz11_nb1	opn1sw1	opn1sw1	q2l6a1_oryla	jun < target this, regulate neurodegeneration - important for ocular hypertension-induced retinal ganglion cell degeneration
  # scaffold_102	2340493	2340494	mz11_nb1	100	+	C	T	scaffold_1	19181229	scaffold_102	2340489	2340501	mz.gene.s102.69	opn1sw1	opn1sw1	q2l6a1_oryla	+	motif_190	mz.gene.s127.35	jun	NA	NA	NA	TCTGACGCTTCCA	1b	0.2	scaffold_102	2340494	.	C	T
  # C/C: 5 (3.704%); C/T: 4 (2.963%); T/T: 126 (93.333%) > (Mz: C|C, Nb: T/T) > most of the lake malawi species have same as N. brichardi however, few have same as M. zebra (including closely related mbuna - tropheops_tropheops; Cynotilapia_axelrodi; Cynotilapia_afra)
  # Is this mbuna specific?!?! The VCF suggests it might be
  # grep -wiF 'opn1sw1'  /tgac/workarea/Research-Groups/RG-cichlids/nb-motifenr_mergedmap2d_motifpromANDgenomeannot.txt | grep 'jun' - Absent in Nb
# mz11_nb1	opn1sw1	opn1sw1	q2l6a1_oryla	runx2
# nb1_mz11	opn1sw1	opn1sw1	q2l6a1_oryla	nkx2-1
# nb1_mz11	opn1sw1	opn1sw1	q2l6a1_oryla	runx2


## rho - focus on mz vs ab, mz vs on
grep -wiF 'rho' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 | cut -f4,15-17,21 | sort -u
# ab1_mz11	rho	rho	rho	atrx - found in all, different sites
# ab1_mz11	rho	rho	rho	irf7 - found in mz and ab, absent in on
# ab1_mz11	rho	rho	rho	kdm5b - found in all
# ab1_mz11	rho	rho	rho	stat1 - unique to ab < target this, retina neural function/mediates apoptosis of retinal pericytes
  # stat1: T/T: 134 (99.259%); T/A: 1 (0.741%); A/A: 0 (0.000%) (Ab: A, Mz: T|T) > this SNP is very specific to the Lake Malawi radiation, potential absence of stat1 regulation of rho
  # scaffold_12	3795981	3795982	ab1_mz11	100	-	T	A	scaffold_134	143380	scaffold_134	143366	143389	ab.gene.s134.8	rho	rho	rho	+	motif_123	ab.gene.s32.25	stat1	7.46518	3.48e-05	0.0398	TTCTTTTTAAATTTTGTTTGCTTA	2a	0.125	scaffold_12	3795982	.	T	A
# mz11_ab1	rho	rho	rho	irf7 - found in mz and ab, absent in on
# mz11_on11	rho	rho	rho	irf7 - found in mz and ab, absent in on
# on11_mz11	rho	rho	rho	gata2 - found in all, different sites though > this one present in On, absent in Mz (on vs mz segregates in Malawi) < target this, migration and differentiation of retinorecipient neurons in the superior colliculus
  # gata2: T/T: 13 (9.630%); T/C: 15 (11.111%); C/C: 107 (79.259%) (On: C, Mz: T/C)
  # scaffold_12	3797654	3797655	on11_mz11	100	-	G	A	LG20	13233110	LG20	13233110	13233132	on.gene.LG20.398	rho	rho	rho	+	motif_29	on.gene.LG5.242	gata2	2.80412	1.97e-05	0.017GACTGATTGTAAGTGATTGTGCT	2a	0.125	scaffold_12	3797655	.	T	C
  # grep -wiF 'rho' /tgac/workarea/Research-Groups/RG-cichlids/mz-motifenr_mergedmap2d_motifpromANDgenomeannot.txt | grep 'gata2'
    # There are four gata2 sites in Mz_sws1 however, none overlap the same site in On
    # scaffold_12	3799401	3799423	150	172	-	mz.gene.s12.95	rho	rho	rho	-	motif_30	mz.gene.s10.82	gata2	11.4505	2.21e-05	0.0413	CACAGACAGCATGAATGACAGAC	2a	0.125
    # scaffold_12	4236871	4236893	406	428	+	mz.gene.s12.113	rho	rho	rho	-	motif_30	mz.gene.s10.82	gata2	13.4685	7.62e-06	0.0248	CTGTGAGAGAGAGAGGGACTCTC	2a	0.125
    # scaffold_12	4236873	4236895	404	426	+	mz.gene.s12.113	rho	rho	rho	-	motif_30	mz.gene.s10.82	gata2	15.5856	2.18e-06	0.0248	GTCTGTGAGAGAGAGAGGGACTC	2a	0.125
    # scaffold_12	4236877	4236899	400	422	+	mz.gene.s12.113	rho	rho	rho	-	motif_30	mz.gene.s10.82	gata2	13.0631	9.53e-06	0.0248	GTGCGTCTGTGAGAGAGAGAGGG	2a	0.125

# rh2
grep -wiF 'rh2' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 | cut -f4,15-17,21 | sort -u
# none of interest

## neurod1 - focus on mz vs any other (but vs ab if you can)
grep -wiF 'neurod1' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.bed2 | cut -f4,15-17,21 | sort -u
# none in there



### In Pn-Victoria, look at opn1sw2, neurod1, dlx5
## sws2 - focus on pn vs ab vs on variants
grep -wiF 'opn1sw2' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed2 | cut -f4,15-17,21 | sort -u
# ab1_pn1	opn1sw2	opn1sw2	q4sr56_tetng	irf7 - X (repetitive)
# ab1_pn1	opn1sw2	opn1sw2	q4sr56_tetng	ncor1 – embryogenesis, myelopoiesis
# ab1_pn1	opn1sw2	opn1sw2	q4sr56_tetng	relb - X (repetitive)
# ab1_pn1	opn1sw2	opn1sw2	q4sr56_tetng	stat1 – Retina neural function/mediates apoptosis of retinal pericytes - X (repetitive)
# ab1_pn1	opn1sw2	opn1sw2	q4sr56_tetng	znf384 - X (repetitive)
# on11_pn1	opn1sw2	opn1sw2	q4sr56_tetng	satb1 – retinal ganglion cell patterning (has eye expr in all species) - would like to target however, there are satb1 sites in pn promoter and could be repetitive
# pn1_ab1	opn1sw2	opn1sw2	q4sr56_tetng	irf7 - X (repetitive)
# pn1_ab1	opn1sw2	opn1sw2	q4sr56_tetng	ncor1 – embryogenesis, myelopoiesis - X (repetitive)
# pn1_ab1	opn1sw2	opn1sw2	q4sr56_tetng	satb1 – retinal ganglion cell patterning
# pn1_ab1	opn1sw2	opn1sw2	q4sr56_tetng	znf384 - X (repetitive)

## neurod1 - focus on pn vs any other (but vs ab if you can)
grep -wiF 'neurod1' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.bed2 | cut -f4,15-17,21 | sort -u
# none of interest





### THIS NEEDS RE-RUNNING
# prep VCF files for IGV - NOTE: THE VCF FILES ARE NOT LOADING INTO IGV, TOO FEW COLS ON LINE 6911?!?!
sort -k1,1 -k2,2n mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_sorted.bed
sed '/#CHROM/ q' /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf > mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_sorted.vcf # get all headers of VCF
sort -k1,1 -k2,2n mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs.bed | awk '{print $1,$2,$4"-"$14"-"$15"-"$19"-"$20"-"$21,$7,$8,$5,"PASS",$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55,$56,$57,$58,$59,$60,$61,$62,$63,$64,$65,$66,$67,$68,$69,$70,$71,$72,$73,$74,$75,$76,$77,$78,$79,$80,$81,$82,$83,$84,$85,$86,$87,$88,$89,$90,$91,$92,$93,$94,$95,$96,$97,$98,$99,$100,$101,$102,$103,$104,$105,$106,$107,$108,$109,$110,$111,$112,$113,$114,$115,$116,$117,$118,$119,$120,$121,$122,$123,$124,$125,$126,$127,$128,$129,$130,$131,$132,$133,$134,$135,$136,$137,$138,$139,$140,$141,$142,$143,$144,$145,$146,$147,$148,$149,$150,$151,$152,$153,$154,$155,$156,$157,$158,$159,$160,$161,$162,$163,$164,$165,$166,$167,$168,$169,$170,$171}' >> mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_sorted.vcf

sort -k1,1 -k2,2n pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs.bed > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_sorted.bed
sed '/#CHROM/ q' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/tarangRegions_Pn-v1liftover.vcf > pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_sorted.vcf # get all headers of VCF
sort -k1,1 -k2,2n pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs.bed | awk '{print $1,$2,$4"-"$14"-"$15"-"$19"-"$20"-"$21,$7,$8,$5,"PASS",$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55,$56,$57,$58,$59,$60,$61,$62,$63,$64,$65,$66,$67,$68,$69,$70,$71,$72,$73,$74,$75,$76,$77,$78,$79,$80,$81,$82,$83,$84,$85,$86,$87,$88,$89,$90,$91,$92,$93,$94,$95,$96,$97,$98,$99,$100,$101,$102,$103,$104,$105,$106,$107,$108,$109,$110,$111,$112,$113,$114,$115,$116,$117,$118,$119,$120,$121,$122,$123,$124,$125,$126,$127,$128,$129,$130,$131,$132,$133,$134,$135,$136,$137,$138,$139}' >> pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_sorted.vcf

# can load *_filteredSNPs_sorted.vcf into IGV with genome + annotation + sortedbed

## search for TFs that you know have retinal function and opsin TG
grep -i 'atrx' mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_sorted.bed | grep 'opn' # this is however between mz and on
awk '$4=="mz11_nb1" || $4=="nb1_mz11"' mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_sorted.bed | grep 'opn1sw1' | cut -f15,21




##########################################################################################################################
#
# Re-annotate some TFs to complete the protein/cds/cdna
#
##########################################################################################################################

# Nb nr2c2 - nb.mrna.s44.59.1
# incomplete start (as compared by aligning to Hs, Mm and Mz)
# cds is scaffold_44:2,361,054-2,367,852
# take scaffold_44:2,330,000-2,373,000 for prediction of new protein
# extract genomic sequence
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
echo -e 'scaffold_44\t2330000\t2373000' > nb_scaff44_2330000-2373000.bed
ml bedtools/2.25.0
ml zlib
bedtools getfasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -bed nb_scaff44_2330000-2373000.bed -fo nb_scaff44_2330000-2373000_nr2c2pred.fasta
cp nb_scaff44_2330000-2373000_nr2c2pred.fasta /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
# using FGENESH+, input the genome sequence and mz_nr2c2 protein to predict the gene

# Mz nr2c2 - mz.mrna.s12.107.1
# get exons
# cds is scaffold_12:4035215-4042366
# take scaffold_12:4035000-4044000 for prediction of new protein
# extract genomic sequence
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
echo -e 'scaffold_12\t4035000\t4044000' > mz_scaff12_4035000-4044000.bed
ml bedtools/2.25.0
ml zlib
bedtools getfasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -bed mz_scaff12_4035000-4044000.bed -fo mz_scaff12_4035000-4044000_nr2c2pred.fasta
cp mz_scaff12_4035000-4044000_nr2c2pred.fasta /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
# using FGENESH+, input the genome sequence and mz_nr2c2 protein to predict the gene
echo -e 'scaffold_12:4035000-4044000\t1052\t1158' > mz_scaff12_4035000-4044000_nr2c2_exon4_5_6.bed
echo -e 'scaffold_12:4035000-4044000\t1292\t1472' >> mz_scaff12_4035000-4044000_nr2c2_exon4_5_6.bed
echo -e 'scaffold_12:4035000-4044000\t1946\t2103' >> mz_scaff12_4035000-4044000_nr2c2_exon4_5_6.bed
bedtools getfasta -fi mz_scaff12_4035000-4044000_nr2c2pred.fasta -bed mz_scaff12_4035000-4044000_nr2c2_exon4_5_6.bed -fo mz_scaff12_4035000-4044000_nr2c2_exon4_5_6.fasta
cp mz_scaff12_4035000-4044000_nr2c2_exon4_5_6.fasta /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/


# Nb rxrb - nb.gene.s79.100
# cds is scaffold_79:2225724-2236786
# take scaffold_79:2224500-2242000 for prediction of new protein
# extract genomic sequence
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
echo -e 'scaffold_79\t2220000\t2242000' > nb_scaff79_2220000-2242000.bed
ml bedtools/2.25.0
ml zlib
bedtools getfasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -bed nb_scaff79_2220000-2242000.bed -fo nb_scaff79_2220000-2242000_rxrbpred.fasta
cp nb_scaff79_2220000-2242000_rxrbpred.fasta /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
# using FGENESH+, input the genome sequence and mz_nr2c2 protein to predict the gene

# Mz rxrb - mz.gene.s4.354
# get exons
# cds is scaffold_4:9166655-9177017
# take scaffold_4:9160000-9185000 for prediction of new protein
# extract genomic sequence
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
echo -e 'scaffold_4\t9160000\t9185000' > mz_scaff4_9160000-9185000.bed
ml bedtools/2.25.0
ml zlib
bedtools getfasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -bed mz_scaff4_9160000-9185000.bed -fo mz_scaff4_9160000-9185000_rxrbpred.fasta
cp mz_scaff4_9160000-9185000_rxrbpred.fasta /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
# using FGENESH+, input the genome sequence and mz_nr2c2 protein to predict the gene

# blast nr2c2 and rxrb primers aginst cds SEQUENCES
# only blast the gene specific region of primer
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
echo '>mz-nb_T7_nr2c2DBD_F1' > nr2c2_rxrbDBDprimers.fasta
echo 'TCAGGAGATTTGAGCCGACCA' >> nr2c2_rxrbDBDprimers.fasta
echo 'mz-nb_T7_nr2c2DBD_R1' >> nr2c2_rxrbDBDprimers.fasta
echo 'GGGCACAATGTCGATGGGTTTCCTCTCACTCTGGACAGACTCAGTCTTCATCCCCAT' >> nr2c2_rxrbDBDprimers.fasta
echo 'mz-nb_T7_rxrbDBD_F1' >> nr2c2_rxrbDBDprimers.fasta
echo 'GCTCACAGCCCGGGAATAATG' >> nr2c2_rxrbDBDprimers.fasta
echo 'mz-nb_T7_rxrbDBD_R1' >> nr2c2_rxrbDBDprimers.fasta
echo 'CTGTCGTTCCTCTTGTACCGCTTCCCTCTTCATTCCCATGGCCAGGCACTTCTGGTA' >> nr2c2_rxrbDBDprimers.fasta

# create blast databases of longest cds
ml blast
cd /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/
makeblastdb -in Mz_alllongestcds.fasta -parse_seqids -dbtype nucl
makeblastdb -in Pn_alllongestcds.fasta -parse_seqids -dbtype nucl
makeblastdb -in Ab_alllongestcds.fasta -parse_seqids -dbtype nucl
makeblastdb -in Nb_alllongestcds.fasta -parse_seqids -dbtype nucl
makeblastdb -in On_alllongestcds.fasta -parse_seqids -dbtype nucl

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/
blastn -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_alllongestcds.fasta -outfmt 6 -evalue 1e-1 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 -out nr2c2_rxrbDBDprimers_blastn-Mzlongestcds.txt -query nr2c2_rxrbDBDprimers.fasta
blastn -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_alllongestcds.fasta -outfmt 6 -evalue 1e-1 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 -out nr2c2_rxrbDBDprimers_blastn-Nblongestcds.txt -query nr2c2_rxrbDBDprimers.fasta
cp nr2c2_rxrbDBDprimers_blastn-*longestcds.txt /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/


# Mz JUN
# no annotated nb protein in our orthogroups
# have a look if exists in main paper orthogrouping - it does as nb.gene.s10.3630001
grep -wiF mz.gene.s127.35 /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/protein_coding/v2/full_orthologs # nb.gene.s10.3630001
less /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/mz/a/mz.gene.s127.35.fa.NbGenome.blast # this is too see whether there is a hit via BLAST of the mz gene on Nb genome


##########################################################################################################################
#
# Scan scrambled negative control oligos for all motifs to ensure no motifs created by chance
#
##########################################################################################################################

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/

echo '>nb_cy5_nr2c2_rxrb-sws1TG-ve_F1' > scrambled_oligos_nb.fa
echo 'TTGGTGAAAAGTGACTGGCTGCAGCAAC' >> scrambled_oligos_nb.fa
echo '>mz_cy5_jun-sws1TG-ve_F1' > scrambled_oligos_mz.fa
echo 'TGGTCCTGTAGTCTAAAACTTCCTTG' >> scrambled_oligos_mz.fa
echo '>ab_cy5_stat1-rhoTG-ve_F1' > scrambled_oligos_ab.fa
echo 'TAATAAAATAAGAAAATACTGAATAAAGGATAAATAC' >> scrambled_oligos_ab.fa
echo '>on_cy5_gata2-rhoTG-ve_F1' > scrambled_oligos_on.fa
echo 'TCAACTGAAACATCCGCAAAACACAAGTTTCGTATT' >> scrambled_oligos_on.fa


nano fimo_scan_oligos.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the (three confidence level) meme files created by Will, at a level of statistical significance (use relaxed parameters of any motif instances of p-value <1e-4)

# N. brichardi
fimo --o scrambled_oligos_nb_CS-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2a_CS_nb.meme scrambled_oligos_nb.fa
fimo --o scrambled_oligos_nb_CW-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2b_CW_nb.meme scrambled_oligos_nb.fa
fimo --o scrambled_oligos_nb_JASPAR-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2c_JASPAR_nb.meme scrambled_oligos_nb.fa

# M. zebra
fimo --o scrambled_oligos_mz_CS-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2a_CS_mz.meme scrambled_oligos_mz.fa
fimo --o scrambled_oligos_mz_CW-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2b_CW_mz.meme scrambled_oligos_mz.fa
fimo --o scrambled_oligos_mz_JASPAR-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2c_JASPAR_mz.meme scrambled_oligos_mz.fa

# A. burtoni
fimo --o scrambled_oligos_ab_CS-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2a_CS_ab.meme scrambled_oligos_ab.fa
fimo --o scrambled_oligos_ab_CW-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2b_CW_ab.meme scrambled_oligos_ab.fa
fimo --o scrambled_oligos_ab_JASPAR-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2c_JASPAR_ab.meme scrambled_oligos_ab.fa

# O. niloticus
fimo --o scrambled_oligos_on_CS-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2a_CS_on.meme scrambled_oligos_on.fa
fimo --o scrambled_oligos_on_CW-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2b_CW_on.meme scrambled_oligos_on.fa
fimo --o scrambled_oligos_on_JASPAR-FIMO --thresh 1e-4 --max-stored-scores 500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/2c_JASPAR_on.meme scrambled_oligos_on.fa

# run the above
qsub -q Test -l select=1:mem=40GB:ncpus=1 fimo_scan_oligos.sh

# You ideally want empty output files, check with this
ls -tlrh *-FIMO/fimo.txt # OK!

##########################################################################################################################
#
# 2.	Focus on all genes
#
##########################################################################################################################

##########################################################################################################################
#
# a.	Pull out all TFBS absence pairwise comparisons from Regulon master tables
#
##########################################################################################################################

regulon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/v1_regulonanalysis)
cd $regulon

for i in Ab Mz Nb On Pn ; do
  awk '$6==1 && $7==1 && $37==1 && $38==1' $i-wholegenome_Regulon_v3.txt > $i-wholegenome_Regulon_v3.txt3 # TF,TG present in genome and modules
done

########################################################################################
# ~~~~~~~ # ai. non-module rewired TFs
########################################################################################

nano nonmodrewiredTF_TFBSabsent.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

abreg=(Ab-wholegenome_Regulon_v3.txt)
mzreg=(Mz-wholegenome_Regulon_v3.txt)
nbreg=(Nb-wholegenome_Regulon_v3.txt)
onreg=(On-wholegenome_Regulon_v3.txt)
pnreg=(Pn-wholegenome_Regulon_v3.txt)

for i in Ab ; do
  awk '$24==0 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzAbTFnoswitch_MzTFBSabsent
  awk '$27==0 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnAbTFnoswitch_PnTFBSabsent
  awk '$30==0 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbAbTFnoswitch_NbTFBSabsent
  awk '$31==0 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnAbTFnoswitch_OnTFBSabsent
done

for i in Mz ; do
  awk '$23==0 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnMzTFnoswitch_PnTFBSabsent
  awk '$24==0 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbMzTFnoswitch_AbTFBSabsent
  awk '$25==0 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbMzTFnoswitch_NbTFBSabsent
  awk '$26==0 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnMzTFnoswitch_OnTFBSabsent
done

for i in Nb ; do
  awk '$25==0 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzNbTFnoswitch_MzTFBSabsent
  awk '$28==0 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnNbTFnoswitch_PnTFBSabsent
  awk '$30==0 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbNbTFnoswitch_AbTFBSabsent
  awk '$32==0 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnNbTFnoswitch_OnTFBSabsent
done

for i in On ; do
  awk '$26==0 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzOnTFnoswitch_MzTFBSabsent
  awk '$29==0 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnOnTFnoswitch_PnTFBSabsent
  awk '$31==0 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbOnTFnoswitch_AbTFBSabsent
  awk '$32==0 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbOnTFnoswitch_NbTFBSabsent
done

for i in Pn ; do
  awk '$23==0 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzPnTFnoswitch_MzTFBSabsent
  awk '$27==0 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbPnTFnoswitch_AbTFBSabsent
  awk '$28==0 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbPnTFnoswitch_NbTFBSabsent
  awk '$29==0 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnPnTFnoswitch_OnTFBSabsent
done

# run the above
qsub -q Test -l select=1:mem=30GB:ncpus=1 nonmodrewiredTF_TFBSabsent.sh

########################################################################################
# ~~~~~~~ # aii. module rewired TFs
########################################################################################

regulon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/v1_regulonanalysis)
cd $regulon

nano modrewiredTF_TFBSabsent.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

abreg=(Ab-wholegenome_Regulon_v3.txt)
mzreg=(Mz-wholegenome_Regulon_v3.txt)
nbreg=(Nb-wholegenome_Regulon_v3.txt)
onreg=(On-wholegenome_Regulon_v3.txt)
pnreg=(Pn-wholegenome_Regulon_v3.txt)

for i in Ab ; do
  awk '$24==1 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzAbTFswitch_MzTFBSabsent
  awk '$27==1 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnAbTFswitch_PnTFBSabsent
  awk '$30==1 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbAbTFswitch_NbTFBSabsent
  awk '$31==1 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnAbTFswitch_OnTFBSabsent
done

for i in Mz ; do
  awk '$23==1 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnMzTFswitch_PnTFBSabsent
  awk '$24==1 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbMzTFswitch_AbTFBSabsent
  awk '$25==1 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbMzTFswitch_NbTFBSabsent
  awk '$26==1 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnMzTFswitch_OnTFBSabsent
done

for i in Nb ; do
  awk '$25==1 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzNbTFswitch_MzTFBSabsent
  awk '$28==1 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnNbTFswitch_PnTFBSabsent
  awk '$30==1 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbNbTFswitch_AbTFBSabsent
  awk '$32==1 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnNbTFswitch_OnTFBSabsent
done

for i in On ; do
  awk '$26==1 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzOnTFswitch_MzTFBSabsent
  awk '$29==1 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnOnTFswitch_PnTFBSabsent
  awk '$31==1 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbOnTFswitch_AbTFBSabsent
  awk '$32==1 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbOnTFswitch_NbTFBSabsent
done

for i in Pn ; do
  awk '$23==1 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzPnTFswitch_MzTFBSabsent
  awk '$27==1 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbPnTFswitch_AbTFBSabsent
  awk '$28==1 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbPnTFswitch_NbTFBSabsent
  awk '$29==1 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnPnTFswitch_OnTFBSabsent
done

# run the above
qsub -q Test -l select=1:mem=30GB:ncpus=1 modrewiredTF_TFBSabsent.sh

########################################################################################
# ~~~~~~~ # aiii. non-module rewired TG
########################################################################################

regulon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/v1_regulonanalysis)
cd $regulon

nano nonmodrewiredTG_TFBSabsent.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

abreg=(Ab-wholegenome_Regulon_v3.txt)
mzreg=(Mz-wholegenome_Regulon_v3.txt)
nbreg=(Nb-wholegenome_Regulon_v3.txt)
onreg=(On-wholegenome_Regulon_v3.txt)
pnreg=(Pn-wholegenome_Regulon_v3.txt)

for i in Ab ; do
  awk '$55==0 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzAbTGnoswitch_MzTFBSabsent
  awk '$58==0 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnAbTGnoswitch_PnTFBSabsent
  awk '$61==0 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbAbTGnoswitch_NbTFBSabsent
  awk '$62==0 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnAbTGnoswitch_OnTFBSabsent
done

for i in Mz ; do
  awk '$54==0 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnMzTGnoswitch_PnTFBSabsent
  awk '$55==0 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbMzTGnoswitch_AbTFBSabsent
  awk '$56==0 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbMzTGnoswitch_NbTFBSabsent
  awk '$57==0 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnMzTGnoswitch_OnTFBSabsent
done

for i in Nb ; do
  awk '$56==0 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzNbTGnoswitch_MzTFBSabsent
  awk '$59==0 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnNbTGnoswitch_PnTFBSabsent
  awk '$61==0 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbNbTGnoswitch_AbTFBSabsent
  awk '$63==0 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnNbTGnoswitch_OnTFBSabsent
done

for i in On ; do
  awk '$57==0 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzOnTGnoswitch_MzTFBSabsent
  awk '$60==0 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnOnTGnoswitch_PnTFBSabsent
  awk '$62==0 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbOnTGnoswitch_AbTFBSabsent
  awk '$63==0 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbOnTGnoswitch_NbTFBSabsent
done

for i in Pn ; do
  awk '$54==0 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzPnTGnoswitch_MzTFBSabsent
  awk '$58==0 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbPnTGnoswitch_AbTFBSabsent
  awk '$59==0 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbPnTGnoswitch_NbTFBSabsent
  awk '$60==0 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnPnTGnoswitch_OnTFBSabsent
done

# run the above
qsub -q Test -l select=1:mem=30GB:ncpus=1 nonmodrewiredTG_TFBSabsent.sh


########################################################################################
# ~~~~~~~ # aiv. module rewired TG
########################################################################################

regulon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/v1_regulonanalysis)
cd $regulon

nano modrewiredTG_TFBSabsent.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

abreg=(Ab-wholegenome_Regulon_v3.txt)
mzreg=(Mz-wholegenome_Regulon_v3.txt)
nbreg=(Nb-wholegenome_Regulon_v3.txt)
onreg=(On-wholegenome_Regulon_v3.txt)
pnreg=(Pn-wholegenome_Regulon_v3.txt)

for i in Ab ; do
  awk '$55==1 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzAbTGswitch_MzTFBSabsent
  awk '$58==1 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnAbTGswitch_PnTFBSabsent
  awk '$61==1 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbAbTGswitch_NbTFBSabsent
  awk '$62==1 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnAbTGswitch_OnTFBSabsent
done

for i in Mz ; do
  awk '$54==1 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnMzTGswitch_PnTFBSabsent
  awk '$55==1 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbMzTGswitch_AbTFBSabsent
  awk '$56==1 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbMzTGswitch_NbTFBSabsent
  awk '$57==1 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnMzTGswitch_OnTFBSabsent
done

for i in Nb ; do
  awk '$56==1 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzNbTGswitch_MzTFBSabsent
  awk '$59==1 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnNbTGswitch_PnTFBSabsent
  awk '$61==1 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbNbTGswitch_AbTFBSabsent
  awk '$63==1 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnNbTGswitch_OnTFBSabsent
done

for i in On ; do
  awk '$57==1 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzOnTGswitch_MzTFBSabsent
  awk '$60==1 && $67==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-PnOnTGswitch_PnTFBSabsent
  awk '$62==1 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbOnTGswitch_AbTFBSabsent
  awk '$63==1 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbOnTGswitch_NbTFBSabsent
done

for i in Pn ; do
  awk '$54==1 && $66==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-MzPnTGswitch_MzTFBSabsent
  awk '$58==1 && $68==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-AbPnTGswitch_AbTFBSabsent
  awk '$59==1 && $69==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-NbPnTGswitch_NbTFBSabsent
  awk '$60==1 && $70==0' $i-wholegenome_Regulon_v3.txt3 > $i-wholegenome_Regulon_v3.txt3-OnPnTGswitch_OnTFBSabsent
done

# run the above
qsub -q Test -l select=1:mem=30GB:ncpus=1 modrewiredTG_TFBSabsent.sh

########################################################################################
# ~~~~~~~ # b. Overlap absent TFBSs with presence of SNP in SNP overlaps of TFBSs
########################################################################################

# For this, you should focus on target genes (col34) only since they have the TFBSs
# col34 in above files and col16 in SNP files

nano TFBSabsent_SNPoverlap.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

ls -1 *TG*TFBSabsent > TFTFBS
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_mz11_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt' > SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_mz11_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_nb1_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_nb1_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_on11_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_on11_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_pn1_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_pn1_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_ab1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_ab1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_nb1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_nb1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_on11_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_on11_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_pn1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_pn1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_ab1_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_ab1_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_mz11_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_mz11_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_on11_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_on11_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_pn1_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_pn1_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_ab1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_ab1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_mz11_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_mz11_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_nb1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_nb1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_pn1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_pn1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_ab1_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_ab1_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_mz11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_mz11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_nb1_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_nb1_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_on11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_on11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt' >> SNPfiles

while read -u 3 -r file1 && read -u 4 -r file2 ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$16]=$0;next}{if(a[$34]){print $0,a[$34];}else{print $0,"REMOVEME";}}' ${file1} ${file2} | grep -v 'REMOVEME' > ${file2}_SNPoverlap
done 3<SNPfiles 4<TFTFBS

qsub -q Test -l select=1:mem=60GB:ncpus=1 TFBSabsent_SNPoverlap.sh

########################################################################################
# ~~~~~~~ # c. Overlap all candidate genes with absent TFBSs with SNP overlaps - focus on switched target genes
########################################################################################

nano grepcandTGswitch.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

cut -f2 ../../Edge_Attributes/Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt | grep -v 'NULL' | xargs -i grep -wiF {} Mz*TGswitch_*TFBSabsent_SNPoverlap > Mz_Candidates_IDs_Fast_Opsins_Hahn_SNPs-TGswitch_TFBSabsent_SNPoverlap.txt
cut -f3 ../../Edge_Attributes/Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt | grep -v 'NULL' | xargs -i grep -wiF {} Pn*TGswitch_*TFBSabsent_SNPoverlap > Pn_Candidates_IDs_Fast_Opsins_Hahn_SNPs-TGswitch_TFBSabsent_SNPoverlap.txt
cut -f4 ../../Edge_Attributes/Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt | grep -v 'NULL' | xargs -i grep -wiF {} Ab*TGswitch_*TFBSabsent_SNPoverlap > Ab_Candidates_IDs_Fast_Opsins_Hahn_SNPs-TGswitch_TFBSabsent_SNPoverlap.txt
cut -f5 ../../Edge_Attributes/Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt | grep -v 'NULL' | xargs -i grep -wiF {} Nb*TGswitch_*TFBSabsent_SNPoverlap > Nb_Candidates_IDs_Fast_Opsins_Hahn_SNPs-TGswitch_TFBSabsent_SNPoverlap.txt
cut -f6 ../../Edge_Attributes/Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt | grep -v 'NULL' | xargs -i grep -wiF {} On*TGswitch_*TFBSabsent_SNPoverlap > On_Candidates_IDs_Fast_Opsins_Hahn_SNPs-TGswitch_TFBSabsent_SNPoverlap.txt

qsub -q Test -l select=1:mem=60GB:ncpus=1 grepcandTGswitch.sh

# within these files, you shoud focus on partiular species pairwise comparisons

# Ab vs rest - any TFs in literature implicated with behavioural phenotype > look for egr1, cfos
# irx1b - forebrain
# cntn4 - neural system and plasticity
# foxp2 - TF; developing neural, gastrointestinal and cardiovascular tissues < egr1 targets
# bmpr1a - morphogen important for neural development, forebrain specification, left-right asymmetry
# gdf10b - brain expression, axonal outgrowth < egr1 targets


# sws1 - mz vs nb
# sws2 - ab vs pn vs on
# rh2
# rho - module variation and under selection; all vs all
