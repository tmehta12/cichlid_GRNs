#!/bin/sh


##########################################################################################################################
# SNP overlap of CNEs, aCNEs, 5kb promoters (TFBSs) and 3' UTRs (miRNA binding sites) - March 2018
#
# 1. Overlap SNPs with different regulatory regions > plot {DONE}
# 2. Test if any pairwise SNPs overlap TFBS motif scans in promoter regions {DONE}
# 3. Run GO enrichment of variant associated regulatory regions {DONE} - AS OF 21/03/18, STILL NEED TO FINISH PLOTTING
# 4. Analyse variants associated with candidate genes
# 5. Plot pairwise SNPs against Euclidian expression divergence {DONE} - AS OF 21/03/18, STILL NEED TO FINISH PLOTTING
#
##########################################################################################################################

########## 1. Overlap SNPs with different regulatory regions > plot {DONE}

# Will used cichlid_5way.maf generated as part of the publication to generate pairwise substitutions of 4 species against O. niloticus
# All vcf files can be found here: /tgac/workarea/group-vh/cichlids/substitutions_VCFs/

# Using bedtools, intersect vcf files with genomic regions of interest
slurm
interactive --mem=8G

cd /tgac/workarea/group-vh/cichlids/substitutions_VCFs/

## have other pairwise variants in another folder, create symbolic links to them here
for i in /tgac/workarea/group-vh/cichlids/MSA/MetZeb1.1prescreen_centred_8wayteleostMultiz/*_1217.vcf ; do ln -s $i ; done
#for i in /tgac/workarea/Research-Groups/RG-cichlids/variants/* ; do ln -s $i ; done

# Amend the chr columns in vcf to have proper annotations then change format into bed (0-based - start-1, start)
# create a shell script to run on all and change in place
echo '#!/bin/bash -e' > sed_me2.sh
echo '#SBATCH -p tgac-medium # partition (queue)' >> sed_me2.sh
echo '#SBATCH -N 1 # number of nodes' >> sed_me2.sh
echo '#SBATCH -n 1 # number of tasks' >> sed_me2.sh
echo '#SBATCH --mem 12000 # memory pool for all cores' >> sed_me2.sh
echo '#SBATCH -t 0-02:59 # time (D-HH:MM)' >> sed_me2.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> sed_me2.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> sed_me2.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> sed_me2.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> sed_me2.sh
printf "\n" >> sed_me2.sh
ls -1 *_1217.vcf | sed 's/_substitutions_1217.vcf//g' | awk -F'_' '{print $1,$2}' | sed "s/^/sed_'s\//g" | sed "s/$/.\/\/g/g" | awk -F' ' '{gsub(/^/,"sed_s/",$2); print $0}' OFS='\t' | awk -F' ' '{gsub(/$/,".//g",$1); print $0}' OFS='\t' | sed "s/sed_/sed /g" | sed "s/sed/sed -i/g" | sed "s/sed -i s/sed -i 's/g" | sed "s/.\/\/g/.\/\/g'/g" > sed_me.sh
ls -1 *_1217.vcf > listfiles ; paste -d'\t' sed_me.sh listfiles | sed "s/ /_/g" | awk '{print $1,$3,$2,$3}' | sed "s/ sed/; sed/g" | sed "s/sed_-i_/sed -i /g" >> sed_me2.sh

rm sed_me.sh
rm listfiles

# run the above
sbatch sed_me2.sh

nano prepme.sh # run as an array

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-55
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-5:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 *_1217.vcf > vcffiles # create a list of all files
mapfile -t vcffiles < vcffiles # assign as elements to $vcffiles variable
awk '{print $1,$2,$2,$3,$4,$5,$6,$7,$8,$9,$10}' OFS='\t' ${vcffiles[${SLURM_ARRAY_TASK_ID}]} | awk '{$2-=1;}1' OFS='\t' | awk '{OFS="\t"} {if ($2=="-1") $2="0"; print $0}' | awk '{OFS="\t"} {if ($3=="0") $3="1"; print $0}' | sort -k1,1 -k2,2n > "$(basename "${vcffiles[${SLURM_ARRAY_TASK_ID}]}" .vcf).bed"

# for i in *_1217.vcf ; do awk '{print $1,$2,$2,$3,$4,$5,$6,$7,$8,$9,$10}' OFS='\t' $i | awk '{$2-=1;}1' OFS='\t' | awk '{OFS="\t"} {if ($2=="-1") $2="0"; print $0}' | awk '{OFS="\t"} {if ($3=="0") $3="1"; print $0}' | sort -k1,1 -k2,2n > "$(basename "$i" .vcf).bed" ; done

# for i in cichlid_Niloticus_*.vcf ; do sed 's/Niloticus://g' $i | awk '{print $1,$2,$2,$3,$4,$5,$6,$7,$8}' OFS='\t' | awk '{$2-=1;}1' OFS='\t' | awk '{OFS="\t"} {if ($2=="-1") $2="0"; print $0}' | awk '{OFS="\t"} {if ($3=="0") $3="1"; print $0}' > "$(basename "$i" .vcf).bed" ; done
# for i in cichlid_Brichardi_*.vcf ; do sed 's/Brichardi://g' $i | awk '{print $1,$2,$2,$3,$4,$5,$6,$7,$8}' OFS='\t' | awk '{$2-=1;}1' OFS='\t' | awk '{OFS="\t"} {if ($2=="-1") $2="0"; print $0}' | awk '{OFS="\t"} {if ($3=="0") $3="1"; print $0}' > "$(basename "$i" .vcf).bed" ; done
# for i in cichlid_Burtoni_*.vcf ; do sed 's/Burtoni://g' $i | awk '{print $1,$2,$2,$3,$4,$5,$6,$7,$8}' OFS='\t' | awk '{$2-=1;}1' OFS='\t' | awk '{OFS="\t"} {if ($2=="-1") $2="0"; print $0}' | awk '{OFS="\t"} {if ($3=="0") $3="1"; print $0}' > "$(basename "$i" .vcf).bed" ; done
# for i in cichlid_Nyererei_*.vcf ; do sed 's/Nyererei://g' $i | awk '{print $1,$2,$2,$3,$4,$5,$6,$7,$8}' OFS='\t' | awk '{$2-=1;}1' OFS='\t' | awk '{OFS="\t"} {if ($2=="-1") $2="0"; print $0}' | awk '{OFS="\t"} {if ($3=="0") $3="1"; print $0}' > "$(basename "$i" .vcf).bed" ; done
# for i in cichlid_Zebra_*.vcf ; do sed 's/Zebra://g' $i | awk '{print $1,$2,$2,$3,$4,$5,$6,$7,$8}' OFS='\t' | awk '{$2-=1;}1' OFS='\t' | awk '{OFS="\t"} {if ($2=="-1") $2="0"; print $0}' | awk '{OFS="\t"} {if ($3=="0") $3="1"; print $0}' > "$(basename "$i" .vcf).bed" ; done

# run the above
sbatch prepme.sh

# Do all pairwise intersects

nano subst_intersect.sh # add below

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 3-15:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

ml bedtools/2.25.0

for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k_* ; do sort -k1,1 -k2,2n $i > /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/"$(basename "$i" .bed)_sorted.bed" ; done

# CNEs
for i in mz11*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k_Mzebra1.1-Onliftover_sorted.bed -wb > ${i}_substCNEoverlap.txt ; done
for i in pn1*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k_Pnyererei-Onliftover_sorted.bed -wb > ${i}_substCNEoverlap.txt ; done
for i in ab1*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k_Aburtoni-Onliftover_sorted.bed -wb > ${i}_substCNEoverlap.txt ; done
for i in nb1*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k_Nbrichardi-Onliftover_sorted.bed -wb > ${i}_substCNEoverlap.txt ; done
for i in on11*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k.amend.bed -wb > ${i}_substCNEoverlap.txt ; done
# bedtools intersect -a cichlid_Niloticus_Zebra_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k.amend.bed -wb > On-Mz_substCNEoverlap.txt # 2380557/47473955
# bedtools intersect -a cichlid_Niloticus_Nyererei_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k.amend.bed -wb > On-Pn_substCNEoverlap.txt  # 2463809/49606035
# bedtools intersect -a cichlid_Niloticus_Burtoni_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k.amend.bed -wb  > On-Ab_substCNEoverlap.txt # 2419585/48920232
# bedtools intersect -a cichlid_Niloticus_Brichardi_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/CNE.528k.amend.bed -wb  > On-Nb_substCNEoverlap.txt # 2638953/52733124

# aCNEs
for i in mz11*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/aCNEs/aCNEs_Mzebra1.1-Onliftover_true_aCNEs.sorted.bed -wb > ${i}_substaCNEoverlap.txt ; done
for i in pn1*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/aCNEs/aCNEs_Pnyererei-Onliftover_true_aCNEs.sorted.bed -wb > ${i}_substaCNEoverlap.txt ; done
for i in ab1*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/aCNEs/aCNEs_Aburtoni-Onliftover_true_aCNEs.sorted.bed -wb > ${i}_substaCNEoverlap.txt ; done
for i in nb1*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/aCNEs/aCNEs_Nbrichardi-Onliftover_true_aCNEs.sorted.bed -wb > ${i}_substaCNEoverlap.txt ; done
for i in on11*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/aCNEs/aCNE.sorted.bed -wb > ${i}_substaCNEoverlap.txt ; done
# bedtools intersect -a cichlid_Niloticus_Zebra_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/aCNEs/aCNE.sorted.bed -wb > On-Mz_substaCNEoverlap.txt  # 1828/47473955
# bedtools intersect -a cichlid_Niloticus_Nyererei_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/aCNEs/aCNE.sorted.bed -wb > On-Pn_substaCNEoverlap.txt # 2052/49606035
# bedtools intersect -a cichlid_Niloticus_Burtoni_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/aCNEs/aCNE.sorted.bed -wb > On-Ab_substaCNEoverlap.txt # 1460/48920232
# bedtools intersect -a cichlid_Niloticus_Brichardi_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/CNEs/aCNEs/aCNE.sorted.bed -wb > On-Nb_substaCNEoverlap.txt # 2194/52733124

# Promoters
# First, pull out switch/non-switched promoters
grep -wiFf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Abvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed > Abvsany_switch.genes.5kb_promoters.stranded.sorted.bed
grep -wiFvf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Abvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed > Abvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed
grep -wiFf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Mzvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed > Mzvsany_switch.genes.5kb_promoters.stranded.sorted.bed
grep -wiFvf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Mzvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed > Mzvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed
grep -wiFf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Nbvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed > Nbvsany_switch.genes.5kb_promoters.stranded.sorted.bed
grep -wiFvf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Nbvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed > Nbvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed
grep -wiFf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Pnvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed > Pnvsany_switch.genes.5kb_promoters.stranded.sorted.bed
grep -wiFvf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Pnvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed > Pnvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed
grep -wiFf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Onvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.bed > Onvsany_switch.genes.5kb_promoters.stranded.sorted.bed
grep -wiFvf /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/Onvsany_switch.genes /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.bed > Onvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed

# switched
for i in mz11*_substitutions_1217.bed ; do bedtools intersect -a $i -b Mzvsany_switch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-swPromoverlap.txt ; done
for i in pn1*_substitutions_1217.bed ; do bedtools intersect -a $i -b Pnvsany_switch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-swPromoverlap.txt ; done
for i in ab1*_substitutions_1217.bed ; do bedtools intersect -a $i -b Abvsany_switch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-swPromoverlap.txt ; done
for i in nb1*_substitutions_1217.bed ; do bedtools intersect -a $i -b Nbvsany_switch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-swPromoverlap.txt ; done
for i in on11*_substitutions_1217.bed ; do bedtools intersect -a $i -b Onvsany_switch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-swPromoverlap.txt ; done
# bedtools intersect -a cichlid_Niloticus_Zebra_substitutions.bed -b Onvsany_switch.genes.5kb_promoters.stranded.sorted.bed -wb > On-Mz_subst-swPromoverlap.txt # 707400/47473955
# bedtools intersect -a cichlid_Niloticus_Nyererei_substitutions.bed -b Onvsany_switch.genes.5kb_promoters.stranded.sorted.bed -wb > On-Pn_subst-swPromoverlap.txt # 725544/49606035
# bedtools intersect -a cichlid_Niloticus_Burtoni_substitutions.bed -b Onvsany_switch.genes.5kb_promoters.stranded.sorted.bed -wb > On-Ab_subst-swPromoverlap.txt # 714706/48920232
# bedtools intersect -a cichlid_Niloticus_Brichardi_substitutions.bed -b Onvsany_switch.genes.5kb_promoters.stranded.sorted.bed -wb > On-Nb_subst-swPromoverlap.txt # 746708/52733124

# no switch
for i in mz11*_substitutions_1217.bed ; do bedtools intersect -a $i -b Mzvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-NoswPromoverlap.txt ; done
for i in pn1*_substitutions_1217.bed ; do bedtools intersect -a $i -b Pnvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-NoswPromoverlap.txt ; done
for i in ab1*_substitutions_1217.bed ; do bedtools intersect -a $i -b Abvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-NoswPromoverlap.txt ; done
for i in nb1*_substitutions_1217.bed ; do bedtools intersect -a $i -b Nbvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-NoswPromoverlap.txt ; done
for i in on11*_substitutions_1217.bed ; do bedtools intersect -a $i -b Onvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed -wb > ${i}_subst-NoswPromoverlap.txt ; done
# bedtools intersect -a cichlid_Niloticus_Zebra_substitutions.bed -b Onvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed -wb > On-Mz_subst-NoswPromoverlap.txt # 4032734/47473955
# bedtools intersect -a cichlid_Niloticus_Nyererei_substitutions.bed -b Onvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed -wb > On-Pn_subst-NoswPromoverlap.txt # 4251931/49606035
# bedtools intersect -a cichlid_Niloticus_Burtoni_substitutions.bed -b Onvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed -wb > On-Ab_subst-NoswPromoverlap.txt # 4216301/48920232
# bedtools intersect -a cichlid_Niloticus_Brichardi_substitutions.bed -b Onvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed -wb > On-Nb_subst-NoswPromoverlap.txt # 4512816/52733124

# 3'UTR
for i in mz11*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/3UTR/Metriaclima_zebra.BROADMZ1.UTR3_longest.bed -wb > ${i}_subst-3UTRoverlap.txt ; done
for i in pn1*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/3UTR/Pundamilia_nyererei.BROADPN1.UTR3_longest.bed -wb > ${i}_subst-3UTRoverlap.txt ; done
for i in ab1*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/3UTR/Astatotilapia_burtoni.BROADAB1.UTR3_longest.bed -wb > ${i}_subst-3UTRoverlap.txt ; done
for i in nb1*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/3UTR/Neolamprologus_brichardi.BROADNB1.UTR3_longest.bed -wb > ${i}_subst-3UTRoverlap.txt ; done
for i in on11*_substitutions_1217.bed ; do bedtools intersect -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/3UTR/Oreochromis_niloticus.BROADON1.UTR3_longest.bed -wb > ${i}_subst-3UTRoverlap.txt ; done
# bedtools intersect -a cichlid_Niloticus_Zebra_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/3UTR/Oreochromis_niloticus.BROADON1.UTR3_longest.bed -wb > On-Mz_subst-3UTRoverlap.txt # 1493200/47473955
# bedtools intersect -a cichlid_Niloticus_Nyererei_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/3UTR/Oreochromis_niloticus.BROADON1.UTR3_longest.bed -wb > On-Pn_subst-3UTRoverlap.txt # 1517366/49606035
# bedtools intersect -a cichlid_Niloticus_Burtoni_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/3UTR/Oreochromis_niloticus.BROADON1.UTR3_longest.bed -wb > On-Ab_subst-3UTRoverlap.txt # 1503244/48920232
# bedtools intersect -a cichlid_Niloticus_Brichardi_substitutions.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/3UTR/Oreochromis_niloticus.BROADON1.UTR3_longest.bed -wb > On-Nb_subst-3UTRoverlap.txt # 1660944/52733124

# Then run the above
sbatch subst_intersect.sh # {DONE}

# Intersect 3' UTR SNPs with miRNA binding sites

# First, we need to prepare a bed file of miRNA binding site locations (since the targetscan output is only the position within the UTR)
# pull out all within UTR cooridnates of bindign sites
cut -f5,6,36,39 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/miRNA_targetpred/targetscan/filtered/combined_families_targets.txt | sed 's/mz/mz.gene/g' | sed 's/pn/pn.gene/g' | sed 's/ab/ab.gene/g' | sed 's/nb/nb.gene/g' | sed 's/on/on.gene/g' | awk '{print $4,$1,$2,$3}' OFS='\t' > miRNA_bindingsite_inUTR.txt
# prepare BED files of UTRs
for i in *coords* ; do sed 's/\//\t/g' $i | sed 's/(-)/(XXX)/g'| sed 's/-/\t/g' | sed 's/(XXX)/\t-/g' | sed 's/(+)/\t+/g' | sed 's/_/\t/g' | awk '{print $1"_"$2,$3,$4,$5,$6}' OFS='\t' > "$(basename "$i" txt)bed" ; done
sed -i 's/mz/mz.gene/g' Metriaclima_zebra.BROADMZ1.UTR3_longest_coords.bed #21672
sed -i 's/pn/pn.gene/g' Pundamilia_nyererei.BROADPN1.UTR3_longest_coords.bed #20609
sed -i 's/ab/ab.gene/g' Astatotilapia_burtoni.BROADAB1.UTR3_longest_coords.bed #23436
sed -i 's/nb/nb.gene/g' Neolamprologus_brichardi.BROADNB1.UTR3_longest_coords.bed #20106
sed -i 's/on/on.gene/g' Oreochromis_niloticus.BROADON1.UTR3_longest_coords.bed #24559
sed -i 's/_/\t/g' Oreochromis_niloticus.BROADON1.UTR3_longest_coords.bed # as the file is different in terms of scaffold we need to change
# prepare files for matching
for i in *_coords.bed ; do awk '{print $5,$1,$2,$3,$4}' OFS='\t' $i > $i.2 ; done
# do the matching then output BED files of binding sites for each species
for i in *_coords.bed.2 ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' $i miRNA_bindingsite_inUTR.txt | grep -v 'REMOVEME' | awk '{if($8 == "-")print $5,$7-$3,$7-$2,$8,$1;else print $5,$6+$2,$6+$3,$8,$1;}' OFS='\t' | sort -u | sort -k1,1 -k2,2n > "$(basename "$i" .bed.2).BindingSites.bed"
done

# Overlap with all pairwise VCFs
ml bedtools/2.25.0
for i in mz11*_subst-3UTRoverlap.txt ; do bedtools intersect -a Metriaclima_zebra.BROADMZ1.UTR3_longest_coords.BindingSites.bed -b $i -wb > ${i}_subst_miRNA-3UTRoverlap.txt ; done
for i in pn1*_subst-3UTRoverlap.txt ; do bedtools intersect -a Pundamilia_nyererei.BROADPN1.UTR3_longest_coords.BindingSites.bed -b $i -wb > ${i}_subst_miRNA-3UTRoverlap.txt ; done
for i in ab1*_subst-3UTRoverlap.txt ; do bedtools intersect -a Astatotilapia_burtoni.BROADAB1.UTR3_longest_coords.BindingSites.bed -b $i -wb > ${i}_subst_miRNA-3UTRoverlap.txt ; done
for i in nb1*_subst-3UTRoverlap.txt ; do bedtools intersect -a Neolamprologus_brichardi.BROADNB1.UTR3_longest_coords.BindingSites.bed -b $i -wb > ${i}_subst_miRNA-3UTRoverlap.txt ; done
for i in on11*_subst-3UTRoverlap.txt ; do bedtools intersect -a Oreochromis_niloticus.BROADON1.UTR3_longest_coords.BindingSites.bed -b $i -wb > ${i}_subst_miRNA-3UTRoverlap.txt ; done
# bedtools intersect -a Oreochromis_niloticus.BROADON1.UTR3_longest_coords.BindingSites.bed -b On-Mz_subst-3UTRoverlap.txt -wb > On-Mz_subst_miRNA-3UTRoverlap.txt #147847
# bedtools intersect -a Oreochromis_niloticus.BROADON1.UTR3_longest_coords.BindingSites.bed -b On-Pn_subst-3UTRoverlap.txt -wb > On-Pn_subst_miRNA-3UTRoverlap.txt #149514
# bedtools intersect -a Oreochromis_niloticus.BROADON1.UTR3_longest_coords.BindingSites.bed -b On-Ab_subst-3UTRoverlap.txt -wb > On-Ab_subst_miRNA-3UTRoverlap.txt #147333
# bedtools intersect -a Oreochromis_niloticus.BROADON1.UTR3_longest_coords.BindingSites.bed -b On-Nb_subst-3UTRoverlap.txt -wb > On-Nb_subst_miRNA-3UTRoverlap.txt #163804

# Put all the numbers into tables to then plot (~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/SNP_overlap.R)

# Prepare files to output the counts which you can then paste directly into Excel (~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/SNP_overlap2.txt)
echo ab1 > mzlistcomp
echo danRer7 >> mzlistcomp
echo gasAcu1 >> mzlistcomp
echo nb1 >> mzlistcomp
echo on11 >> mzlistcomp
echo oryLat2 >> mzlistcomp
echo pn1 >> mzlistcomp

echo ab1 > pnlistcomp
echo danRer7 >> pnlistcomp
echo gasAcu1 >> pnlistcomp
echo nb1 >> pnlistcomp
echo on11 >> pnlistcomp
echo oryLat2 >> pnlistcomp
echo mz11 >> pnlistcomp

echo danRer7 >> ablistcomp
echo gasAcu1 >> ablistcomp
echo nb1 >> ablistcomp
echo on11 >> ablistcomp
echo oryLat2 >> ablistcomp
echo pn1 >> ablistcomp
echo mz11 >> ablistcomp

echo ab1 > nblistcomp
echo danRer7 >> nblistcomp
echo gasAcu1 >> nblistcomp
echo on11 >> nblistcomp
echo oryLat2 >> nblistcomp
echo pn1 >> nblistcomp
echo mz11 >> nblistcomp

echo ab1 > onlistcomp
echo danRer7 >> onlistcomp
echo gasAcu1 >> onlistcomp
echo nb1 >> onlistcomp
echo oryLat2 >> onlistcomp
echo pn1 >> onlistcomp
echo mz11 >> onlistcomp


while read F ; do
	for i in mz11_${F}*_subst*.txt ; do wc -l $i >> mz11_${F}.counts
	for j in {1..2} ; do cut -d' ' -f"$j" mz11_${F}.counts | paste -s | sed "s/mz11_${F}_substitutions_1217.bed_//g" | sed 's/.txt//g' >> mz11_${F}.counts2
	tail -2 mz11_${F}.counts2 > mz11_${F}.counts3
done
done
done < mzlistcomp
while read F; do
	awk '$7=FILENAME' OFS='\t' mz11_${F}.counts3 |
	awk '{print $7,$1,$2,$3,$4,$5,$6}' OFS='\t' |
	tac |
	sed '0,/mz11_ab1.counts3/{s/mz11_ab1.counts3/RefSpecies_CompSpecies/}' > mz11_${F}.counts4
	head -1 mz11_${F}.counts4 | grep -v mz11 >> colheads
	tail -1 mz11_${F}.counts4 > mz11_${F}.counts5
	cat colheads mz11*.counts5 | sed 's/.counts3//g' > mz11.counts6
done < mzlistcomp

while read F ; do
	for i in pn1_${F}*_subst*.txt ; do wc -l $i >> pn1_${F}.counts
	for j in {1..2} ; do cut -d' ' -f"$j" pn1_${F}.counts | paste -s | sed "s/pn1_${F}_substitutions_1217.bed_//g" | sed 's/.txt//g' >> pn1_${F}.counts2
	tail -2 pn1_${F}.counts2 > pn1_${F}.counts3
done
done
done < pnlistcomp
while read F; do
	awk '$7=FILENAME' OFS='\t' pn1_${F}.counts3 |
	awk '{print $7,$1,$2,$3,$4,$5,$6}' OFS='\t' |
	tac |
	sed '0,/pn1_ab1.counts3/{s/pn1_ab1.counts3/RefSpecies_CompSpecies/}' > pn1_${F}.counts4
	head -1 pn1_${F}.counts4 | grep -v pn1 >> pn1colheads
	tail -1 pn1_${F}.counts4 > pn1_${F}.counts5
	cat colheads pn1*.counts5 | sed 's/.counts3//g' > pn1.counts6
done < pnlistcomp


while read F ; do
	for i in ab1_${F}*_subst*.txt ; do wc -l $i >> ab1_${F}.counts
	for j in {1..2} ; do cut -d' ' -f"$j" ab1_${F}.counts | paste -s | sed "s/ab1_${F}_substitutions_1217.bed_//g" | sed 's/.txt//g' >> ab1_${F}.counts2
	tail -2 ab1_${F}.counts2 > ab1_${F}.counts3
done
done
done < ablistcomp
while read F; do
	awk '$7=FILENAME' OFS='\t' ab1_${F}.counts3 |
	awk '{print $7,$1,$2,$3,$4,$5,$6}' OFS='\t' |
	tac |
	sed '0,/ab1_danRer7.counts3/{s/ab1_danRer7.counts3/RefSpecies_CompSpecies/}' > ab1_${F}.counts4
	head -1 ab1_${F}.counts4 | grep -v ab1 >> ab1colheads
	tail -1 ab1_${F}.counts4 > ab1_${F}.counts5
	cat colheads ab1*.counts5 | sed 's/.counts3//g' > ab1.counts6
done < ablistcomp

while read F ; do
	for i in nb1_${F}*_subst*.txt ; do wc -l $i >> nb1_${F}.counts
	for j in {1..2} ; do cut -d' ' -f"$j" nb1_${F}.counts | paste -s | sed "s/nb1_${F}_substitutions_1217.bed_//g" | sed 's/.txt//g' >> nb1_${F}.counts2
	tail -2 nb1_${F}.counts2 > nb1_${F}.counts3
done
done
done < nblistcomp
while read F; do
	awk '$7=FILENAME' OFS='\t' nb1_${F}.counts3 |
	awk '{print $7,$1,$2,$3,$4,$5,$6}' OFS='\t' |
	tac |
	sed '0,/nb1_ab1.counts3/{s/nb1_ab1.counts3/RefSpecies_CompSpecies/}' > nb1_${F}.counts4
	head -1 nb1_${F}.counts4 | grep -v nb1 >> nb1colheads
	tail -1 nb1_${F}.counts4 > nb1_${F}.counts5
	cat colheads nb1*.counts5 | sed 's/.counts3//g' > nb1.counts6
done < nblistcomp

while read F ; do
	for i in on11_${F}*_subst*.txt ; do wc -l $i >> on11_${F}.counts
	for j in {1..2} ; do cut -d' ' -f"$j" on11_${F}.counts | paste -s | sed "s/on11_${F}_substitutions_1217.bed_//g" | sed 's/.txt//g' >> on11_${F}.counts2
	tail -2 on11_${F}.counts2 > on11_${F}.counts3
done
done
done < onlistcomp
while read F; do
	awk '$7=FILENAME' OFS='\t' on11_${F}.counts3 |
	awk '{print $7,$1,$2,$3,$4,$5,$6}' OFS='\t' |
	tac |
	sed '0,/on11_ab1.counts3/{s/on11_ab1.counts3/RefSpecies_CompSpecies/}' > on11_${F}.counts4
	head -1 on11_${F}.counts4 | grep -v on11 >> on11colheads
	tail -1 on11_${F}.counts4 > on11_${F}.counts5
	cat colheads on11*.counts5 | sed 's/.counts3//g' > on11.counts6
done < onlistcomp

rm *.counts
rm *.counts2
rm *.counts3
rm *.counts4
rm *.counts5

# get total SNPs
for i in mz11 pn1 ab1 nb1 on11 ; do
	for j in ${i}_*.bed ; do
	wc -l $j | awk -F' ' '{print $1}'
done
done

for i in pn1 ab1 nb1 on11 ; do
	for j in ${i}_*.bed ; do
	wc -l $j | awk -F' ' '{print $1}'
done
done


##########################################################################################################################

########## 2. Test if any pairwise SNPs overlap TFBS motif scans in promoter regions

WDpromtfbs=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2) #out to here
cd $WDpromtfbs

# a. Prep species specific, reordered, trimmed files based on final TF-TG edge table (these are FDR<0.05 or extrapolated, conf >70% Hs/Mm conservation)
awk '{print $6,$7,$8,$9,$3,$4,$1,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' motifenr_merged-TFBSs_map2d.txt > motifenr_merged-TFBSs_map2d.simplified.txt
for i in mz pn ab nb on ; do
  grep $i.gene motifenr_merged-TFBSs_map2d.simplified.txt > $i-motifenr_mergedmap2d.txt
done
rm motifenr_merged-TFBSs_map2d.simplified.txt

# b. Match gene promoter with bed file location of promoter and then calculate the whole genome region of the TFBS based on fimo hit (and strand!)
mkdir SNP_overlap
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap

for i in ../*-motifenr_mergedmap2d.txt ; do
  ln -s $i
done

for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_14032018.5kb_promoters.stranded.bed ; do ln -s $i ; done
for i in *.bed; do
 awk '{print $5,$1,$2,$3,$4}' OFS='\t' $i > "$(basename "$i" .bed)_reorder.bed"
done

ls -1 *-motifenr_mergedmap2d.txt > FIMOres
ls -1 *stranded_reorder.bed > PROMbed
while read -u 3 -r file1 && read -u 4 -r file2
do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' ${file1} ${file2} |
 grep -v 'REMOVEME' |
 awk '{if($21 == "+")print $18,$19+$8,$19+$9,$10,$1,$2,$3,$4,$21,$7,$5,$6,$11,$12,$13,$14,$15,$16;else print $18,$20-$9,$20-$8,$10,$1,$2,$3,$4,$21,$7,$5,$6,$11,$12,$13,$14,$15,$16;}' OFS='\t' > "$(basename "${file2}" .txt)_motifgenomeannot.txt"
done 3<PROMbed 4<FIMOres

# # Checking that the correct region is taken for -ve strands!
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' Mz_14032018.5kb_promoters.stranded_reorder.bed mz-motifenr_mergedmap2d.txt | grep -v 'REMOVEME' |
# awk '{if($21 == "+")print $18,$19+$8,$19+$9,$10,$1,$2,$3,$4,$21,$7,$5,$6,$11,$12,$13,$14,$15,$16;else print $18,$20-$9,$20-$8,$10,$1,$2,$3,$4,$21,$7,$5,$6,$11,$12,$13,$14,$15,$16;}' OFS='\t' | head

# Cols of above:- 1. motif_chr; 2. motif_start; 3. motif_end; 4. motif_strand; 5. gene; 6. motif_genesymbolDr; 7. motif_genesymbolGa; 8. motif_genesymbolSp; 9. promoter_strand; 7. motif_ID; 8. motif_cichlidID; 9. motif_gene_symbol; 10. fimo_score; 11. fimo_pval; 12. fimo_qval; 13. motif_seq; 14. conf_level; 15. conf_score

# c. Intersect the complete TFBS region with all pairwise substitution VCFs (except the teleost reference files)
ml bedtools/2.25.0
ml zlib
for i in /tgac/workarea/group-vh/cichlids/substitutions_VCFs/*_substitutions_1217.bed ; do ln -s $i ; done
ls -1 *_substitutions_1217.bed | grep -v '^danRer7' | grep -v '^gasAcu1' | grep -v '^oryLat2' > substvcf_files
printf 'ab-motifenr_mergedmap2d_motifgenomeannot.txt\n%.0s' {1..7} > motifGenomeannot_files
printf 'mz-motifenr_mergedmap2d_motifgenomeannot.txt\n%.0s' {1..7} >> motifGenomeannot_files
printf 'nb-motifenr_mergedmap2d_motifgenomeannot.txt\n%.0s' {1..7} >> motifGenomeannot_files
printf 'on-motifenr_mergedmap2d_motifgenomeannot.txt\n%.0s' {1..7} >> motifGenomeannot_files
printf 'pn-motifenr_mergedmap2d_motifgenomeannot.txt\n%.0s' {1..7} >> motifGenomeannot_files
while read -u 3 -r file1 && read -u 4 -r file2
do
  bedtools intersect -a ${file1} -b ${file2} -wb > ${file1}_${file2}
done 3<substvcf_files 4<motifGenomeannot_files

rm *_1217.vcf
rm slurm*

##########################################################################################################################

# For Will, do the following

# a. Intersect the Mz motifs in promoter regions with the Lake Malawi VCF
# b. Intersect the Mz promoters with the Lake Malawi VCF

cd /tgac/workarea/Research-Groups/RG-cichlids/

nano mzPromandTFBS-MalawiVCFoverlap.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

ml bedtools/2.25.0
ml zlib

# this does the promoter TFBS overlap with VCF
bedtools intersect -a mz-motifenr_mergedmap2d_motifpromANDgenomeannot.txt -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> mz-motifenr_mergedmap2d_motifpromANDgenomeannot_MalawiVCFoverlap.bed
# this does the promoter overlap with the VCF
bedtools intersect -a /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_14032018.5kb_sorted.promoters.stranded.bed -b /tgac/workarea/Research-Groups/RG-cichlids/Malinsky_et_al_2017_filteredSNPs.vcf -wa -wb >> Mz_14032018.5kb_sorted.promoters.stranded_MalawiVCFoverlap.bed

# run the above
qsub -q Test -l select=1:mem=250GB:ncpus=1 mzPromandTFBS-MalawiVCFoverlap.sh # {DONE}

##########################################################################################################################

# d. Characterise how many of the variant-TFBSs intersects are from switched or non-switched promoters

# total overlaps
for i in mz11 pn1 ab1 nb1 on11 ; do
 for j in ${i}*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 wc -l $j | awk -F' ' '{print $1}'
done
done

for i in /tgac/workarea/group-vh/cichlids/substitutions_VCFs/*switch* ; do ln -s $i ; done

nano SNPTFBSoverlap.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 0-00:45 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

echo -n > SNPs_overlap_TFBSs_switch-noswitch.counts.txt

# switched
for i in mz11*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tState-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Mzvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done

for i in pn1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tState-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Pnvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done

for i in ab1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tState-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Abvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done

for i in nb1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tState-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Nbvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done

for i in on11*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tState-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Onvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done


# no switch
for i in mz11*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tNon_State-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Mzvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done

for i in pn1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tNon_State-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Pnvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done

for i in ab1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tNon_State-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Abvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done

for i in nb1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tNon_State-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Nbvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done

for i in on11*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
	printf $i >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	printf '\tNon_State-changed_Promoters_TFBSs\t' >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"REMOVEME";}}' Onvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v 'REMOVEME' | wc -l >> SNPs_overlap_TFBSs_switch-noswitch.counts.txt
done


# run the above
sbatch SNPTFBSoverlap.sh




# d - i. TFBS mutations enrichment in promoters

### Reviewer may ask to:
## 1.	Categorize TFBS mutations as breaking/gaining
# separate SNPs overlapping TFBSs into
# a. gain - whether the motif (that the SNP overlaps, say in Mz) is absent in an ancestral species (like Ab), or
# b. loss - whether the motif (that the SNP overlaps, say in Nb) is absent in a descendent species (like Mz).
# then rank each of the above according to the motif's fimo P-value - SNPs that change a PWM are marked as having 'regulatory potential' if P-value < 3 × 10−6 (similar to PMC5210548)
## 2.	Then, run an enrichment to test whether those that break/gain associate with switched/non-switched gene promoters


## 1.	Categorize TFBS mutations as breaking/gaining
# prepping a python script for this motifgainbreak.py (currently incomplete)


## 2.	Then, run an enrichment to test whether those that break/gain associate with switched/non-switched gene promoters

## if you run GAT then A set of intervals S with segments of interest to test should be annotation of all SNPs overlapping TFBSs (not all SNPs!)


######### Another (better) option is to:
# 1. replace all SNPs overlapping motifs in the promoter with the ALT genotype, then
# 2. re-scan for motifs (this tests whether the SNP breaks the site as flanking is unchanged) and compare the p-value for the original and mutated hit, then
# 3. Categorise according to motif gain (from ancestor REF SNP), loss (from ancestor from REF SNP)


# strategy
# first file are your fasta files with >gene IDs followed by promoter sequence
# second file is a tabulated index file with all SNP details
# then:
	# a. read fasta and index file: match ID
	# b. read along fasta seq and stop when get to index col1 position in promoter
	# c. read col2 and check that position nt in promoter seq matches REF
	# d. read col3 and change that position nt in promoter seq to ALT
		# c and d. if fasta_pos=index_col2, fasta_pos=index_col3


# # Files to use
# /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta
# /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta

# /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/*_*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt



## Test files:
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed

# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed

# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/mz11_pn1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/pn1_mz11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt

cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/

cat -n mz11_pn1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt | gshuf -n 500000 > mz11_pn1_subst_SAMPLE.txt # randomly sample 500,000 lines
head -100 Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta > Mz_Prom_SAMPLE.fasta
head -100 Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Mz_PromBED.txt

# prepped a python script for this /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/motifgainbreak.py
# run the above on the cluster

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/
mkdir Motif_gain_and_loss/
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/

## for the script to run, you need to number each line in the VCF SNP file
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/
for i in *_*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt; do
	cat -n $i | sed 's/ //g' > "$(basename "$i" .txt)_Num.txt"
done

## split the numbered SNP vcf file into 5,000 chunks and then run as an array so that it runs quicker
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/
for f in /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot_Num.txt; do
	split -d -a4 -l5000 --additional-suffix=.txt "$f" "${f%.txt}-"
done

## run as an array by mapping the commands to five species specific files (as there is a limit of 6000 on the number of submission jobs)
nano makearraycommands.sh

#!/bin/sh

echo -n > motifgainbreak_commands2_Ab # create an empty file to append to
echo -n > motifgainbreak_commands2_Mz # create an empty file to append to
echo -n > motifgainbreak_commands2_Pn # create an empty file to append to
echo -n > motifgainbreak_commands2_Nb # create an empty file to append to
echo -n > motifgainbreak_commands2_On # create an empty file to append to

# Ab
for i in {0000..0508} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_mz11_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Ab
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed AbProm_mzSNP-'${i} >> motifgainbreak_commands2_Ab
	printf '.fasta\n' >> motifgainbreak_commands2_Ab
done
for i in {0000..0482} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_pn1_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Ab
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed AbProm_pnSNP-'${i} >> motifgainbreak_commands2_Ab
	printf '.fasta\n' >> motifgainbreak_commands2_Ab
done
for i in {0000..0809} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_nb1_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Ab
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed AbProm_nbSNP-'${i} >> motifgainbreak_commands2_Ab
	printf '.fasta\n' >> motifgainbreak_commands2_Ab
done
for i in {0000..1275} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_on11_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Ab
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed AbProm_onSNP-'${i} >> motifgainbreak_commands2_Ab
	printf '.fasta\n' >> motifgainbreak_commands2_Ab
done
# Mz
for i in {0000..0379} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_pn1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Mz
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed MzProm_pnSNP-'${i} >> motifgainbreak_commands2_Mz
	printf '.fasta\n' >> motifgainbreak_commands2_Mz
done
for i in {0000..0458} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_ab1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Mz
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed MzProm_abSNP-'${i} >> motifgainbreak_commands2_Mz
	printf '.fasta\n' >> motifgainbreak_commands2_Mz
done
for i in {0000..0764} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_nb1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Mz
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed MzProm_nbSNP-'${i} >> motifgainbreak_commands2_Mz
	printf '.fasta\n' >> motifgainbreak_commands2_Mz
done
for i in {0000..1212} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_on11_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Mz
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed MzProm_onSNP-'${i} >> motifgainbreak_commands2_Mz
	printf '.fasta\n' >> motifgainbreak_commands2_Mz
done
# Pn
for i in {0000..0415} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_mz11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Pn
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed PnProm_mzSNP-'${i} >> motifgainbreak_commands2_Pn
	printf '.fasta\n' >> motifgainbreak_commands2_Pn
done
for i in {0000..0461} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_ab1_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Pn
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed PnProm_abSNP-'${i} >> motifgainbreak_commands2_Pn
	printf '.fasta\n' >> motifgainbreak_commands2_Pn
done
for i in {0000..0760} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_nb1_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Pn
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed PnProm_nbSNP-'${i} >> motifgainbreak_commands2_Pn
	printf '.fasta\n' >> motifgainbreak_commands2_Pn
done
for i in {0000..1191} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_on11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Pn
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed PnProm_onSNP-'${i} >> motifgainbreak_commands2_Pn
	printf '.fasta\n' >> motifgainbreak_commands2_Pn
done
# Nb
for i in {0000..0557} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_mz11_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Nb
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed NbProm_mzSNP-'${i} >> motifgainbreak_commands2_Nb
	printf '.fasta\n' >> motifgainbreak_commands2_Nb
done
for i in {0000..0531} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_pn1_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Nb
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed NbProm_pnSNP-'${i} >> motifgainbreak_commands2_Nb
	printf '.fasta\n' >> motifgainbreak_commands2_Nb
done
for i in {0000..0534} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_ab1_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Nb
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed NbProm_abSNP-'${i} >> motifgainbreak_commands2_Nb
	printf '.fasta\n' >> motifgainbreak_commands2_Nb
done
for i in {0000..0799} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_on11_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_Nb
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed NbProm_onSNP-'${i} >> motifgainbreak_commands2_Nb
	printf '.fasta\n' >> motifgainbreak_commands2_Nb
done
# On
for i in {0000..1313} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_mz11_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_On
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed OnProm_mzSNP-'${i} >> motifgainbreak_commands2_On
	printf '.fasta\n' >> motifgainbreak_commands2_On
done
for i in {0000..1239} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_pn1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_On
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed OnProm_pnSNP-'${i} >> motifgainbreak_commands2_On
	printf '.fasta\n' >> motifgainbreak_commands2_On
done
for i in {0000..1239} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_ab1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_On
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed OnProm_abSNP-'${i} >> motifgainbreak_commands2_On
	printf '.fasta\n' >> motifgainbreak_commands2_On
done
for i in {0000..1207} ; do
	printf 'python motifgainbreak.py /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_nb1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot_Num-'${i} >> motifgainbreak_commands2_On
	printf '.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed OnProm_nbSNP-'${i} >> motifgainbreak_commands2_On
	printf '.fasta\n' >> motifgainbreak_commands2_On
done

# run the above
sh makearraycommands.sh

wc -l motifgainbreak_commands2* # detemine number of array items
# 3078 motifgainbreak_commands2_Ab
# 2831 motifgainbreak_commands2_Pn
# 2817 motifgainbreak_commands2_Mz
# 2425 motifgainbreak_commands2_Nb
# 5002 motifgainbreak_commands2_On
# 16153 total

nano motifgainbreak_array2_Ab.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-3077
#SBATCH --mem-per-cpu 18000
#SBATCH -t 0-0:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t motifgainbreak_commands2_Ab < motifgainbreak_commands2_Ab # assign as elements to $motifgainbreak_commands variable
${motifgainbreak_commands2_Ab[${SLURM_ARRAY_TASK_ID}]}

nano motifgainbreak_array2_Pn.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-2830
#SBATCH --mem-per-cpu 18000
#SBATCH -t 0-0:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t motifgainbreak_commands2_Pn < motifgainbreak_commands2_Pn # assign as elements to $motifgainbreak_commands variable
${motifgainbreak_commands2_Pn[${SLURM_ARRAY_TASK_ID}]}

nano motifgainbreak_array2_Mz.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-2816
#SBATCH --mem-per-cpu 18000
#SBATCH -t 0-0:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t motifgainbreak_commands2_Mz < motifgainbreak_commands2_Mz # assign as elements to $motifgainbreak_commands variable
${motifgainbreak_commands2_Mz[${SLURM_ARRAY_TASK_ID}]}

nano motifgainbreak_array2_Nb.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-2424
#SBATCH --mem-per-cpu 24000
#SBATCH -t 0-0:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t motifgainbreak_commands2_Nb < motifgainbreak_commands2_Nb # assign as elements to $motifgainbreak_commands variable
${motifgainbreak_commands2_Nb[${SLURM_ARRAY_TASK_ID}]}

nano motifgainbreak_array2_On.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-5001
#SBATCH --mem-per-cpu 24000
#SBATCH -t 0-0:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t motifgainbreak_commands2_On < motifgainbreak_commands2_On # assign as elements to $motifgainbreak_commands variable
${motifgainbreak_commands2_On[${SLURM_ARRAY_TASK_ID}]}


# run the above ensuring not to exceed 6000 jobs
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/
sbatch motifgainbreak_array2_Ab.sh # DONE - 3078 jobs - ALL COMPLETED
sbatch motifgainbreak_array2_Pn.sh # DONE - 2831 jobs - ALL COMPLETED
sbatch motifgainbreak_array2_Mz.sh # DONE - 2817 jobs - ALL COMPLETED
sbatch motifgainbreak_array2_Nb.sh # DONE - 2425 jobs - ALL COMPLETED
sbatch motifgainbreak_array2_On.sh # DONE - 5002 jobs - ALL COMPLETED

# Once completed, remove all /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/*_*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot_Num-*.txt
rm /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot_Num-*.txt # DONE
rm /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot_Num-*.txt # DONE
rm /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot_Num-*.txt # DONE
rm /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot_Num-*.txt # DONE
rm /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot_Num-*.txt # DONE
rm /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/slurm.t* # DONE
##### NOTE: All Num files have been deleted and hence to correlate to the SNP in the SNP file, count to the nth line of the SNPID in the fasta e.g. _SNPid:5 means the fifth SNP line

## then join the files
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/

# assign in files to a file to map as elements in a variable
echo 'AbProm_mzSNP-*.fasta' > join_commands
echo 'AbProm_pnSNP-*.fasta' >> join_commands
echo 'AbProm_nbSNP-*.fasta' >> join_commands
echo 'AbProm_onSNP-*.fasta' >> join_commands
echo 'MzProm_pnSNP-*.fasta' >> join_commands
echo 'MzProm_abSNP-*.fasta' >> join_commands
echo 'MzProm_nbSNP-*.fasta' >> join_commands
echo 'MzProm_onSNP-*.fasta' >> join_commands
echo 'PnProm_mzSNP-*.fasta' >> join_commands
echo 'PnProm_abSNP-*.fasta' >> join_commands
echo 'PnProm_nbSNP-*.fasta' >> join_commands
echo 'PnProm_onSNP-*.fasta' >> join_commands
echo 'NbProm_mzSNP-*.fasta' >> join_commands
echo 'NbProm_pnSNP-*.fasta' >> join_commands
echo 'NbProm_abSNP-*.fasta' >> join_commands
echo 'NbProm_onSNP-*.fasta' >> join_commands
echo 'OnProm_mzSNP-*.fasta' >> join_commands
echo 'OnProm_pnSNP-*.fasta' >> join_commands
echo 'OnProm_abSNP-*.fasta' >> join_commands
echo 'OnProm_nbSNP-*.fasta' >> join_commands

# assign out files to a file to map as elements in a variable
echo 'AbProm_mzSNP.5kbProm.fasta' > join_commands2
echo 'AbProm_pnSNP.5kbProm.fasta' >> join_commands2
echo 'AbProm_nbSNP.5kbProm.fasta' >> join_commands2
echo 'AbProm_onSNP.5kbProm.fasta' >> join_commands2
echo 'MzProm_pnSNP.5kbProm.fasta' >> join_commands2
echo 'MzProm_abSNP.5kbProm.fasta' >> join_commands2
echo 'MzProm_nbSNP.5kbProm.fasta' >> join_commands2
echo 'MzProm_onSNP.5kbProm.fasta' >> join_commands2
echo 'PnProm_mzSNP.5kbProm.fasta' >> join_commands2
echo 'PnProm_abSNP.5kbProm.fasta' >> join_commands2
echo 'PnProm_nbSNP.5kbProm.fasta' >> join_commands2
echo 'PnProm_onSNP.5kbProm.fasta' >> join_commands2
echo 'NbProm_mzSNP.5kbProm.fasta' >> join_commands2
echo 'NbProm_pnSNP.5kbProm.fasta' >> join_commands2
echo 'NbProm_abSNP.5kbProm.fasta' >> join_commands2
echo 'NbProm_onSNP.5kbProm.fasta' >> join_commands2
echo 'OnProm_mzSNP.5kbProm.fasta' >> join_commands2
echo 'OnProm_pnSNP.5kbProm.fasta' >> join_commands2
echo 'OnProm_abSNP.5kbProm.fasta' >> join_commands2
echo 'OnProm_nbSNP.5kbProm.fasta' >> join_commands2

nano joinfasta.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-19
#SBATCH --mem-per-cpu 78000
#SBATCH -t 0-1:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t join_commands < join_commands # assign as elements to variable $join_commands
mapfile -t join_commands2 < join_commands2 # assign as elements to variable $join_commands2

cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}
cat ${join_commands[${SLURM_ARRAY_TASK_ID}]} > ${join_commands2[${SLURM_ARRAY_TASK_ID}]}

# run the above
sbatch joinfasta.sh

### After files have been joined, remove all split fasta files and scan all joined fasta sequences for motifs
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/
rm slurm.*
rm *Prom_*SNP-*.fasta

## Create 0-order backgrounds of the mutant promoter sequences for fimo runs
# Providing an appropriate background model is the single most important step you can take to improve your FIMO results.
# The FIMO match score is the log-odds ratio of the probability of observing a sequence given the motif PWM to the probability of observing the same sequence given the background model.
# Ideally you'd create a background model from sequence data that is biologically similar to the sequences you are analyzing with FIMO, but that doesn't contain any instances of the motifs you are trying to match. If you are scanning a full genome then using the nucleotide frequencies for that genome as a background model is a reasonable approach.

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/

# a. create array commands for creating 0-order nbackgrounds

echo 'fasta-get-markov -dna AbProm_mzSNP.5kbProm.fasta AbProm_mzSNP.5kbProm.0orderMarkovBackgrnd' > getmarkovBG_commands
echo 'fasta-get-markov -dna AbProm_pnSNP.5kbProm.fasta AbProm_pnSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna AbProm_nbSNP.5kbProm.fasta AbProm_nbSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna AbProm_onSNP.5kbProm.fasta AbProm_onSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna MzProm_pnSNP.5kbProm.fasta MzProm_pnSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna MzProm_abSNP.5kbProm.fasta MzProm_abSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna MzProm_nbSNP.5kbProm.fasta MzProm_nbSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna MzProm_onSNP.5kbProm.fasta MzProm_onSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna PnProm_mzSNP.5kbProm.fasta PnProm_mzSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna PnProm_abSNP.5kbProm.fasta PnProm_abSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna PnProm_nbSNP.5kbProm.fasta PnProm_nbSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna PnProm_onSNP.5kbProm.fasta PnProm_onSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna NbProm_mzSNP.5kbProm.fasta NbProm_mzSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna NbProm_pnSNP.5kbProm.fasta NbProm_pnSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna NbProm_abSNP.5kbProm.fasta NbProm_abSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna NbProm_onSNP.5kbProm.fasta NbProm_onSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna OnProm_mzSNP.5kbProm.fasta OnProm_mzSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna OnProm_pnSNP.5kbProm.fasta OnProm_pnSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna OnProm_abSNP.5kbProm.fasta OnProm_abSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands
echo 'fasta-get-markov -dna OnProm_nbSNP.5kbProm.fasta OnProm_nbSNP.5kbProm.0orderMarkovBackgrnd' >> getmarkovBG_commands

nano getmarkovBG.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-19
#SBATCH --mem-per-cpu 68000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ml meme/4.11.1 # note: the newest version is meme-5.0.1 but we want to replicate method for motif calling on native promoters

mapfile -t getmarkovBG_commands < getmarkovBG_commands # assign as elements to $motifgainbreak_commands variable
${getmarkovBG_commands[${SLURM_ARRAY_TASK_ID}]}

# run above
sbatch getmarkovBG.sh

rm slurm.*

# b. then run fimo to scan
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/

# motif IDs, orthogrouping and cichlid gene IDs here:
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/merged_mat_qual_pvals_ALL.out1.TF.* # e.g. merged_mat_qual_pvals_ALL.out1.TF.ab - These are the
# MEME files of all motifs for motif scanning according to confidence levels: 2a - cichlid-specific matrices; 2b - non-species specific matrices; 2c - JASPAR matrices
#~ Derived from Human - * is lowercase short e.g. mz pn ab nb on
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_*.meme
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_*.meme
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_*.meme
# all meta data on above motifs, including motifID, genesymbol, EnsemblID, orthoGroup, cichlid_geneIDs, collated Optimal p-values for scanning: /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/cichlid_allHs_motifs_meta_info.txt
#~ Derived from Mouse - * is lowercase short e.g. mz pn ab nb on
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_*.meme
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_*.meme
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_*.meme
# all meta data on above motifs, including motifID, genesymbol, EnsemblID, orthoGroup, cichlid_geneIDs, collated Optimal p-values for scanning: /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/cichlid_allMm_motifs_meta_info.txt



# to run fimo much quicker, split each meme file into separate motifs with a python script: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/split_meme.py (copied to /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/split_meme.py)
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/
mkdir -p /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/{CS,JP,CW}/{mz,pn,ab,nb,on}
mkdir -p /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/{CS,JP,CW}/{mz,pn,ab,nb,on}
# to run: python3 split_meme.py -i input_file -o output_dir -t file_type -s species

for h in HumanDerived MouseDerived ; do
	for i in mz pn ab nb on ; do
		python3 split_meme.py -i ${h}/2a_CS_${i}.meme -o /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/${h}/CS/${i} -t CS -s ${i}
	done
done

for h in HumanDerived MouseDerived ; do
	for i in mz pn ab nb on ; do
		python3 split_meme.py -i ${h}/2b_CW_${i}.meme -o /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/${h}/CW/${i} -t CW -s ${i}
	done
done

for h in HumanDerived MouseDerived ; do
	for i in mz pn ab nb on ; do
		python3 split_meme.py -i ${h}/2c_JASPAR_${i}.meme -o /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/${h}/JP/${i} -t JP -s ${i}
	done
done

# remove random additional files created
rm HumanDerived/{CS,JP,CW}/{mz,pn,ab,nb,on}/MEME*.meme
rm MouseDerived/{CS,JP,CW}/{mz,pn,ab,nb,on}/MEME*.meme


cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/
mkdir -p /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/{HumanDerived,MouseDerived}/{CS_out,CW_out,JP_out}/{mz,pn,ab,nb,on}

nano createfimocommands.sh

#!/bin/sh

# Define the arrays
Ab=("mz" "pn" "nb" "on")
Mz=("pn" "ab" "nb" "on")
Pn=("mz" "ab" "nb" "on")
Nb=("mz" "pn" "ab" "on")
On=("mz" "pn" "ab" "nb")

length=${#Ab[@]} # get the length of the array - applicable to all

### NOTE: to keep output files and size down, immediately delete *.xml, *.gff and *.html

# the following looped commands will create all array commands for running species and SNP specific fimo on gene promoters
for h in HumanDerived MouseDerived ; do
	for g in CS CW JP ; do # NOTE: will remove MouseDerived JASPAR commands as these are same as HumanDerived - only need to run one
		for i in /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/${h}/${g}/ab/*.meme ; do
			for ((j=0;j<=$length-1;j++)); do
				echo fimo --bgfile AbProm_${Ab[$j]}SNP.5kbProm.0orderMarkovBackgrnd --o ${h}/${g}_out/ab/${i##*/}_${Ab[$j]}SNP --thresh 1e-4 --max-stored-scores 100000 ${i} AbProm_${Ab[$j]}SNP.5kbProm.fasta >> fimocommands
			done
		done
	done
done

for h in HumanDerived MouseDerived ; do
	for g in CS CW JP ; do # NOTE: will remove MouseDerived JASPAR commands as these are same as HumanDerived - only need to run one
		for i in /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/${h}/${g}/mz/*.meme ; do
			for ((j=0;j<=$length-1;j++)); do
				echo fimo --bgfile MzProm_${Mz[$j]}SNP.5kbProm.0orderMarkovBackgrnd --o ${h}/${g}_out/mz/${i##*/}_${Mz[$j]}SNP --thresh 1e-4 --max-stored-scores 100000 ${i} MzProm_${Mz[$j]}SNP.5kbProm.fasta >> fimocommands
			done
		done
	done
done

for h in HumanDerived MouseDerived ; do
	for g in CS CW JP ; do # NOTE: will remove MouseDerived JASPAR commands as these are same as HumanDerived - only need to run one
		for i in /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/${h}/${g}/pn/*.meme ; do
			for ((j=0;j<=$length-1;j++)); do
				echo fimo --bgfile PnProm_${Pn[$j]}SNP.5kbProm.0orderMarkovBackgrnd --o ${h}/${g}_out/pn/${i##*/}_${Pn[$j]}SNP --thresh 1e-4 --max-stored-scores 100000 ${i} PnProm_${Pn[$j]}SNP.5kbProm.fasta >> fimocommands
			done
		done
	done
done

for h in HumanDerived MouseDerived ; do
	for g in CS CW JP ; do # NOTE: will remove MouseDerived JASPAR commands as these are same as HumanDerived - only need to run one
		for i in /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/${h}/${g}/nb/*.meme ; do
			for ((j=0;j<=$length-1;j++)); do
				echo fimo --bgfile NbProm_${Nb[$j]}SNP.5kbProm.0orderMarkovBackgrnd --o ${h}/${g}_out/nb/${i##*/}_${Nb[$j]}SNP --thresh 1e-4 --max-stored-scores 100000 ${i} NbProm_${Nb[$j]}SNP.5kbProm.fasta >> fimocommands
			done
		done
	done
done

for h in HumanDerived MouseDerived ; do
	for g in CS CW JP ; do # NOTE: will remove MouseDerived JASPAR commands as these are same as HumanDerived - only need to run one
		for i in /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/${h}/${g}/on/*.meme ; do
			for ((j=0;j<=$length-1;j++)); do
				echo fimo --bgfile OnProm_${On[$j]}SNP.5kbProm.0orderMarkovBackgrnd --o ${h}/${g}_out/on/${i##*/}_${On[$j]}SNP --thresh 1e-4 --max-stored-scores 100000 ${i} OnProm_${On[$j]}SNP.5kbProm.fasta >> fimocommands
			done
		done
	done
done

grep -v 'MouseDerived/JP_out/' fimocommands > fimocommands2 # remove the mouse derived JP

split -l 5988 fimocommands2 fimocommandssplit

# run the above
sh createfimocommands.sh

wc -l fimocommandssplit*
5988 fimocommandssplitaa
5988 fimocommandssplitab
5988 fimocommandssplitac
868 fimocommandssplitad

# create sbatch fimo scripts

for i in fimocommandssplitaa fimocommandssplitab fimocommandssplitac ; do
	echo '#!/bin/bash -e' > ${i}.sh
	echo '#SBATCH -p ei-largemem # partition (queue)' >> ${i}.sh
	echo '#SBATCH -N 1 # number of nodes' >> ${i}.sh
	echo '#SBATCH -n 1 # number of tasks' >> ${i}.sh
	echo '#SBATCH --array=0-5987' >> ${i}.sh
	echo '#SBATCH --mem-per-cpu 200GB' >> ${i}.sh
	echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> ${i}.sh
	echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> ${i}.sh
	echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> ${i}.sh
	echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> ${i}.sh
	printf '\n' >> ${i}.sh
	echo 'ml meme/4.11.1' >> ${i}.sh
	printf '\n' >> ${i}.sh
	echo mapfile -t ${i} '<' ${i} >> ${i}.sh
	echo '${'${i}'[${SLURM_ARRAY_TASK_ID}]}' >> ${i}.sh
done

echo '#!/bin/bash -e' > fimocommandssplitad.sh
echo '#SBATCH -p ei-largemem # partition (queue)' >> fimocommandssplitad.sh
echo '#SBATCH -N 1 # number of nodes' >> fimocommandssplitad.sh
echo '#SBATCH -n 1 # number of tasks' >> fimocommandssplitad.sh
echo '#SBATCH --array=0-867' >> fimocommandssplitad.sh
echo '#SBATCH --mem-per-cpu 200GB' >> fimocommandssplitad.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> fimocommandssplitad.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlhad.ac.uk # send-to address' >> fimocommandssplitad.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> fimocommandssplitad.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> fimocommandssplitad.sh
printf '\n' >> fimocommandssplitad.sh
echo 'ml meme/4.11.1' >> fimocommandssplitad.sh
printf '\n' >> fimocommandssplitad.sh
echo mapfile -t fimocommandssplitad '<' fimocommandssplitad >> fimocommandssplitad.sh
echo '${fimocommandssplitad[${SLURM_ARRAY_TASK_ID}]}' >> fimocommandssplitad.sh

# run the above
sbatch fimocommandssplitaa.sh # DONE
sbatch fimocommandssplitab.sh # DONE
sbatch fimocommandssplitac.sh # RUNNING
sbatch fimocommandssplitad.sh # DONE

# periodically keep running the following to free up space for new outputs:
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/
rm {HumanDerived,MouseDerived}/{CS_out,CW_out,JP_out}/{mz,pn,ab,nb,on}/*/{*.xml,*.gff,*.html} # remove the output files you do not need to free up space
# another option is to temporarily move files to home directory - as of 18/02/2019, /hpc-home/mehtat/ is using 250Gb/1Tb


# ### Below is for running without the splitting the MEME files
# # create array commands for running fimo
#
# echo 'fimo --bgfile AbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_mzSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_ab.meme AbProm_mzSNP.5kbProm.fasta' > runfimo_commands
# echo 'fimo --bgfile AbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_mzSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_ab.meme AbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_mzSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_ab.meme AbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_mzSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_ab.meme AbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_mzSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_ab.meme AbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_mzSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_ab.meme AbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile AbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_pnSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_ab.meme AbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_pnSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_ab.meme AbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_pnSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_ab.meme AbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_pnSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_ab.meme AbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_pnSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_ab.meme AbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_pnSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_ab.meme AbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile AbProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_nbSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_ab.meme AbProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_nbSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_ab.meme AbProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_nbSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_ab.meme AbProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_nbSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_ab.meme AbProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_nbSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_ab.meme AbProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_nbSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_ab.meme AbProm_nbSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile AbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_onSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_ab.meme AbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_onSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_ab.meme AbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_onSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_ab.meme AbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_onSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_ab.meme AbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_onSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_ab.meme AbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile AbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o AbProm_onSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_ab.meme AbProm_onSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile MzProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_pnSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_mz.meme MzProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_pnSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_mz.meme MzProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_pnSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_mz.meme MzProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_pnSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_mz.meme MzProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_pnSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_mz.meme MzProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_pnSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_mz.meme MzProm_pnSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile MzProm_abSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_abSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_mz.meme MzProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_abSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_abSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_mz.meme MzProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_abSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_abSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_mz.meme MzProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_abSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_abSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_mz.meme MzProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_abSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_abSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_mz.meme MzProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_abSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_abSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_mz.meme MzProm_abSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile MzProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_nbSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_mz.meme MzProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_nbSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_mz.meme MzProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_nbSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_mz.meme MzProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_nbSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_mz.meme MzProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_nbSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_mz.meme MzProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_nbSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_mz.meme MzProm_nbSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile MzProm_onSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_onSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_mz.meme MzProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_onSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_onSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_mz.meme MzProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_onSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_onSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_mz.meme MzProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_onSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_onSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_mz.meme MzProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_onSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_onSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_mz.meme MzProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile MzProm_onSNP.5kbProm.0orderMarkovBackgrnd --o MzProm_onSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_mz.meme MzProm_onSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile PnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_mzSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_pn.meme PnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_mzSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_pn.meme PnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_mzSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_pn.meme PnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_mzSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_pn.meme PnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_mzSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_pn.meme PnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_mzSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_pn.meme PnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile PnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_abSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_pn.meme PnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_abSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_pn.meme PnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_abSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_pn.meme PnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_abSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_pn.meme PnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_abSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_pn.meme PnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_abSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_pn.meme PnProm_abSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile PnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_nbSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_pn.meme PnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_nbSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_pn.meme PnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_nbSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_pn.meme PnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_nbSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_pn.meme PnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_nbSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_pn.meme PnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_nbSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_pn.meme PnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile PnProm_onSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_onSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_pn.meme PnProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_onSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_onSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_pn.meme PnProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_onSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_onSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_pn.meme PnProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_onSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_onSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_pn.meme PnProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_onSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_onSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_pn.meme PnProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile PnProm_onSNP.5kbProm.0orderMarkovBackgrnd --o PnProm_onSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_pn.meme PnProm_onSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile NbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_mzSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_nb.meme NbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_mzSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_nb.meme NbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_mzSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_nb.meme NbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_mzSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_nb.meme NbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_mzSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_nb.meme NbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_mzSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_nb.meme NbProm_mzSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile NbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_pnSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_nb.meme NbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_pnSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_nb.meme NbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_pnSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_nb.meme NbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_pnSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_nb.meme NbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_pnSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_nb.meme NbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_pnSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_nb.meme NbProm_pnSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile NbProm_abSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_abSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_nb.meme NbProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_abSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_abSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_nb.meme NbProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_abSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_abSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_nb.meme NbProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_abSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_abSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_nb.meme NbProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_abSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_abSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_nb.meme NbProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_abSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_abSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_nb.meme NbProm_abSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile NbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_onSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_nb.meme NbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_onSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_nb.meme NbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_onSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_nb.meme NbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_onSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_nb.meme NbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_onSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_nb.meme NbProm_onSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile NbProm_onSNP.5kbProm.0orderMarkovBackgrnd --o NbProm_onSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_nb.meme NbProm_onSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile OnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_mzSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_on.meme OnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_mzSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_on.meme OnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_mzSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_on.meme OnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_mzSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_on.meme OnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_mzSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_on.meme OnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_mzSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_mzSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_on.meme OnProm_mzSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile OnProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_pnSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_on.meme OnProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_pnSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_on.meme OnProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_pnSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_on.meme OnProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_pnSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_on.meme OnProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_pnSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_on.meme OnProm_pnSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_pnSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_pnSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_on.meme OnProm_pnSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile OnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_abSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_on.meme OnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_abSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_on.meme OnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_abSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_on.meme OnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_abSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_on.meme OnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_abSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_on.meme OnProm_abSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_abSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_abSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_on.meme OnProm_abSNP.5kbProm.fasta' >> runfimo_commands
#
# echo 'fimo --bgfile OnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_nbSNP.5kbProm_BG_HsCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2a_CS_on.meme OnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_nbSNP.5kbProm_BG_MmCS --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2a_CS_on.meme OnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_nbSNP.5kbProm_BG_HsCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2b_CW_on.meme OnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_nbSNP.5kbProm_BG_MmCW --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2b_CW_on.meme OnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_nbSNP.5kbProm_BG_HsJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/HumanDerived/2c_JASPAR_on.meme OnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
# echo 'fimo --bgfile OnProm_nbSNP.5kbProm.0orderMarkovBackgrnd --o OnProm_nbSNP.5kbProm_BG_MmJP --thresh 1e-4 --max-stored-scores 2500000 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/2c_JASPAR_on.meme OnProm_nbSNP.5kbProm.fasta' >> runfimo_commands
#
# nano runfimo.sh # run on large memory (512GB) with default 90 day time limit
#
# #!/bin/bash -e
# #SBATCH -p ei-largemem # partition (queue)
# #SBATCH -N 1 # number of nodes
# #SBATCH -n 1 # number of tasks
# #SBATCH --array=0-119
# #SBATCH --mem-per-cpu 512GB
# #SBATCH --mail-type=ALL # notifications for job done & fail
# #SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
# #SBATCH -o slurm.%N.%j.out # STDOUT
# #SBATCH -e slurm.%N.%j.err # STDERR
#
# ml meme/4.11.1
#
# mapfile -t runfimo_commands < runfimo_commands # assign as elements to $motifgainbreak_commands variable
# ${runfimo_commands[${SLURM_ARRAY_TASK_ID}]}
#
# # run above
# sbatch runfimo.sh # Using 512Gb per run, CW files roughly took around 12-15 hours, CS files around 4-7 days and JASPAR files around XX days.

################### TEST FILES FOR PYTHON SCRIPT - motifgainbreak_pt2.py ############################################################################
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss/

nano mk_testMutant_promoter_fimo.sh

#!/bin/sh

cat HumanDerived/CS_out/mz/mz.gene.s10.105_CS_mz.meme_pnSNP/fimo.txt MouseDerived/CS_out/mz/RG-cich_GTRDdata_mouse_sites_TF_ig_mz.gene.s10.105.ig_1_CS_mz.meme_pnSNP/fimo.txt |
awk '($8 < 0.05 )' | # multiple test correction <0.05
sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
sort -k4,4 -rn | uniq -f 1 |
sed $'s/;/\t/g' |
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="2a",$11="0.125"}' OFS='\t' >> MzProm_pnSNP.5kbProm-CS_fimo.txt

TFhsmmpvalmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.mz)
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpvalmz MzProm_pnSNP.5kbProm-CS_fimo.txt | awk '{print $2,$18,$20,$1,$17,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > MzProm_pnSNP.5kbProm-CS_fimo.txt1

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz MzProm_pnSNP.5kbProm-CS_fimo.txt1 |
grep -v 'REMOVEME' |
awk '{print $1,$2,$3,$4,$5,$17,$24,$26,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' | awk '$7=tolower($7)' OFS='\t' | awk '$8=tolower($8)' OFS='\t' > MzProm_pnSNP.5kbProm-CS_fimo.txt2

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/geneNamesMapping2.txt MzProm_pnSNP.5kbProm-CS_fimo.txt2 |
grep -v 'REMOVEME' |
awk '$19=tolower($19)' OFS='\t' |
awk '{print $2,$3,$4,$5,$9,$6,$1,$8,$7,$19,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' |
awk '{OFS="\t"} {if ($8=="null"||$8=="none") $8=$9; print $0}' |
awk '{OFS="\t"} {if ($9=="null"||$9=="none") $9=$10; print $0}' |
awk '{OFS="\t"} {if ($8=="null"||$8=="none") $8=$9; print $0}' |
awk '{OFS="\t"} {if ($10=="null"||$10=="none") $10=$8; print $0}' |
awk '{OFS="\t"} {if ($9=="null"||$9=="none") $9=$10; print $0}' > MzProm_pnSNP.5kbProm-CS_fimo.txt3

# Final Table format to output
# SNP line no. - 2
# motif_pattern - 3
# GeneA_OGID - 4
# GeneA - 5
# GeneA_Symbol - 9
# GeneB_OGID - 6
# GeneB - 1
# GeneB_SymbolDr - 8
# GeneB_SymbolGa - 7
# GeneB_SymbolSp - 19
# start - 10
# stop - 11
# strand - 12
# score - 13
# p-value - 14
# q-value - 15
# sequence - 16
# conf_level - 17
# conf_score - 18
# {Sp}Prom_{Sp}SNP - added below
# interaction_type - added below
# effect - added below
# interaction - added below
# directness - added below
# direction - added below
# layer - added below
# source - added below

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$4]{print $0}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Mz-speciesspecnames_clusterassign.txt MzProm_pnSNP.5kbProm-CS_fimo.txt3 > MzProm_pnSNP.5kbProm-CS_fimo.txt3.2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$7]{print $0}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Mz-speciesspecnames_clusterassign.txt MzProm_pnSNP.5kbProm-CS_fimo.txt3.2 > MzProm_pnSNP.5kbProm-CS_fimo.txt3.3

# At this point, merge same species files but add column for {Sp}Prom_{Sp}SNP from the filename
awk '{print $0, FILENAME}' OFS='\t' MzProm_pnSNP.5kbProm-CS_fimo.txt3.3 | sed 's/.5kbProm-CS_fimo.txt3//g' | sed 's/.5kbProm-CW_fimo.txt3//g' | sed 's/.5kbProm-JP_fimo.txt3//g' | sed 's/SNP.3/SNP/g' > MzProm_pnSNP.5kbProm-CS_fimo.txt3.4


# add colheaders
printf 'SNP_line\tmotif_pattern\tGeneA_OGID\tGeneA\tGeneA_Symbol\tGeneB_OGID\tGeneB\tGeneB_SymbolDr\tGeneB_SymbolGa\tGeneB_SymbolSp\tstart\tstop\tstrand\tscore\tp-value\tq-value\tsequence\tconf_level\tconf_score\tSpProm_spSNP\n' > colheads_tfbs
cat colheads_tfbs MzProm_pnSNP.5kbProm-CS_fimo.txt3.4 > MzProm_pnSNP.5kbProm-CS_fimo.txt3.5

# add interaction_type column, effect column, interaction column, directness column, direction column, layer column, and source column
awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"PD_directed"}' MzProm_pnSNP.5kbProm-CS_fimo.txt3.5 |
awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}' |
awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}' |
awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"direct"}' |
awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"directed"}' |
awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Transcriptional_regulation"}' |
awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"promoter_motif"}' > MzProm_pnSNP.5kbProm-CS_fimo.txt3.6 ### this is the test Mutant_promoter_fimo(mainsp_vs_ancspSNP) you can use

rm MzProm_pnSNP.5kbProm-CS_fimo.txt
rm MzProm_pnSNP.5kbProm-CS_fimo.txt1
rm MzProm_pnSNP.5kbProm-CS_fimo.txt2
rm MzProm_pnSNP.5kbProm-CS_fimo.txt3
rm MzProm_pnSNP.5kbProm-CS_fimo.txt3.2
rm MzProm_pnSNP.5kbProm-CS_fimo.txt3.3
rm MzProm_pnSNP.5kbProm-CS_fimo.txt3.4
rm MzProm_pnSNP.5kbProm-CS_fimo.txt3.5


############################################################################################################################



## join the individual scan files based on per species-SNP and dump into separate folders with nomenclature AbProm_mzSNP.5kbProm_BG_HsCS/fimo.txt

nano joinfimo.sh

#!/bin/bash -e
#SBATCH -p ei-largemem # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 512GB # memory pool for all cores
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

WD=('/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss')

Ab=("mz" "pn" "nb" "on")
Mz=("pn" "ab" "nb" "on")
Pn=("mz" "ab" "nb" "on")
Nb=("mz" "pn" "ab" "on")
On=("mz" "pn" "ab" "nb")

length=${#Ab[@]} # get the length of the array - applicable to all

echo '1. arrays defined'

# do the loop
for ((i=0;i<=$length-1;i++)); do
	mkdir -p $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_MmCS
	mkdir -p $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_MmCS
	mkdir -p $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_MmCS
	mkdir -p $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_MmCS
	mkdir -p $WD/OnProm_${On[$i]}SNP.5kbProm_BG_MmCS
	mkdir -p $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsCS
	mkdir -p $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_HsCS
	mkdir -p $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_HsCS
	mkdir -p $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_HsCS
	mkdir -p $WD/OnProm_${On[$i]}SNP.5kbProm_BG_HsCS
	mkdir -p $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_MmCW
	mkdir -p $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_MmCW
	mkdir -p $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_MmCW
	mkdir -p $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_MmCW
	mkdir -p $WD/OnProm_${On[$i]}SNP.5kbProm_BG_MmCW
	mkdir -p $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsCW
	mkdir -p $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_HsCW
	mkdir -p $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_HsCW
	mkdir -p $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_HsCW
	mkdir -p $WD/OnProm_${On[$i]}SNP.5kbProm_BG_HsCW
	mkdir -p $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsJP
	mkdir -p $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_HsJP
	mkdir -p $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_HsJP
	mkdir -p $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_HsJP
	mkdir -p $WD/OnProm_${On[$i]}SNP.5kbProm_BG_HsJP
	cat HumanDerived/CS_out/ab/*ab.meme_${Ab[$i]}SNP/fimo.txt >> $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsCS/fimo.txt
	cat HumanDerived/CS_out/mz/*mz.meme_${Mz[$i]}SNP/fimo.txt >> $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_HsCS/fimo.txt
	cat HumanDerived/CS_out/pn/*pn.meme_${Pn[$i]}SNP/fimo.txt >> $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_HsCS/fimo.txt
	cat HumanDerived/CS_out/nb/*nb.meme_${Nb[$i]}SNP/fimo.txt >> $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_HsCS/fimo.txt
	cat HumanDerived/CS_out/on/*on.meme_${On[$i]}SNP/fimo.txt >> $WD/OnProm_${On[$i]}SNP.5kbProm_BG_HsCS/fimo.txt
	cat MouseDerived/CS_out/ab/*ab.meme_${Ab[$i]}SNP/fimo.txt >> $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_MmCS/fimo.txt
	cat MouseDerived/CS_out/mz/*mz.meme_${Mz[$i]}SNP/fimo.txt >> $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_MmCS/fimo.txt
	cat MouseDerived/CS_out/pn/*pn.meme_${Pn[$i]}SNP/fimo.txt >> $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_MmCS/fimo.txt
	cat MouseDerived/CS_out/nb/*nb.meme_${Nb[$i]}SNP/fimo.txt >> $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_MmCS/fimo.txt
	cat MouseDerived/CS_out/on/*on.meme_${On[$i]}SNP/fimo.txt >> $WD/OnProm_${On[$i]}SNP.5kbProm_BG_MmCS/fimo.txt
	cat HumanDerived/CW_out/ab/*${Ab[$i]}SNP/fimo.txt >> $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsCW/fimo.txt
	cat HumanDerived/CW_out/mz/*${Mz[$i]}SNP/fimo.txt >> $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_HsCW/fimo.txt
	cat HumanDerived/CW_out/pn/*${Pn[$i]}SNP/fimo.txt >> $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_HsCW/fimo.txt
	cat HumanDerived/CW_out/nb/*${Nb[$i]}SNP/fimo.txt >> $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_HsCW/fimo.txt
	cat HumanDerived/CW_out/on/*${On[$i]}SNP/fimo.txt >> $WD/OnProm_${On[$i]}SNP.5kbProm_BG_HsCW/fimo.txt
	cat MouseDerived/CW_out/ab/*${Ab[$i]}SNP/fimo.txt >> $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_MmCW/fimo.txt
	cat MouseDerived/CW_out/mz/*${Mz[$i]}SNP/fimo.txt >> $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_MmCW/fimo.txt
	cat MouseDerived/CW_out/pn/*${Pn[$i]}SNP/fimo.txt >> $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_MmCW/fimo.txt
	cat MouseDerived/CW_out/nb/*${Nb[$i]}SNP/fimo.txt >> $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_MmCW/fimo.txt
	cat MouseDerived/CW_out/on/*${On[$i]}SNP/fimo.txt >> $WD/OnProm_${On[$i]}SNP.5kbProm_BG_MmCW/fimo.txt
	cat HumanDerived/JP_out/ab/*${Ab[$i]}SNP/fimo.txt >> $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsJP/fimo.txt
	cat HumanDerived/JP_out/mz/*${Mz[$i]}SNP/fimo.txt >> $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_HsJP/fimo.txt
	cat HumanDerived/JP_out/pn/*${Pn[$i]}SNP/fimo.txt >> $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_HsJP/fimo.txt
	cat HumanDerived/JP_out/nb/*${Nb[$i]}SNP/fimo.txt >> $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_HsJP/fimo.txt
	cat HumanDerived/JP_out/on/*${On[$i]}SNP/fimo.txt >> $WD/OnProm_${On[$i]}SNP.5kbProm_BG_HsJP/fimo.txt
done

# run the above
sbatch joinfimo.sh


# 1ai. create unified (where relevant) mouse and human sets
# 1aii. Filter outputs (where relevant) for significant q-val (0.05)
# 1aiii. Create confidence levels for each scan - NOTE: will not have extrapolated so only conf. level 2

# Confidence level 2 # {DONE}
# Confidence level 2a - FDR corrected interactions resulting from the FIMO scans using the extrapolated cichlid species specific matrices: 0.125 #{DONE}


nano mergeCSfimo.sh

#!/bin/bash -e
#SBATCH -p ei-largemem # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of cores
#SBATCH --mem 512GB # memory pool for all cores
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

hspval=(/tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mat_qual_pvals_ALL.out2)
mmpval=(/tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mm10_mat_qual_pvals_ALL.out1)
OGIDab=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab)
OGIDmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz)
OGIDpn=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn)
OGIDnb=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb)
OGIDon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on)
TFhsmmpval=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF)
ENShsmmpval=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.ENS)

WD=('/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss')

# Define the arrays
Ab=("mz" "pn" "nb" "on")
Mz=("pn" "ab" "nb" "on")
Pn=("mz" "ab" "nb" "on")
Nb=("mz" "pn" "ab" "on")
On=("mz" "pn" "ab" "nb")

length=${#Ab[@]} # get the length of the array - applicable to all

echo '1. arrays defined'

# do the loop
for ((i=0;i<=$length-1;i++)); do
	# echo -e "$WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsCS/fimo.txt $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_MmCS/fimo.txt"
	cat $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsCS/fimo.txt $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_MmCS/fimo.txt | # cat human and mouse Ab CS files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="2a",$11="0.125"}' OFS='\t' >> $WD/AbProm_${Ab[$i]}SNP.5kbProm-CS_fimo.txt ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	echo "Ab files done"
	cat $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_HsCS/fimo.txt $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_MmCS/fimo.txt | # cat human and mouse Mz CS files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="2a",$11="0.125"}' OFS='\t' >> $WD/MzProm_${Mz[$i]}SNP.5kbProm-CS_fimo.txt ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	echo "Mz files done"
	cat $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_HsCS/fimo.txt $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_MmCS/fimo.txt | # cat human and mouse Pn CS files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="2a",$11="0.125"}' OFS='\t' >> $WD/PnProm_${Pn[$i]}SNP.5kbProm-CS_fimo.txt ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	echo "Pn files done"
	cat $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_HsCS/fimo.txt $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_MmCS/fimo.txt | # cat human and mouse Nb CS files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="2a",$11="0.125"}' OFS='\t' >> $WD/NbProm_${Nb[$i]}SNP.5kbProm-CS_fimo.txt ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	echo "Nb files done"
	cat $WD/OnProm_${On[$i]}SNP.5kbProm_BG_HsCS/fimo.txt $WD/OnProm_${On[$i]}SNP.5kbProm_BG_MmCS/fimo.txt | # cat human and mouse On CS files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="2a",$11="0.125"}' OFS='\t' >> $WD/OnProm_${On[$i]}SNP.5kbProm-CS_fimo.txt ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	echo "On files done"
done

echo '2. concatenation of mouse + human CS files, filtering, sorting and reporting unique hits done'

# for i in mz pn nb on ; do
# 	for j in $WD/AbProm_${i}SNP.5kbProm_BG_HsCS/fimo.txt ; do
# 		for k in $WD/AbProm_${i}SNP.5kbProm_BG_MmCS/fimo.txt ; do
# 			cat $j $k | # cat human and mouse files
# 			awk '($8 < 0.05 )' | # multiple test correction <0.05
# 			sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
# 			sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
# 			awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
# 			sort -k4,4 -rn | uniq -f 1 |
# 			sed $'s/;/\t/g' |
# 			awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="2a",$11="0.125"}' OFS='\t' >> $WD/AbProm_${i}SNP.5kbProm-CS_fimo.txt ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
# 			# cut -f1 $i-CS_default_fimo.out | sort -u > $i-CS_default_fimo.CUT.out # cut out merged first col and create no duplicates file
# 			# use cut file to match only the first hits in the original file which should take the top hit q-val
# 		done
# 	done
# done

# Run above
sbatch mergeCSfimo.sh


# Confidence level 2b - FDR corrected interactions resulting from the FIMO scans using the extrapolated non-species specific matrices: 0.110 # {DONE}

nano mergeCWfimo.sh

#!/bin/bash -e
#SBATCH -p ei-largemem # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of cores
#SBATCH --mem 512GB # memory pool for all cores
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

WD=('/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss')
ENShsmmpval=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.ENS)

# Define the arrays
Ab=("mz" "pn" "nb" "on")
Mz=("pn" "ab" "nb" "on")
Pn=("mz" "ab" "nb" "on")
Nb=("mz" "pn" "ab" "on")
On=("mz" "pn" "ab" "nb")

length=${#Ab[@]} # get the length of the array - applicable to all

echo '1. arrays defined'

# do the loop
for ((i=0;i<=$length-1;i++)); do
	# echo -e "$WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsCW/fimo.txt $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_MmCW/fimo.txt"
	cat $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsCW/fimo.txt $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_MmCW/fimo.txt | # cat human and mouse Ab CW files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level.ig_1//g' |
	sed 's/RG-cich_GTRDhuman_human_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level//g' >> $WD/AbProm_${Ab[$i]}SNP.5kbProm-CW_fimo.txt1 ;
	cat $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_HsCW/fimo.txt $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_MmCW/fimo.txt | # cat human and mouse Ab CW files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level.ig_1//g' |
	sed 's/RG-cich_GTRDhuman_human_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level//g' >> $WD/MzProm_${Mz[$i]}SNP.5kbProm-CW_fimo.txt1 ;
	cat $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_HsCW/fimo.txt $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_MmCW/fimo.txt | # cat human and mouse Ab CW files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level.ig_1//g' |
	sed 's/RG-cich_GTRDhuman_human_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level//g' >> $WD/PnProm_${Pn[$i]}SNP.5kbProm-CW_fimo.txt1 ;
	cat $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_HsCW/fimo.txt $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_MmCW/fimo.txt | # cat human and mouse Ab CW files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level.ig_1//g' |
	sed 's/RG-cich_GTRDhuman_human_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level//g' >> $WD/NbProm_${Nb[$i]}SNP.5kbProm-CW_fimo.txt1 ;
	cat $WD/OnProm_${On[$i]}SNP.5kbProm_BG_HsCW/fimo.txt $WD/OnProm_${On[$i]}SNP.5kbProm_BG_MmCW/fimo.txt | # cat human and mouse Ab CW files
	awk '($8 < 0.05 )' | # multiple test correction <0.05
	sed 's/RG-cich_GTRDdata_mouse_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level.ig_1//g' |
	sed 's/RG-cich_GTRDhuman_human_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level//g' >> $WD/OnProm_${On[$i]}SNP.5kbProm-CW_fimo.txt1 ;
done

echo '2. All CW files merged by species, filtered for q-val and formatted > '

for ((i=0;i<=$length-1;i++)); do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $ENShsmmpval $WD/AbProm_${Ab[$i]}SNP.5kbProm-CW_fimo.txt1 |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19="2b",$20="0.110"}' OFS='\t' >> $WD/AbProm_${Ab[$i]}SNP.5kbProm-CW_fimo.txt2 ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $ENShsmmpval $WD/MzProm_${Mz[$i]}SNP.5kbProm-CW_fimo.txt1 |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19="2b",$20="0.110"}' OFS='\t' >> $WD/MzProm_${Mz[$i]}SNP.5kbProm-CW_fimo.txt2 ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $ENShsmmpval $WD/PnProm_${Pn[$i]}SNP.5kbProm-CW_fimo.txt1 |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19="2b",$20="0.110"}' OFS='\t' >> $WD/PnProm_${Pn[$i]}SNP.5kbProm-CW_fimo.txt2 ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $ENShsmmpval $WD/NbProm_${Nb[$i]}SNP.5kbProm-CW_fimo.txt1 |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19="2b",$20="0.110"}' OFS='\t' >> $WD/NbProm_${Nb[$i]}SNP.5kbProm-CW_fimo.txt2 ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $ENShsmmpval $WD/OnProm_${On[$i]}SNP.5kbProm-CW_fimo.txt1 |
	sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
	awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' | # merge the 5 cols, separated by colon
	sort -k4,4 -rn | uniq -f 1 |
	sed $'s/;/\t/g' |
	awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19="2b",$20="0.110"}' OFS='\t' >> $WD/OnProm_${On[$i]}SNP.5kbProm-CW_fimo.txt2 ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
done

# run the above
sbatch mergeCWfimo.sh
rm *SNP.5kbProm-CW_fimo.txt1 # remove intermediate files

nano TFmap_mergeJPfimo.sh

#!/bin/bash -e
#SBATCH -p ei-largemem # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of cores
#SBATCH --mem 512GB # memory pool for all cores
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

WD=('/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss')

# Confidence level 2c - FDR corrected interactions resulting from the FIMO scans using the Jaspar matrices: 0.115 # {DONE}
# this is done with TF mapping below

# Map the first column (TF) to gene symbol/cichlid gene ID (where relevant - for 1a, 1b, 1c, 2a, 2b)
# Also, for TG mapping below, will need to split the geneID and the SNPID into their own columns - col header SNPID

# species-specific file - cichlid gene first column
TFhsmmpvalmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.mz)
TFhsmmpvalpn=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.pn)
TFhsmmpvalab=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.ab)
TFhsmmpvalnb=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.nb)
TFhsmmpvalon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.on)

# CS TFs
for i in MzProm_*SNP.5kbProm-CS_fimo.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpvalmz $i | awk '{print $2,$18,$20,$1,$17,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > ${i}1
done
for i in PnProm_*SNP.5kbProm-CS_fimo.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpvalpn $i | awk '{print $2,$18,$20,$1,$17,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > ${i}1
done
for i in AbProm_*SNP.5kbProm-CS_fimo.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpvalab $i | awk '{print $2,$18,$20,$1,$17,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > ${i}1
done
for i in NbProm_*SNP.5kbProm-CS_fimo.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpvalnb $i | awk '{print $2,$18,$20,$1,$17,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > ${i}1
done
for i in OnProm_*SNP.5kbProm-CS_fimo.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpvalon $i | awk '{print $2,$18,$20,$1,$17,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > ${i}1
done
# rm *Prom_*SNP.5kbProm-CS_fimo.txt

# CW TFs
# just select required columns from files created above, then remove old files
# these files are ready for TG mapping (as first column)
for i in MzProm_*SNP.5kbProm-CW_fimo.txt2; do
	awk '{print $2,$16,$18,$11,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' $i | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > "$(basename "${i}" 2)3"
done
for i in PnProm_*SNP.5kbProm-CW_fimo.txt2; do
	awk '{print $2,$16,$18,$12,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' $i | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > "$(basename "${i}" 2)3"
done
for i in AbProm_*SNP.5kbProm-CW_fimo.txt2; do
	awk '{print $2,$16,$18,$13,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' $i | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > "$(basename "${i}" 2)3"
done
for i in NbProm_*SNP.5kbProm-CW_fimo.txt2; do
	awk '{print $2,$16,$18,$14,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' $i | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > "$(basename "${i}" 2)3"
done
for i in OnProm_*SNP.5kbProm-CW_fimo.txt2; do
	awk '{print $2,$16,$18,$15,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' $i | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > "$(basename "${i}" 2)3"
done
rm *Prom_*SNP.5kbProm-CW_fimo.txt1
for i in *Prom_*SNP.5kbProm-CW_fimo.txt3 ; do
	mv $i "$(basename "$i" 3)1"
done
rm *Prom_*SNP.5kbProm-CW_fimo.txt2

# JASPAR TFs
TFhsmmpval=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF)

# Define the arrays
Ab=("mz" "pn" "nb" "on")
Mz=("pn" "ab" "nb" "on")
Pn=("mz" "ab" "nb" "on")
Nb=("mz" "pn" "ab" "on")
On=("mz" "pn" "ab" "nb")

length=${#Ab[@]} # get the length of the array - applicable to all

# do the loop
for ((i=0;i<=$length-1;i++)); do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' $TFhsmmpval $WD/AbProm_${Ab[$i]}SNP.5kbProm_BG_HsJP/fimo.txt | grep -v 'REMOVEME' |
	awk '($8 < 0.05 )' | awk '{print $2,$16,$18,$13,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' | awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > $WD/AbProm_${Ab[$i]}SNP.5kbProm-JP_fimo.txt1 # this file is ready for TG mapping (as first column)
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' $TFhsmmpval $WD/MzProm_${Mz[$i]}SNP.5kbProm_BG_HsJP/fimo.txt | grep -v 'REMOVEME' |
	awk '($8 < 0.05 )' | awk '{print $2,$16,$18,$11,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' | awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > $WD/MzProm_${Mz[$i]}SNP.5kbProm-JP_fimo.txt1 # this file is ready for TG mapping (as first column)
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' $TFhsmmpval $WD/PnProm_${Pn[$i]}SNP.5kbProm_BG_HsJP/fimo.txt | grep -v 'REMOVEME' |
	awk '($8 < 0.05 )' | awk '{print $2,$16,$18,$12,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' | awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > $WD/PnProm_${Pn[$i]}SNP.5kbProm-JP_fimo.txt1 # this file is ready for TG mapping (as first column)
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' $TFhsmmpval $WD/NbProm_${Nb[$i]}SNP.5kbProm_BG_HsJP/fimo.txt | grep -v 'REMOVEME' |
	awk '($8 < 0.05 )' | awk '{print $2,$16,$18,$14,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' | awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > $WD/NbProm_${Nb[$i]}SNP.5kbProm-JP_fimo.txt1 # this file is ready for TG mapping (as first column)
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' $TFhsmmpval $WD/OnProm_${On[$i]}SNP.5kbProm_BG_HsJP/fimo.txt | grep -v 'REMOVEME' |
	awk '($8 < 0.05 )' | awk '{print $2,$16,$18,$15,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' | awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' | awk -F'_SNPid:' '{print $1,$2}' OFS='\t' > $WD/OnProm_${On[$i]}SNP.5kbProm-JP_fimo.txt1 # this file is ready for TG mapping (as first column)
done

# At this point, all TF mapped files are *_fimo.txt1
# 4. Map the first column (TG) to gene symbol (run for all levels)- creating a file with the required columns of motif_ID,TF_OGID, TF_cichlid_ID, TF_gene_symbol,TG_OGID, TG_cichlid_ID, TG_Dr_symbol, TG_Ga_symbol

# Create a filelist of species-specific OGID files that are repeated for the while loop
for i in {1..3} ; do
	for j in mz pn ab nb on ; do
		for k in {1..4} ; do
			echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-$j >> OGIDs_list
		done
	done
done

# Create filelist of processed fimo files in order of Mz Pn Ab Nb On
for i in Mz Pn Ab Nb On ; do
	ls -1 ${i}Prom_*SNP.5kbProm-CS_fimo.txt1 >> conf1a1b1c2a2bTG
done
for i in Mz Pn Ab Nb On ; do
	ls -1 ${i}Prom_*SNP.5kbProm-CW_fimo.txt1 >> conf1a1b1c2a2bTG
done
for i in Mz Pn Ab Nb On ; do
	ls -1 ${i}Prom_*SNP.5kbProm-JP_fimo.txt1 >> conf1a1b1c2a2bTG
done

# Final Table format to output (below we are outputting the TG as first column for mapping in next step)
# SNP line no. - 2
# motif_pattern - 3
# GeneA_OGID - 4
# GeneA - 5
# GeneA_Symbol - 9
# GeneB_OGID - 6
# GeneB - 1
# GeneB_SymbolDr - 8
# GeneB_SymbolGa - 7
# GeneB_SymbolSp - 19
# start - 10
# stop - 11
# strand - 12
# score - 13
# p-value - 14
# q-value - 15
# sequence - 16
# conf_level - 17
# conf_score - 18
# {Sp}Prom_{Sp}SNP - added below
# interaction_type - added below
# effect - added below
# interaction - added below
# directness - added below
# direction - added below
# layer - added below
# source - added below

while read -u 3 -r file1 && read -u 4 -r file2
do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' ${file1} ${file2} |
	grep -v 'REMOVEME' |
	awk '{print $1,$2,$3,$4,$5,$17,$24,$26,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' | awk '$7=tolower($7)' OFS='\t' | awk '$8=tolower($8)' OFS='\t' > "$(basename "${file2}" 1)2"
done 3<OGIDs_list 4<conf1a1b1c2a2bTG

# mapping to geneNames file here to create an extra column
for i in Mz Pn Ab Nb On ; do
	for j in ${i}Prom_*SNP.5kbProm-*.txt2 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/geneNamesMapping2.txt $j |
		grep -v 'REMOVEME' |
		awk '$19=tolower($19)' OFS='\t' |
		awk '{print $2,$3,$4,$5,$9,$6,$1,$8,$7,$19,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' |
		awk '{OFS="\t"} {if ($8=="null"||$8=="none") $8=$9; print $0}' |
		awk '{OFS="\t"} {if ($9=="null"||$9=="none") $9=$10; print $0}' |
		awk '{OFS="\t"} {if ($8=="null"||$8=="none") $8=$9; print $0}' |
		awk '{OFS="\t"} {if ($10=="null"||$10=="none") $10=$8; print $0}' |
		awk '{OFS="\t"} {if ($9=="null"||$9=="none") $9=$10; print $0}' > "$(basename "$j" 2)3"
	done
done

# filter, for the presence of genes in modules only - need to do for both GeneA and GeneB in modules
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Mz-speciesspecnames_clusterassign.txt .
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Pn-speciesspecnames_clusterassign.txt .
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Ab-speciesspecnames_clusterassign.txt .
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Nb-speciesspecnames_clusterassign.txt .
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/On-speciesspecnames_clusterassign.txt .

for i in Mz Pn Ab Nb On ; do
	for j in ${i}Prom_*SNP.5kbProm-*.txt3 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$4]{print $0}' $i-speciesspecnames_clusterassign.txt $j > $j.2
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$7]{print $0}' $i-speciesspecnames_clusterassign.txt $j.2 > $j.3
	done
done

rm *.txt3.2

# At this point, merge same species files but add column for {Sp}Prom_{Sp}SNP from the filename
for i in Mz Pn Ab Nb On ; do
	for j in ${i}Prom_*SNP.5kbProm-*.txt3.3 ; do
		awk '{print $0, FILENAME}' OFS='\t' $j | sed 's/.5kbProm-CS_fimo.txt3.3//g' | sed 's/.5kbProm-CW_fimo.txt3.3//g' | sed 's/.5kbProm-JP_fimo.txt3.3//g' > "$(basename "$j" .3).4"
	done
done

cat Mz*.txt3.4 > MzProm_allSNP.5kbProm-CS_CW_JP.txt
cat Pn*.txt3.4 > PnProm_allSNP.5kbProm-CS_CW_JP.txt
cat Ab*.txt3.4 > AbProm_allSNP.5kbProm-CS_CW_JP.txt
cat Nb*.txt3.4 > NbProm_allSNP.5kbProm-CS_CW_JP.txt
cat On*.txt3.4 > OnProm_allSNP.5kbProm-CS_CW_JP.txt

rm {*txt1,*txt2,*txt3,*txt3.4}
# keep *txt3.3 files that are final species and confidence specific files - perfect for checking number of rows etc.


# add colheaders, interaction_type column, effect column, interaction column, directness column, direction column, layer column, and source column
printf 'SNP_line\tmotif_pattern\tGeneA_OGID\tGeneA\tGeneA_Symbol\tGeneB_OGID\tGeneB\tGeneB_SymbolDr\tGeneB_SymbolGa\tGeneB_SymbolSp\tstart\tstop\tstrand\tscore\tp-value\tq-value\tsequence\tconf_level\tconf_score\tSpProm_spSNP\n' > colheads_tfbs
for i in *Prom_allSNP.5kbProm-CS_CW_JP.txt ; do
	cat colheads_tfbs $i |
	awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"PD_directed"}' |
	awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}' |
	awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}' |
	awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"direct"}' |
	awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"directed"}' |
	awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Transcriptional_regulation"}' |
	awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"promoter_motif"}' > ${i}1
done

rm *Prom_allSNP.5kbProm-CS_CW_JP.txt
# these are the collated Mutant_promoter_fimo(mainsp_vs_ancspSNP) files
mv MzProm_allSNP.5kbProm-CS_CW_JP.txt1 MzProm_allSNP.5kbProm-CS_CW_JP.txt
mv PnProm_allSNP.5kbProm-CS_CW_JP.txt1 PnProm_allSNP.5kbProm-CS_CW_JP.txt
mv AbProm_allSNP.5kbProm-CS_CW_JP.txt1 AbProm_allSNP.5kbProm-CS_CW_JP.txt
mv NbProm_allSNP.5kbProm-CS_CW_JP.txt1 NbProm_allSNP.5kbProm-CS_CW_JP.txt
mv OnProm_allSNP.5kbProm-CS_CW_JP.txt1 OnProm_allSNP.5kbProm-CS_CW_JP.txt

# Create species specific files for the SNP - these are the final Mutant_promoter_fimo(mainsp_vs_ancspSNP) files to input into motifgainbreak_pt2.py
awk 'NR==1 || /MzProm_pnSNP/' MzProm_allSNP.5kbProm-CS_CW_JP.txt > MzProm_pnSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /MzProm_abSNP/' MzProm_allSNP.5kbProm-CS_CW_JP.txt > MzProm_abSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /MzProm_nbSNP/' MzProm_allSNP.5kbProm-CS_CW_JP.txt > MzProm_nbSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /MzProm_onSNP/' MzProm_allSNP.5kbProm-CS_CW_JP.txt > MzProm_onSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /PnProm_mzSNP/' PnProm_allSNP.5kbProm-CS_CW_JP.txt > PnProm_mzSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /PnProm_abSNP/' PnProm_allSNP.5kbProm-CS_CW_JP.txt > PnProm_abSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /PnProm_nbSNP/' PnProm_allSNP.5kbProm-CS_CW_JP.txt > PnProm_nbSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /PnProm_onSNP/' PnProm_allSNP.5kbProm-CS_CW_JP.txt > PnProm_onSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /AbProm_mzSNP/' AbProm_allSNP.5kbProm-CS_CW_JP.txt > AbProm_mzSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /AbProm_pnSNP/' AbProm_allSNP.5kbProm-CS_CW_JP.txt > AbProm_pnSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /AbProm_nbSNP/' AbProm_allSNP.5kbProm-CS_CW_JP.txt > AbProm_nbSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /AbProm_onSNP/' AbProm_allSNP.5kbProm-CS_CW_JP.txt > AbProm_onSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /NbProm_mzSNP/' NbProm_allSNP.5kbProm-CS_CW_JP.txt > NbProm_mzSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /NbProm_pnSNP/' NbProm_allSNP.5kbProm-CS_CW_JP.txt > NbProm_pnSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /NbProm_abSNP/' NbProm_allSNP.5kbProm-CS_CW_JP.txt > NbProm_abSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /NbProm_onSNP/' NbProm_allSNP.5kbProm-CS_CW_JP.txt > NbProm_onSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /OnProm_mzSNP/' OnProm_allSNP.5kbProm-CS_CW_JP.txt > OnProm_mzSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /OnProm_pnSNP/' OnProm_allSNP.5kbProm-CS_CW_JP.txt > OnProm_pnSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /OnProm_abSNP/' OnProm_allSNP.5kbProm-CS_CW_JP.txt > OnProm_abSNP.5kbProm-CS_CW_JP.txt
awk 'NR==1 || /OnProm_nbSNP/' OnProm_allSNP.5kbProm-CS_CW_JP.txt > OnProm_nbSNP.5kbProm-CS_CW_JP.txt

# check all rows have same number of columns
for i in *Prom_allSNP.5kbProm-CS_CW_JP.txt ; do
	awk '{print NF}' $i | sort -nu | head -n 1
	awk '{print NF}' $i | sort -nu | tail -n 1
done

# remove orginal individual fimo scans
rm -r {HumanDerived,MouseDerived}

# tar up the merged FIMO runs
tar -czvf SpProm_SpSNP.5kbProm_BG_HsMm_CSCWJP_FIMOcollated_NOFILT.tar.gz *Prom_*SNP.5kbProm_BG_{Hs,Mm}*
rm -r *Prom_*SNP.5kbProm_BG_{Hs,Mm}*


# run the above
sbatch TFmap_mergeJPfimo.sh

# tar up the mutant promoter FASTA files too

nano tarFASTA.sh

#!/bin/bash -e
#SBATCH -p ei-largemem # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of cores
#SBATCH --mem 512GB # memory pool for all cores
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

tar -czvf SpProm_SpSNP_5kbProm.fasta.tar.gz *Prom_*SNP.5kbProm.fasta
rm *Prom_*SNP.5kbProm.fasta

# run the above
sbatch tarFASTA.sh

# then, prepare a python script (/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Motif_gain_and_loss/motifgainbreak_pt2.py)
# to parse the filtered FIMO output for whether the motif is gained or lost as a result of the SNP
# for this, add OGIDs to the files: /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt

for i in /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/ab1_{nb1,on11,pn1,mz11}_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$4]=$1;next}{if(a[$16]){print $0,a[$16];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i > $i.OGID2
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$4]=$1;next}{if(a[$22]){print $0,a[$22];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i.OGID2 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$30,$17,$18,$19,$20,$21,$22,$31,$23,$24,$25,$26,$27,$28,$29}' OFS='\t' > $i.OGID
done

for i in /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/pn1_{ab1,nb1,on11,mz11}_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$3]=$1;next}{if(a[$16]){print $0,a[$16];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i > $i.OGID2
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$3]=$1;next}{if(a[$22]){print $0,a[$22];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i.OGID2 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$30,$17,$18,$19,$20,$21,$22,$31,$23,$24,$25,$26,$27,$28,$29}' OFS='\t' > $i.OGID
done

for i in /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/mz11_{ab1,nb1,on11,pn1}_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$2]=$1;next}{if(a[$16]){print $0,a[$16];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i > $i.OGID2
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$2]=$1;next}{if(a[$22]){print $0,a[$22];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i.OGID2 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$30,$17,$18,$19,$20,$21,$22,$31,$23,$24,$25,$26,$27,$28,$29}' OFS='\t' > $i.OGID
done

for i in /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/nb1_{ab1,on11,pn1,mz11}_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$1;next}{if(a[$16]){print $0,a[$16];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i > $i.OGID2
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$1;next}{if(a[$22]){print $0,a[$22];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i.OGID2 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$30,$17,$18,$19,$20,$21,$22,$31,$23,$24,$25,$26,$27,$28,$29}' OFS='\t' > $i.OGID
done

for i in /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/on11_{ab1,nb1,pn1,mz11}_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$6]=$1;next}{if(a[$16]){print $0,a[$16];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i > $i.OGID2
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$6]=$1;next}{if(a[$22]){print $0,a[$22];}else{print $0,"NULL";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 $i.OGID2 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$30,$17,$18,$19,$20,$21,$22,$31,$23,$24,$25,$26,$27,$28,$29}' OFS='\t' > $i.OGID
done

rm /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/*.OGID2


# run the motif gain and loss script (motifgainbreak_pt2.py) by creatinf commands to run in an array

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/Motif_gain_and_loss

echo 'python -u motifgainbreak_pt2.py ../mz11_pn1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../pn1_mz11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.OGID MzProm_pnSNP.5kbProm-CS_CW_JP.txt Mz Pn MzProm_pnSNP_Motif_gain_loss.txt' > run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../mz11_ab1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../ab1_mz11_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.OGID MzProm_abSNP.5kbProm-CS_CW_JP.txt Mz Ab MzProm_abSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../mz11_nb1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../nb1_mz11_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.OGID MzProm_nbSNP.5kbProm-CS_CW_JP.txt Mz Nb MzProm_nbSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../mz11_on11_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../on11_mz11_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.OGID MzProm_onSNP.5kbProm-CS_CW_JP.txt Mz On MzProm_onSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../pn1_mz11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../mz11_pn1_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.OGID PnProm_mzSNP.5kbProm-CS_CW_JP.txt Pn Mz PnProm_mzSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../pn1_ab1_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../ab1_pn1_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.OGID PnProm_abSNP.5kbProm-CS_CW_JP.txt Pn Ab PnProm_abSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../pn1_nb1_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../nb1_pn1_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.OGID PnProm_nbSNP.5kbProm-CS_CW_JP.txt Pn Nb PnProm_nbSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../pn1_on11_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../on11_pn1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.OGID PnProm_onSNP.5kbProm-CS_CW_JP.txt Pn On PnProm_onSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../ab1_nb1_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../nb1_ab1_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.OGID AbProm_nbSNP.5kbProm-CS_CW_JP.txt Ab Nb AbProm_nbSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../ab1_on11_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../on11_ab1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.OGID AbProm_onSNP.5kbProm-CS_CW_JP.txt Ab On AbProm_onSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands
echo 'python -u motifgainbreak_pt2.py ../nb1_on11_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.OGID ../on11_nb1_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.OGID NbProm_onSNP.5kbProm-CS_CW_JP.txt Nb On NbProm_onSNP_Motif_gain_loss.txt' >> run_motifgainbreak_pt2_commands

nano run_motifgainbreak_pt2.sh # run on large memory (200GB) with default 90 day time limit

#!/bin/bash -e
#SBATCH -p ei-largemem # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-10
#SBATCH --mem-per-cpu 200GB
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t run_motifgainbreak_pt2_commands < run_motifgainbreak_pt2_commands # assign as elements to $motifgainbreak_commands variable
${run_motifgainbreak_pt2_commands[${SLURM_ARRAY_TASK_ID}]}

# run above
sbatch run_motifgainbreak_pt2.sh










# e. GO enrichment of TFBSs associated variants only

nano runGOINPUT_and_enrichment.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

ml gcc
ml zlib

# Promoters TFBSs {DONE}
for i in mz11*_substitutions_1217*_motifgenomeannot.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 mz11*_substitutions_1217*_motifgenomeannot.txt.IDs > mz1TFBSslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1TFBSslist
for i in mz11*_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> mz11_PromTFBS_collated-GOINPUT ; done

for i in pn1*_substitutions_1217*_motifgenomeannot.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 pn1*_substitutions_1217*_motifgenomeannot.txt.IDs > pn1TFBSslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < pn1TFBSslist
for i in pn1*_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> pn1_PromTFBS_collated-GOINPUT ; done

for i in ab1*_substitutions_1217*_motifgenomeannot.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 ab1*_substitutions_1217*_motifgenomeannot.txt.IDs > ab1TFBSslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < ab1TFBSslist
for i in ab1*_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> ab1_PromTFBS_collated-GOINPUT ; done

for i in nb1*_substitutions_1217*_motifgenomeannot.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 nb1*_substitutions_1217*_motifgenomeannot.txt.IDs > nb1TFBSslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < nb1TFBSslist
for i in nb1*_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> nb1_PromTFBS_collated-GOINPUT ; done

for i in on11*_substitutions_1217*_motifgenomeannot.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 on11*_substitutions_1217*_motifgenomeannot.txt.IDs > on1TFBSslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < on1TFBSslist
for i in on11*_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> on11_PromTFBS_collated-GOINPUT ; done

rm *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs
rm *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt
rm *TFBSslist


# Promoters_TFBSs_switched {DONE}
for i in mz11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Mzvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 mz11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > mz1promTFBSswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1promTFBSswlist
for i in mz11*_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> mz11_promTFBSsw_collated-GOINPUT ; done

for i in pn1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Pnvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 pn1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > pn1promTFBSswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < pn1promTFBSswlist
for i in pn1*_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> pn1_promTFBSsw_collated-GOINPUT ; done

for i in ab1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Abvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 ab1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > ab1promTFBSswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < ab1promTFBSswlist
for i in ab1*_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> ab1_promTFBSsw_collated-GOINPUT ; done

for i in nb1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Nbvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 nb1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > nb1promTFBSswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < nb1promTFBSswlist
for i in nb1*_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> nb1_promTFBSsw_collated-GOINPUT ; done

for i in on11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Onvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 on11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > on11promTFBSswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < on11promTFBSswlist
for i in on11*_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> on11_promTFBSsw_collated-GOINPUT ; done

rm *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs
rm *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt
rm *promTFBSswlist

# Promoters_TFBSs_noswitch {DONE}
for i in mz11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Mzvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 mz11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > mz1promTFBSnoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1promTFBSnoswlist
for i in mz11*_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> mz11_promTFBSnosw_collated-GOINPUT ; done

for i in pn1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Pnvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 pn1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > pn1promTFBSnoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < pn1promTFBSnoswlist
for i in pn1*_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> pn1_promTFBSnosw_collated-GOINPUT ; done

for i in ab1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Abvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 ab1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > ab1promTFBSnoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < ab1promTFBSnoswlist
for i in ab1*_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> ab1_promTFBSnosw_collated-GOINPUT ; done

for i in nb1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Nbvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 nb1*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > nb1promTFBSnoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < nb1promTFBSnoswlist
for i in nb1*_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> nb1_promTFBSnosw_collated-GOINPUT ; done

for i in on11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Onvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | cut -f16 | sort -u > $i.IDs
done
ls -1 on11*_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs > on11promTFBSnoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.IDs//g' > ${F}-GOINPUT.txt ; done < on11promTFBSnoswlist
for i in on11*_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> on11_promTFBSnosw_collated-GOINPUT ; done

rm *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs
rm *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt.IDs-GOINPUT.txt
rm *promTFBSnoswlist

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib2/

#create variables for each GO Input file
GOINPUT_ab1_PromTFBS_collated_GOINPUT=(ab1_PromTFBS_collated-GOINPUT)
GOINPUT_ab1_promTFBSnosw_collated_GOINPUT=(ab1_promTFBSnosw_collated-GOINPUT)
GOINPUT_ab1_promTFBSsw_collated_GOINPUT=(ab1_promTFBSsw_collated-GOINPUT)
GOINPUT_mz11_PromTFBS_collated_GOINPUT=(mz11_PromTFBS_collated-GOINPUT)
GOINPUT_mz11_promTFBSnosw_collated_GOINPUT=(mz11_promTFBSnosw_collated-GOINPUT)
GOINPUT_mz11_promTFBSsw_collated_GOINPUT=(mz11_promTFBSsw_collated-GOINPUT)
GOINPUT_nb1_PromTFBS_collated_GOINPUT=(nb1_PromTFBS_collated-GOINPUT)
GOINPUT_nb1_promTFBSnosw_collated_GOINPUT=(nb1_promTFBSnosw_collated-GOINPUT)
GOINPUT_nb1_promTFBSsw_collated_GOINPUT=(nb1_promTFBSsw_collated-GOINPUT)
GOINPUT_on11_PromTFBS_collated_GOINPUT=(on11_PromTFBS_collated-GOINPUT)
GOINPUT_on11_promTFBSnosw_collated_GOINPUT=(on11_promTFBSnosw_collated-GOINPUT)
GOINPUT_on11_promTFBSsw_collated_GOINPUT=(on11_promTFBSsw_collated-GOINPUT)
GOINPUT_pn1_PromTFBS_collated_GOINPUT=(pn1_PromTFBS_collated-GOINPUT)
GOINPUT_pn1_promTFBSnosw_collated_GOINPUT=(pn1_promTFBSnosw_collated-GOINPUT)
GOINPUT_pn1_promTFBSsw_collated_GOINPUT=(pn1_promTFBSsw_collated-GOINPUT)

BG_abgene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab.gene-genelist)
BG_mzgene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz.gene-genelist)
BG_nbgene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb.gene-genelist)
BG_ongene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on.gene-genelist)
BG_pngene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn.gene-genelist)

mkdir GOOUTPUT

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_PromTFBS_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_promTFBSnosw_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_promTFBSsw_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_PromTFBS_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_promTFBSnosw_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_promTFBSsw_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_PromTFBS_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_promTFBSnosw_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_promTFBSsw_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_PromTFBS_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_promTFBSnosw_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_promTFBSsw_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_PromTFBS_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_promTFBSnosw_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_promTFBSsw_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

# Select those with a q-value (adjusted p-value) of < 0.05 and sort ascending
for i in GOOUTPUT/* ; do awk -F'\t' '($4 < 0.05)' $i | sort -k4,4 -n -t $'\t' > GOOUTPUT/"$(basename "$i" .txt)_filtered.txt" ; done


# for plotting the GO (just doing based on log10 fold enrichment (FDR < 0.05), sorting based on p-value here) prepare the OUTPUT file so that matches will show the GO output on each line for each gene
for i in GOOUTPUT/*details_filtered.txt ; do
 sed 's/$/;/g' $i |
 awk '{ gsub(";", ";"$1";"$2";"$3";"$4";"$5";"$6";"$7";"$8";"$9"\n") } 1' |
 cut -f10 |
 awk '{ gsub(";", "\t") } 1' |
 sed '/^$/d' > GOOUTPUT/"$(basename "$i" )2" ;
done

# remove first column then only display unique entries as final column displays no. of genes assigned to the term if you need it!
for i in GOOUTPUT/*details_filtered.txt2 ; do
 cut -f2-10 $i |
 sort -u | sort -k3,3 > GOOUTPUT/"$(basename "$i" 2)3" ;
done

# then combine same genic region files
cat /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/GOOUTPUT/*_PromTFBS_collated-GOOUTPUT_details_filtered.txt3 > /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/GOOUTPUT/PromTFBS_collated-GOOUTPUT_details_filtered.txt3
cat /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/GOOUTPUT/*_promTFBSnosw_collated-GOOUTPUT_details_filtered.txt3 > /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/GOOUTPUT/promTFBSnosw_collated-GOOUTPUT_details_filtered.txt3
cat /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/GOOUTPUT/*_promTFBSsw_collated-GOOUTPUT_details_filtered.txt3 > /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/GOOUTPUT/promTFBSsw_collated-GOOUTPUT_details_filtered.txt3

# Copied all cat files above to local here: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis
### then plotting /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Multiz_PairwiseVCF_GenomicReg_GOenrichmentplots.R

# Run on uv2k2
qsub -q Test -l select=1:mem=100GB:ncpus=1 runGOINPUT_and_enrichment.sh



##########################################################################################################################

########## 3. GO enrichment of pairwise variants in each gene region to get a good idea of associations {DONE - done TFBSs above}

cd /tgac/workarea/group-vh/cichlids/substitutions_VCFs/

# 1. Prepare GO Input files

uv2k2
nano runGOINPUT.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

# CNEs {DONE}
# run bedtools closest to get the closest promoter and use that as gene ID for GO enrichment
ml bedtools/2.25.0
ml GCC
ml zlib

for i in mz11*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 mz11*_substitutions_1217.bed_substCNEoverlap.txt.IDs > mz11CNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz11CNEslist
for i in mz11*_substitutions_1217.bed_substCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> mz11_CNE_collated-GOINPUT ; done

for i in pn1*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 pn1*_substitutions_1217.bed_substCNEoverlap.txt.IDs > pn1CNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < pn1CNEslist
for i in pn1*_substitutions_1217.bed_substCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> pn1_CNE_collated-GOINPUT ; done

for i in ab1*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 ab1*_substitutions_1217.bed_substCNEoverlap.txt.IDs > ab1CNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < ab1CNEslist
for i in ab1*_substitutions_1217.bed_substCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> ab1_CNE_collated-GOINPUT ; done

for i in nb1*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 nb1*_substitutions_1217.bed_substCNEoverlap.txt.IDs > nb1CNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < nb1CNEslist
for i in nb1*_substitutions_1217.bed_substCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> nb1_CNE_collated-GOINPUT ; done

for i in on11*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 on11*_substitutions_1217.bed_substCNEoverlap.txt.IDs > on11CNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < on11CNEslist
for i in on11*_substitutions_1217.bed_substCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> on11_CNE_collated-GOINPUT ; done

rm *_substitutions_1217.bed_substCNEoverlap.txt.IDs-GOINPUT.txt
rm *_substitutions_1217.bed_substCNEoverlap.txt.IDs
rm *CNEslist


# aCNEs {DONE}
# run bedtools closest to get the closest promoter and use that as gene ID for GO enrichment
ml bedtools/2.25.0
ml GCC
ml zlib

for i in mz11*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 mz11*_substitutions_1217.bed_substaCNEoverlap.txt.IDs > mz11aCNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz11aCNEslist
for i in mz11*_substitutions_1217.bed_substaCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> mz11_aCNE_collated-GOINPUT ; done

for i in pn1*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 pn1*_substitutions_1217.bed_substaCNEoverlap.txt.IDs > pn1aCNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < pn1aCNEslist
for i in pn1*_substitutions_1217.bed_substaCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> pn1_aCNE_collated-GOINPUT ; done

for i in ab1*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 ab1*_substitutions_1217.bed_substaCNEoverlap.txt.IDs > ab1aCNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < ab1aCNEslist
for i in ab1*_substitutions_1217.bed_substaCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> ab1_aCNE_collated-GOINPUT ; done

for i in nb1*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 nb1*_substitutions_1217.bed_substaCNEoverlap.txt.IDs > nb1aCNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < nb1aCNEslist
for i in nb1*_substitutions_1217.bed_substaCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> nb1_aCNE_collated-GOINPUT ; done

for i in on11*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed | cut -f22 | sort -u > $i.IDs ; done
ls -1 on11*_substitutions_1217.bed_substaCNEoverlap.txt.IDs > on11aCNEslist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < on11aCNEslist
for i in on11*_substitutions_1217.bed_substaCNEoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> on11_aCNE_collated-GOINPUT ; done

rm *_substitutions_1217.bed_substaCNEoverlap.txt.IDs-GOINPUT.txt
rm *_substitutions_1217.bed_substaCNEoverlap.txt.IDs
rm *aCNEslist


# Promoters
# switched {DONE}
for i in mz11*_subst-swPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 mz11*_subst-swPromoverlap.txt.IDs > mz1Promswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1Promswlist
for i in mz11*_substitutions_1217.bed_subst-swPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> mz11_Promsw_collated-GOINPUT ; done

for i in pn1*_subst-swPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 pn1*_subst-swPromoverlap.txt.IDs > mz1Promswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1Promswlist
for i in pn1*_substitutions_1217.bed_subst-swPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> pn1_Promsw_collated-GOINPUT ; done

for i in ab1*_subst-swPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 ab1*_subst-swPromoverlap.txt.IDs > mz1Promswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1Promswlist
for i in ab1*_substitutions_1217.bed_subst-swPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> ab1_Promsw_collated-GOINPUT ; done

for i in nb1*_subst-swPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 nb1*_subst-swPromoverlap.txt.IDs > mz1Promswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1Promswlist
for i in nb1*_substitutions_1217.bed_subst-swPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> nb1_Promsw_collated-GOINPUT ; done

for i in on11*_subst-swPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 on11*_subst-swPromoverlap.txt.IDs > mz1Promswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1Promswlist
for i in on11*_substitutions_1217.bed_subst-swPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> on11_Promsw_collated-GOINPUT ; done

rm *_substitutions_1217.bed_subst-swPromoverlap.txt.IDs
rm *_substitutions_1217.bed_subst-swPromoverlap.txt.IDs-GOINPUT.txt
rm *Promswlist

# no switch {DONE}
for i in mz11*_subst-NoswPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 mz11*_subst-NoswPromoverlap.txt.IDs > mz1PromNoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1PromNoswlist
for i in mz11*_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> mz11_PromNosw_collated-GOINPUT ; done

for i in pn1*_subst-NoswPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 pn1*_subst-NoswPromoverlap.txt.IDs > mz1PromNoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1PromNoswlist
for i in pn1*_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> pn1_PromNosw_collated-GOINPUT ; done

for i in ab1*_subst-NoswPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 ab1*_subst-NoswPromoverlap.txt.IDs > mz1PromNoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1PromNoswlist
for i in ab1*_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> ab1_PromNosw_collated-GOINPUT ; done

for i in nb1*_subst-NoswPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 nb1*_subst-NoswPromoverlap.txt.IDs > mz1PromNoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1PromNoswlist
for i in nb1*_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> nb1_PromNosw_collated-GOINPUT ; done

for i in on11*_subst-NoswPromoverlap.txt ; do cut -f16 $i | sort -u > $i.IDs ; done
ls -1 on11*_subst-NoswPromoverlap.txt.IDs > mz1PromNoswlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz1PromNoswlist
for i in on11*_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> on11_PromNosw_collated-GOINPUT ; done

rm *_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs
rm *_substitutions_1217.bed_subst-NoswPromoverlap.txt.IDs-GOINPUT.txt
rm *PromNoswlist

# 3'UTR {DONE}
for i in mz11*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' | sort -u > $i.IDs ; done
ls -1 mz11*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs > mz113UTRlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz113UTRlist
for i in mz11*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> mz11_3UTR_collated-GOINPUT ; done

for i in pn1*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' | sort -u > $i.IDs ; done
ls -1 pn1*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs > pn13UTRlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < pn13UTRlist
for i in pn1*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> pn1_3UTR_collated-GOINPUT ; done

for i in ab1*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' | sort -u > $i.IDs ; done
ls -1 ab1*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs > ab13UTRlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < ab13UTRlist
for i in ab1*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> ab1_3UTR_collated-GOINPUT ; done

for i in nb1*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' | sort -u > $i.IDs ; done
ls -1 nb1*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs > nb13UTRlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < nb13UTRlist
for i in nb1*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i | sed 's/bn./nb./g' >> nb1_3UTR_collated-GOINPUT ; done

for i in on11*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' | sort -u > $i.IDs ; done
ls -1 on11*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs > on113UTRlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < on113UTRlist
for i in on11*_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> on11_3UTR_collated-GOINPUT ; done

rm *_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs
rm *_substitutions_1217.bed_subst-3UTRoverlap.txt.IDs-GOINPUT.txt
rm *3UTRlist

# 3' UTR miRNA binding sites {DONE}
for i in mz11*_subst_miRNA-3UTRoverlap.txt ; do cut -f5 $i | sort -u > $i.IDs ; done
ls -1 mz11*_subst_miRNA-3UTRoverlap.txt.IDs > mz13UTRmiRNAlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz13UTRmiRNAlist
for i in mz11*_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> mz11_3UTRmiRNA_collated-GOINPUT ; done

for i in pn1*_subst_miRNA-3UTRoverlap.txt ; do cut -f5 $i | sort -u > $i.IDs ; done
ls -1 pn1*_subst_miRNA-3UTRoverlap.txt.IDs > mz13UTRmiRNAlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz13UTRmiRNAlist
for i in pn1*_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> pn1_3UTRmiRNA_collated-GOINPUT ; done

for i in ab1*_subst_miRNA-3UTRoverlap.txt ; do cut -f5 $i | sort -u > $i.IDs ; done
ls -1 ab1*_subst_miRNA-3UTRoverlap.txt.IDs > mz13UTRmiRNAlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz13UTRmiRNAlist
for i in ab1*_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> ab1_3UTRmiRNA_collated-GOINPUT ; done

for i in nb1*_subst_miRNA-3UTRoverlap.txt ; do cut -f5 $i | sort -u > $i.IDs ; done
ls -1 nb1*_subst_miRNA-3UTRoverlap.txt.IDs > mz13UTRmiRNAlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz13UTRmiRNAlist
for i in nb1*_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> nb1_3UTRmiRNA_collated-GOINPUT ; done

for i in on11*_subst_miRNA-3UTRoverlap.txt ; do cut -f5 $i | sort -u > $i.IDs ; done
ls -1 on11*_subst_miRNA-3UTRoverlap.txt.IDs > mz13UTRmiRNAlist
while read F ; do sed "1i ${F}" $F | tr '\n' '#' | sed 's/#/\t/' | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs//g' > ${F}-GOINPUT.txt ; done < mz13UTRmiRNAlist
for i in on11*_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs-GOINPUT.txt ; do awk '{print}' $i >> on11_3UTRmiRNA_collated-GOINPUT ; done

rm *_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs
rm *_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.IDs-GOINPUT.txt
rm *3UTRmiRNAlist


# run the above to generate GO Input files
qsub -q Test -l select=1:mem=50GB:ncpus=1 runGOINPUT.sh


# 2. Use unix to create a script to run the GO enrichment

# Background files should be all genes in genome - use the geneNamesTree file for this
for i in mz.gene pn.gene ab.gene nb.gene on.gene ; do grep -wiF $i geneNamesTree | cut -f1 > $i-genelist; done

# create script to run GO enrichment - this only creates some of it, need to amend some!
echo '#!/bin/bash' > runGOenrichment.sh
printf "\n" >> runGOenrichment.sh
echo 'cd $PBS_O_WORKDIR;' >> runGOenrichment.sh
echo 'ml gcc' >> runGOenrichment.sh
echo 'ml zlib' >> runGOenrichment.sh
printf "\n" >> runGOenrichment.sh
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib2/' >> runGOenrichment.sh
printf "\n" >> runGOenrichment.sh
echo '#create variables for each GO Input file' >> runGOenrichment.sh
ls -1 *_collated-GOINPUT | awk '{print "GOINPUT_"$1"=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/"$1")"}' | sed 's/-GOINPUT=/_GOINPUT=/g' >> runGOenrichment.sh
printf "\n" >> runGOenrichment.sh
ls -1 *-genelist | awk '{print "BG_"$1"=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/"$1")"}' | sed 's/-genelist=/_genelist=/g' >> runGOenrichment.sh
printf "\n" >> runGOenrichment.sh
ls -1 *_collated-GOINPUT | awk '{print "/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_"$1" /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/ 1 $(basename -GOINPUT)-GOOUTPUT persg"}' | sed 's/_collated-GOINPUT/_collated_GOINPUT/g' >> runGOenrichment.sh


nano runGOenrichment.sh

#!/bin/bash

cd $PBS_O_WORKDIR;
ml gcc
ml zlib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib2/

#create variables for each GO Input file
GOINPUT_ab1_3UTR_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab1_3UTR_collated-GOINPUT)
GOINPUT_ab1_3UTRmiRNA_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab1_3UTRmiRNA_collated-GOINPUT)
GOINPUT_ab1_aCNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab1_aCNE_collated-GOINPUT)
GOINPUT_ab1_CNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab1_CNE_collated-GOINPUT)
GOINPUT_ab1_PromNosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab1_PromNosw_collated-GOINPUT)
GOINPUT_ab1_Promsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab1_Promsw_collated-GOINPUT)
GOINPUT_ab1_PromTFBS_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab1_PromTFBS_collated-GOINPUT)
GOINPUT_ab1_promTFBSnosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab1_promTFBSnosw_collated-GOINPUT)
GOINPUT_ab1_promTFBSsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab1_promTFBSsw_collated-GOINPUT)
GOINPUT_mz11_3UTR_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz11_3UTR_collated-GOINPUT)
GOINPUT_mz11_3UTRmiRNA_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz11_3UTRmiRNA_collated-GOINPUT)
GOINPUT_mz11_aCNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz11_aCNE_collated-GOINPUT)
GOINPUT_mz11_CNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz11_CNE_collated-GOINPUT)
GOINPUT_mz11_PromNosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz11_PromNosw_collated-GOINPUT)
GOINPUT_mz11_Promsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz11_Promsw_collated-GOINPUT)
GOINPUT_mz11_PromTFBS_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz11_PromTFBS_collated-GOINPUT)
GOINPUT_mz11_promTFBSnosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz11_promTFBSnosw_collated-GOINPUT)
GOINPUT_mz11_promTFBSsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz11_promTFBSsw_collated-GOINPUT)
GOINPUT_nb1_3UTR_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb1_3UTR_collated-GOINPUT)
GOINPUT_nb1_3UTRmiRNA_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb1_3UTRmiRNA_collated-GOINPUT)
GOINPUT_nb1_aCNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb1_aCNE_collated-GOINPUT)
GOINPUT_nb1_CNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb1_CNE_collated-GOINPUT)
GOINPUT_nb1_PromNosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb1_PromNosw_collated-GOINPUT)
GOINPUT_nb1_Promsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb1_Promsw_collated-GOINPUT)
GOINPUT_nb1_PromTFBS_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb1_PromTFBS_collated-GOINPUT)
GOINPUT_nb1_promTFBSnosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb1_promTFBSnosw_collated-GOINPUT)
GOINPUT_nb1_promTFBSsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb1_promTFBSsw_collated-GOINPUT)
GOINPUT_on11_3UTR_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on11_3UTR_collated-GOINPUT)
GOINPUT_on11_3UTRmiRNA_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on11_3UTRmiRNA_collated-GOINPUT)
GOINPUT_on11_aCNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on11_aCNE_collated-GOINPUT)
GOINPUT_on11_CNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on11_CNE_collated-GOINPUT)
GOINPUT_on11_PromNosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on11_PromNosw_collated-GOINPUT)
GOINPUT_on11_Promsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on11_Promsw_collated-GOINPUT)
GOINPUT_on11_PromTFBS_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on11_PromTFBS_collated-GOINPUT)
GOINPUT_on11_promTFBSnosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on11_promTFBSnosw_collated-GOINPUT)
GOINPUT_on11_promTFBSsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on11_promTFBSsw_collated-GOINPUT)
GOINPUT_pn1_3UTR_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn1_3UTR_collated-GOINPUT)
GOINPUT_pn1_3UTRmiRNA_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn1_3UTRmiRNA_collated-GOINPUT)
GOINPUT_pn1_aCNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn1_aCNE_collated-GOINPUT)
GOINPUT_pn1_CNE_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn1_CNE_collated-GOINPUT)
GOINPUT_pn1_PromNosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn1_PromNosw_collated-GOINPUT)
GOINPUT_pn1_Promsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn1_Promsw_collated-GOINPUT)
GOINPUT_pn1_PromTFBS_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn1_PromTFBS_collated-GOINPUT)
GOINPUT_pn1_promTFBSnosw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn1_promTFBSnosw_collated-GOINPUT)
GOINPUT_pn1_promTFBSsw_collated_GOINPUT=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn1_promTFBSsw_collated-GOINPUT)

BG_abgene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/ab.gene-genelist)
BG_mzgene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/mz.gene-genelist)
BG_nbgene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/nb.gene-genelist)
BG_ongene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/on.gene-genelist)
BG_pngene_genelist=(/tgac/workarea/group-vh/cichlids/substitutions_VCFs/pn.gene-genelist)

mkdir GOOUTPUT

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_3UTR_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_3UTR_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_3UTRmiRNA_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_3UTRmiRNA_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_aCNE_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_aCNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_CNE_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_CNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_PromNosw_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_PromNosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_Promsw_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_Promsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_PromTFBS_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_promTFBSnosw_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_ab1_promTFBSsw_collated_GOINPUT $BG_abgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Asbu_ 1 GOOUTPUT/"$(basename "${GOINPUT_ab1_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_3UTR_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_3UTR_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_3UTRmiRNA_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_3UTRmiRNA_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_aCNE_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_aCNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_CNE_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_CNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_PromNosw_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_PromNosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_Promsw_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_Promsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_PromTFBS_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_promTFBSnosw_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_mz11_promTFBSsw_collated_GOINPUT $BG_mzgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Meze_ 1 GOOUTPUT/"$(basename "${GOINPUT_mz11_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_3UTR_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_3UTR_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_3UTRmiRNA_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_3UTRmiRNA_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_aCNE_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_aCNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_CNE_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_CNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_PromNosw_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_PromNosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_Promsw_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_Promsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_PromTFBS_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_promTFBSnosw_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_nb1_promTFBSsw_collated_GOINPUT $BG_nbgene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Nebr_ 1 GOOUTPUT/"$(basename "${GOINPUT_nb1_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_3UTR_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_3UTR_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_3UTRmiRNA_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_3UTRmiRNA_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_aCNE_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_aCNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_CNE_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_CNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_PromNosw_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_PromNosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_Promsw_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_Promsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_PromTFBS_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_promTFBSnosw_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_on11_promTFBSsw_collated_GOINPUT $BG_ongene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Orni_ 1 GOOUTPUT/"$(basename "${GOINPUT_on11_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_3UTR_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_3UTR_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_3UTRmiRNA_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_3UTRmiRNA_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_aCNE_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_aCNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_CNE_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_CNE_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_PromNosw_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_PromNosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_Promsw_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_Promsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_PromTFBS_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_PromTFBS_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_promTFBSnosw_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_promTFBSnosw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer $GOINPUT_pn1_promTFBSsw_collated_GOINPUT $BG_pngene_genelist /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/Puny_ 1 GOOUTPUT/"$(basename "${GOINPUT_pn1_promTFBSsw_collated_GOINPUT}" -GOINPUT)-GOOUTPUT" persg

# Run on uv2k2
qsub -q Test -l select=1:mem=60GB:ncpus=1 runGOenrichment.sh #{DONE}

# Select those with a q-value (adjusted p-value) of < 0.05 and sort ascending

for i in GOOUTPUT/* ; do awk -F'\t' '$4 < 0.05' $i | sort -k4,4 -n -t $'\t' > GOOUTPUT/"$(basename "$i" .txt)_filtered.txt" ; done


# for plotting the GO (just doing based on log10 fold enrichment (FDR < 0.05), sorting based on p-value here) prepare the OUTPUT file so that matches will show the GO output on each line for each gene
for i in GOOUTPUT/*details_filtered.txt ; do
 sed 's/$/;/g' $i |
 awk '{ gsub(";", ";"$1";"$2";"$3";"$4";"$5";"$6";"$7";"$8";"$9"\n") } 1' |
 cut -f10 |
 awk '{ gsub(";", "\t") } 1' |
 sed '/^$/d' > GOOUTPUT/"$(basename "$i" )2" ;
done

# remove first column then only display unique entries as final column displays no. of genes assigned to the term if you need it!
for i in GOOUTPUT/*details_filtered.txt2 ; do
 cut -f2-10 $i |
 sort -u | sort -k3,3 > GOOUTPUT/"$(basename "$i" 2)3" ;
done

# then combine same genic region files
cat GOOUTPUT/*_3UTR_collated-GOOUTPUT_details_filtered.txt3 > GOOUTPUT/3UTR_collated-GOOUTPUT_details_filtered.txt3
cat GOOUTPUT/*_3UTRmiRNA_collated-GOOUTPUT_details_filtered.txt3 > GOOUTPUT/3UTRmiRNA_collated-GOOUTPUT_details_filtered.txt3
cat GOOUTPUT/*_aCNE_collated-GOOUTPUT_details_filtered.txt3 > GOOUTPUT/aCNE_collated-GOOUTPUT_details_filtered.txt3
cat GOOUTPUT/*_CNE_collated-GOOUTPUT_details_filtered.txt3 > GOOUTPUT/CNE_collated-GOOUTPUT_details_filtered.txt3
cat GOOUTPUT/*_PromNosw_collated-GOOUTPUT_details_filtered.txt3 > GOOUTPUT/PromNosw_collated-GOOUTPUT_details_filtered.txt3
cat GOOUTPUT/*_Promsw_collated-GOOUTPUT_details_filtered.txt3 > GOOUTPUT/Promsw_collated-GOOUTPUT_details_filtered.txt3

# Copied all cat files above to local here: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis
### then plotting /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Multiz_PairwiseVCF_GenomicReg_GOenrichmentplots.R # Uses files created below


##########################################################################################################################

########## 4. Analyse TFBS variants associated with candidate genes




# first, create a new file that retains the motif position within the predicted promoter with all the original details
# col4 and 5 will be motif position within the predicted promoter

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap

### Since i filtered the Candidates to create a new file > Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt2
# check how many TFBS SNPs overlap all candidate genes in M. zebra and P. nyererei
cut -f2 Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt2 | grep -v NULL | xargs -i grep -wiF {} mz-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed | wc -l # 60537 - M. zebra
cut -f3 Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt2 | grep -v NULL | xargs -i grep -wiF {} pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.bed | wc -l # 61079 - P. nyererei


nano candSNP.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

while read -u 3 -r file1 && read -u 4 -r file2
do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' ${file1} ${file2} |
 grep -v 'REMOVEME' |
 awk '{if($21 == "+")print $18,$19+$8,$19+$9,$8,$9,$10,$1,$2,$3,$4,$21,$7,$5,$6,$11,$12,$13,$14,$15,$16;else print $18,$20-$9,$20-$8,$8,$9,$10,$1,$2,$3,$4,$21,$7,$5,$6,$11,$12,$13,$14,$15,$16;}' OFS='\t' > "$(basename "${file2}" .txt)_motifpromANDgenomeannot.txt"
done 3<PROMbed 4<FIMOres

# Pull out SWS1 TFBSs
for i in mz.gene.s102.69 ab.gene.s2279.1 nb.gene.s1.386 on.gene.LG17.282 ; do grep -wiF $i *motifenr_mergedmap2d_motifpromANDgenomeannot.txt >> SWS1_TFBSs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done

# Pull out SWS2a TFBSs
for i in pn.gene.s177.2 ab.gene.s9.15 on.gene.LG5.187 ; do grep -wiF $i *motifenr_mergedmap2d_motifpromANDgenomeannot.txt >> SWS2a_TFBSs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done

# Pull out Rh2 (opn1mw1) TFBSs
for i in mz.gene.s18.149 pn.gene.s11.51 ab.gene.s86.69 nb.gene.s4.225 on.gene.LG5.742 ; do grep -wiF $i *motifenr_mergedmap2d_motifpromANDgenomeannot.txt >> Rh2_TFBSs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done

# Pull out Rho TFBSs
for i in mz.gene.s12.95 ab.gene.s134.8 on.gene.LG20.398 ; do grep -wiF $i *motifenr_mergedmap2d_motifpromANDgenomeannot.txt >> Rho_TFBSs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done

# # get promoter sequence lengths
# grep -A 999999 mz.gene.s102.69 /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta | awk 'NR>1 && /^>/{exit} 1' | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
# grep -A 999999 ab.gene.s2279.1 /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta | awk 'NR>1 && /^>/{exit} 1' | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
# grep -A 999999 nb.gene.s1.386 /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta | awk 'NR>1 && /^>/{exit} 1' | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
#
# ## SWS1 rxra/b/g; nr2f6 and nr2c2 binding sites in Nb and Mz
# # N. brichardi
# scaffold_1:19181348-19181361 # rxra/b/g/nr2f6
# scaffold_1:19181348-19181362 # nr2c2
# scaffold_1:19181359 # genotype: C/C - polymorphic site as compared to Mz
#
# # M. zebra
# scaffold_102:2340611-2340624 # lack of binding sites prediction mapped to N. brichardi
# scaffold_102:2340622 # genotype: T|T - polymorphism to Nb
# # P. nyererei: scaffold_56:1,598,722
# # A. burtoni: scaffold_2279:1,085
# # O. niloticus (on11): LG17:11,755,299

# Pull out SWS1 TFBSs SNPs - as per wavelength palette, compare Mz and Nb here and SNPs between them
for i in mz.gene.s102.69 ab.gene.s2279.1 nb.gene.s1.386 on.gene.LG17.282 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> SWS1_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done
grep 'mz11_nb1_' SWS1_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt
grep 'nb1_mz11_' SWS1_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt

# Pull out SWS2a TFBSs SNPs - as per wavelength palette, compare On, Ab and Pn here and SNPs between them
for i in pn.gene.s177.2 ab.gene.s9.15 on.gene.LG5.187 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> SWS2a_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done

# Pull out Rh2 (opn1mw1) TFBSs SNPs - as per wavelength palette, can compare all
for i in mz.gene.s18.149 pn.gene.s11.51 ab.gene.s86.69 nb.gene.s4.225 on.gene.LG5.742 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Rh2_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done

# Pull out Rho TFBSs SNPs
for i in mz.gene.s12.95 ab.gene.s134.8 on.gene.LG20.398 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Rho_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done

## Also, use a list of candidate genes to pull out a comprehensive set of SNPs overlapping TFBSs in candidates (including the opsins) - this is to share with Ole and Walter to check for segregating SNPs in other species radiations
# irx1b - forebrain
for i in mz.gene.s92.14	pn.gene.s158.15	ab.gene.s221.8	nb.gene.s167.12	on.gene.LG22.514 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt > Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done
# cntn4 - neural system and plasticity
for i in mz.gene.s36.117	pn.gene.s48.53	ab.gene.s607.12	nb.gene.s8.206	on.gene.LG20.843 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done
# foxp2 - TF; developing neural, gastrointestinal and cardiovascular tissues
for i in mz.gene.s69.29	pn.gene.s13.57	ab.gene.s10.20	nb.gene.s1.62	on.gene.LG17.648 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done
# bmpr1a - morphogen important for neural development, forebrain specification, left-right asymmetry
for i in mz.gene.s378.5	pn.gene.s132.19	ab.gene.s467.18	nb.gene.s172.33	on.gene.UNK31.47 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done
# gdf10b - brain expression, axonal outgrowth
for i in mz.gene.s51.93	pn.gene.s4.110	ab.gene.s237.13	nb.gene.s185.10	on.gene.LG8-24.192 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done
## Opsins (sws1, sws2, rh2 and rho) - module variation and under selection
for i in mz.gene.s102.69 ab.gene.s2279.1 nb.gene.s1.386 on.gene.LG17.282 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done #sws1
for i in pn.gene.s177.2 ab.gene.s9.15 on.gene.LG5.187 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done #sws2a
for i in mz.gene.s18.149 pn.gene.s11.51 ab.gene.s86.69 nb.gene.s4.225 on.gene.LG5.742 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done #rh2
for i in mz.gene.s12.95 ab.gene.s134.8 on.gene.LG20.398 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done #rho

# These are from Hahn et al. 2017
# •	visual perception and signal transduction
for i in mz.gene.s45.9 pn.gene.s284.14 ab.gene.s343.2 nb.gene.s3.105 on.gene.LG7.996 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # actr1
for i in mz.gene.s197.7 pn.gene.s237.11 nb.gene.s346.2 on.gene.LG12.944 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # ppiaa
for i in mz.gene.s39.32 pn.gene.s28.25 ab.gene.s688.3 nb.gene.s9.75 on.gene.LG19.762 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # sgpp1

# •	craniofacial and eye morphogenesis
for i in mz.gene.s45.12 pn.gene.s284.11 ab.gene.s343.5 nb.gene.s3.107 on.gene.LG7.998 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # fgfr1
for i in mz.gene.s72.74 pn.gene.s68.37 ab.gene.s30.67 nb.gene.s19.136 on.gene.LG8-24.519 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # zmiz1
for i in mz.gene.s111.51 on.gene.LG1.693 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # foxb1
for i in mz.gene.s193.28 pn.gene.s153.40 nb.gene.s422.3 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # rpgrip1l
for i in mz.gene.s162.49 pn.gene.s323.16 ab.gene.s330.22 on.gene.LG19.581 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # pxdn
for i in mz.gene.s15.43 pn.gene.s33.20 ab.gene.s20.33 nb.gene.s4.127 on.gene.LG5.615 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # draxin
for i in mz.gene.s81.8 pn.gene.s4.96 ab.gene.s259.39 nb.gene.s76.8 on.gene.LG8-24.209 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # foxj1b
for i in mz.gene.s79.13 pn.gene.s71.5 ab.gene.s59.5 on.gene.LG17.3 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # sox2
for i in mz.gene.s36.169 pn.gene.s48.4 ab.gene.s7.142 nb.gene.s863.2 on.gene.LG20.789 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # alx3
for i in pn.gene.s158.9 ab.gene.s221.1 nb.gene.s167.7 on.gene.LG22.505 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # sox4a
for i in mz.gene.s197.6 pn.gene.s237.10 nb.gene.s346.3 ; do grep -wiF $i *_substitutions_1217.bed_*-motifenr_mergedmap2d_motifgenomeannot.txt >> Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done # zmiz2

# Separate the candidate TFBSs SNPs according to species
for i in mz pn ab nb on ; do grep $i.gene Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt > $i-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; done

# Filter for repetitive TFBSs
for i in *-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.txt ; do awk '$20!~"ZNF|RREB1|KLF9|CTCF"' $i > "$(basename "$i" .txt).znf-rreb1-klf9-ctcfFILT.txt" ; done

# Column headers in the files are:
# 1.	PairwiseSNP_scaffold-chr
# 2.	SNP_start
# 3.	SNP_end
# 4.	ID
# 5.	Reference_base
# 6.	Alternate_allele
# 7.	Quality
# 8.	Filter
# 9.	Info
# 10.	Alternate_allele_scaff-chr
# 11.	Alternate_allele_SNP_Pos
# 12.	TFBS_scaff-chr
# 13.	TFBS_start
# 14.	TFBS_end
# 15.	TFBS_strand
# 16.	geneID
# 17.	gene_symbol
# 18.	gene_strand
# 19.	motifID
# 20.	motif_symbol
# 21.	TFBS_score
# 22.	TFBS_pval
# 23.	TFBS_qval
# 24.	TFBS_sequence

qsub -q Test -l select=1:mem=50GB:ncpus=1 candSNP.sh

## change the above file into regions of interest to share with Ole instead
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap
cut -f1-3 pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT.txt | # cut the required columns of chr, start, end
awk -F':' '{print $2,$3}' OFS='\t' | # separate out of the chr location
sort -k1,1 -k2,2n | # sort the scaffolds start and end
awk -F'\t' 'NF>1{a[$1] = a[$1]$2"\t"$3"\t"};END{for(i in a)print i"\t"a[i]}' | # if the first column string is the same, then join every second and third column entry into one row, retaining first column entry, then tab then columns of the second and third line entries
awk '{print $1,$2,$NF}' OFS='\t' | # retain the chr, col2 (start of region in scaffold with SNPs), last column (end of region in scaffold with SNPs)
sort -V -k1,1 > pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT_RegionsOfInterest.txt # then sort and output file
cp pn-Candidates_TFBSsSNPs_motifenr_mergedmap2d_motifpromANDgenomeannot.znf-rreb1-klf9-ctcfFILT_RegionsOfInterest.txt /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap # copy to workarea to share

### The above files have been shared with Ole and Walter to check for genotypes in other radiation species


## Can load into IGV to view
# Nb - load against ATAC peaks
N_brichardi_v1.assembly.fasta # Genome
Neolamprologus_brichardi.BROADNB2.sorted.gtf # annotation
aCNEs_Nbrichardi-Onliftover_true_aCNEs.bed # aCNEs
CNE.528k_Nbrichardi-Onliftover.bed #CNEs
awk '{print $1,$2,$3,$5"-"$6,$7,$4,$8,$9,$10,$11}' OFS='\t' nb1_mz11_substitutions_1217.bed > nb1_mz11_substitutions_1217-IGV.bed # Nb-Mz pairwise variants
awk '{print $1,$2,$3,$10,"0",$4,$5,$6,$7,$8,$9,$11,$12,$13,$14}' OFS='\t' N_brichardi_1e-4_JASPAR2018-allfimo_qval0.05_motifpromANDgenomeannot.txt.5 > N_brichardi_1e-4_JASPAR2018-allfimo_qval0.05_motifpromANDgenomeannot.txt.5-IGV.bed # all TFBSs
# Ab
# Mz
# Pn
P_nyererei_v1.assembly.fasta # Genome
Pundamilia_nyererei.BROADPN2.sorted.gtf # annotation
aCNEs_Pnyererei-Onliftover_true_aCNEs.bed # aCNEs
CNE.528k_Pnyererei-Onliftover.bed #CNEs
awk '{print $1,$2,$3,$5"-"$6,$7,$4,$8,$9,$10,$11}' OFS='\t' pn1_ab1_substitutions_1217.bed > pn1_ab1_substitutions_1217-IGV.bed # Pn-Ab pairwise variants > for sws2a comparisons
awk '{print $1,$2,$3,$5"-"$6,$7,$4,$8,$9,$10,$11}' OFS='\t' pn1_nb1_substitutions_1217.bed > pn1_nb1_substitutions_1217-IGV.bed # Pn-Nb pairwise variants
awk '{print $1,$2,$3,$10,"0",$4,$5,$6,$7,$8,$9,$11,$12,$13,$14}' OFS='\t' P_nyererei_1e-4_JASPAR2018-allfimo_qval0.05_motifpromANDgenomeannot.txt.5 > P_nyererei_1e-4_JASPAR2018-allfimo_qval0.05_motifpromANDgenomeannot.txt.5-IGV.bed # all TFBSs

# On


## Ole and Joana Meier shared a VCF file using an anchored unpublished P. nyererei assembly

## A. USE THESE - this is the vcf they provided for the new scans
# zipped vcf found here: /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/tarangRegions.secondSet.vcf.gz
cd /tgac/scratch/mehtat/Cichlid_GRNs
mkdir functional_variant_analysis/
cd functional_variant_analysis/
cp -r /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/ .
cd /tgac/scratch/mehtat/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/2.fullTFBSNPs_overlap_May2018/
# they provided a liftover of coordinates to the old assembly: secondSet_RegionsOfInterest_forTarang_withNewPos.txt
# need to amend the vcf file so that the coordinates match the old assembly
# Prepped a python script to do this - /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/Pnv1.1_vcf_liftover.py


# ### B. this is for the old TFBS scans (using JASPAR)
# # they provided a liftover of coordinates to the old assembly
# # need to amend the vcf file so that the coordinates match the old assembly and can be loaded into IGV
# cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018 # local
# awk '{OFS="\t"} {if ($1=="chr15" && $2<6163000) $2=$2-3872491; print $0}' tarangRegions.vcf |
# awk '{OFS="\t"} {if ($1=="chr15" && $2>2260000) $1="scaffold_4"; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr13" && $2<9299000) $2=$2-8449235; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr13" && $2<850000) $1="scaffold_33"; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr17") $2=$2-28204324; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr17") $1="scaffold_48"; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr15" && $2>13837000) $2=$2-12707138; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr15" && $2<1157000) $1="scaffold_68"; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr12") $2=$2-2809604; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr12") $1="scaffold_71"; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr2" && $2>58885000) $2=$2-57594939; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr2" && $2>1290000) $1="scaffold_153"; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr3") $2=$2-2099411; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr3") $1="scaffold_158"; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr13" && $2>26335000) $2=$2-26306570; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr13" && $2<54000) $1="scaffold_177"; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr2" && $2<26026000) $2=$2-25753991; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr2" && $2<272000) $1="scaffold_284"; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr8") $2=$2-19537502; print $0}' |
# awk '{OFS="\t"} {if ($1=="chr8") $1="scaffold_323"; print $0}' > tarangRegions_Pn-v1liftover.vcf
# # then sorted the above file to load into IGV
# # Candidate TGs to view
# # scaffold_153	pn.gene.s153.40	NOX5
# # scaffold_158	pn.gene.s158.9	SOX4
# # scaffold_177	pn.gene.s177.2	sws2
# # scaffold_237	pn.gene.s237.10	ZMIZ2
# # scaffold_237	pn.gene.s237.11	ppiaa
# # scaffold_284	pn.gene.s284.11	FGFR1
# # scaffold_323	pn.gene.s323.16	PXDNL
# # scaffold_33	pn.gene.s33.20	draxin
# # scaffold_48	pn.gene.s48.4	ALX3
# # scaffold_4	pn.gene.s4.96	foxj1b
# # scaffold_68	pn.gene.s68.37	ZMIZ1
# # scaffold_71	pn.gene.s71.5	SOX2
#
#
# ## To run IGV with higher memory
# java -Xmx4000m -jar /Applications/IGV_2.3.91.app/Contents/Java/igv.jar

##########################################################################################################################

########## 5. Plot pairwise SNPs against Euclidian expression divergence {DONE}

# for each region, create a file for each speces; col1=geneID, col2=pairwise_species, col3=refNT-variantSNP, col4=SNPchr, col5=SNPstart, col6=SNPend, col7=aCNE/CNE/geneID > col1=OGID, col2=geneID, col3=pairwise_species, col4=refNT-variantSNP, col5=SNPchr, col6=SNPstart, col7=SNPend, col8=aCNE/CNE/geneID

cd /tgac/workarea/group-vh/cichlids/substitutions_VCFs
mkdir SNP_ExprDiv

nano prepSNP_ExprDiv.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
ml bedtools/2.25.0
ml GCC
ml zlib

# CNEs {DONE}
# run bedtools closest to get the closest promoter and use that as gene ID
for i in mz11*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done
for i in pn1*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done
for i in ab1*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done
for i in nb1*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done
for i in on11*_substitutions_1217.bed_substCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done

# aCNEs {DONE}
# run bedtools closest to get the closest promoter and use that as gene ID for GO enrichment
for i in mz11*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done
for i in pn1*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done
for i in ab1*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done
for i in nb1*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done
for i in on11*_substitutions_1217.bed_substaCNEoverlap.txt ; do bedtools closest -d -a $i -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk -v i=$i '{print $22,i,$5"-"$6,$1,$2,$3,$15}' OFS='\t' | sed 's/_substitutions_1217.bed_substaCNEoverlap.txt//g' | awk '$1!="."' > SNP_ExprDiv/$i.ED1 ; done

# Promoters
# switched {DONE}
for i in mz11*_subst-swPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in pn1*_subst-swPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in ab1*_subst-swPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in nb1*_subst-swPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in on11*_subst-swPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-swPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done

# no switch {DONE}
for i in mz11*_subst-NoswPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in pn1*_subst-NoswPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in ab1*_subst-NoswPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in nb1*_subst-NoswPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in on11*_subst-NoswPromoverlap.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-NoswPromoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done

# 3'UTR {DONE}

for i in mz11*_subst-3UTRoverlap.txt ; do awk -v i=$i '{print i,$5"-"$6,$1,$2,$3}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1a ; done
for i in mz11*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' > SNP_ExprDiv/$i.ED1b ; done

for i in pn1*_subst-3UTRoverlap.txt ; do awk -v i=$i '{print i,$5"-"$6,$1,$2,$3}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1a ; done
for i in pn1*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' > SNP_ExprDiv/$i.ED1b ; done

for i in ab1*_subst-3UTRoverlap.txt ; do awk -v i=$i '{print i,$5"-"$6,$1,$2,$3}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1a ; done
for i in ab1*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' > SNP_ExprDiv/$i.ED1b ; done

for i in nb1*_subst-3UTRoverlap.txt ; do awk -v i=$i '{print i,$5"-"$6,$1,$2,$3}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1a ; done
for i in nb1*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' > SNP_ExprDiv/$i.ED1b ; done

for i in on11*_subst-3UTRoverlap.txt ; do awk -v i=$i '{print i,$5"-"$6,$1,$2,$3}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1a ; done
for i in on11*_subst-3UTRoverlap.txt ; do cut -f15 $i | awk -F';' '{print $2}'| sed 's/Parent=//g' | sed 's/mrna/gene/g' | awk -F'.' '{print $1,$2,$3,$4}' OFS='.' > SNP_ExprDiv/$i.ED1b ; done

ls -1 SNP_ExprDiv/*_subst-3UTRoverlap.txt.ED1a > SNP_ExprDiv/3UTRED1a
ls -1 SNP_ExprDiv/*_subst-3UTRoverlap.txt.ED1b > SNP_ExprDiv/3UTRED1b
while read -u 3 -r file1 && read -u 4 -r file2
do
 paste -d'\t' ${file2} ${file1} ${file2} > SNP_ExprDiv/"$(basename "${file1}" .ED1a).ED1"
done 3<SNP_ExprDiv/3UTRED1a 4<SNP_ExprDiv/3UTRED1b
for i in SNP_ExprDiv/*_subst-3UTRoverlap.txt.ED1 ; do sed -i 's/bn.gene/nb.gene/g' $i ; done
rm SNP_ExprDiv/*.ED1a
rm SNP_ExprDiv/*.ED1b
rm SNP_ExprDiv/3UTRED1a
rm SNP_ExprDiv/3UTRED1b

# 3' UTR miRNA binding sites {DONE}
for i in mz11*_subst_miRNA-3UTRoverlap.txt ; do awk -v i=$i '{print $5,i,$10"-"$11,$1,$2,$3,$5}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in pn1*_subst_miRNA-3UTRoverlap.txt ; do awk -v i=$i '{print $5,i,$10"-"$11,$1,$2,$3,$5}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in ab1*_subst_miRNA-3UTRoverlap.txt ; do awk -v i=$i '{print $5,i,$10"-"$11,$1,$2,$3,$5}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in nb1*_subst_miRNA-3UTRoverlap.txt ; do awk -v i=$i '{print $5,i,$10"-"$11,$1,$2,$3,$5}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done
for i in on11*_subst_miRNA-3UTRoverlap.txt ; do awk -v i=$i '{print $5,i,$10"-"$11,$1,$2,$3,$5}' OFS='\t' $i | sed 's/_substitutions_1217.bed_subst-3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt//g' > SNP_ExprDiv/$i.ED1 ; done

# Promoters TFBSs {DONE}
for i in /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/SNP_overlap/*_substitutions_1217*_motifgenomeannot.txt ; do ln -s $i ; done # create symbolic links to files created another folder
for i in mz11*_substitutions_1217*_motifgenomeannot.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_mz-motifenr_mergedmap2d_motifgenomeannot.txt//g' > SNP_ExprDiv/$i.PromTFBS.ED1 ; done
for i in pn1*_substitutions_1217*_motifgenomeannot.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_pn-motifenr_mergedmap2d_motifgenomeannot.txt//g' > SNP_ExprDiv/$i.PromTFBS.ED1 ; done
for i in ab1*_substitutions_1217*_motifgenomeannot.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_ab-motifenr_mergedmap2d_motifgenomeannot.txt//g' > SNP_ExprDiv/$i.PromTFBS.ED1 ; done
for i in nb1*_substitutions_1217*_motifgenomeannot.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_nb-motifenr_mergedmap2d_motifgenomeannot.txt//g' > SNP_ExprDiv/$i.PromTFBS.ED1 ; done
for i in on11*_substitutions_1217*_motifgenomeannot.txt ; do awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' $i | sed 's/_substitutions_1217.bed_on-motifenr_mergedmap2d_motifgenomeannot.txt//g' > SNP_ExprDiv/$i.PromTFBS.ED1 ; done

# Promoters_TFBSs_switched {DONE}
for i in mz11*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Mzvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_M_zebra_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSsw.ED1
done

for i in pn1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Pnvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_P_nyererei_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSsw.ED1
done

for i in ab1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Abvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_A_burtoni_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSsw.ED1
done

for i in nb1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Nbvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_N_brichardi_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSsw.ED1
done

for i in on11*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Onvsany_switch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_O_niloticus_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSsw.ED1
done


# Promoters_TFBSs_noswitch {DONE}

for i in mz11*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Mzvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_M_zebra_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSnosw.ED1
done

for i in pn1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Pnvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_P_nyererei_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSnosw.ED1
done

for i in ab1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Abvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_A_burtoni_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSnosw.ED1
done

for i in nb1*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Nbvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_N_brichardi_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSnosw.ED1
done

for i in on11*_substitutions_1217.bed_*_motifgenomeannot.txt ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$5;};NR>FNR{if($16 in a){print $0,a[$16];}else{print $0,"NA";}}' Onvsany_noswitch.genes.5kb_promoters.stranded.sorted.bed $i | grep -v NA | awk -v i=$i '{print $16,i,$5"-"$6,$1,$2,$3,$16}' OFS='\t' |
 sed 's/_substitutions_1217.bed_O_niloticus_1e-4_JASPAR2018-allfimo_qval0.05_motifgenomeannot.txt.6//g' > SNP_ExprDiv/$i.PromTFBSnosw.ED1
done


# Map all geneIDs to OGIDS > col1=OGID, col2=geneID, col3=pairwise_species, col4=refNT-variantSNP, col5=SNPchr, col6=SNPstart, col7=SNPend, col8=aCNE/CNE/geneID

for i in SNP_ExprDiv/mz11*ED1 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz $i | awk '{print $8,$1,$2,$3,$4,$5,$6,$7}' OFS='\t' | grep -v NA > SNP_ExprDiv/"$(basename "$i" .ED1).ED2"
done

for i in SNP_ExprDiv/pn1*ED1 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn $i | awk '{print $8,$1,$2,$3,$4,$5,$6,$7}' OFS='\t' | grep -v NA > SNP_ExprDiv/"$(basename "$i" .ED1).ED2"
done

for i in SNP_ExprDiv/ab1*ED1 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab $i | awk '{print $8,$1,$2,$3,$4,$5,$6,$7}' OFS='\t' | grep -v NA > SNP_ExprDiv/"$(basename "$i" .ED1).ED2"
done

for i in SNP_ExprDiv/nb1*ED1 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb $i | awk '{print $8,$1,$2,$3,$4,$5,$6,$7}' OFS='\t' | grep -v NA > SNP_ExprDiv/"$(basename "$i" .ED1).ED2"
done

for i in SNP_ExprDiv/on11*ED1 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on $i | awk '{print $8,$1,$2,$3,$4,$5,$6,$7}' OFS='\t' | grep -v NA > SNP_ExprDiv/"$(basename "$i" .ED1).ED2"
done

rm SNP_ExprDiv/*.ED1

# Map OGIDs to pairwise expression divergence for each tissue and evol rate
# copied euclidian files ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/*_Euclidian_evolRate_1to1OGIDs.txt to /tgac/workarea/group-vh/cichlids/substitutions_VCFs/SNP_ExprDiv

# add colhead for OGIDs column
cd /tgac/workarea/group-vh/cichlids/substitutions_VCFs/SNP_ExprDiv
for i in *_Euclidian_evolRate_1to1OGIDs.txt ; do sed -i 's/Mz/OGID\tMz/' $i ; done

# 9:OGID	10:Mz-Pn	11:Mz-Ab	12:Mz-Nb	13:Mz-On	14:Pn-Ab	15:Pn-Nb	16:Pn-On	17:Ab-Nb	18:Ab-On	19:Nb-On	20:ED.rowsum	21:4fold	22:prom	23:Ga_gene	24:Dr_gene	25:Hs_gene



for i in Br Ey Ht Kd Ms Ts ; do
 for j in mz11_ab1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,11,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in mz11_pn1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,10,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in mz11_nb1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,12,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in mz11_on11*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,13,21-25 > ${i}_${j}
done
done


for i in Br Ey Ht Kd Ms Ts ; do
 for j in pn1_ab1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,14,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in pn1_mz11*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,10,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in pn1_nb1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,15,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in pn1_on11*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,16,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in ab1_pn1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,14,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in ab1_mz11*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,11,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in ab1_nb1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,17,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in ab1_on11*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,18,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in nb1_pn1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,15,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in nb1_mz11*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,12,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in nb1_ab1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,17,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in nb1_on11*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,19,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in on11_pn1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,16,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in on11_mz11*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,13,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in on11_ab1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,18,21-25 > ${i}_${j}
done
done

for i in Br Ey Ht Kd Ms Ts ; do
 for j in on11_nb1*.ED2 ; do
 awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ${i}_Euclidian_evolRate_1to1OGIDs.txt $j | grep -v NA | cut -f1-8,19,21-25 > ${i}_${j}
done
done


for i in Br Ey Ht Kd Ms Ts ; do
 cat ${i}_*PromTFBS.ED2 > ${i}_PromTFBS.collated.ED2
 cat ${i}_*PromTFBSnosw.ED2 > ${i}_PromTFBSnosw.collated.ED2
 cat ${i}_*PromTFBSsw.ED2 > ${i}_PromTFBSsw.collated.ED2
 cat ${i}_*subst-3UTRoverlap.txt.ED2 > ${i}_3UTR.collated.ED2
 cat ${i}_*3UTRoverlap.txt_subst_miRNA-3UTRoverlap.txt.ED2 > ${i}_3UTRmiRNA.collated.ED2
 cat ${i}_*substaCNEoverlap.txt.ED2 > ${i}_aCNE.collated.ED2
 cat ${i}_*substCNEoverlap.txt.ED2 > ${i}_CNE.collated.ED2
 cat ${i}_*subst-NoswPromoverlap.txt.ED2 > ${i}_Promnosw.collated.ED2
 cat ${i}_*subst-swPromoverlap.txt.ED2 > ${i}_Promsw.collated.ED2
done

# Run on uv2k2
qsub -q Test -l select=1:mem=70GB:ncpus=1 prepSNP_ExprDiv.sh #{DONE}

cd SNP_ExprDiv

nano Collate.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

# Create files that have the counts for SNPs in each gene
# Resulting files *.ED4 have the following columns
# col1=OGID, col2=geneID, col3=pairwise_species, col4=SNPchr, col5=pairwise_ED, col6=4foldevolrate, col7=promevolrate, col8=Ga_symbol, col9=Dr_symbol, col10=Hs_symbol, col11=SNPcountsingene, col12=1/0 for col5>1 and col11>32
for i in *.collated.ED2 ; do
 cut -f1-3,5,9-14 $i | sort | uniq -c | awk -F' ' '{print $1,$2,$0}' OFS='\t' | cut -f1,2,4-12 | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$1}' OFS='\t' > "$(basename "$i" .ED2).ED3"
done

for i in *.collated.ED3 ; do
 awk '{if($5>1 && $11>32){print$0,$12="1";}else{print$0,$12="0";}}' OFS='\t' $i | sort -k5,5 -k11,11 -rn > "$(basename "$i" .ED3).ED4"
done

# tar up the intermediate files
tar -zcvf ED2_ED3-files.tar.gz *.ED2 *.ED3
rm *.ED2
rm *.ED3

qsub -q Test -l select=1:mem=60GB:ncpus=1 Collate.sh

# copied above cat files *.ED4 to local directory /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis
# plotted in ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/SNP_overlap.R
