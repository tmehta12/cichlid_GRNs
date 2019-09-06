#!/bin/sh

## prep Marton's rewiring csv files to look for overlap of candidate genes
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton

for i in /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/Outputs/Result_files/*.csv ; do
  sed $'s/,/\t/g' $i | sed 's/"//g' > $i.tsv
done

cand=('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt')
haplorew=('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/Outputs/Result_files/Haplochromines_Mz_Pn_Ab_node.csv.tsv')
fivesprew=('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/Outputs/Result_files/All_five_species_node_noExpression.csv.tsv')

awk 'BEGIN{OFS="\t"}NR==FNR{a[$7]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $haplorew $cand | grep -v 'REMOVEME' | sort -k21,21 -rn > /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/Outputs/Result_files/Haplochromines_Mz_Pn_Ab_node_CAND.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$7]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $fivesprew $cand | grep -v 'REMOVEME' | sort -k21,21 -rn > /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/Outputs/Result_files/All_five_species_node_noExpression_CAND.tsv
# 73 out of 132 candidate genes are found in 1:1 rewired edges - sorted the above based on rewiring score
# the above files are sorted based on the DyNet Rewiring (Dn-score) (non-corrected) score - ideally you should order based on the correct score (-f20)

# OG1652_0 - rh2 - Green-sensitive opsin - Opsin – associated trait under selection
# OG2764_0 - foxp2 - TF; developing neural, gastrointestinal and cardiovascular tissues - Fast evolving in cichlids (Brawand, D. et al. 2014)
# OG346_0 - cntn4 - Neural system and plasticity - Fast evolving in cichlids (Brawand, D. et al. 2014)
# OG2764_0 - foxp2 - TF; developing neural, gastrointestinal and cardiovascular tissues - Fast evolving in cichlids (Brawand, D. et al. 2014)
# OG2800_0 - bmpr1aa - morphogen important for neural development, forebrain specification, left-right asymmetry - Fast evolving in cichlids (Brawand, D. et al. 2014)
# OG8558_0 - gdf10b - brain expression, axonal outgrowth - Fast evolving in cichlids (Brawand, D. et al. 2014)

## what is the average DN-score (degree corrected) and DN-rewiring score of module switched TGs?
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only
# the following file describes any switches between pairs in 1to1 orthologs of TGs
awk '$5==1 && $6==1 && $7==1 && $8==1 && $9==1' OFS='\t' Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | awk '$21==1 || $23==1 || $25==1 || $27==1 || $29==1 || $31==1 || $33==1 || $35==1 || $37==1 || $39==1' OFS='\t' > Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID-OnNbAbPnMz_AllTGsw.txt
# then on this file, you should pull out the TGs and sort -u
cut -f3 Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID-OnNbAbPnMz_AllTGsw.txt | sort -u > TG_1to1switchedgenes
# copy this file to local:
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton
# do an awk match on the rewiring results files (both 5 species and 5speciesCAND) and work out on average rewiring score
fivesprew=('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/Outputs/Result_files/All_five_species_node_noExpression.csv.tsv')
fivesprewCAND=('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/Outputs/Result_files/All_five_species_node_noExpression_CAND.tsv')
awk 'BEGIN{OFS="\t"}NR==FNR{a[$7]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $fivesprew TG_1to1switchedgenes | grep -v 'REMOVEME' | awk '{ sum += $5; n++ } END { if (n > 0) print sum / n; }' #3458 genes = 11.3894 (DyNet rewiring score)
awk 'BEGIN{OFS="\t"}NR==FNR{a[$7]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $fivesprew TG_1to1switchedgenes | grep -v 'REMOVEME' | awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' #3458 genes = 0.159497 (Dn-score, degree corrected)
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $fivesprewCAND TG_1to1switchedgenes | grep -v 'REMOVEME' | awk '{ sum += $22; n++ } END { if (n > 0) print sum / n; }' # 36 genes = 6.09722 (DyNet rewiring score)
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $fivesprewCAND TG_1to1switchedgenes | grep -v 'REMOVEME' | awk '{ sum += $21; n++ } END { if (n > 0) print sum / n; }' # 36 genes = 0.141111 (Dn-score, degree corrected)
# look at their scores
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $fivesprewCAND TG_1to1switchedgenes | grep -v 'REMOVEME' | cut -f1,11,13,21,22 | sort -k5,5 -rn

### Prep new alternative files for Marton to run DyNet
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/
# 0. only1to1 nodes
/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt #215810
# 1. All TF-TG edges
# remove second row that is redundant
sed -i '2d' Edge_Attributes_Collated4c.coexpr_promONLY.matrix3.txt
cp Edge_Attributes_Collated4c.coexpr_promONLY.matrix3.txt /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/ # 1053165 (not including colheader)
# 2. TF nodes that are 1:1 and all TG nodes (1:1 and non 1:1)
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}{if(a[$2]){print $0,a[$2];}else{print $0,"REMOVEME";}}' 1to1only/OGIDS.1to1.txt5 Edge_Attributes_Collated4c.coexpr_promONLY.matrix3.txt | grep -v 'REMOVEME' > Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneA1to1OGID.txt #533271
# add colheaders
head -1 Edge_Attributes_Collated4c.coexpr_promONLY.matrix3.txt > Edge_Attributes_Collated4c.coexpr_promONLY.matrix3.colheads
cat Edge_Attributes_Collated4c.coexpr_promONLY.matrix3.colheads Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneA1to1OGID.txt > Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneA1to1OGID.txt2
rm Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneA1to1OGID.txt
mv Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneA1to1OGID.txt2 Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneA1to1OGID.txt
cp Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneA1to1OGID.txt /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/ # 533270 (not including colheader)

### Created boxplots to assess DyNet scores from the three different sources - all edges, 1to1TF, all1to1 (then, as per below, also assessed module constrained)
/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/DyNet_scores_boxplots.R

### Using the original files that you shared with Marton, you should constrain those files by module assignment
# geneA and geneB for each species need to be in the same module
# Marton will then use these new files to re-calculate DyNet scores
# We will then use these to compare when module assignment is not used (can compare by the boxplots)
# It is likely that it will make sense to use the scores when constrained by module assignment since these are the networks we display plus highlight expression based rewiring.

cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton

# 0. only1to1 nodes - module constrained in each species
awk -F'\t' 'NR==1 {print};$10==$15 && $11==$16 && $12==$17 && $13==$18 && $14==$19' 1.1to1_files_shared/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt > 1.1to1_files_shared/1_all1to1_ModuleConstrained.txt # 6546/215810 (not incl. colheader)

# 1. All TF-TG edges - module constrained in each species (here, you need to input the conditional of unless == "NULL")
awk -F'\t' 'NR==1 {print};($10==$15 || $10=="NULL" || $15=="NULL") && ($11==$16 || $11=="NULL" || $16=="NULL") && ($12==$17 || $12=="NULL" || $17=="NULL") && ($13==$18 || $13=="NULL" || $18=="NULL") && ($14==$19 || $14=="NULL" || $19=="NULL")' 3.alledges_no1to1/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3.txt > 3.alledges_no1to1/3_allTFTG_ModuleConstrained.txt # 57620/1053165 (not incl. colheader)

# 2. TF nodes that are 1:1 and all TG nodes (1:1 and non 1:1) - module constrained in each species (here, you need to input the conditional of unless == "NULL")
awk -F'\t' 'NR==1 {print};($10==$15 || $10=="NULL" || $15=="NULL") && ($11==$16 || $11=="NULL" || $16=="NULL") && ($12==$17 || $12=="NULL" || $17=="NULL") && ($13==$18 || $13=="NULL" || $18=="NULL") && ($14==$19 || $14=="NULL" || $19=="NULL")' 4.TF1to1_TGallnodes/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneA1to1OGID.txt > 4.TF1to1_TGallnodes/4_TF1to1_allTG_ModuleConstrained.txt # 21723/533270 (not incl. colheader)


####### After several tests, we decided on the following - 09/07/18

# 1. In cases where the network is large and can be linked to a tissue-specific context e.g. neurod1, you can module constrain to find fine-grained tissue-specific regulators.
# 2. In cases where the networks are small e.g. sws1 and rho, do not module constrain and instead, show the network with the related tissue expression pattern e.g. eye.
# 3. In all cases, DyNet needs to be ran using 1to1TF and non1to1TG edges that can be confirmed as not being present in the genome.
  # For each of the candidate genes, define in which species the comparison is made, based on the presence of TG. For example if TG found in Ab and Nb only, then TFs only need to be 1-to-1 in those species for that edge and networks created accordingly.
  # This is the bit where I need to find the best way to prepare files for you. It would be best to still run DyNet for all TGs (with 1to1 TFs) and not just candidate genes so this will require some thought (on how we deal with all the genes that are not found in each species, and run DyNet just on those species comparisons). I’ll have a think about it over the weekend and maybe we can meet on Monday if you’re available.


## A. Create a new edge matrices for the coexpr_prom files that contain more of the interaction information
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix


nano NEWcoexpr_prom_matrices.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

# Pull out the relevant columns, including all additional information for those edges
for i in /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/*-Edge_Attributes_Collated4c.txt ; do
  awk '$18=="co-expression_TF-TG" || $18=="promoter_motif"' $i | cut -f1-3,5,7-9,15,17,18,27-36,49 > "$(basename "$i" .txt).coexpr_promONLY.NEW.txt"
done

# Comparisons for presence/absence

# 0. re-structure the column tables so that edges can be matched
for i in *ONLY.NEW.txt ; do
  awk '{print $3":"$7,$0}' OFS='\t' $i > $i.tmp1
done

# 1. collate all possible node interactions from each species, removing redundant edges
cat *.coexpr_promONLY.NEW.txt.tmp1 | sort -u > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp1

# 2. Then, do the species presence/absence
# Mz - presence/absence of edges
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,"1";}else{print $0,"0";}}' Mz-Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp1 Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp1 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp2

# Pn - presence/absence of edges
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,"1";}else{print $0,"0";}}' Pn-Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp1 Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp2 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp3

# Ab - presence/absence of edges
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,"1";}else{print $0,"0";}}' Ab-Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp1 Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp3 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp4

# Nb - presence/absence of edges
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,"1";}else{print $0,"0";}}' Nb-Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp1 Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp4 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp5

# On - presence/absence of edges
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,"1";}else{print $0,"0";}}' On-Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp1 Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp5 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp6

# 3. create colheaders and add
printf 'interaction\tCichlid_GeneA\tCichlid_GeneA_symbol\tGeneA\tCichlid_GeneB\tCichlid_GeneB_symbolGa\tCichlid_GeneB_symbolDr\tGeneB\tdirectness\tlayer\tsource\tmotif_pattern[pm_cm]\tstart[pm_cm]\tstop[pm_cm]\tstrand[pm_cm]\tscore[pm_cm]\tp.value[pm_cm]\tq.value[pm_cm]\tsequence[pm_cm]\tconf_level[pm_cm]\tconf_score[pm_cm]\tconfidence[coTFTG]\tMz.P-1_A-0\tPn.P-1_A-0\tAb.P-1_A-0\tNb.P-1_A-0\tOn.P-1_A-0\n' > matrix_colheads.NEW
cat matrix_colheads.NEW Edge_Attributes_Collated4c.coexpr_promONLY.NEW.txt.tmp6 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.txt

# remove all temp files
rm *.tmp*

# run the above
qsub -q Test -l select=1:mem=60GB:ncpus=1 NEWcoexpr_prom_matrices.sh # DONE


nano NEWcoexpr_prom_matrices-pt2.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

# Add module assignment

# iaa. create module mapping files - already done in previous scripts

# iab. map to add module assignment
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$4]){print $0,a[$4];}else{print $0,"NULL","NULL";}}' Mz-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.txt | cut -f1-27,29 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$4]){print $0,a[$4];}else{print $0,"NULL","NULL";}}' Pn-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp1 | cut -f1-28,30 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$4]){print $0,a[$4];}else{print $0,"NULL","NULL";}}' Ab-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp2 | cut -f1-29,31 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp3
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$4]){print $0,a[$4];}else{print $0,"NULL","NULL";}}' Nb-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp3 | cut -f1-30,32 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp4
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$4]){print $0,a[$4];}else{print $0,"NULL","NULL";}}' On-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp4 | cut -f1-31,33 | sed '1d' > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp5
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$8]){print $0,a[$8];}else{print $0,"NULL","NULL";}}' Mz-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp5 | cut -f1-32,34 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp6
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$8]){print $0,a[$8];}else{print $0,"NULL","NULL";}}' Pn-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp6 | cut -f1-33,35 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp7
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$8]){print $0,a[$8];}else{print $0,"NULL","NULL";}}' Ab-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp7 | cut -f1-34,36 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp8
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$8]){print $0,a[$8];}else{print $0,"NULL","NULL";}}' Nb-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp8 | cut -f1-35,37 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp9
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$8]){print $0,a[$8];}else{print $0,"NULL","NULL";}}' On-speciesspecnames_clusterassign.OGID.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp9 | cut -f1-36,38 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp10

printf 'interaction\tCichlid_GeneA\tCichlid_GeneA_symbol\tGeneA\tCichlid_GeneB\tCichlid_GeneB_symbolGa\tCichlid_GeneB_symbolDr\tGeneB\tdirectness\tlayer\tsource\tmotif_pattern[pm_cm]\tstart[pm_cm]\tstop[pm_cm]\tstrand[pm_cm]\tscore[pm_cm]\tp.value[pm_cm]\tq.value[pm_cm]\tsequence[pm_cm]\tconf_level[pm_cm]\tconf_score[pm_cm]\tconfidence[coTFTG]\tMz.P-1_A-0\tPn.P-1_A-0\tAb.P-1_A-0\tNb.P-1_A-0\tOn.P-1_A-0\tMz_GeneAModule\tPn_GeneAModule\tAb_GeneAModule\tNb_GeneAModule\tOn_GeneAModule\tMz_GeneBModule\tPn_GeneBModule\tAb_GeneBModule\tNb_GeneBModule\tOn_GeneBModule\n' > matrix_colheads.NEW2
cat matrix_colheads.NEW2 Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp10 > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix2.txt
rm Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.txt # DONE
rm Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix.tmp* # DONE

# iac. check all column numbers for each row
# awk '{print NF}' Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix2.txt | sort -nu | head -n 1 # 37
# awk '{print NF}' Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix2.txt | sort -nu | tail -n 1 # 37

# ib. then compare the GeneA modules per species, and same for GeneB to highlight module switches - add 1 for switch, 0 for noswitch, and "NULL" if no comparison can be made for all pairwise comparisons
# Mz-Pn - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($28=="NULL" || $29=="NULL")print $0,$38="NULL";else if($28!=$29)print $0,$38=1;else print $0,$38=0;}' OFS='\t' $i > $i.1
  awk '{if($33=="NULL" || $34=="NULL")print $0,$39="NULL";else if($33!=$34)print $0,$39=1;else print $0,$39=0;}' OFS='\t' $i.1 > $i.2
done

# Mz-Ab - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($28=="NULL" || $30=="NULL")print $0,$40="NULL";else if($28!=$30)print $0,$40=1;else print $0,$40=0;}' OFS='\t' $i.2 > $i.3
  awk '{if($33=="NULL" || $35=="NULL")print $0,$41="NULL";else if($33!=$35)print $0,$41=1;else print $0,$41=0;}' OFS='\t' $i.3 > $i.4
done
# Mz-Nb - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($28=="NULL" || $31=="NULL")print $0,$42="NULL";else if($28!=$31)print $0,$42=1;else print $0,$42=0;}' OFS='\t' $i.4 > $i.5
  awk '{if($33=="NULL" || $36=="NULL")print $0,$43="NULL";else if($33!=$36)print $0,$43=1;else print $0,$43=0;}' OFS='\t' $i.5 > $i.6
done
# Mz-On - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($28=="NULL" || $32=="NULL")print $0,$44="NULL";else if($28!=$32)print $0,$44=1;else print $0,$44=0;}' OFS='\t' $i.6 > $i.7
  awk '{if($33=="NULL" || $37=="NULL")print $0,$45="NULL";else if($33!=$37)print $0,$45=1;else print $0,$45=0;}' OFS='\t' $i.7 > $i.8
done
# Pn-Ab - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($29=="NULL" || $30=="NULL")print $0,$46="NULL";else if($29!=$30)print $0,$46=1;else print $0,$46=0;}' OFS='\t' $i.8 > $i.9
  awk '{if($34=="NULL" || $35=="NULL")print $0,$47="NULL";else if($34!=$35)print $0,$47=1;else print $0,$47=0;}' OFS='\t' $i.9 > $i.10
done
# Pn-Nb - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($29=="NULL" || $31=="NULL")print $0,$48="NULL";else if($29!=$31)print $0,$48=1;else print $0,$48=0;}' OFS='\t' $i.10 > $i.11
  awk '{if($34=="NULL" || $36=="NULL")print $0,$49="NULL";else if($34!=$36)print $0,$49=1;else print $0,$49=0;}' OFS='\t' $i.11 > $i.12
done
# Pn-On - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($29=="NULL" || $32=="NULL")print $0,$50="NULL";else if($29!=$32)print $0,$50=1;else print $0,$50=0;}' OFS='\t' $i.12 > $i.13
  awk '{if($34=="NULL" || $37=="NULL")print $0,$51="NULL";else if($34!=$37)print $0,$51=1;else print $0,$51=0;}' OFS='\t' $i.13 > $i.14
done
# Ab-Nb - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($30=="NULL" || $31=="NULL")print $0,$52="NULL";else if($30!=$31)print $0,$52=1;else print $0,$52=0;}' OFS='\t' $i.14 > $i.15
  awk '{if($35=="NULL" || $36=="NULL")print $0,$53="NULL";else if($35!=$36)print $0,$53=1;else print $0,$53=0;}' OFS='\t' $i.15 > $i.16
done
# Ab-On - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($30=="NULL" || $32=="NULL")print $0,$54="NULL";else if($30!=$32)print $0,$54=1;else print $0,$54=0;}' OFS='\t' $i.16 > $i.17
  awk '{if($35=="NULL" || $37=="NULL")print $0,$55="NULL";else if($35!=$37)print $0,$55=1;else print $0,$55=0;}' OFS='\t' $i.17 > $i.18
done
# Nb-On - GeneA; GeneB
for i in *NEW.matrix2.txt ; do
  awk '{if($31=="NULL" || $32=="NULL")print $0,$56="NULL";else if($31!=$32)print $0,$56=1;else print $0,$56=0;}' OFS='\t' $i.18 > $i.19
  awk '{if($35=="NULL" || $37=="NULL")print $0,$57="NULL";else if($36!=$37)print $0,$57=1;else print $0,$57=0;}' OFS='\t' $i.19 > $i.20
done

printf 'interaction\tCichlid_GeneA\tCichlid_GeneA_symbol\tGeneA\tCichlid_GeneB\tCichlid_GeneB_symbolGa\tCichlid_GeneB_symbolDr\tGeneB\tdirectness\tlayer\tsource\tmotif_pattern[pm_cm]\tstart[pm_cm]\tstop[pm_cm]\tstrand[pm_cm]\tscore[pm_cm]\tp.value[pm_cm]\tq.value[pm_cm]\tsequence[pm_cm]\tconf_level[pm_cm]\tconf_score[pm_cm]\tconfidence[coTFTG]\tMz.P-1_A-0\tPn.P-1_A-0\tAb.P-1_A-0\tNb.P-1_A-0\tOn.P-1_A-0\tMz_GeneAModule\tPn_GeneAModule\tAb_GeneAModule\tNb_GeneAModule\tOn_GeneAModule\tMz_GeneBModule\tPn_GeneBModule\tAb_GeneBModule\tNb_GeneBModule\tOn_GeneBModule\tMz-Pn_GeneA\tMz-Pn_GeneB\tMz-Ab_GeneA\tMz-Ab_GeneB\tMz-Nb_GeneA\tMz-Nb_GeneB\tMz-On_GeneA\tMz-On_GeneB\tPn-Ab_GeneA\tPn-Ab_GeneB\tPn-Nb_GeneA\tPn-Nb_GeneB\tPn-On_GeneA\tPn-On_GeneB\tAb-Nb_GeneA\tAb-Nb_GeneB\tAb-On_GeneA\tAb-On_GeneB\tNb-On_GeneA\tNb-On_GeneB\n' > matrix_colheads.NEW3
cat matrix_colheads.NEW3 Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix2.txt.20 | sed '2d' > Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix3.txt

rm *.matrix2.txt.*
rm Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix2.txt

# awk '{print NF}' Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix3.txt | sort -nu | head -n 1 # 39
# awk '{print NF}' Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix3.txt | sort -nu | tail -n 1 # 39


qsub -q Test -l select=1:mem=60GB:ncpus=1 NEWcoexpr_prom_matrices-pt2.sh

# Filter your edge matrix by removing the non unique columns like cichlid gene id
cut -f1,3,4,6-12,20-57 Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix3.txt | sort -u > Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.txt # 1131811 edges

## final edge matrices
/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEW.matrix3.txt # these are per edge type for each species - elongated table with replicated edges across the species (28948929 edges)
/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.txt # these are the simplified tables of only one edge per line (1131811 edges)




## B. Confirm how many of the 18,799 OGs as used for modules of co-expressed that are non-1-to-1 orthologous are truly not present in the genome
# 11,955 OGs are non-1-to-1 orthologous

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/
mkdir OGIDtblastx
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx

# 1. Do a tblastx of all longest cds from each species to other genomes

# 1a. amend the bed files into correct format
cp /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.bed .
cp /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.bed On_GeneAnnotation_11092017_FINALcorrected.bed
for i in *.bed ; do awk '{print $1,$2,$3,$5,"1000",$4}' OFS='\t' $i > "$(basename "$i" .bed)format.bed" ; done
rm *corrected.bed

# 1b. pull out the sequences of the longest cds in each species
nano getlongestcdsfasta.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml bedtools/2.25.0
ml zlib

# get fasta of the longest cds using the name column to name each sequence and force strandedness (if in -ve strand, sequence will be rev comp)
bedtools getfasta -name -s -fo Ab_alllongestcds.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -bed Ab_GeneAnnotation_11092017_FINALcorrectedformat.bed
bedtools getfasta -name -s -fo Mz_alllongestcds.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -bed Mz_GeneAnnotation_11092017_FINALcorrectedformat.bed
bedtools getfasta -name -s -fo Nb_alllongestcds.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -bed Nb_GeneAnnotation_11092017_FINALcorrectedformat.bed
bedtools getfasta -name -s -fo On_alllongestcds.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -bed On_GeneAnnotation_11092017_FINALcorrectedformat.bed
bedtools getfasta -name -s -fo Pn_alllongestcds.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -bed Pn_GeneAnnotation_11092017_FINALcorrectedformat.bed

# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=1 getlongestcdsfasta.sh #{DONE}

# copy all files to workarea for any future use - updated in ref locations
cp *correctedformat.bed /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/
cp *_alllongestcds.fasta /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/

# 1c. run the tblastx of ALL predicted genes against all other genomes including the PacBio genomes (no filtering on OGID) - run with e-val <1e-5

# best to split up the fasta of each species file
mkdir fastasplit
mkdir -p fastasplit/{mz,pn,ab,nb,on}
mkdir -p fastasplit/{mz,pn,ab,nb,on}/{a,b,c,d,e} # create folders so that split files can go in there

nano fastasplits.sh

#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 28000 # memory pool for all cores
#SBATCH -t 0-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/mz
awk '/^>mz/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Mz_alllongestcds.fasta # split the fasta by each gene
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #1673

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/pn
awk '/^>pn/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Pn_alllongestcds.fasta # split the fasta by each gene
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #611

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/ab
awk '/^>ab/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Ab_alllongestcds.fasta # split the fasta by each gene
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #3436

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/nb
awk '/^>nb/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Nb_alllongestcds.fasta # split the fasta by each gene
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #119

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/on
awk '/^>on/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../On_alllongestcds.fasta # split the fasta by each gene
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #4559

## create one job script for the fixed 5000 jobs, and separate scripts for the smaller folder e's

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx

nano tblastx-a_d-Mz-array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-4999
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

cd fastasplit/mz

ls -1 a/*.fa > mzreads_a # create a list of all fastq files
mapfile -t mzreads_a < mzreads_a # assign as elements to $mzreads_a variable
ls -1 b/*.fa > mzreads_b # create a list of all fastq files
mapfile -t mzreads_b < mzreads_b # assign as elements to $mzreads_b variable
ls -1 c/*.fa > mzreads_c # create a list of all fastq files
mapfile -t mzreads_c < mzreads_c # assign as elements to $mzreads_c variable
ls -1 d/*.fa > mzreads_d # create a list of all fastq files
mapfile -t mzreads_d < mzreads_d # assign as elements to $mzreads_d variable

#load the latest blast module
ml blast/2.3.0

srun tblastx -query ${mzreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${mzreads_a[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${mzreads_a[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${mzreads_a[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${mzreads_a[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${mzreads_a[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${mzreads_a[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${mzreads_a[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${mzreads_b[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${mzreads_b[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${mzreads_b[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${mzreads_b[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${mzreads_b[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${mzreads_b[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${mzreads_b[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${mzreads_c[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${mzreads_c[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${mzreads_c[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${mzreads_c[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${mzreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${mzreads_c[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${mzreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${mzreads_d[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${mzreads_d[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${mzreads_d[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${mzreads_d[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${mzreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${mzreads_d[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${mzreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5


nano tblastx-a_d-Pn-array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-4999
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

cd fastasplit/pn

ls -1 a/*.fa > pnreads_a # create a list of all fastq files
mapfile -t pnreads_a < pnreads_a # assign as elements to $pnreads_a variable
ls -1 b/*.fa > pnreads_b # create a list of all fastq files
mapfile -t pnreads_b < pnreads_b # assign as elements to $pnreads_b variable
ls -1 c/*.fa > pnreads_c # create a list of all fastq files
mapfile -t pnreads_c < pnreads_c # assign as elements to $pnreads_c variable
ls -1 d/*.fa > pnreads_d # create a list of all fastq files
mapfile -t pnreads_d < pnreads_d # assign as elements to $pnreads_d variable

#load the latest blast module
ml blast/2.3.0

srun tblastx -query ${pnreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${pnreads_a[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${pnreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${pnreads_a[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${pnreads_a[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${pnreads_a[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${pnreads_a[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${pnreads_a[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${pnreads_a[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${pnreads_b[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${pnreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${pnreads_b[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${pnreads_b[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${pnreads_b[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${pnreads_b[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${pnreads_b[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${pnreads_b[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${pnreads_c[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${pnreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${pnreads_c[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${pnreads_c[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${pnreads_c[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${pnreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${pnreads_c[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${pnreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${pnreads_d[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${pnreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${pnreads_d[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${pnreads_d[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${pnreads_d[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${pnreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${pnreads_d[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${pnreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5


nano tblastx-a_d-Ab-array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-4999
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

cd fastasplit/ab

ls -1 a/*.fa > abreads_a # create a list of all fastq files
mapfile -t abreads_a < abreads_a # assign as elements to $abreads_a variable
ls -1 b/*.fa > abreads_b # create a list of all fastq files
mapfile -t abreads_b < abreads_b # assign as elements to $abreads_b variable
ls -1 c/*.fa > abreads_c # create a list of all fastq files
mapfile -t abreads_c < abreads_c # assign as elements to $abreads_c variable
ls -1 d/*.fa > abreads_d # create a list of all fastq files
mapfile -t abreads_d < abreads_d # assign as elements to $abreads_d variable

#load the latest blast module
ml blast/2.3.0

srun tblastx -query ${abreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${abreads_a[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${abreads_a[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${abreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${abreads_a[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${abreads_a[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${abreads_a[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${abreads_a[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${abreads_a[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${abreads_b[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${abreads_b[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${abreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${abreads_b[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${abreads_b[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${abreads_b[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${abreads_b[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${abreads_b[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${abreads_c[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${abreads_c[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${abreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${abreads_c[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${abreads_c[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${abreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${abreads_c[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${abreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${abreads_d[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${abreads_d[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${abreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${abreads_d[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${abreads_d[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${abreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${abreads_d[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${abreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5


nano tblastx-a_d-Nb-array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-4999
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

cd fastasplit/nb

ls -1 a/*.fa > nbreads_a # create a list of all fastq files
mapfile -t nbreads_a < nbreads_a # assign as elements to $nbreads_a variable
ls -1 b/*.fa > nbreads_b # create a list of all fastq files
mapfile -t nbreads_b < nbreads_b # assign as elements to $nbreads_b variable
ls -1 c/*.fa > nbreads_c # create a list of all fastq files
mapfile -t nbreads_c < nbreads_c # assign as elements to $nbreads_c variable
ls -1 d/*.fa > nbreads_d # create a list of all fastq files
mapfile -t nbreads_d < nbreads_d # assign as elements to $nbreads_d variable

#load the latest blast module
ml blast/2.3.0

srun tblastx -query ${nbreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${nbreads_a[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${nbreads_a[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${nbreads_a[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${nbreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${nbreads_a[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${nbreads_a[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${nbreads_a[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${nbreads_a[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${nbreads_b[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${nbreads_b[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${nbreads_b[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${nbreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${nbreads_b[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${nbreads_b[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${nbreads_b[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${nbreads_b[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${nbreads_c[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${nbreads_c[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${nbreads_c[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${nbreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${nbreads_c[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${nbreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${nbreads_c[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${nbreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${nbreads_d[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${nbreads_d[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${nbreads_d[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${nbreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${nbreads_d[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${nbreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${nbreads_d[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${nbreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5


nano tblastx-a_d-On-array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-4999
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

cd fastasplit/on

ls -1 a/*.fa > onreads_a # create a list of all fastq files
mapfile -t onreads_a < onreads_a # assign as elements to $onreads_a variable
ls -1 b/*.fa > onreads_b # create a list of all fastq files
mapfile -t onreads_b < onreads_b # assign as elements to $onreads_b variable
ls -1 c/*.fa > onreads_c # create a list of all fastq files
mapfile -t onreads_c < onreads_c # assign as elements to $onreads_c variable
ls -1 d/*.fa > onreads_d # create a list of all fastq files
mapfile -t onreads_d < onreads_d # assign as elements to $onreads_d variable

#load the latest blast module
ml blast/2.3.0

srun tblastx -query ${onreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${onreads_a[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${onreads_a[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${onreads_a[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${onreads_a[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${onreads_a[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${onreads_a[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_a[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${onreads_a[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${onreads_b[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${onreads_b[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${onreads_b[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${onreads_b[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${onreads_b[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${onreads_b[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_b[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${onreads_b[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${onreads_c[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${onreads_c[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${onreads_c[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${onreads_c[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${onreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${onreads_c[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_c[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${onreads_c[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${onreads_d[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${onreads_d[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${onreads_d[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${onreads_d[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${onreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${onreads_d[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_d[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${onreads_d[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5

## run all of the above
sbatch tblastx-a_d-Mz-array.sh # {DONE}
sbatch tblastx-a_d-Pn-array.sh # {DONE}
sbatch tblastx-a_d-Ab-array.sh # {DONE}
sbatch tblastx-a_d-Nb-array.sh # {DONE}
sbatch tblastx-a_d-On-array.sh # {DONE}

# create other scripts for the remaining smaller number of files
# Mz -e
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/mz

nano tblastx-e-Mz-array.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-1672
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-13:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > mzreads_e # create a list of all fastq files
mapfile -t mzreads_e < mzreads_e # assign as elements to $mzreads_e variable

#load the latest blast module
ml blast/2.3.0

# srun tblastx -query ${mzreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${mzreads_e[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${mzreads_e[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${mzreads_e[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${mzreads_e[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${mzreads_e[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${mzreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${mzreads_e[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${mzreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${mzreads_e[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5


# Pn -e
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/pn

nano tblastx-e-Pn-array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-610
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > pnreads_e # create a list of all fastq files
mapfile -t pnreads_e < pnreads_e # assign as elements to $pnreads_e variable

#load the latest blast module
ml blast/2.3.0

srun tblastx -query ${pnreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${pnreads_e[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${pnreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${pnreads_e[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${pnreads_e[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${pnreads_e[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${pnreads_e[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${pnreads_e[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${pnreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${pnreads_e[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5


# Ab -e
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/ab

nano tblastx-e-Ab-array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-3435
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > abreads_e # create a list of all fastq files
mapfile -t abreads_e < abreads_e # assign as elements to $abreads_e variable

#load the latest blast module
ml blast/2.3.0

srun tblastx -query ${abreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${abreads_e[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${abreads_e[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${abreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${abreads_e[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${abreads_e[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${abreads_e[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${abreads_e[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${abreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${abreads_e[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5


# Nb -e
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/nb

nano tblastx-e-Nb-array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-118
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > nbreads_e # create a list of all fastq files
mapfile -t nbreads_e < nbreads_e # assign as elements to $nbreads_e variable

#load the latest blast module
ml blast/2.3.0

srun tblastx -query ${nbreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${nbreads_e[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${nbreads_e[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${nbreads_e[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${nbreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${nbreads_e[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${nbreads_e[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${nbreads_e[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${nbreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${nbreads_e[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5


# On -e
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/on

nano tblastx-e-On-array.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --array=0-4558
#SBATCH --mem-per-cpu 12000
#SBATCH -t 0-00:45
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > onreads_e # create a list of all fastq files
mapfile -t onreads_e < onreads_e # assign as elements to $onreads_e variable

#load the latest blast module
ml blast/2.3.0

srun tblastx -query ${onreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -out ${onreads_e[${SLURM_ARRAY_TASK_ID}]}.MzGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -out ${onreads_e[${SLURM_ARRAY_TASK_ID}]}.PnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -out ${onreads_e[${SLURM_ARRAY_TASK_ID}]}.AbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -out ${onreads_e[${SLURM_ARRAY_TASK_ID}]}.NbGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -out ${onreads_e[${SLURM_ARRAY_TASK_ID}]}.OnGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
srun tblastx -query ${onreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -out ${onreads_e[${SLURM_ARRAY_TASK_ID}]}.MzPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5
# srun tblastx -query ${onreads_e[${SLURM_ARRAY_TASK_ID}]} -db /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta -out ${onreads_e[${SLURM_ARRAY_TASK_ID}]}.OnPBGenome.blast -outfmt 6 -num_threads 2 -evalue 1e-5


## run all of the above
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/mz
sbatch tblastx-e-Mz-array.sh #1673  # {DONE}
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/pn
sbatch tblastx-e-Pn-array.sh #611 # {DONE}
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/ab
sbatch tblastx-e-Ab-array.sh #3436 # {DONE}
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/nb
sbatch tblastx-e-Nb-array.sh #119 # {DONE}
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/on
sbatch tblastx-e-On-array.sh #4559 #{DONE}


# 2. Will prepared a portable python script (tblastx_parser2.py - originally used for Arabidopsis) that
# a. take raw blast output - ensure that original output to be used as input here is *.blast
# b. script will filter the blast of results of hits lower than pid of 30 and e value higher than 1e-10
# c. The parsed remaining results are returned along with an additional column which works out the percent coverage of the CDS that the hit is made up of. This can be useful for filtering later.
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/tblastx_parser2.py .

nano create_tblastx_scripts.sh

#!/bin/sh

echo '#!/bin/bash -e' > runtblastx_parser-mz_a.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-mz_a.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-mz_a.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-mz_a.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-mz_a.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-mz_a.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-mz_a.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-mz_a.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-mz_a.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-mz_a.sh
echo -e '\n' >> runtblastx_parser-mz_a.sh
echo 'for a in fastasplit/mz/a/*.blast; do' >> runtblastx_parser-mz_a.sh
echo -e '\tpython tblastx_parser2.py $a Mz_alllongestcds.fasta || true ;' >> runtblastx_parser-mz_a.sh
echo 'done' >> runtblastx_parser-mz_a.sh

echo '#!/bin/bash -e' > runtblastx_parser-mz_b.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-mz_b.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-mz_b.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-mz_b.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-mz_b.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-mz_b.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-mz_b.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-mz_b.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-mz_b.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-mz_b.sh
echo -e '\n' >> runtblastx_parser-mz_b.sh
echo 'for b in fastasplit/mz/b/*.blast; do' >> runtblastx_parser-mz_b.sh
echo -e '\tpython tblastx_parser2.py $b Mz_alllongestcds.fasta || true ;' >> runtblastx_parser-mz_b.sh
echo 'done' >> runtblastx_parser-mz_b.sh

echo '#!/bin/bash -e' > runtblastx_parser-mz_c.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-mz_c.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-mz_c.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-mz_c.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-mz_c.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-mz_c.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-mz_c.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-mz_c.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-mz_c.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-mz_c.sh
echo -e '\n' >> runtblastx_parser-mz_c.sh
echo 'for c in fastasplit/mz/c/*.blast; do' >> runtblastx_parser-mz_c.sh
echo -e '\tpython tblastx_parser2.py $c Mz_alllongestcds.fasta || true ;' >> runtblastx_parser-mz_c.sh
echo 'done' >> runtblastx_parser-mz_c.sh

echo '#!/bin/bash -e' > runtblastx_parser-mz_d.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-mz_d.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-mz_d.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-mz_d.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-mz_d.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-mz_d.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-mz_d.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-mz_d.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-mz_d.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-mz_d.sh
echo -e '\n' >> runtblastx_parser-mz_d.sh
echo 'for d in fastasplit/mz/d/*.blast; do' >> runtblastx_parser-mz_d.sh
echo -e '\tpython tblastx_parser2.py $d Mz_alllongestcds.fasta || true ;' >> runtblastx_parser-mz_d.sh
echo 'done' >> runtblastx_parser-mz_d.sh

echo '#!/bin/bash -e' > runtblastx_parser-mz_e.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-mz_e.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-mz_e.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-mz_e.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-mz_e.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-mz_e.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-mz_e.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-mz_e.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-mz_e.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-mz_e.sh
echo -e '\n' >> runtblastx_parser-mz_e.sh
echo 'for e in fastasplit/mz/e/*.blast; do' >> runtblastx_parser-mz_e.sh
echo -e '\tpython tblastx_parser2.py $e Mz_alllongestcds.fasta || true ;' >> runtblastx_parser-mz_e.sh
echo 'done' >> runtblastx_parser-mz_e.sh


echo '#!/bin/bash -e' > runtblastx_parser-pn_a.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-pn_a.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-pn_a.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-pn_a.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-pn_a.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-pn_a.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-pn_a.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-pn_a.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-pn_a.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-pn_a.sh
echo -e '\n' >> runtblastx_parser-pn_a.sh
echo 'for a in fastasplit/pn/a/*.blast; do' >> runtblastx_parser-pn_a.sh
echo -e '\tpython tblastx_parser2.py $a Pn_alllongestcds.fasta || true ;' >> runtblastx_parser-pn_a.sh
echo 'done' >> runtblastx_parser-pn_a.sh

echo '#!/bin/bash -e' > runtblastx_parser-pn_b.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-pn_b.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-pn_b.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-pn_b.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-pn_b.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-pn_b.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-pn_b.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-pn_b.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-pn_b.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-pn_b.sh
echo -e '\n' >> runtblastx_parser-pn_b.sh
echo 'for b in fastasplit/pn/b/*.blast; do' >> runtblastx_parser-pn_b.sh
echo -e '\tpython tblastx_parser2.py $b Pn_alllongestcds.fasta || true ;' >> runtblastx_parser-pn_b.sh
echo 'done' >> runtblastx_parser-pn_b.sh

echo '#!/bin/bash -e' > runtblastx_parser-pn_c.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-pn_c.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-pn_c.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-pn_c.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-pn_c.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-pn_c.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-pn_c.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-pn_c.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-pn_c.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-pn_c.sh
echo -e '\n' >> runtblastx_parser-pn_c.sh
echo 'for c in fastasplit/pn/c/*.blast; do' >> runtblastx_parser-pn_c.sh
echo -e '\tpython tblastx_parser2.py $c Pn_alllongestcds.fasta || true ;' >> runtblastx_parser-pn_c.sh
echo 'done' >> runtblastx_parser-pn_c.sh

echo '#!/bin/bash -e' > runtblastx_parser-pn_d.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-pn_d.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-pn_d.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-pn_d.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-pn_d.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-pn_d.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-pn_d.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-pn_d.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-pn_d.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-pn_d.sh
echo -e '\n' >> runtblastx_parser-pn_d.sh
echo 'for d in fastasplit/pn/d/*.blast; do' >> runtblastx_parser-pn_d.sh
echo -e '\tpython tblastx_parser2.py $d Pn_alllongestcds.fasta || true ;' >> runtblastx_parser-pn_d.sh
echo 'done' >> runtblastx_parser-pn_d.sh

echo '#!/bin/bash -e' > runtblastx_parser-pn_e.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-pn_e.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-pn_e.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-pn_e.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-pn_e.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-pn_e.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-pn_e.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-pn_e.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-pn_e.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-pn_e.sh
echo -e '\n' >> runtblastx_parser-pn_e.sh
echo 'for e in fastasplit/pn/e/*.blast; do' >> runtblastx_parser-pn_e.sh
echo -e '\tpython tblastx_parser2.py $e Pn_alllongestcds.fasta || true ;' >> runtblastx_parser-pn_e.sh
echo 'done' >> runtblastx_parser-pn_e.sh

echo '#!/bin/bash -e' > runtblastx_parser-ab_a.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-ab_a.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-ab_a.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-ab_a.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-ab_a.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-ab_a.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-ab_a.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-ab_a.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-ab_a.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-ab_a.sh
echo -e '\n' >> runtblastx_parser-ab_a.sh
echo 'for a in fastasplit/ab/a/*.blast; do' >> runtblastx_parser-ab_a.sh
echo -e '\tpython tblastx_parser2.py $a Ab_alllongestcds.fasta || true ;' >> runtblastx_parser-ab_a.sh
echo 'done' >> runtblastx_parser-ab_a.sh

echo '#!/bin/bash -e' > runtblastx_parser-ab_b.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-ab_b.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-ab_b.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-ab_b.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-ab_b.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-ab_b.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-ab_b.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-ab_b.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-ab_b.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-ab_b.sh
echo -e '\n' >> runtblastx_parser-ab_b.sh
echo 'for b in fastasplit/ab/b/*.blast; do' >> runtblastx_parser-ab_b.sh
echo -e '\tpython tblastx_parser2.py $b Ab_alllongestcds.fasta || true ;' >> runtblastx_parser-ab_b.sh
echo 'done' >> runtblastx_parser-ab_b.sh

echo '#!/bin/bash -e' > runtblastx_parser-ab_c.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-ab_c.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-ab_c.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-ab_c.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-ab_c.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-ab_c.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-ab_c.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-ab_c.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-ab_c.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-ab_c.sh
echo -e '\n' >> runtblastx_parser-ab_c.sh
echo 'for c in fastasplit/ab/c/*.blast; do' >> runtblastx_parser-ab_c.sh
echo -e '\tpython tblastx_parser2.py $c Ab_alllongestcds.fasta || true ;' >> runtblastx_parser-ab_c.sh
echo 'done' >> runtblastx_parser-ab_c.sh

echo '#!/bin/bash -e' > runtblastx_parser-ab_d.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-ab_d.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-ab_d.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-ab_d.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-ab_d.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-ab_d.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-ab_d.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-ab_d.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-ab_d.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-ab_d.sh
echo -e '\n' >> runtblastx_parser-ab_d.sh
echo 'for d in fastasplit/ab/d/*.blast; do' >> runtblastx_parser-ab_d.sh
echo -e '\tpython tblastx_parser2.py $d Ab_alllongestcds.fasta || true ;' >> runtblastx_parser-ab_d.sh
echo 'done' >> runtblastx_parser-ab_d.sh

echo '#!/bin/bash -e' > runtblastx_parser-ab_e.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-ab_e.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-ab_e.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-ab_e.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-ab_e.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-ab_e.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-ab_e.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-ab_e.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-ab_e.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-ab_e.sh
echo -e '\n' >> runtblastx_parser-ab_e.sh
echo 'for e in fastasplit/ab/e/*.blast; do' >> runtblastx_parser-ab_e.sh
echo -e '\tpython tblastx_parser2.py $e Ab_alllongestcds.fasta || true ;' >> runtblastx_parser-ab_e.sh
echo 'done' >> runtblastx_parser-ab_e.sh

echo '#!/bin/bash -e' > runtblastx_parser-nb_a.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-nb_a.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-nb_a.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-nb_a.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-nb_a.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-nb_a.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-nb_a.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-nb_a.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-nb_a.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-nb_a.sh
echo -e '\n' >> runtblastx_parser-nb_a.sh
echo 'for a in fastasplit/nb/a/*.blast; do' >> runtblastx_parser-nb_a.sh
echo -e '\tpython tblastx_parser2.py $a Nb_alllongestcds.fasta || true ;' >> runtblastx_parser-nb_a.sh
echo 'done' >> runtblastx_parser-nb_a.sh

echo '#!/bin/bash -e' > runtblastx_parser-nb_b.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-nb_b.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-nb_b.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-nb_b.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-nb_b.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-nb_b.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-nb_b.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-nb_b.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-nb_b.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-nb_b.sh
echo -e '\n' >> runtblastx_parser-nb_b.sh
echo 'for b in fastasplit/nb/b/*.blast; do' >> runtblastx_parser-nb_b.sh
echo -e '\tpython tblastx_parser2.py $b Nb_alllongestcds.fasta || true ;' >> runtblastx_parser-nb_b.sh
echo 'done' >> runtblastx_parser-nb_b.sh

echo '#!/bin/bash -e' > runtblastx_parser-nb_c.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-nb_c.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-nb_c.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-nb_c.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-nb_c.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-nb_c.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-nb_c.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-nb_c.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-nb_c.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-nb_c.sh
echo -e '\n' >> runtblastx_parser-nb_c.sh
echo 'for c in fastasplit/nb/c/*.blast; do' >> runtblastx_parser-nb_c.sh
echo -e '\tpython tblastx_parser2.py $c Nb_alllongestcds.fasta || true ;' >> runtblastx_parser-nb_c.sh
echo 'done' >> runtblastx_parser-nb_c.sh

echo '#!/bin/bash -e' > runtblastx_parser-nb_d.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-nb_d.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-nb_d.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-nb_d.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-nb_d.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-nb_d.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-nb_d.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-nb_d.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-nb_d.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-nb_d.sh
echo -e '\n' >> runtblastx_parser-nb_d.sh
echo 'for d in fastasplit/nb/d/*.blast; do' >> runtblastx_parser-nb_d.sh
echo -e '\tpython tblastx_parser2.py $d Nb_alllongestcds.fasta || true ;' >> runtblastx_parser-nb_d.sh
echo 'done' >> runtblastx_parser-nb_d.sh

echo '#!/bin/bash -e' > runtblastx_parser-nb_e.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-nb_e.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-nb_e.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-nb_e.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-nb_e.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-nb_e.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-nb_e.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-nb_e.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-nb_e.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-nb_e.sh
echo -e '\n' >> runtblastx_parser-nb_e.sh
echo 'for e in fastasplit/nb/e/*.blast; do' >> runtblastx_parser-nb_e.sh
echo -e '\tpython tblastx_parser2.py $e Nb_alllongestcds.fasta || true ;' >> runtblastx_parser-nb_e.sh
echo 'done' >> runtblastx_parser-nb_e.sh

echo '#!/bin/bash -e' > runtblastx_parser-on_a.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-on_a.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-on_a.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-on_a.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-on_a.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-on_a.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-on_a.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-on_a.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-on_a.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-on_a.sh
echo -e '\n' >> runtblastx_parser-on_a.sh
echo 'for a in fastasplit/on/a/*.blast; do' >> runtblastx_parser-on_a.sh
echo -e '\tpython tblastx_parser2.py $a On_alllongestcds.fasta || true ;' >> runtblastx_parser-on_a.sh
echo 'done' >> runtblastx_parser-on_a.sh

echo '#!/bin/bash -e' > runtblastx_parser-on_b.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-on_b.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-on_b.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-on_b.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-on_b.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-on_b.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-on_b.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-on_b.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-on_b.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-on_b.sh
echo -e '\n' >> runtblastx_parser-on_b.sh
echo 'for b in fastasplit/on/b/*.blast; do' >> runtblastx_parser-on_b.sh
echo -e '\tpython tblastx_parser2.py $b On_alllongestcds.fasta || true ;' >> runtblastx_parser-on_b.sh
echo 'done' >> runtblastx_parser-on_b.sh

echo '#!/bin/bash -e' > runtblastx_parser-on_c.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-on_c.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-on_c.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-on_c.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-on_c.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-on_c.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-on_c.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-on_c.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-on_c.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-on_c.sh
echo -e '\n' >> runtblastx_parser-on_c.sh
echo 'for c in fastasplit/on/c/*.blast; do' >> runtblastx_parser-on_c.sh
echo -e '\tpython tblastx_parser2.py $c On_alllongestcds.fasta || true ;' >> runtblastx_parser-on_c.sh
echo 'done' >> runtblastx_parser-on_c.sh

echo '#!/bin/bash -e' > runtblastx_parser-on_d.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-on_d.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-on_d.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-on_d.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-on_d.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-on_d.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-on_d.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-on_d.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-on_d.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-on_d.sh
echo -e '\n' >> runtblastx_parser-on_d.sh
echo 'for d in fastasplit/on/d/*.blast; do' >> runtblastx_parser-on_d.sh
echo -e '\tpython tblastx_parser2.py $d On_alllongestcds.fasta || true ;' >> runtblastx_parser-on_d.sh
echo 'done' >> runtblastx_parser-on_d.sh

echo '#!/bin/bash -e' > runtblastx_parser-on_e.sh
echo '#SBATCH -p tgac-long # partition (queue)' >> runtblastx_parser-on_e.sh
echo '#SBATCH -N 1 # number of nodes' >> runtblastx_parser-on_e.sh
echo '#SBATCH -n 1 # number of tasks' >> runtblastx_parser-on_e.sh
echo '#SBATCH --mem 18000 # memory pool for all cores' >> runtblastx_parser-on_e.sh
echo '#SBATCH -t 2-15:59 # time (D-HH:MM)' >> runtblastx_parser-on_e.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> runtblastx_parser-on_e.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> runtblastx_parser-on_e.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> runtblastx_parser-on_e.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> runtblastx_parser-on_e.sh
echo -e '\n' >> runtblastx_parser-on_e.sh
echo 'for e in fastasplit/on/e/*.blast; do' >> runtblastx_parser-on_e.sh
echo -e '\tpython tblastx_parser2.py $e On_alllongestcds.fasta || true ;' >> runtblastx_parser-on_e.sh
echo 'done' >> runtblastx_parser-on_e.sh

# create the scripts
sh create_tblastx_scripts.sh

# run all of the above
for i in mz pn ab nb on ; do
  sbatch runtblastx_parser-${i}_a.sh
  sbatch runtblastx_parser-${i}_b.sh
  sbatch runtblastx_parser-${i}_c.sh
  sbatch runtblastx_parser-${i}_d.sh
  sbatch runtblastx_parser-${i}_e.sh
done

# check that input files are the same number as output files - if ok, then proceed
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/
for i in mz pn ab nb on ; do
  ls -1 fastasplit/${i}/a/*Genome.blast | wc -l # input files in folder a
  ls -1 fastasplit/${i}/b/*Genome.blast | wc -l # input files in folder b
  ls -1 fastasplit/${i}/c/*Genome.blast | wc -l # input files in folder c
  ls -1 fastasplit/${i}/d/*Genome.blast | wc -l # input files in folder d
  ls -1 fastasplit/${i}/e/*Genome.blast | wc -l # input files in folder e
  ls -1 fastasplit/${i}/a/*.filtered.blast | wc -l # output files in folder a
  ls -1 fastasplit/${i}/b/*.filtered.blast | wc -l # output files in folder b
  ls -1 fastasplit/${i}/c/*.filtered.blast | wc -l # output files in folder c
  ls -1 fastasplit/${i}/d/*.filtered.blast | wc -l # output files in folder d
  ls -1 fastasplit/${i}/e/*.filtered.blast | wc -l # output files in folder e
done

# Mz
# 32753 - input files in folder a
# 26419 - input files in folder b
# 18214 - input files in folder c
# 14365 - input files in folder d
# 8143 - input files in folder e
# 32753 - output files in folder a
# 26419 - output files in folder a
# 18214 - output files in folder a
# 14365 - output files in folder a
# 8143 - output files in folder a

# Pn
# 22658
# 8174
# 1892
# 358
# 2702
# 22658
# 8174
# 1892
# 358
# 2702

# Ab
# 22770
# 8563
# 2268
# 519
# 16059
# 22770
# 8563
# 2268
# 519
# 16059

# Nb
# 23089
# 9126
# 2443
# 590
# 458
# 23089
# 9126
# 2443
# 590
# 458

# On
# 19049
# 9812
# 3835
# 1167
# 19568
# 19049
# 9812
# 3835
# 1167
# 19568

# 2d. Merge all individual gene blast outputs from above into single species-genome files using array
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/

# create files to run slurm arrays
nano createcatcommands.sh

#!/bin/sh

echo 'cat a/*PnGenome.filtered.blast b/*PnGenome.filtered.blast c/*PnGenome.filtered.blast d/*PnGenome.filtered.blast e/*PnGenome.filtered.blast' > mzcat
echo 'cat a/*AbGenome.filtered.blast b/*AbGenome.filtered.blast c/*AbGenome.filtered.blast d/*AbGenome.filtered.blast e/*AbGenome.filtered.blast' >> mzcat
echo 'cat a/*NbGenome.filtered.blast b/*NbGenome.filtered.blast c/*NbGenome.filtered.blast d/*NbGenome.filtered.blast e/*NbGenome.filtered.blast' >> mzcat
echo 'cat a/*OnGenome.filtered.blast b/*OnGenome.filtered.blast c/*OnGenome.filtered.blast d/*OnGenome.filtered.blast e/*OnGenome.filtered.blast' >> mzcat
echo 'cat a/*OnPBGenome.filtered.blast b/*OnPBGenome.filtered.blast c/*OnPBGenome.filtered.blast d/*OnPBGenome.filtered.blast e/*OnPBGenome.filtered.blast' >> mzcat

echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Mzalllongestcds-Pn_Genome.filtered.blast' > mzcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Mzalllongestcds-Ab_Genome.filtered.blast' >> mzcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Mzalllongestcds-Nb_Genome.filtered.blast' >> mzcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Mzalllongestcds-On_Genome.filtered.blast' >> mzcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Mz_alllongestcds-OnPacBio_Genome.filtered.blast' >> mzcatout

echo 'cat a/*MzGenome.filtered.blast b/*MzGenome.filtered.blast c/*MzGenome.filtered.blast d/*MzGenome.filtered.blast e/*MzGenome.filtered.blast' > pncat
echo 'cat a/*AbGenome.filtered.blast b/*AbGenome.filtered.blast c/*AbGenome.filtered.blast d/*AbGenome.filtered.blast e/*AbGenome.filtered.blast' >> pncat
echo 'cat a/*NbGenome.filtered.blast b/*NbGenome.filtered.blast c/*NbGenome.filtered.blast d/*NbGenome.filtered.blast e/*NbGenome.filtered.blast' >> pncat
echo 'cat a/*OnGenome.filtered.blast b/*OnGenome.filtered.blast c/*OnGenome.filtered.blast d/*OnGenome.filtered.blast e/*OnGenome.filtered.blast' >> pncat
echo 'cat a/*MzPBGenome.filtered.blast b/*MzPBGenome.filtered.blast c/*MzPBGenome.filtered.blast d/*MzPBGenome.filtered.blast e/*MzPBGenome.filtered.blast' >> pncat
echo 'cat a/*OnPBGenome.filtered.blast b/*OnPBGenome.filtered.blast c/*OnPBGenome.filtered.blast d/*OnPBGenome.filtered.blast e/*OnPBGenome.filtered.blast' >> pncat

echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Pnalllongestcds-Mz_Genome.filtered.blast' > pncatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Pnalllongestcds-Ab_Genome.filtered.blast' >> pncatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Pnalllongestcds-Nb_Genome.filtered.blast' >> pncatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Pnalllongestcds-On_Genome.filtered.blast' >> pncatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Pn_alllongestcds-MzPacBio_Genome.filtered.blast' >> pncatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Pn_alllongestcds-OnPacBio_Genome.filtered.blast' >> pncatout

echo 'cat a/*MzGenome.filtered.blast b/*MzGenome.filtered.blast c/*MzGenome.filtered.blast d/*MzGenome.filtered.blast e/*MzGenome.filtered.blast' > abcat
echo 'cat a/*PnGenome.filtered.blast b/*PnGenome.filtered.blast c/*PnGenome.filtered.blast d/*PnGenome.filtered.blast e/*PnGenome.filtered.blast' >> abcat
echo 'cat a/*NbGenome.filtered.blast b/*NbGenome.filtered.blast c/*NbGenome.filtered.blast d/*NbGenome.filtered.blast e/*NbGenome.filtered.blast' >> abcat
echo 'cat a/*OnGenome.filtered.blast b/*OnGenome.filtered.blast c/*OnGenome.filtered.blast d/*OnGenome.filtered.blast e/*OnGenome.filtered.blast' >> abcat
echo 'cat a/*MzPBGenome.filtered.blast b/*MzPBGenome.filtered.blast c/*MzPBGenome.filtered.blast d/*MzPBGenome.filtered.blast e/*MzPBGenome.filtered.blast' >> abcat
echo 'cat a/*OnPBGenome.filtered.blast b/*OnPBGenome.filtered.blast c/*OnPBGenome.filtered.blast d/*OnPBGenome.filtered.blast e/*OnPBGenome.filtered.blast' >> abcat

echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Aballlongestcds-Mz_Genome.filtered.blast' > abcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Aballlongestcds-Pn_Genome.filtered.blast' >> abcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Aballlongestcds-Nb_Genome.filtered.blast' >> abcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Aballlongestcds-On_Genome.filtered.blast' >> abcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Ab_alllongestcds-MzPacBio_Genome.filtered.blast' >> abcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Ab_alllongestcds-OnPacBio_Genome.filtered.blast' >> abcatout

echo 'cat a/*MzGenome.filtered.blast b/*MzGenome.filtered.blast c/*MzGenome.filtered.blast d/*MzGenome.filtered.blast e/*MzGenome.filtered.blast' > nbcat
echo 'cat a/*PnGenome.filtered.blast b/*PnGenome.filtered.blast c/*PnGenome.filtered.blast d/*PnGenome.filtered.blast e/*PnGenome.filtered.blast' >> nbcat
echo 'cat a/*AbGenome.filtered.blast b/*AbGenome.filtered.blast c/*AbGenome.filtered.blast d/*AbGenome.filtered.blast e/*AbGenome.filtered.blast' >> nbcat
echo 'cat a/*OnGenome.filtered.blast b/*OnGenome.filtered.blast c/*OnGenome.filtered.blast d/*OnGenome.filtered.blast e/*OnGenome.filtered.blast' >> nbcat
echo 'cat a/*MzPBGenome.filtered.blast b/*MzPBGenome.filtered.blast c/*MzPBGenome.filtered.blast d/*MzPBGenome.filtered.blast e/*MzPBGenome.filtered.blast' >> nbcat
echo 'cat a/*OnPBGenome.filtered.blast b/*OnPBGenome.filtered.blast c/*OnPBGenome.filtered.blast d/*OnPBGenome.filtered.blast e/*OnPBGenome.filtered.blast' >> nbcat

echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Nballlongestcds-Mz_Genome.filtered.blast' > nbcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Nballlongestcds-Pn_Genome.filtered.blast' >> nbcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Nballlongestcds-Ab_Genome.filtered.blast' >> nbcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Nballlongestcds-On_Genome.filtered.blast' >> nbcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Nb_alllongestcds-MzPacBio_Genome.filtered.blast' >> nbcatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Nb_alllongestcds-OnPacBio_Genome.filtered.blast' >> nbcatout

echo 'cat a/*MzGenome.filtered.blast b/*MzGenome.filtered.blast c/*MzGenome.filtered.blast d/*MzGenome.filtered.blast e/*MzGenome.filtered.blast' > oncat
echo 'cat a/*PnGenome.filtered.blast b/*PnGenome.filtered.blast c/*PnGenome.filtered.blast d/*PnGenome.filtered.blast e/*PnGenome.filtered.blast' >> oncat
echo 'cat a/*AbGenome.filtered.blast b/*AbGenome.filtered.blast c/*AbGenome.filtered.blast d/*AbGenome.filtered.blast e/*AbGenome.filtered.blast' >> oncat
echo 'cat a/*NbGenome.filtered.blast b/*NbGenome.filtered.blast c/*NbGenome.filtered.blast d/*NbGenome.filtered.blast e/*NbGenome.filtered.blast' >> oncat
echo 'cat a/*MzPBGenome.filtered.blast b/*MzPBGenome.filtered.blast c/*MzPBGenome.filtered.blast d/*MzPBGenome.filtered.blast e/*MzPBGenome.filtered.blast' >> oncat

echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Onalllongestcds-Mz_Genome.filtered.blast' > oncatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Onalllongestcds-Pn_Genome.filtered.blast' >> oncatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Onalllongestcds-Ab_Genome.filtered.blast' >> oncatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/Onalllongestcds-Nb_Genome.filtered.blast' >> oncatout
echo '/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/On_alllongestcds-MzPacBio_Genome.filtered.blast' >> oncatout

sh createcatcommands.sh

nano mergeblastfiles-mz.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-4
#SBATCH --mem-per-cpu 48000
#SBATCH -t 0-5:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t mzcat < mzcat # assign cat commands as elements to variable
mapfile -t mzcatout < mzcatout # assing out files as variable

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/mz
${mzcat[${SLURM_ARRAY_TASK_ID}]} || true >> ${mzcatout[${SLURM_ARRAY_TASK_ID}]}


nano mergeblastfiles-pn.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-5
#SBATCH --mem-per-cpu 48000
#SBATCH -t 0-5:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t pncat < pncat # assign cat commands as elements to variable
mapfile -t pncatout < pncatout # assing out files as variable

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/pn
${pncat[${SLURM_ARRAY_TASK_ID}]} || true  >> ${pncatout[${SLURM_ARRAY_TASK_ID}]}


nano mergeblastfiles-ab.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-5
#SBATCH --mem-per-cpu 48000
#SBATCH -t 0-5:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t abcat < abcat # assign cat commands as elements to variable
mapfile -t abcatout < abcatout # assing out files as variable

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/ab
${abcat[${SLURM_ARRAY_TASK_ID}]} || true  >> ${abcatout[${SLURM_ARRAY_TASK_ID}]}


nano mergeblastfiles-nb.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-5
#SBATCH --mem-per-cpu 48000
#SBATCH -t 0-5:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t nbcat < nbcat # assign cat commands as elements to variable
mapfile -t nbcatout < nbcatout # assing out files as variable

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/nb
${nbcat[${SLURM_ARRAY_TASK_ID}]} || true  >> ${nbcatout[${SLURM_ARRAY_TASK_ID}]}


nano mergeblastfiles-on.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-4
#SBATCH --mem-per-cpu 48000
#SBATCH -t 0-5:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

mapfile -t oncat < oncat # assign cat commands as elements to variable
mapfile -t oncatout < oncatout # assing out files as variable

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/fastasplit/on
${oncat[${SLURM_ARRAY_TASK_ID}]} >> ${oncatout[${SLURM_ARRAY_TASK_ID}]}

#run all the above
sbatch mergeblastfiles-mz.sh
sbatch mergeblastfiles-pn.sh
sbatch mergeblastfiles-ab.sh
sbatch mergeblastfiles-nb.sh
sbatch mergeblastfiles-on.sh

rm slurm.t*

# 3. then filter the above to only focus on OGID-based hits where OGID is non1to1

# 3a. create an OGID file that only has the non-1-to-1 orthologs and species specific
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/
cp ../Module_genesandexpr/OGIDS.txt5 .
awk '($2 == "NULL" || $3 == "NULL" || $4 == "NULL" || $5 == "NULL" || $6 == "NULL")' OGIDS.txt5 > OGIDS.no1to1.txt5
awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' OGIDS.no1to1.txt5 > OGIDS.no1to1.txt5-mz
awk '{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' OGIDS.no1to1.txt5 > OGIDS.no1to1.txt5-pn
awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' OGIDS.no1to1.txt5 > OGIDS.no1to1.txt5-ab
awk '{print $5,$1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' OGIDS.no1to1.txt5 > OGIDS.no1to1.txt5-nb
awk '{print $6,$1,$2,$3,$4,$5,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' OGIDS.no1to1.txt5 > OGIDS.no1to1.txt5-on

# 3b. filter the filtered blast outputs for genes that are non1to1 only and add OGID as last column
for i in Mz*-*_Genome.filtered.blast ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' OGIDS.no1to1.txt5-mz $i | grep -v 'REMOVEME' | cut -f1-7,9 > "$(basename "$i" .filtered.blast).non1to1.filtered.blast" ;
done

for i in Pn*-*_Genome.filtered.blast ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' OGIDS.no1to1.txt5-pn $i | grep -v 'REMOVEME' | cut -f1-7,9 > "$(basename "$i" .filtered.blast).non1to1.filtered.blast" ;
done

for i in Ab*-*_Genome.filtered.blast ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' OGIDS.no1to1.txt5-ab $i | grep -v 'REMOVEME' | cut -f1-7,9 > "$(basename "$i" .filtered.blast).non1to1.filtered.blast" ;
done

for i in Nb*-*_Genome.filtered.blast ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' OGIDS.no1to1.txt5-nb $i | grep -v 'REMOVEME' | cut -f1-7,9 > "$(basename "$i" .filtered.blast).non1to1.filtered.blast" ;
done

for i in On*-*_Genome.filtered.blast ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' OGIDS.no1to1.txt5-on $i | grep -v 'REMOVEME' | cut -f1-7,9 > "$(basename "$i" .filtered.blast).non1to1.filtered.blast" ;
done

# Output is:
# 0. gene (query)
# 1. chr (hit)
# 2. number of BLAST hits
# 3. query length
# 4. total nt coverage of hits
# 5. query covered
# 6. Proportion of CDS (query) covered by BLAST hits
# 7. OGID

# 3c. Then, output the OGIDS where a gene looks like it is present in the genome where originally thought absent based on blasting with orthologs from all other species
# the first check for a missing orthogroup should be with the closest species, then, if not found, the next closest and so on
# absent in mz, check vs pn, then ab, nb, on
# absent in pn, check vs mz, then ab, nb, on
# absent in ab, check vs pn, then mz, nb, on
# absent in nb, check vs on, then ab, pn, mz
# absent in on, check vs nb, then ab, pn, mz

# doing a test
# python tblastx_parser2.py Mz_alllongestcds-Pn_Genome.blast Mz_alllongestcds.fasta
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' OGIDS.no1to1.txt5-mz Mz_alllongestcds-Pn_Genome.filtered.blast | grep -v 'REMOVEME' | cut -f1-7,9 > test.blastfilt

# NOTE: created a python script to do the above matching and retaining of OGIDS that are actually present in each genome: Arboretum_GT_v3/OGIDtblastx/tblastx_genePA-OGID.py
# run the python script:
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/

nano tblastx_genePA-OGID.sh

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

python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-mz Pnalllongestcds-Mz_Genome.non1to1.filtered.blast MzGenome-PnBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-mz Aballlongestcds-Mz_Genome.non1to1.filtered.blast MzGenome-AbBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-mz Nballlongestcds-Mz_Genome.non1to1.filtered.blast MzGenome-NbBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-mz Onalllongestcds-Mz_Genome.non1to1.filtered.blast MzGenome-OnBLAST_PresentNULLOGIDS.txt

python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-pn Mzalllongestcds-Pn_Genome.non1to1.filtered.blast PnGenome-MzBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-pn Aballlongestcds-Pn_Genome.non1to1.filtered.blast PnGenome-AbBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-pn Nballlongestcds-Pn_Genome.non1to1.filtered.blast PnGenome-NbBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-pn Onalllongestcds-Pn_Genome.non1to1.filtered.blast PnGenome-OnBLAST_PresentNULLOGIDS.txt

python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-ab Pnalllongestcds-Ab_Genome.non1to1.filtered.blast AbGenome-PnBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-ab Mzalllongestcds-Ab_Genome.non1to1.filtered.blast AbGenome-MzBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-ab Nballlongestcds-Ab_Genome.non1to1.filtered.blast AbGenome-NbBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-ab Onalllongestcds-Ab_Genome.non1to1.filtered.blast AbGenome-OnBLAST_PresentNULLOGIDS.txt

python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-nb Onalllongestcds-Nb_Genome.non1to1.filtered.blast NbGenome-OnBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-nb Aballlongestcds-Nb_Genome.non1to1.filtered.blast NbGenome-AbBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-nb Pnalllongestcds-Nb_Genome.non1to1.filtered.blast NbGenome-PnBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-nb Mzalllongestcds-Nb_Genome.non1to1.filtered.blast NbGenome-MzBLAST_PresentNULLOGIDS.txt

python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-on Nballlongestcds-On_Genome.non1to1.filtered.blast OnGenome-NbBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-on Aballlongestcds-On_Genome.non1to1.filtered.blast OnGenome-AbBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-on Pnalllongestcds-On_Genome.non1to1.filtered.blast OnGenome-PnBLAST_PresentNULLOGIDS.txt
python tblastx_genePA-OGID.py OGIDS.no1to1.txt5-on Mzalllongestcds-On_Genome.non1to1.filtered.blast OnGenome-MzBLAST_PresentNULLOGIDS.txt

cat MzGenome-*BLAST_PresentNULLOGIDS.txt | sort -u > MzGenome-PnAbNbOnBLAST_PresentNULLOGIDS.txt
cat PnGenome-*BLAST_PresentNULLOGIDS.txt | sort -u > PnGenome-MzAbNbOnBLAST_PresentNULLOGIDS.txt
cat AbGenome-*BLAST_PresentNULLOGIDS.txt | sort -u > AbGenome-MzPnNbOnBLAST_PresentNULLOGIDS.txt
cat NbGenome-*BLAST_PresentNULLOGIDS.txt | sort -u > NbGenome-MzPnAbOnBLAST_PresentNULLOGIDS.txt
cat OnGenome-*BLAST_PresentNULLOGIDS.txt | sort -u > OnGenome-MzPnAbNbBLAST_PresentNULLOGIDS.txt

#run the above
sbatch tblastx_genePA-OGID.sh


## TOTAL non1to1OGIDS TO FILTER
# Out of the 18799 Orthogroups, 11955 are non1to1
# Out of the 11955 non1to1 OGs, the following number are actually present in each genome:
# Mz
wc -l MzGenome-PnAbNbOnBLAST_PresentNULLOGIDS.txt # 876
# Pn
wc -l PnGenome-MzAbNbOnBLAST_PresentNULLOGIDS.txt # 1397
# Ab
wc -l AbGenome-MzPnNbOnBLAST_PresentNULLOGIDS.txt # 1185
# Nb
wc -l NbGenome-MzPnAbOnBLAST_PresentNULLOGIDS.txt # 1730
# On
wc -l OnGenome-MzPnAbNbBLAST_PresentNULLOGIDS.txt # 1084

cat MzGenome-PnAbNbOnBLAST_PresentNULLOGIDS.txt PnGenome-MzAbNbOnBLAST_PresentNULLOGIDS.txt AbGenome-MzPnNbOnBLAST_PresentNULLOGIDS.txt NbGenome-MzPnAbOnBLAST_PresentNULLOGIDS.txt OnGenome-MzPnAbNbBLAST_PresentNULLOGIDS.txt | sort -u | wc -l # 4387
# In total, there are 4387 OGIDs where in one genome, there is a NULL, but the gene is actually present in the genome
cat MzGenome-PnAbNbOnBLAST_PresentNULLOGIDS.txt PnGenome-MzAbNbOnBLAST_PresentNULLOGIDS.txt AbGenome-MzPnNbOnBLAST_PresentNULLOGIDS.txt NbGenome-MzPnAbOnBLAST_PresentNULLOGIDS.txt OnGenome-MzPnAbNbBLAST_PresentNULLOGIDS.txt | sort -u > MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS.txt # 4387



## C. Create matrices for the 1to1TF and non1to1TG edges.
# FIRST, filter the simplified edge table for all NULLOGIDs present in the genome
# For each of the candidate genes, define in which species the comparison is made, based on the presence of TG.
# For example if TG found in Ab and Nb only, then TFs only need to be 1-to-1 in those species for that edge and networks created accordingly.
# create a python script to do this

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/

# 0a. filter the simplified edge table for NULL OGIDs present in the genome
# this will be for all NULL OGIDs that ARE NOT TFs - so here, create a unified list but only remove from the simplified edge table if it is TG (col6) ONLY (not a TF, col3)

cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/OGIDtblastx/MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS.txt .
wc -l Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.txt # in the original files there are 1131812 lines
# check, how many of the NULL_OGIDs_PresentInGenome are candidates
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}a[$1]{print $0,a[$1]}' MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS.txt ../Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt
# 27 in total including rho, sws1 and sws2 - many of these genes are those that we are interested and hence, we will not remove them
# check how many TFs are removed
awk '{print $9,$6}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.mz | sort -u | awk '$1!="NA"' > TF_OGIDs
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}a[$1]{print $0,a[$1]}' MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS.txt TF_OGIDs
# 163 in the TF set - retain those as TFs DO NOT NEED TO BE REMOVED
# filter your MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS.txt for OGIDS that are candidates and TFs
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"KEEPME";}}' ../Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS.txt | grep 'KEEPME' | cut -f1 > MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS-noCand.txt1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"KEEPME";}}' TF_OGIDs MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS-noCand.txt1 | grep 'KEEPME' | cut -f1 > MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS-noCand-noTF.txt2
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS* /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/
# AFTER FILTERING FOR CANDIDATES AND TFs, THIS LEAVES 4209 PRESENT NULL OGIDS
# Finally, filter the simplified edge matrix using the new PresentNULLOGIDS filtered for candidates and TFs
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$6]){print $0,a[$6];}else{print $0,"KEEPME";}}' MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS-noCand-noTF.txt2 Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.txt | grep 'KEEPME' | cut -f1-48 > Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt # 907418/1131812 lines retained
# final file to use - Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt


# 0b. create an OGID presence/absence as 1/0 table in bash
awk '{if($2 == "NULL")print $0,$17=0;else print $0,$17=1;}' OFS='\t' OGIDS.txt5 | awk '{if($3 == "NULL")print $0,$18=0;else print $0,$18=1;}' OFS='\t' | awk '{if($4 == "NULL")print $0,$19=0;else print $0,$19=1;}' OFS='\t' | awk '{if($5 == "NULL")print $0,$20=0;else print $0,$20=1;}' OFS='\t' | awk '{if($6 == "NULL")print $0,$21=0;else print $0,$21=1;}' OFS='\t' | cut -f1,17-21 > OGIDS.txt5.PA

# 1. read in the OGIDs file - store each species OGID and geneID presence/absence as 1/0 in species dictionaries
# 2. read each line of the edge file, and:
  # a. assess which geneB (TG) have 1/0, then
  # b. for every 1 in geneB (TG), check if the corresponding species geneB (TG) is present (1) in OGID, if so
  # c. output the line

nano createtftglines-sbatch.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 0-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

python createtftglines.py


nano createtftglines.py

#!/usr/bin/env python

from __future__ import division
import time

# 1. Read in yours OGIDs presence/absence file and create a dictionary e.g. OGIDid 1 0 1 1 0
OGID = []
Mz = []
Pn = []
Ab = []
Nb = []
On = []
# for line in open("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5"):
for line in open("OGIDS.txt5.PA"):
    fields = line.rstrip("\n").split("\t")
    OGID.append(fields[0])
    Mz.append(fields[1])
    Pn.append(fields[2])
    Ab.append(fields[3])
    Nb.append(fields[4])
    On.append(fields[5])

# 1d. Create a merged dictionary of all species OGID presence/absence (key = OGID, values = presence/absence of OGID in species e.g. 1 0 1 1 0)
OGID_PA=dict(zip(OGID, zip(Mz,Pn,Ab,Nb,On)))


# # 1. Turn your original OGIDs file into a presence/absence matrix dictionary e.g. OGIDid 1 0 1 1 0
# OGID = []
# Mz = []
# Pn = []
# Ab = []
# Nb = []
# On = []
# # for line in open("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5"):
# for line in open("OGIDS.txt5"):
#     fields = line.rstrip("\n").split("\t")
#     OGID.append(fields[0])
#     Mz.append(fields[1])
#     Pn.append(fields[2])
#     Ab.append(fields[3])
#     Nb.append(fields[4])
#     On.append(fields[5])
#
# # 1a. define a function for creating the lists of presence/absence
# def create_list(sp,list): # provide the species list from above and an empty list name
#     list = [] # create an empty list in each case
#     for line in sp: # read line by line from the species list
#         if line == "NULL":
#             # if no gene present, then append '0' value for absence
#             list.append(0)
#         else:
#             # otherwise, append a '1' for the presence of the gene as an array in dict
#             list.append(1)
#     return list
#
# # 1b. call the function for creating the lists for each species
# Mz2 = create_list(Mz,list)
# Pn2 = create_list(Pn,list)
# Ab2 = create_list(Ab,list)
# Nb2 = create_list(Nb,list)
# On2 = create_list(On,list)
#
# # 1c. Create dictionaries from 3 lists (key = OGID, value1 = gene, value2 = presence/absence) for each species - probably won't need any of these
# MzOGID=dict(zip(OGID, zip(Mz, Mz2)))
# PnOGID=dict(zip(OGID, zip(Pn, Pn2)))
# AbOGID=dict(zip(OGID, zip(Ab, Ab2)))
# NbOGID=dict(zip(OGID, zip(Nb, Nb2)))
# OnOGID=dict(zip(OGID, zip(On, On2)))
#
# # 1d. Create a merged dictionary of all species OGID presence/absence (key = OGID, values = presence/absence of OGID in species e.g. 1 0 1 1 0)
# OGID_PA=dict(zip(OGID, zip(Mz2,Pn2,Ab2,Nb2,On2)))


# 2. For each line of TG in your edge table, pull out the presence/absence columns of that edge, have tf:tg as key and 1/0 as five values - col13-17 in edge table > store as dict
tftg = []
tftg_mz_pa = []
tftg_pn_pa = []
tftg_ab_pa = []
tftg_nb_pa = []
tftg_on_pa = []
with open("Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt") as f:
    next(f) # skip the first row that is a header
    for line in f:
        edges = line.rstrip("\n").split("\t")
        tftg.append(edges[0])
        tftg_mz_pa.append(edges[13])
        tftg_pn_pa.append(edges[14])
        tftg_ab_pa.append(edges[15])
        tftg_nb_pa.append(edges[16])
        tftg_on_pa.append(edges[17])

# 2a. Create a merged dictionary of all tf:tg presence/absence (key = OGID, values = presence/absence of OGID in species e.g. 1 0 1 1 0)
tftg_PA=dict(zip(tftg, zip(tftg_mz_pa,tftg_pn_pa,tftg_ab_pa,tftg_nb_pa,tftg_on_pa)))

# 3. Then, map the TF in each line to OGID presence/absence, again with tf:tg as key and 1/0 as five values - col1-5 in OGID > store as dict

# 3a. create a list of the TFs of each edge
TF_tftg = [] # this will store the TF OGID
TF_tftg2 = [] # this will store the tftg OGID
with open("Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt") as f:
    next(f) # skip the first row that is a header
    for line in f:
        edges = line.rstrip("\n").split("\t")
        TF_tftg.append(edges[2])
        TF_tftg2.append(edges[0])


# 3b. map the TF of edge list to the OGID presence/absence dict and output another dict
TF_tftgOGID = [] # create an empty list to store the matched TF OGID
TF_tftgOGID_PA = [] # create an empty list to store the matched TF OGID presence/absence values
for line in TF_tftg: # read the TFs from tf:tg line by line
    for OGID, values in OGID_PA.items(): # for each OG (key) and its values (presence/absence) in the OGID presence/absence dict
        if line == OGID: # if the TF OGID matches the OGID from presence/absence
            TF_tftgOGID.append(OGID) # append the OGID to this list
            TF_tftgOGID_PA.append(values) # append the presence/absence values to another list
            # print(OGID_PA[OGID])

# ensure to apply the zip in this way so that this dict has the same format as 'tftg_PA'
TF_tftgOGID_PAdict=dict(zip(TF_tftg2,TF_tftgOGID_PA)) # this creates a dict of key: tftg and values : Presence/Absence of TF in tftg in 5 species

# 4. Then check for coherence of 1's as you iterate between each value in both dicts using the keys as common:
    # For this, iterate through each value of the tf:tg of tg dict and for each 1 in here, there needs to be a 1 in the corresponding position in the tf:tg of tf OGID dict
    # e.g. tf:tg of tf OGID: 1 0 1 1 0 AND tf:tg of tg dict: 0 1 1 1 0 > Dropped as no TF for Pn when tf:tg exists
    # e.g. tf:tg of tf OGID: 1 1 1 1 0 AND tf:tg of tg dict: 0 1 1 1 0 > Kept as TF exists for all species where tf:tg exists

# compare the two dicts: tftg_PA vs TF_tftgOGID_PAdict - can compare directly since the order is the same

valelement = [0,1,2,3,4] # create a list of the elements in the dict values to iterate over
finaltftg = [] # 5. retain the tftg lines based on the logic above in this list
for key, value in tftg_PA.items():
    for x in valelement:
        if value[x] == '1' and value[x] == TF_tftgOGID_PAdict[key][x]: # If the tftg present and TF present
            print(key, value,"pass-PT1",TF_tftgOGID_PAdict[key][x])
            finaltftg.append(key)
        elif value[x] == '0' and value[x] == TF_tftgOGID_PAdict[key][x]: # If the tftg absent and TF absent
            print(key, value,"pass-PT2",TF_tftgOGID_PAdict[key][x])
            finaltftg.append(key)
        elif value[x] == '0' and value[x] != TF_tftgOGID_PAdict[key][x]: # If the tftg absent but TF present
            print(key, value,"pass-PT3",TF_tftgOGID_PAdict[key][x])
            finaltftg.append(key)
        else:
            print(key, value,"fail",TF_tftgOGID_PAdict[key][x]) # all other cases fail e.g. If the tftg present and TF absent

# keep only unique lines
finaltftgunique = set(finaltftg)
# write out to file
finaltftg_lines = open("finaltftg_lines.txt", "w")
finaltftg_lines.write(str('\n'.join(finaltftgunique)))
finaltftg_lines.close() # close the file

# then, use the above file to select lines from your main edge table file in awk script
wc -l finaltftg_lines.txt # 843042
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' finaltftg_lines.txt Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | grep -v 'REMOVEME' | wc -l # 907418/907418 - this keeps all lines so just share the original file
# Final File - Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt
cp Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix

###### Run species-specific GO enrichment of geneA (TF) and geneB (TG) separately of 1to1 and non1to1 edges

cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/
mkdir GOenrichment_BPonly_alledges_and1to1
cd GOenrichment_BPonly_alledges_and1to1

nano mkGOINPUT.sh

#!/bin/sh

## 1to1 edges only
# geneA
awk '{if($5 == "1")print $2;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > mz
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 mz | cut -f2 | sed "1i Mz_TF_1to1edges" | tr '\n' '#' | sed 's/#/\t/' > Mz_GOINPUT # mz
echo "" >> Mz_GOINPUT
awk '{if($6 == "1")print $2;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > pn
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 pn | cut -f2 | sed "1i Pn_TF_1to1edges" | tr '\n' '#' | sed 's/#/\t/' > Pn_GOINPUT # pn
echo "" >> Pn_GOINPUT
awk '{if($7 == "1")print $2;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > ab
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$4;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 ab | cut -f2 | sed "1i Ab_TF_1to1edges" | tr '\n' '#' | sed 's/#/\t/' > Ab_GOINPUT # ab
echo "" >> Ab_GOINPUT
awk '{if($8 == "1")print $2;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > nb
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$5;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 nb | cut -f2 | sed "1i Nb_TF_1to1edges" | tr '\n' '#' | sed 's/#/\t/' > Nb_GOINPUT # nb
echo "" >> Nb_GOINPUT
awk '{if($9 == "1")print $2;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > on
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$6;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 on | cut -f2 | sed "1i On_TF_1to1edges" | tr '\n' '#' | sed 's/#/\t/' > On_GOINPUT # on
echo "" >> On_GOINPUT

#geneB
awk '{if($5 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > mz
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 mz | cut -f2 | sed "1i Mz_TG_1to1edges" | tr '\n' '#' | sed 's/#/\t/' >> Mz_GOINPUT # mz
echo "" >> Mz_GOINPUT
awk '{if($6 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > pn
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 pn | cut -f2 | sed "1i Pn_TG_1to1edges" | tr '\n' '#' | sed 's/#/\t/' >> Pn_GOINPUT # pn
echo "" >> Pn_GOINPUT
awk '{if($7 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > ab
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$4;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 ab | cut -f2 | sed "1i Ab_TG_1to1edges" | tr '\n' '#' | sed 's/#/\t/' >> Ab_GOINPUT # ab
echo "" >> Ab_GOINPUT
awk '{if($8 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > nb
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$5;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 nb | cut -f2 | sed "1i Nb_TG_1to1edges" | tr '\n' '#' | sed 's/#/\t/' >> Nb_GOINPUT # nb
echo "" >> Nb_GOINPUT
awk '{if($9 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/Edge_Attributes_Collated4c.coexpr_promONLY.matrix3-geneAB1to1OGID.txt | sort -u > on
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$6;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 on | cut -f2 | sed "1i On_TG_1to1edges" | tr '\n' '#' | sed 's/#/\t/' >> On_GOINPUT # on
echo "" >> On_GOINPUT

## all edges (1to1 TF in species comparisons where edge present and non1to1TG and confirmed absent from each genome)
# geneA
awk '{if($14 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > mz
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 mz | cut -f2 | sed "1i Mz_TF_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> Mz_GOINPUT # mz
echo "" >> Mz_GOINPUT
awk '{if($15 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > pn
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 pn | cut -f2 | sed "1i Pn_TF_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> Pn_GOINPUT # pn
echo "" >> Pn_GOINPUT
awk '{if($16 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > ab
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$4;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 ab | cut -f2 | sed "1i Ab_TF_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> Ab_GOINPUT # ab
echo "" >> Ab_GOINPUT
awk '{if($17 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > nb
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$5;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 nb | cut -f2 | sed "1i Nb_TF_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> Nb_GOINPUT # nb
echo "" >> Nb_GOINPUT
awk '{if($18 == "1")print $3;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > on
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$6;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 on | cut -f2 | sed "1i On_TF_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> On_GOINPUT # on
echo "" >> On_GOINPUT

#geneB
awk '{if($14 == "1")print $6;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > mz
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 mz | cut -f2 | sed "1i Mz_TG_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> Mz_GOINPUT # mz
echo "" >> Mz_GOINPUT
awk '{if($15 == "1")print $6;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > pn
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 pn | cut -f2 | sed "1i Pn_TG_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> Pn_GOINPUT # pn
echo "" >> Pn_GOINPUT
awk '{if($16 == "1")print $6;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > ab
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$4;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 ab | cut -f2 | sed "1i Ab_TG_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> Ab_GOINPUT # ab
echo "" >> Ab_GOINPUT
awk '{if($17 == "1")print $6;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > nb
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$5;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 nb | cut -f2 | sed "1i Nb_TG_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> Nb_GOINPUT # nb
echo "" >> Nb_GOINPUT
awk '{if($18 == "1")print $6;}' OFS='\t' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt | sort -u > on
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$6;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' ../OGIDS.txt5 on | cut -f2 | sed "1i On_TG_TF1to1TGnon1to1PresentNULLOGIDedges" | tr '\n' '#' | sed 's/#/\t/' >> On_GOINPUT # on
echo "" >> On_GOINPUT

# run the above
sh mkGOINPUT.sh

# Run GO enrichment
mkdir GOOUTPUT_BPonly

ml gcc
ml zlib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib2/

/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer Ab_GOINPUT /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Ab-genelist.txt /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Asbu_ 1 GOOUTPUT_BPonly/Ab_GOOUTPUT persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer Mz_GOINPUT /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Mz-genelist.txt /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Meze_ 1 GOOUTPUT_BPonly/Mz_GOOUTPUT persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer Pn_GOINPUT /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Pn-genelist.txt /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Puny_ 1 GOOUTPUT_BPonly/Pn_GOOUTPUT persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer Nb_GOINPUT /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Nb-genelist.txt /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Nebr_ 1 GOOUTPUT_BPonly/Nb_GOOUTPUT persg
/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer On_GOINPUT /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/On-genelist.txt /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Orni_ 1 GOOUTPUT_BPonly/On_GOOUTPUT persg

# filter output (qval < 0.05)
for i in GOOUTPUT_BPonly/*details.txt ; do awk '($4 < 0.05 )' $i > GOOUTPUT_BPonly/"$(basename "$i" .txt)_filtered.txt" ; done

# for plotting the GO (just doing based on log10 fold enrichment (FDR < 0.05), sorting based on p-value here) prepare the OUTPUT file so that matches will show the GO output on each line for each gene
for i in GOOUTPUT_BPonly/*_filtered.txt ; do
 sed 's/$/;/g' $i |
 awk '{ gsub(";", ";"$1";"$2";"$3";"$4";"$5";"$6";"$7";"$8";"$9"\n") } 1' |
 cut -f10 |
 awk '{ gsub(";", "\t") } 1' |
 sed '/^$/d' > GOOUTPUT_BPonly/"$(basename "$i" )2" ;
done

# remove first column then only display unique entries as final column displays no. of genes assigned to the term if you need it!
for i in GOOUTPUT_BPonly/*_filtered.txt2 ; do
 cut -f2-10 $i |
 sort -u | sort -k3,3 > GOOUTPUT_BPonly/"$(basename "$i" 2)3" ;
done

# collate same type edges for each species
grep -h _TF_1to1edges GOOUTPUT_BPonly/*_filtered.txt3 | sed 's/_TF_1to1edges//g' >> GOOUTPUT_BPonly/collated_TF_1to1edges_GOOUTPUT_details_filtered.txt3
grep -h _TG_1to1edges GOOUTPUT_BPonly/*_filtered.txt3 | sed 's/_TG_1to1edges//g' >> GOOUTPUT_BPonly/collated_TG_1to1edges_GOOUTPUT_details_filtered.txt3
grep -h _TF_TF1to1TGnon1to1PresentNULLOGIDedges GOOUTPUT_BPonly/*_filtered.txt3 | sed 's/_TF_TF1to1TGnon1to1PresentNULLOGIDedges//g' >> GOOUTPUT_BPonly/collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT_details_filtered.txt3
grep -h _TG_TF1to1TGnon1to1PresentNULLOGIDedges GOOUTPUT_BPonly/*_filtered.txt3 | sed 's/_TG_TF1to1TGnon1to1PresentNULLOGIDedges//g' >> GOOUTPUT_BPonly/collated_TG_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT_details_filtered.txt3


# mkdir /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix
# mkdir /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only
cd ../
cp -r GOenrichment_BPonly_alledges_and1to1 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/

# Copied all cat files above to local here: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Edge_Attributes/matrix/GOenrichment_BPonly_alledges_and1to1 #DONE
### then plotting /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Edge_Attributes/matrix/GOenrichment_BPonly_alledges_and1to1/TF-TGedges_GOenrichmentplots.R
### Created GO enrichment plots for separated comparisons (all 5 species, haplo and pairwise - geneA/B SW) > plots here: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/TFTGconsensusedges



#####################################################################################
# Prep finalised DyNet tables for manuscript:
#
# Extended Data Table S-R3A - All_five_species_node_noExpression.csv > Rename: Extended_Data_Table_S-R3A_TFTG1to1_DyNet_rewiring_allfive_species.tsv
# Extended Data Table S-R3B - DyNet_rewiring_allfive_species.csv > Rename: Extended_Data_Table_S-R3B_TFTGalledges_DyNet_rewiring_allfive_species.tsv
# Extended Data Table S-R3C - All_five_species_node_noExpression_CAND.tsv > Rename: Extended_Data_Table_S-R3C_TFTG1to1_DyNet_rewiring_allfive_species_candidates.tsv
# Extended Data Table S-R3D - DyNet_rewiring_allfive_species_candidates.csv > Rename: Extended_Data_Table_S-R3D_TFTGalledges_DyNet_rewiring_allfive_species_candidates.tsv
#
#####################################################################################

cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton

cand=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Candidates_IDs_Fast_Opsins_Hahn_SNPs.txt2)
ogids=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5)
geneN=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/geneNamesMapping.txt)
onetoonerew=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/1.1to1_files_shared/Result_files/All_five_species_node_noExpression.csv)
allrew=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/7.FINALedgetables_TFno1to1_TG1to1PresenceNULLOGIDS/Dynet_rewiring_08082018_cichlid_dataset/DyNet_rewiring_allfive_species.csv)
onetoonerewout=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/1.1to1_files_shared/Result_files/All_five_species_node_noExpression.tsv)
allrewout=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/7.FINALedgetables_TFno1to1_TG1to1PresenceNULLOGIDS/Dynet_rewiring_08082018_cichlid_dataset/DyNet_rewiring_allfive_species.tsv)

onetoonerewout2=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Extended_Data_Table_S-R3A_TFTG1to1_DyNet_rewiring_allfive_species.tsv)
allrewout2=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Extended_Data_Table_S-R3C_TFTGalledges_DyNet_rewiring_allfive_species.tsv)
onetoonerewoutcand=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Extended_Data_Table_S-R3B_TFTG1to1_DyNet_rewiring_allfive_species_candidates.tsv)
allrewoutcand=(/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Extended_Data_Table_S-R3D_TFTGalledges_DyNet_rewiring_allfive_species_candidates.tsv)


sed $'s/,/\t/g' $onetoonerew | sed 's/"//g' | sed $'s/ /_/g' | awk '{print $7,$1,$3,$4,$5,$6,$10,$2,$8,$9}' OFS='\t' > $onetoonerewout
sed $'s/,/\t/g' $allrew | sed 's/"//g' | sed $'s/ /_/g' | cut -f1,4-7,10,11,14,17,20,21,22 | awk '{print $7,$1,$3,$4,$5,$6,$10,$2,$8,$9}' OFS='\t' > $allrewout

head -1 $onetoonerewout | sed $'s/$/\tGene_Symbol/g' | sed 's/.csv//g' > out_colheader

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $ogids $onetoonerewout | cut -f1-10,14,20,22,25 | sed 's/_(1_of_many)//g' | awk '{OFS="\t"} {if ($12=="NULL" && $13!="NULL") $12=$13; print $0}' | awk '{OFS="\t"} {if ($12=="NULL" && $14!="NULL") $12=$14; print $0}' > XX
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;next}{if(a[$11]){print $0,a[$11];}else{print $0,"NULL";}}' $geneN XX | awk '{OFS="\t"} {if ($12=="NULL" && $15!="NULL") $12=$15; print $0}' | awk '{OFS="\t"} {if ($12~"si:dkey" && $15!="NULL") $12=$15; print $0}' | awk '$12=tolower($12)' OFS='\t' | awk '{OFS="\t"} {if ($12=="null") $12=$1; print $0}' | cut -f1-10,12 > XX2
cat out_colheader XX2 > $onetoonerewout2
rm XX XX2

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $ogids $allrewout | cut -f1-10,14,20,22,25 | sed 's/_(1_of_many)//g' | awk '{OFS="\t"} {if ($12=="NULL" && $13!="NULL") $12=$13; print $0}' | awk '{OFS="\t"} {if ($12=="NULL" && $14!="NULL") $12=$14; print $0}' > XX
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;next}{if(a[$11]){print $0,a[$11];}else{print $0,"NULL";}}' $geneN XX | awk '{OFS="\t"} {if ($12=="NULL" && $15!="NULL") $12=$15; print $0}' | awk '{OFS="\t"} {if ($12~"si:dkey" && $15!="NULL") $12=$15; print $0}' | awk '$12=tolower($12)' OFS='\t' | awk '{OFS="\t"} {if ($12=="null") $12=$1; print $0}' | cut -f1-10,12 > XX2
cat out_colheader XX2 > $allrewout2
rm XX XX2

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $cand $onetoonerewout2 | grep -v 'REMOVEME' | cut -f1-11 > XX
cat out_colheader XX > $onetoonerewoutcand
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}{if(a[$1]){print $0,a[$1];}else{print $0,"REMOVEME";}}' $cand $allrewout2 | grep -v 'REMOVEME' | cut -f1-11 > XX
cat out_colheader XX > $allrewoutcand
rm XX
