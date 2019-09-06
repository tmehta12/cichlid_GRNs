#!/bin/sh

##########################################################################################################################
#
# Network Reconstruction constrained by Arboretum output - 5 cichlids (O. niloticus, P. nyererei, A. burtoni, M .zebra, N. brichardi)
#
# Interaction inputs - miRNA targets, TFBSs, CNEs
#
# By Tarang Mehta, Earlham, UK
# Version 2.0a 2018 - Simplified
#
##########################################################################################################################

# No need to add species or module column to edge files, this information can come from node attributes
# After each stage input into MS Excel to check table structure and also data assignment to column headers are correct

##########################################################################################################################

# COLUMN DETAILS TO ADD TO EACH INTERACTION SET

#interaction_type
# Post-transcriptional_directed #miRNAs
# coexpression #TF-TG coexpression
# PD_directed #TFBSs, CNEs
# PD_undirected #TFBSs, CNEs

#effect
# unknown
# stimulation
# inhibition

#interaction - make sure this is not mixed up with other 'interaction' columns
# unknown
# stimulation
# inhibition

#directness
# direct
# indirect

#direction
# directed
# undirected

#layer
# Post-transcription_regulation #miRNAs
# Co-expression #TF-TG: PGG algorithm
# Transcriptional_regulation #TFBSs, CNEs

#source
# targetscan7[ts7]
# co-expression_TF-TG[coTFTG]
# promoter_motif[pm]
# CNE_motif[cm]
# CNE_Proximal[cp]
# aCNE_motif[cm]
# aCNE_Proximal[cp]

slurm
homeWD=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/)
cd $homeWD

##########################################################################################################################

# 0. Genes and expression in each module >--{DONE - RAN ON SCRATCH > COPIED LOCAL}--<
# cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr
mkdir $homeWD/Module_genesandexpr
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr

# No Pattern>Cluster conversion is required

# Extra files required:
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v2/full_orthologs_map3c.txt # . contains the published orthologs gene IDs for mapping used previously
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v2/Module_genesandexpr/list_k10 . # list of k10 module numbers
# /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v2/Module_genesandexpr/list_k12 - list of k12 module numbers
# /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v2/Module_genesandexpr/list_singletissue - list of single tissue module numbers
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v2/Module_genesandexpr/list_species . # list of species
#ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v2/Module_genesandexpr/list_species-longer # list of species slightly longer length
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt . # the most recent OGids file for the new prediction mode modules


# You will initially need to map module genes for each species separately to zebrafish orthologs (five lots for each module) - need to map Dr to OGIDS file using full_orthology file
# cut the relevant columns and run an awk to match and output columns plus NA when not available, assign variable to each as these will be used throughout the script
# output file - col1 cichlidgeneID, col2 module number, col3 DrEnsemblID

# First, download the module assignment and expression files
ssh software
mkdir /usr/users/ga004/mehtat/Cichlid_GRNs/Arboretum_GT_v3/
cd /usr/users/ga004/mehtat/Cichlid_GRNs/Arboretum_GT_v3/
mkdir Module_genesandexpr
cd Module_genesandexpr
wget  -r -nH --cut-dirs=4 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~ckoch/cichlids/prediction_mode_results/results/
cd ../
wget  -r -nH --cut-dirs=3 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~ckoch/cichlids/prediction_mode_results/relevant_ogs_new.txt
wget  -r -nH --cut-dirs=3 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~ckoch/cichlids/prediction_mode_results/relevant_ogs_old.txt
wget  -r -nH --cut-dirs=3 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~ckoch/cichlids/prediction_mode_results/OGIDS.txt
# all d/l files can be found in /usr/users/ga004/mehtat/Cichlid_GRNs/Arboretum_GT_v3/ so just create symbolic links to files in this folder
exit

# Just copy all files to the scratch folder as you will manipulate filenames
cp /usr/users/ga004/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/* .
OGIDS=OGIDS.txt
WD=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/)

# Map the OGIDS file to Dr symbols using full_orthology file
# already have below files where species gene placed first for awk mapping
awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' full_orthologs_map3c.txt > full_orthologs_map3c_mz.txt
fullorthologsmz=full_orthologs_map3c_mz.txt
awk '{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' full_orthologs_map3c.txt > full_orthologs_map3c_pn.txt
fullorthologspn=full_orthologs_map3c_pn.txt
awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' full_orthologs_map3c.txt > full_orthologs_map3c_ab.txt
fullorthologsab=full_orthologs_map3c_ab.txt
awk '{print $5,$1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' full_orthologs_map3c.txt > full_orthologs_map3c_nb.txt
fullorthologsnb=full_orthologs_map3c_nb.txt
awk '{print $6,$1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' full_orthologs_map3c.txt > full_orthologs_map3c_on.txt
fullorthologson=full_orthologs_map3c_on.txt
# process the OGIDs file - first make tab delimited, then change order for species specific awk match
OGIDS2=$WD$(echo $(basename $OGIDS) | sed -e 's/.txt/.txt2/')
On_OGIDS2=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-On.txt2/')
Ab_OGIDS2=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Ab.txt2/')
Mz_OGIDS2=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Mz.txt2/')
Pn_OGIDS2=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Pn.txt2/')
Nb_OGIDS2=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Nb.txt2/')
sed 's/,/\t/g' $OGIDS > $OGIDS2 ; awk '{print $2,$3,$4,$5,$6,$1}' OFS='\t' $OGIDS2 > $On_OGIDS2 ; awk '{print $3,$2,$4,$5,$6,$1}' OFS='\t' $OGIDS2 > $Ab_OGIDS2 ; awk '{print $4,$2,$3,$5,$6,$1}' OFS='\t' $OGIDS2 > $Mz_OGIDS2 ; awk '{print $5,$2,$3,$4,$6,$1}' OFS='\t' $OGIDS2 > $Pn_OGIDS2 ; awk '{print $6,$2,$3,$4,$5,$1}' OFS='\t' $OGIDS2 > $Nb_OGIDS2
# awk match and pull out Dr ensembl IDs for all species - note though, that each species orthology will be different to the full_orthologs file, so each species will have a shared/unique ID
# Note that 97 rows/18799 rows of Ensembl IDS are not equal
# Also, as missing genes are NULL, they will not match to the NA in full_orthologs file
Mz_OGIDS3=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Mz.txt3/')
Pn_OGIDS3=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Pn.txt3/')
Ab_OGIDS3=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Ab.txt3/')
Nb_OGIDS3=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Nb.txt3/')
On_OGIDS3=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-On.txt3/')
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$10;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL";}}' $fullorthologsmz $Mz_OGIDS2 | awk '{print $6,$1,$7}' OFS='\t' > $Mz_OGIDS3 # for the very first file that we paste we will take the orthogroup
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$10;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL";}}' $fullorthologspn $Pn_OGIDS2 | cut -f1,7 > $Pn_OGIDS3
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$10;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL";}}' $fullorthologsab $Ab_OGIDS2 | cut -f1,7 > $Ab_OGIDS3
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$10;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL";}}' $fullorthologsnb $Nb_OGIDS2 | cut -f1,7 > $Nb_OGIDS3
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$10;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL";}}' $fullorthologson $On_OGIDS2 | cut -f1,7 > $On_OGIDS3
# paste the files together
OGIDS3=$WD$(echo $(basename $OGIDS) | sed -e 's/.txt/.txt3/')
paste -d'\t' $Mz_OGIDS3 $Pn_OGIDS3 $Ab_OGIDS3 $Nb_OGIDS3 $On_OGIDS3 > $OGIDS3 ## THIS IS OGIDS FILE TO USE FOR BELOW col1-OGID,col2-mzID,col3-MzEns etc. - WILL CREATE A COMPLETER ONE FURTHER DOWN THAT WILL BE THE FINAL FILE

# these are generic files that will be used to prepare module genes and expression
cut -f2,3 $OGIDS3 > Mz-Dr_orth
cut -f4,5 $OGIDS3 > Pn-Dr_orth
cut -f6,7 $OGIDS3 > Ab-Dr_orth
cut -f8,9 $OGIDS3 > Nb-Dr_orth
cut -f10,11 $OGIDS3 > On-Dr_orth

# convert filenames to run while loops below - for cluster
rename Asbu_ Ab- Asbu_*
rename Meze_ Mz- Meze_*
rename Nebr_ Nb- Nebr_*
rename Orni_ On- Orni_*
rename Puny_ Pn- Puny_*

# convert filenames to run while loops - this is locally
#rename 's/Asbu_/Ab-/' k10/Asbu_*
#rename 's/Meze_/Mz-/' k10/Meze_*
#rename 's/Nebr_/Nb-/' k10/Nebr_*
#rename 's/Orni_/On-/' k10/Orni_*
#rename 's/Puny_/Pn-/' k10/Puny_*

# while loop to generate the mapped files that will be specific to the modules
while read F ; do awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $F-Dr_orth $F-speciesspecnames_clusterassign.txt > $F-speciesspecnames_clusterassign_Drmapped.txt ; done < list_species

# prepare files according to modules for other scripts e.g. miRNAs - this is not species-specific but module-specific
while read F ; do awk -v module="$F" '{if ($2 == module) print $1;}' *-speciesspecnames_clusterassign.txt > $F-modulegenes.txt ; done < list_k10 # need to add -v to pass the variable as module to awk

# Process the expression files for input into MRNETb
# delete last line from each file in each of the three folders as they have a 'dummy' row
for file in *-exprtab.txt; do sed '$d' "$file" > "$(basename "$file" .txt)_a.txt"; done

# add the module numbers to the expression files by mapping using awk
while read F ; do awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $F-speciesspecnames_clusterassign.txt $F-exprtab_a.txt > $F-expr-mod.txt ; done < list_species

# seperate the files by modules for each species
while read F ; do awk -v module="$F" '{if ($8 == module) print $0;}' Ab-expr-mod.txt > $F-Ab-expr-mod.txt ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($8 == module) print $0;}' Nb-expr-mod.txt > $F-Nb-expr-mod.txt ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($8 == module) print $0;}' Mz-expr-mod.txt > $F-Mz-expr-mod.txt ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($8 == module) print $0;}' Pn-expr-mod.txt > $F-Pn-expr-mod.txt ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($8 == module) print $0;}' On-expr-mod.txt > $F-On-expr-mod.txt ; done < list_k10 # need to add -v to pass the variable as module to awk

# remove the last column (module number) from each file
for file in *-*-expr-mod.txt; do cut -f1-7 "$file" > "$(basename "$file" .txt)_a.txt"; done

# Finally, create a more complete OGIDs file, like the full_orthology file
Mz_OGIDS4=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Mz.txt4/')
Pn_OGIDS4=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Pn.txt4/')
Ab_OGIDS4=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Ab.txt4/')
Nb_OGIDS4=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-Nb.txt4/')
On_OGIDS4=$WD$(echo $(basename $OGIDS2) | sed -e 's/.txt2/-On.txt4/')
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' $fullorthologsmz $Mz_OGIDS2 | awk '{print $6,$1,$13,$14,$15,$16,$17,$18,$19,$20,$21}' OFS='\t' > $Mz_OGIDS4 # take OGIDs ID for first file when pasting
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' $fullorthologspn $Pn_OGIDS2 | cut -f1,13-21 > $Pn_OGIDS4
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' $fullorthologsab $Ab_OGIDS2 | cut -f1,13-21 > $Ab_OGIDS4
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' $fullorthologsnb $Nb_OGIDS2 | cut -f1,13-21 > $Nb_OGIDS4
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' $fullorthologson $On_OGIDS2 | cut -f1,13-21 > $On_OGIDS4
# paste the files together
OGIDS4=$WD$(echo $(basename $OGIDS) | sed -e 's/.txt/.txt4/')
paste -d'\t' $Mz_OGIDS4 $Pn_OGIDS4 $Ab_OGIDS4 $Nb_OGIDS4 $On_OGIDS4 | sed 's/NA/NULL/g' > $OGIDS4

# Here you will be merging columns by cutting the comparison columns, removing all NULL from columns, then removing empty columns but retaining rows (even if empty), then adding NULL to empty rows and then taking the last column in the file (could even be the first) as this will be the evolutionary closest to the ortholog in teleost/human and then add to the main table.
cut -f3 $OGIDS4 > EnOl
cut -f13 $OGIDS4 > EnOl2
cut -f23 $OGIDS4 > EnOl3
cut -f33 $OGIDS4 > EnOl4
cut -f43 $OGIDS4 > EnOl5
paste -d'\t' EnOl EnOl2 EnOl3 EnOl4 EnOl5 | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3,$4,$5}' | awk '!NF{$0="NULL"}1' | awk '{print $NF}' > EnOlOGID

cut -f4 $OGIDS4 > EnTn
cut -f14 $OGIDS4 > EnTn2
cut -f24 $OGIDS4 > EnTn3
cut -f34 $OGIDS4 > EnTn4
cut -f44 $OGIDS4 > EnTn5
paste -d'\t' EnTn EnTn2 EnTn3 EnTn4 EnTn5 | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3,$4,$5}' | awk '!NF{$0="NULL"}1' | awk '{print $NF}' > EnTnOGID

cut -f5 $OGIDS4 > EnGa
cut -f15 $OGIDS4 > EnGa2
cut -f25 $OGIDS4 > EnGa3
cut -f35 $OGIDS4 > EnGa4
cut -f45 $OGIDS4 > EnGa5
paste -d'\t' EnGa EnGa2 EnGa3 EnGa4 EnGa5 | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3,$4,$5}' | awk '!NF{$0="NULL"}1' | awk '{print $NF}' > EnGaOGID

cut -f6 $OGIDS4 > EnDr
cut -f16 $OGIDS4 > EnDr2
cut -f26 $OGIDS4 > EnDr3
cut -f36 $OGIDS4 > EnDr4
cut -f46 $OGIDS4 > EnDr5
paste -d'\t' EnDr EnDr2 EnDr3 EnDr4 EnDr5 | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3,$4,$5}' | awk '!NF{$0="NULL"}1' | awk '{print $NF}' > EnDrOGID

cut -f7 $OGIDS4 > Ga
cut -f17 $OGIDS4 > Ga2
cut -f27 $OGIDS4 > Ga3
cut -f37 $OGIDS4 > Ga4
cut -f47 $OGIDS4 > Ga5
paste -d'\t' Ga Ga2 Ga3 Ga4 Ga5 | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3,$4,$5}' | awk '!NF{$0="NULL"}1' | awk '{print $NF}' > GaOGID

cut -f8 $OGIDS4 > Dr
cut -f18 $OGIDS4 > Dr2
cut -f28 $OGIDS4 > Dr3
cut -f38 $OGIDS4 > Dr4
cut -f48 $OGIDS4 > Dr5
paste -d'\t' Dr Dr2 Dr3 Dr4 Dr5 | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3,$4,$5}' | awk '!NF{$0="NULL"}1' | awk '{print $NF}' > DrOGID

cut -f9 $OGIDS4 > EnOn
cut -f19 $OGIDS4 > EnOn2
cut -f29 $OGIDS4 > EnOn3
cut -f39 $OGIDS4 > EnOn4
cut -f49 $OGIDS4 > EnOn5
paste -d'\t' EnOn EnOn2 EnOn3 EnOn4 EnOn5 | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3,$4,$5}' | awk '!NF{$0="NULL"}1' | awk '{print $NF}' > EnOnOGID

cut -f10 $OGIDS4 > EnHs
cut -f20 $OGIDS4 > EnHs2
cut -f30 $OGIDS4 > EnHs3
cut -f40 $OGIDS4 > EnHs4
cut -f50 $OGIDS4 > EnHs5
paste -d'\t' EnHs EnHs2 EnHs3 EnHs4 EnHs5 | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3,$4,$5}' | awk '!NF{$0="NULL"}1' | awk '{print $NF}' > EnHsOGID

cut -f11 $OGIDS4 > Hs
cut -f21 $OGIDS4 > Hs2
cut -f31 $OGIDS4 > Hs3
cut -f41 $OGIDS4 > Hs4
cut -f51 $OGIDS4 > Hs5
paste -d'\t' Hs Hs2 Hs3 Hs4 Hs5 | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3,$4,$5}' | awk '!NF{$0="NULL"}1' | awk '{print $NF}' > HsOGID

## Add mouse orthologs by using human ortholog IDs
# d/l all Human to Mouse Ensembl IDs from biomart > saved /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Hs-Mm_EnsemblOrthology.txt
perl -pe 's/\t/\tNULL/g' Hs-Mm_EnsemblOrthology.txt | sed 's/NULLENSMU/ENSMU/g' > Hs-Mm_EnsemblOrthology.txt2 # replace all empty spaces with NULL
# copied to scratch dir
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL";}}' Hs-Mm_EnsemblOrthology.txt2 EnHsOGID | cut -f2 > EnMmOGID
paste -d'\t' OGIDS.txt2 EnOlOGID EnTnOGID EnGaOGID GaOGID EnDrOGID DrOGID EnOnOGID EnHsOGID HsOGID EnMmOGID | awk '{print $1,$4,$5,$3,$6,$2,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' > OGIDS.txt5
OGIDS5=$WD$(echo $(basename $OGIDS) | sed -e 's/.txt/.txt5/')
paste -d'\t' $OGIDS2 EnOlOGID EnTnOGID EnGaOGID GaOGID EnDrOGID DrOGID EnOnOGID EnHsOGID HsOGID EnMmOGID | awk '{print $1,$4,$5,$3,$6,$2,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' > $OGIDS5 ##### THIS IS YOUR FINAL ORTHOLOGY FILE - OGIDS.txt5 #####

## Finally, leverage the gene names from the Nature publication for the ab gene (as it has most genes annotated) to have a comprehensive set
awk '{OFS="\t"} {if ($2=="NONE") $3="NONE"; print $0}' /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/protein_coding/v2/meta/geneNamesTree | sed 's/ /_/g' | sed 's/Q4SR56_TETNG/sws2/g' | sed 's/OPSG_ORYLA/rh2/g' | sed 's/Q2L6A1_ORYLA/sws1/g' > geneNamesTree2
# a. split OGIDs file into per species
OGIDS5ab=$WD$(echo $(basename $OGIDS) | sed -e 's/.txt/.abtxt5/')
OGIDS6=$WD$(echo $(basename $OGIDS) | sed -e 's/.txt/.txt6/')
awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' $OGIDS5 > $OGIDS5ab
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' geneNamesTree2 $OGIDS5ab | sed 's/NA/NULL/g' | awk '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' OFS='\t' > $OGIDS6

# $OGIDS6 structure
# 1. OGID
# 2. mz.gene
# 3. pn.gene
# 4. ab.gene
# 5. nb.gene
# 6. on.gene
# 7. Medaka Ensembl
# 8. Tetraodon Ensembl
# 9. Stickleback Ensembl
# 10. Stickleback gene symbol
# 11. Zebrafish Ensembl
# 12. Zebrafish gene symbol
# 13. Tilapia Ensembl
# 14. Human Ensembl
# 15. Human gene symbol
# 16. Mouse Ensembl
# 17. A. burtoni gene symbol

### { DONE } - RAN ON SCRATCH 21/03/17-22/03/17, added OGIDS.txt6 on 8/12/17

##########################################################################################################################

# 1. TFBSs - Motif Discovery {DONE}
# >--{Originally had several iterations of just using JASPAR scanned promoters - see NetworkReconstruction_v5.sh however, now we have a whole pipeline for cichlid-specific PSSMs}--<


mkdir /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2


# Will ran all motif scanning for various sources (if you need scripts, get off him)
# Species-specific promoter sequences used as 0-order markov backgrounds
# FIMO ran with default p-value (1e-4)
# All files are here: /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse AND /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/human
# NOTE: JASPAR runs in both folders are exactly the same - only process for one


# 1. Create variables for relevant folders and mapping files
WDpromtfbs=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2) #out to here
WDpromtfbs_fimo=(/tgac/workarea/Research-Groups/RG-cichlids/matQualOutput)
cd $WDpromtfbs
hspval=(/tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mat_qual_pvals_ALL.out2)
mmpval=(/tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mm10_mat_qual_pvals_ALL.out1)
OGIDab=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab)
OGIDmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz)
OGIDpn=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn)
OGIDnb=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb)
OGIDon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on)

# 1a. amend the p-vals files so that they are species-specific for awk matching # {DONE}
for i in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mat_qual_pvals_ALL.out2 /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mm10_mat_qual_pvals_ALL.out1; do
	awk '{print $4,$1,$2,$3,$5,$6,$7,$8}' OFS='\t' $i | sed 's/;/\t/' | sed 's/;/\t/' | sed 's/;/\t/' | sed 's/;/\t/' > $i.mz
	awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' $i.mz > $i.pn
	awk '{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' $i.mz > $i.ab
	awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' $i.mz > $i.nb
	awk '{print $5,$1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' $i.mz > $i.on
done
for i in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/*out1.* ; do mv $i . ; done
for i in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/*out2.* ; do mv $i . ; done

hspvalmz=(mat_qual_pvals_ALL.out2.mz)
mmpvalmz=(mm10_mat_qual_pvals_ALL.out1.mz)
hspvalpn=(mat_qual_pvals_ALL.out2.pn)
mmpvalpn=(mm10_mat_qual_pvals_ALL.out1.pn)
hspvalab=(mat_qual_pvals_ALL.out2.ab)
mmpvalab=(mm10_mat_qual_pvals_ALL.out1.ab)
hspvalnb=(mat_qual_pvals_ALL.out2.nb)
mmpvalnb=(mm10_mat_qual_pvals_ALL.out1.nb)
hspvalon=(mat_qual_pvals_ALL.out2.on)
mmpvalon=(mm10_mat_qual_pvals_ALL.out1.on)

# 1b. create merged Hs and Mm mapping files
hsmmpvalmz=(merged_mat_qual_pvals_ALL.out1.mz)
cat $hspvalmz $mmpvalmz | sort -u > $hsmmpvalmz
hsmmpvalpn=(merged_mat_qual_pvals_ALL.out1.pn)
cat $hspvalpn $mmpvalpn | sort -u > $hsmmpvalpn
hsmmpvalab=(merged_mat_qual_pvals_ALL.out1.ab)
cat $hspvalab $mmpvalab | sort -u > $hsmmpvalab
hsmmpvalnb=(merged_mat_qual_pvals_ALL.out1.nb)
cat $hspvalnb $mmpvalnb | sort -u > $hsmmpvalnb
hsmmpvalon=(merged_mat_qual_pvals_ALL.out1.on)
cat $hspvalon $mmpvalon | sort -u > $hsmmpvalon

# 1c. create simplified merged mapping files (TF_gene_symbol,cichlidIDs,motif_ID,Hs/Mm_EnsemblID,TF_OGID) - one unified and others species specific to use later

TFhsmmpval=(merged_mat_qual_pvals_ALL.out1.TF)
sed 's/;/\t/' merged_mat_qual_pvals_ALL.out1.mz | sed 's/;/\t/' | awk '{print $8,$1,$2,$3,$4,$5,$6,$9,$10}' OFS='\t' > $TFhsmmpval # unified file, TF_genesymbol first column
ENShsmmpval=(merged_mat_qual_pvals_ALL.out1.ENS)
sed 's/;/\t/' merged_mat_qual_pvals_ALL.out1.mz | sed 's/;/\t/' | awk '{print $9,$1,$2,$3,$4,$5,$6,$8,$10}' OFS='\t' > $ENShsmmpval # unified file, TF_genesymbol first column
# species-specific file - cichlid gene first column
TFhsmmpvalmz=(merged_mat_qual_pvals_ALL.out1.TF.mz)
sed 's/;/\t/' merged_mat_qual_pvals_ALL.out1.mz | sed 's/;/\t/' | awk '{print $1,$2,$3,$4,$5,$8,$6,$9,$10}' OFS='\t' > $TFhsmmpvalmz
TFhsmmpvalpn=(merged_mat_qual_pvals_ALL.out1.TF.pn)
sed 's/;/\t/' merged_mat_qual_pvals_ALL.out1.pn | sed 's/;/\t/' | awk '{print $1,$2,$3,$4,$5,$8,$6,$9,$10}' OFS='\t' > $TFhsmmpvalpn
TFhsmmpvalab=(merged_mat_qual_pvals_ALL.out1.TF.ab)
sed 's/;/\t/' merged_mat_qual_pvals_ALL.out1.ab | sed 's/;/\t/' | awk '{print $1,$2,$3,$4,$5,$8,$6,$9,$10}' OFS='\t' > $TFhsmmpvalab
TFhsmmpvalnb=(merged_mat_qual_pvals_ALL.out1.TF.nb)
sed 's/;/\t/' merged_mat_qual_pvals_ALL.out1.nb | sed 's/;/\t/' | awk '{print $1,$2,$3,$4,$5,$8,$6,$9,$10}' OFS='\t' > $TFhsmmpvalnb
TFhsmmpvalon=(merged_mat_qual_pvals_ALL.out1.TF.on)
sed 's/;/\t/' merged_mat_qual_pvals_ALL.out1.on | sed 's/;/\t/' | awk '{print $1,$2,$3,$4,$5,$8,$6,$9,$10}' OFS='\t' > $TFhsmmpvalon

# 2a. Create the edge tables by doing the following in order (reflected in loops)

# Table format to output:
# motif_pattern
# GeneA_OGID
# GeneA
# GeneA_Symbol
# GeneB_OGID
# GeneB
# GeneB_SymbolDr
# GeneB_SymbolGa
# GeneB_SymbolSp
# start
# stop
# strand
# score
# p-value
# q-value
# sequence
# conf_level
# conf_score

# 2ai. create unified (where relevant) mouse and human sets
# 2aii. Filter outputs (where relevant) for significant q-val (0.05)
# 2aiii. Create confidence levels for each scan
	# There are six confidence levels so scale this from 0 to 1 - however, do not give proportional waiting to level 2b as not all species will have this representation

# Confidence level 1
# Confidence level 1a - interactions in cichlids extrapolated from mouse that OVERLAP with interactions extrapolated from human: 0.3
# interactions in cichlids extrapolated from mouse that OVERLAP with interactions extrapolated from human

for i in /tgac/workarea/Research-Groups/RG-cichlids/GTRDdata/mouse_sites/mm10_tftgtfbs_*_07.extrap.annotations.blast_out ; do cp $i . ; done # create symbolic links to the mouse extrapolated files with coordinates here
for i in /tgac/workarea/Research-Groups/RG-cichlids/GTRDhuman/human_sites/hg38_tftgtfbs_*_07.extrap.annotations.blast_out ; do cp $i . ; done # create symbolic links to the human extrapolated files with coordinates here

# rename files to a unified format {DONE}
mv mm10_tftgtfbs_Ab_GeneAnnotation_11092017_07.extrap.annotations.blast_out ab-mm10_07.extrap.annotations.blast_out
mv mm10_tftgtfbs_Mz_GeneAnnotation_11092017_07.extrap.annotations.blast_out mz-mm10_07.extrap.annotations.blast_out
mv mm10_tftgtfbs_Nb_GeneAnnotation_11092017_07.extrap.annotations.blast_out nb-mm10_07.extrap.annotations.blast_out
mv mm10_tftgtfbs_Oreochromis_niloticus_07.extrap.annotations.blast_out on-mm10_07.extrap.annotations.blast_out
mv mm10_tftgtfbs_Pn_GeneAnnotation_11092017_07.extrap.annotations.blast_out pn-mm10_07.extrap.annotations.blast_out
mv hg38_tftgtfbs_Ab_GeneAnnotation_11092017_07.extrap.annotations.blast_out ab-hg38_07.extrap.annotations.blast_out
mv hg38_tftgtfbs_Mz_GeneAnnotation_11092017_07.extrap.annotations.blast_out mz-hg38_07.extrap.annotations.blast_out
mv hg38_tftgtfbs_Nb_GeneAnnotation_11092017_07.extrap.annotations.blast_out nb-hg38_07.extrap.annotations.blast_out
mv hg38_tftgtfbs_Oreochromis_niloticus_07.extrap.annotations.blast_out on-hg38_07.extrap.annotations.blast_out
mv hg38_tftgtfbs_Pn_GeneAnnotation_11092017_07.extrap.annotations.blast_out pn-hg38_07.extrap.annotations.blast_out

# find the intersect by catting, sorting and retaining duplicates {DONE}
for i in mz pn ab nb on ; do
	for j in ${i}-mm10_07.extrap.annotations.blast_out ; do
		for k in ${i}-hg38_07.extrap.annotations.blast_out ; do
			sort -u $j > $i-mm10_07.extrap.annotations.blast_out.nodup
			sort -u $k > $i-hg38_07.extrap.annotations.blast_out.nodup
			cat $i-mm10_07.extrap.annotations.blast_out.nodup $i-hg38_07.extrap.annotations.blast_out.nodup | sort | uniq -d |
			awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="1a",$11="0.3"}' OFS='\t' >> $i-mm10hg38_07.extrap.annotations.blast_out.intersect # find the duplicated lines
		done
	done
done

# Total number of intersected edges
# 891 ab-mm10hg38_07.extrap.annotations.blast_out.intersect
# 923 mz-mm10hg38_07.extrap.annotations.blast_out.intersect
# 461 nb-mm10hg38_07.extrap.annotations.blast_out.intersect
# 955 on-mm10hg38_07.extrap.annotations.blast_out.intersect
# 866 pn-mm10hg38_07.extrap.annotations.blast_out.intersect
# 4096 total

# Total number of intersected TFs
# 16 Ab
# 19 Mz
# 12 Nb
# 18 On
# 17 Pn

# Confidence level 1b - interactions in cichlids extrapolated from mouse ONLY: 0.2 {DONE}

for i in mz pn ab nb on ; do
	for j in ${i}-mm10hg38_07.extrap.annotations.blast_out.intersect ; do
		for k in ${i}-mm10_07.extrap.annotations.blast_out ; do
			cat $j $k | sort | uniq -u | # remove the intersected lines from the original mouse file
			awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="1b",$11="0.2"}' OFS='\t' > $i-mm10_07.extrap.annotations.blast_out.unique
		done
	done
done

# Confidence level 1c - interactions in cichlids extrapolated from human ONLY: 0.15 {DONE}

for i in mz pn ab nb on ; do
	for j in ${i}-mm10hg38_07.extrap.annotations.blast_out.intersect ; do
		for k in ${i}-hg38_07.extrap.annotations.blast_out ; do
			cat $j $k | sort | uniq -u | # remove the intersected lines from the original human file
			awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="1c",$11="0.15"}' OFS='\t' > $i-hg38_07.extrap.annotations.blast_out.unique
		done
	done
done

# merge the extrapolated files {DONE}
for i in mz pn ab nb on ; do
	cat $i-mm10hg38_07.extrap.annotations.blast_out.intersect $i-mm10_07.extrap.annotations.blast_out.unique $i-hg38_07.extrap.annotations.blast_out.unique >> $i-mm10hg38_07.extrap.annotations.blast_out.merged
done

# Confidence level 2 # {DONE}
# Confidence level 2a - FDR corrected interactions resulting from the FIMO scans using the extrapolated cichlid species specific matrices: 0.125 #{DONE}
nano mergeCSfimo.sh # {DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;

for i in mz pn ab nb on ; do
	for j in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse/CS_default_${i}_results.out ; do
		for k in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/human/CS_default_${i}_results.out ; do
			cat $j $k | # cat human and mosue files
			awk '($8 < 0.05 )' | # multiple test correction <0.05
			sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
			sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
			awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
			sort -k4,4 -rn | uniq -f 1 |
			sed $'s/;/\t/g' |
			awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="2a",$11="0.125"}' OFS='\t' >> $i-CS_default_fimo.out ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
			# cut -f1 $i-CS_default_fimo.out | sort -u > $i-CS_default_fimo.CUT.out # cut out merged first col and create no duplicates file
			# use cut file to match only the first hits in the original file which should take the top hit q-val
		done
	done
done

# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=1 mergeCSfimo.sh # {DONE}


# Confidence level 2b - FDR corrected interactions resulting from the FIMO scans using the extrapolated non-species specific matrices: 0.110 # {DONE}

for i in mz pn ab nb on ; do
	for j in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse/CW_default_${i}_results.out ; do
		for k in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/human/CW_default_${i}_results.out ; do
			cat $j $k | # cat human and mouse files
			awk '($8 < 0.05 )' | # multiple test correction <0.05
			sed 's/RG-cich_GTRDdata_mouse_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level.ig_1//g' |
			sed 's/RG-cich_GTRDhuman_human_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level//g' >> $i-CW_default_fimo.out1 ;
		done
	done
done

for i in mz pn ab nb on ; do
	for j in $i-CW_default_fimo.out1 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $ENShsmmpval $j |
		sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
		awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' | # merge the 5 cols, separated by colon
		sort -k4,4 -rn | uniq -f 1 |
		sed $'s/;/\t/g' |
		awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19="2b",$20="0.110"}' OFS='\t' >> $i-CW_default_fimo.out2 ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	done
done


# Confidence level 2c - FDR corrected interactions resulting from the FIMO scans using the Jaspar matrices: 0.115 # {DONE}
# this is done with TF mapping below


# 3. Map the first column (TF) to gene symbol/cichlid gene ID (where relevant - for 1a, 1b, 1c, 2a, 2b)

# Extrapolated TFs #{DONE}
for i in mz pn ab nb on ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' merged_mat_qual_pvals_ALL.out1.TF.$i $i-mm10hg38_07.extrap.annotations.blast_out.merged |
	grep -v 'REMOVEME' |
	awk '{print $2,$18,$20,$1,$17,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' > $i-mm10hg38_07.extrap.out1
done
rm *-mm10hg38_07.extrap.annotations.blast_out.intersect
rm *-mm10_07.extrap.annotations.blast_out.unique
rm *-hg38_07.extrap.annotations.blast_out.unique
rm *.nodup
rm *-mm10hg38_07.extrap.annotations.blast_out.merged

# CS TFs #{DONE}
for i in mz pn ab nb on ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' merged_mat_qual_pvals_ALL.out1.TF.$i $i-CS_default_fimo.out | awk '{print $2,$18,$20,$1,$17,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' > $i-CS_default_fimo.out1
done
rm *-CS_default_fimo.out

# CW TFs #{DONE}
# just select required columns from files created above, then remove old files
# these files are ready for TG mapping (as first column)
awk '{print $2,$16,$18,$11,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' mz-CW_default_fimo.out2 | awk '$5=tolower($5)' OFS='\t' > mz-CW_default_fimo.out3
awk '{print $2,$16,$18,$12,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' pn-CW_default_fimo.out2 | awk '$5=tolower($5)' OFS='\t' > pn-CW_default_fimo.out3
awk '{print $2,$16,$18,$13,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' ab-CW_default_fimo.out2 | awk '$5=tolower($5)' OFS='\t' > ab-CW_default_fimo.out3
awk '{print $2,$16,$18,$14,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' nb-CW_default_fimo.out2 | awk '$5=tolower($5)' OFS='\t' > nb-CW_default_fimo.out3
awk '{print $2,$16,$18,$15,$17,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' on-CW_default_fimo.out2 | awk '$5=tolower($5)' OFS='\t' > on-CW_default_fimo.out3
rm *-CW_default_fimo.out1
for i in mz pn ab nb on ; do
 mv $i-CW_default_fimo.out3 $i-CW_default_fimo.out1
done
rm *-CW_default_fimo.out2

# JASPAR TFs {DONE}
for i in mz pn ab nb on ; do
	for j in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/human/JASPAR_default_${i}_results.out ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' $TFhsmmpval $j |
		grep -v 'REMOVEME' |
		awk '($8 < 0.05 )' > $i-JASPARTF_interm.out1
	done
done

awk '{print $2,$16,$18,$11,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' mz-JASPARTF_interm.out1 |awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' > mz-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $2,$16,$18,$12,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' pn-JASPARTF_interm.out1 |awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' > pn-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $2,$16,$18,$13,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' ab-JASPARTF_interm.out1 |awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' > ab-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $2,$16,$18,$14,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' nb-JASPARTF_interm.out1 |awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' > nb-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $2,$16,$18,$15,$1,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' on-JASPARTF_interm.out1 |awk '{print $0,$13="2c",$14=0.115}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' > on-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)

rm *-JASPARTF_interm.out1

# and (at this point, all TF mapped files are *.out1)
# 4. Map the second column (TG) to gene symbol (run for all levels)- creating a file with the required columns of motif_ID,TF_OGID, TF_cichlid_ID, TF_gene_symbol,TG_OGID, TG_cichlid_ID, TG_Dr_symbol, TG_Ga_symbol

echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz > OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb >> OGIDs_list
echo /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on >> OGIDs_list

echo mz-mm10hg38_07.extrap.out1 > conf1a1b1c2a2bTG
echo pn-mm10hg38_07.extrap.out1 >> conf1a1b1c2a2bTG
echo ab-mm10hg38_07.extrap.out1 >> conf1a1b1c2a2bTG
echo nb-mm10hg38_07.extrap.out1 >> conf1a1b1c2a2bTG
echo on-mm10hg38_07.extrap.out1 >> conf1a1b1c2a2bTG
echo mz-CS_default_fimo.out1 >> conf1a1b1c2a2bTG
echo pn-CS_default_fimo.out1 >> conf1a1b1c2a2bTG
echo ab-CS_default_fimo.out1 >> conf1a1b1c2a2bTG
echo nb-CS_default_fimo.out1 >> conf1a1b1c2a2bTG
echo on-CS_default_fimo.out1 >> conf1a1b1c2a2bTG
echo mz-CW_default_fimo.out1 >> conf1a1b1c2a2bTG
echo pn-CW_default_fimo.out1 >> conf1a1b1c2a2bTG
echo ab-CW_default_fimo.out1 >> conf1a1b1c2a2bTG
echo nb-CW_default_fimo.out1 >> conf1a1b1c2a2bTG
echo on-CW_default_fimo.out1 >> conf1a1b1c2a2bTG
echo mz-JASPAR_default_fimo.out1 >> conf1a1b1c2a2bTG
echo pn-JASPAR_default_fimo.out1 >> conf1a1b1c2a2bTG
echo ab-JASPAR_default_fimo.out1 >> conf1a1b1c2a2bTG
echo nb-JASPAR_default_fimo.out1 >> conf1a1b1c2a2bTG
echo on-JASPAR_default_fimo.out1 >> conf1a1b1c2a2bTG

# Final Table format to output (below we are outputting the TG as first column for mapping in next step)
# motif_pattern
# GeneA_OGID
# GeneA
# GeneA_Symbol
# GeneB_OGID
# GeneB
# GeneB_SymbolDr
# GeneB_SymbolGa
# GeneB_SymbolSp
# start
# stop
# strand
# score
# p-value
# q-value
# sequence
# conf_level
# conf_score

while read -u 3 -r file1 && read -u 4 -r file2
do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' ${file1} ${file2} |
	grep -v 'REMOVEME' |
	awk '{print $1,$2,$3,$4,$5,$16,$24,$26,$6,$7,$8,$9,$10,$11,$12,$13,$14}' OFS='\t' | awk '$7=tolower($7)' OFS='\t' | awk '$8=tolower($8)' OFS='\t' > "$(basename "${file2}" 1)2"
done 3<OGIDs_list 4<conf1a1b1c2a2bTG

# mapping to geneNames file here to create an extra column
for i in mz pn ab nb on ; do
	for j in $i-*out2 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVEME";}}' ../Edge_Attributes/geneNamesMapping2.txt $j |
		grep -v 'REMOVEME' |
		awk '$18=tolower($18)' OFS='\t' |
		awk '{print $2,$3,$4,$5,$6,$1,$7,$8,$18,$9,$10,$11,$12,$13,$14,$15,$16,$17}' OFS='\t' |
		awk '{OFS="\t"} {if ($7=="null"||$7=="none") $7=$8; print $0}' |
		awk '{OFS="\t"} {if ($8=="null"||$8=="none") $8=$9; print $0}' |
		awk '{OFS="\t"} {if ($7=="null"||$7=="none") $7=$8; print $0}' |
		awk '{OFS="\t"} {if ($9=="null"||$9=="none") $9=$7; print $0}' |
		awk '{OFS="\t"} {if ($8=="null"||$8=="none") $8=$9; print $0}' > "$(basename "$j" 2)3"
	done
done

# filter, for the presence of genes in modules only - need to do for both GeneA and GeneB in modules
cp ../Module_genesandexpr/Mz-speciesspecnames_clusterassign.txt mz-speciesspecnames_clusterassign.txt
cp ../Module_genesandexpr/Pn-speciesspecnames_clusterassign.txt pn-speciesspecnames_clusterassign.txt
cp ../Module_genesandexpr/Ab-speciesspecnames_clusterassign.txt ab-speciesspecnames_clusterassign.txt
cp ../Module_genesandexpr/Nb-speciesspecnames_clusterassign.txt nb-speciesspecnames_clusterassign.txt
cp ../Module_genesandexpr/On-speciesspecnames_clusterassign.txt on-speciesspecnames_clusterassign.txt

for i in mz pn ab nb on ; do
	for j in $i-*out3 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$3]{print $0}' $i-speciesspecnames_clusterassign.txt $j > $j.2
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$6]{print $0}' $i-speciesspecnames_clusterassign.txt $j.2 > $j.3
	done
done

rm *.out3.2

# merge
cat *out3.3 > motifenr_merged-TFBSs_map2.txt

rm *out1
rm *out2
rm *out3

# keep *out3.3 files that are final species and confidence specific files - perfect for checking number of rows etc.

# add colheaders
printf 'motif_pattern\tGeneA_OGID\tGeneA\tGeneA_Symbol\tGeneB_OGID\tGeneB\tGeneB_SymbolDr\tGeneB_SymbolGa\tGeneB_SymbolSp\tstart\tstop\tstrand\tscore\tp-value\tq-value\tsequence\tconf_level\tconf_score\n' > colheads_tfbs
cat colheads_tfbs motifenr_merged-TFBSs_map2.txt > motifenr_merged-TFBSs_map2a.txt

# add interaction_type column, effect column, interaction column, directness column, direction column, layer column, and source column
awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"PD_directed"}' motifenr_merged-TFBSs_map2a.txt |
awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}' |
awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}' |
awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"direct"}' |
awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"directed"}' |
awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Transcriptional_regulation"}' |
awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"promoter_motif"}' > motifenr_merged-TFBSs_map2d.txt ### This is the final TF-Promoter interaction file ## {{DONE 18/03/2018}}

# check all rows have same number of columns
awk '{print NF}' motifenr_merged-TFBSs_map2d.txt | sort -nu | head -n 1 #25
awk '{print NF}' motifenr_merged-TFBSs_map2d.txt | sort -nu | tail -n 1 #25

# remove intermediate files
rm motifenr_merged-TFBSs_map2.txt
rm motifenr_merged-TFBSs_map2a.txt

##########################################################################################################################

# 2. CNEs - JASPAR Motif Discovery >--{DONE 21/03/18}--<
# For this you have to check which TF's are found in all modules and include only them where geneA will be the TF and geneB will be the 'CNE' and then another (hypothetical) cis-interaction will be the CNE (geneA) predicted to act as a enhancer for proximal gene (geneB).

# Run motif matching on all CNEs using new motifs - Will to do this, get scripts from him
# Done liftover of CNE coordinates in ATAC scripts - done using http://em-x1.gurdon.cam.ac.uk/cgi-bin/hgLiftOver
# All liftover files saved in /Users/mehtat/Desktop/cichlid_genomes/CNEs and /tgac/scratch/mehtat/Cichlid_GRNs/CNEs

mkdir /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs_v2
WDcnes=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs_v2)
cd $WDcnes

mkdir aCNEs
cd aCNEs

# First, extract aCNE sequences using BED annotation files
slurm
interative # can run this interactively
for i in /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/aCNEs/*.bed ; do ln -s $i ; done
#load the latest module
ml bedtools/2.25.0
ml GCC
ml zlib
# bedtools getfasta -fi test.fa -bed test.bed -name -fo test.fa.out
# the -name parameter will use the fourth column in your bed file for the fasta header
bedtools getfasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -bed aCNEs_Oniloticus.bed -name -fo On.aCNE.fa
bedtools getfasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -bed aCNEs_Mzebra1.1-Onliftover_true_aCNEs.bed -name -fo Mz.aCNE.fa
bedtools getfasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -bed aCNEs_Aburtoni-Onliftover_true_aCNEs.bed -name -fo Ab.aCNE.fa
bedtools getfasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -bed aCNEs_Pnyererei-Onliftover_true_aCNEs.bed -name -fo Pn.aCNE.fa
bedtools getfasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -bed aCNEs_Nbrichardi-Onliftover_true_aCNEs.bed -name -fo Nb.aCNE.fa

for i in *.fa ; do cp $i /tgac/workarea/Research-Groups/RG-cichlids/CNEs/ ; done # Copy the aCNEs fasta files to RG-cichlids folder so that Will can run the new motif scanning on CNEs and aCNEs

cd ../
# A-B. TF-CNE and TF-aCNE interaction tables

# 1. Create variables for all the mapping files required - these have already been created in TFTG files above
hspval=(/tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mat_qual_pvals_ALL.out2)
mmpval=(/tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mm10_mat_qual_pvals_ALL.out1)
OGIDab=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab)
OGIDmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz)
OGIDpn=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn)
OGIDnb=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb)
OGIDon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on)

hspvalmz=(../1.TFBSs_v2/mat_qual_pvals_ALL.out2.mz)
mmpvalmz=(../1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.mz)
hspvalpn=(../1.TFBSs_v2/mat_qual_pvals_ALL.out2.pn)
mmpvalpn=(../1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.pn)
hspvalab=(../1.TFBSs_v2/mat_qual_pvals_ALL.out2.ab)
mmpvalab=(../1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.ab)
hspvalnb=(../1.TFBSs_v2/mat_qual_pvals_ALL.out2.nb)
mmpvalnb=(../1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.nb)
hspvalon=(../1.TFBSs_v2/mat_qual_pvals_ALL.out2.on)
mmpvalon=(../1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.on)

hsmmpvalmz=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.mz)
hsmmpvalpn=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.pn)
hsmmpvalab=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.ab)
hsmmpvalnb=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.nb)
hsmmpvalon=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.on)

TFhsmmpval=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF)
ENShsmmpval=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.ENS)
# species-specific file - cichlid gene first column
TFhsmmpvalmz=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.mz)
TFhsmmpvalpn=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.pn)
TFhsmmpvalab=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.ab)
TFhsmmpvalnb=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.nb)
TFhsmmpvalon=(../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.on)

# A. TF-CNE

# 2a. Create the edge tables by doing the following in order (reflected in loops)

# Table format to output:
# motif_pattern
# GeneA_OGID
# GeneA
# GeneA_Symbol
# GeneB
# start
# stop
# strand
# score
# p-value
# q-value
# sequence
# conf_level
# conf_score

# 2ai. create unified (where relevant) mouse and human sets
# 2aii. Filter outputs (where relevant) for significant q-val (0.05)
# 2aiii. Create confidence levels for each scan
	# There are three confidence levels so scale this from 0 to 1 - no extrapolated set

# Confidence level 1
# Confidence level 1a - FDR corrected interactions resulting from the FIMO scans using the extrapolated cichlid species specific matrices: 0.5
nano mergeCSfimo.sh # {DONE}

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 32000 # memory pool for all cores
#SBATCH -t 0-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

for i in mz pn ab nb on ; do
	for j in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_CNEs/mouse/CS_default_${i}_results.out ; do
		for k in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_CNEs/human/CS_default_${i}_results.out ; do
			cat $j $k | # cat human and mouse files
			awk '($8 < 0.05 )' | # multiple test correction <0.05
			sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
			sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
			awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
			sort -k4,4 -rn | uniq -f 1 |
			sed $'s/;/\t/g' |
			awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="1a",$11="0.5"}' OFS='\t' >> $i-CS_default_fimo.out ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
			# cut -f1 $i-CS_default_fimo.out | sort -u > $i-CS_default_fimo.CUT.out # cut out merged first col and create no duplicates file
			# use cut file to match only the first hits in the original file which should take the top hit q-val
		done
	done
done

# Run all on uv2
sbatch mergeCSfimo.sh # {DONE - running}


# Confidence level 1b - FDR corrected interactions resulting from the FIMO scans using the extrapolated non-species specific matrices: 0.3 # {DONE}

for i in mz pn ab nb on ; do
	for j in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_CNEs/mouse/CW_default_${i}_results.out ; do
		for k in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_CNEs/human/CW_default_${i}_results.out ; do
			cat $j $k | # cat human and mouse files
			awk '($8 < 0.05 )' | # multiple test correction <0.05
			sed 's/RG-cich_GTRDdata_mouse_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level.ig_1//g' |
			sed 's/RG-cich_GTRDhuman_human_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level//g' >> $i-CW_default_fimo.out1 ;
		done
	done
done # {DONE}

for i in mz pn ab nb on ; do
	for j in $i-CW_default_fimo.out1 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $ENShsmmpval $j |
		sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
		awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' | # merge the 5 cols, separated by colon
		sort -k4,4 -rn | uniq -f 1 |
		sed $'s/;/\t/g' |
		awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19="1b",$20="0.3"}' OFS='\t' >> $i-CW_default_fimo.out2 ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	done
done # {DONE}


# Confidence level 1c - FDR corrected interactions resulting from the FIMO scans using the Jaspar matrices: 0.2
# this is done with TF mapping below

# Table format to output:
# motif_pattern
# GeneA_OGID
# GeneA
# GeneA_Symbol
# GeneB
# start
# stop
# strand
# score
# p-value
# q-value
# sequence
# conf_level
# conf_score

# 3. Map the first column (TF) to gene symbol/cichlid gene ID (where relevant - for 1a, 1b to TF symbol and 1c to cichlid ID)

# CS TFs # {DONE}
for i in mz pn ab nb on ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.$i $i-CS_default_fimo.out | awk '{print $18,$20,$1,$17,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > $i-CS_default_fimo.out1
done
rm *-CS_default_fimo.out # {DONE}

# CW TFs # {DONE}
# just select required columns from files created above, then remove old files
awk '{print $16,$18,$11,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' mz-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > mz-CW_default_fimo.out3
awk '{print $16,$18,$12,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' pn-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > pn-CW_default_fimo.out3
awk '{print $16,$18,$13,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' ab-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > ab-CW_default_fimo.out3
awk '{print $16,$18,$14,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' nb-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > nb-CW_default_fimo.out3
awk '{print $16,$18,$15,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' on-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > on-CW_default_fimo.out3
rm *-CW_default_fimo.out1
for i in mz pn ab nb on ; do
	mv $i-CW_default_fimo.out3 $i-CW_default_fimo.out1
done
rm *-CW_default_fimo.out2 # {DONE}

# JASPAR TFs # {DONE}
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_CNEs/human/JASPAR_default_mz_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > mz-JASPARTF_interm.out1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_CNEs/human/JASPAR_default_pn_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > pn-JASPARTF_interm.out1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_CNEs/human/JASPAR_default_ab_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > ab-JASPARTF_interm.out1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_CNEs/human/JASPAR_default_nb_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > nb-JASPARTF_interm.out1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_CNEs/human/JASPAR_default_on_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > on-JASPARTF_interm.out1

awk '{print $16,$18,$11,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' mz-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > mz-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $16,$18,$12,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' pn-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > pn-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $16,$18,$13,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' ab-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > ab-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $16,$18,$14,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' nb-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > nb-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $16,$18,$15,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' on-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > on-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)

rm *-JASPARTF_interm.out1

# filter, for the presence of genes in modules only - need to do for GeneA only
cp ../Module_genesandexpr/Mz-speciesspecnames_clusterassign.txt mz-speciesspecnames_clusterassign.txt
cp ../Module_genesandexpr/Pn-speciesspecnames_clusterassign.txt pn-speciesspecnames_clusterassign.txt
cp ../Module_genesandexpr/Ab-speciesspecnames_clusterassign.txt ab-speciesspecnames_clusterassign.txt
cp ../Module_genesandexpr/Nb-speciesspecnames_clusterassign.txt nb-speciesspecnames_clusterassign.txt
cp ../Module_genesandexpr/On-speciesspecnames_clusterassign.txt on-speciesspecnames_clusterassign.txt

for i in mz pn ab nb on ; do
	for j in $i-*_fimo.out1 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$3]{print $0}' $i-speciesspecnames_clusterassign.txt $j > "$(basename "$j" 1)2"
	done
done

# at this point, all TF mapped files are *.out2, then merge
cat *-*_fimo.out2 > All_528k_map.txt # {DONE}

# add colheaders
printf 'motif_pattern\tGeneA_OGID\tGeneA\tGeneA_Symbol\tGeneB\tstart\tstop\tstrand\tscore\tp-value\tq-value\tsequence\tconf_level\tconf_score\n' > colheads_528kcnes # prepare colheaders file
cat colheads_528kcnes All_528k_map.txt > All_528k_map9.txt # add column headers - ### This is the final file for TF-CNE interactions ##
# add all other columns of interaction, effect etc. - add all of these into one line, test
awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}' All_528k_map9.txt | # add effect column
awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"PD_directed"}' |
awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}' | # add interaction column
awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"unknown"}' | # add directness column
awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"unknown"}' | # add direction column
awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Transcriptional_regulation"}' | # add layer column
awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"CNE_motif"}' > All_528k_map10.txt # add source column # {DONE}

rm All_528k_map.txt

# B. TF-aCNE

cd aCNEs

nano generate_aCNEs.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 32000 # memory pool for all cores
#SBATCH -t 0-23:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address


# B. TF-aCNE

hspval=(/tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mat_qual_pvals_ALL.out2)
mmpval=(/tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mm10_mat_qual_pvals_ALL.out1)
OGIDab=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab)
OGIDmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz)
OGIDpn=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn)
OGIDnb=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb)
OGIDon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on)

hspvalmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mat_qual_pvals_ALL.out2.mz)
mmpvalmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.mz)
hspvalpn=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mat_qual_pvals_ALL.out2.pn)
mmpvalpn=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.pn)
hspvalab=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mat_qual_pvals_ALL.out2.ab)
mmpvalab=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.ab)
hspvalnb=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mat_qual_pvals_ALL.out2.nb)
mmpvalnb=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.nb)
hspvalon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mat_qual_pvals_ALL.out2.on)
mmpvalon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/mm10_mat_qual_pvals_ALL.out1.on)

hsmmpvalmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.mz)
hsmmpvalpn=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.pn)
hsmmpvalab=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.ab)
hsmmpvalnb=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.nb)
hsmmpvalon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.on)

TFhsmmpval=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF)
ENShsmmpval=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.ENS)
# species-specific file - cichlid gene first column
TFhsmmpvalmz=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.mz)
TFhsmmpvalpn=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.pn)
TFhsmmpvalab=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.ab)
TFhsmmpvalnb=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.nb)
TFhsmmpvalon=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.on)

# 2a. Create the edge tables by doing the following in order (reflected in loops)

# Table format to output:
# motif_pattern
# GeneA_OGID
# GeneA
# GeneA_Symbol
# GeneB
# start
# stop
# strand
# score
# p-value
# q-value
# sequence
# conf_level
# conf_score

# 2ai. create unified (where relevant) mouse and human sets
# 2aii. Filter outputs (where relevant) for significant q-val (0.05)
# 2aiii. Create confidence levels for each scan
	# There are three confidence levels so scale this from 0 to 1 - no extrapolated set

# Confidence level 1
# Confidence level 1a - FDR corrected interactions resulting from the FIMO scans using the extrapolated cichlid species specific matrices: 0.5

for i in mz pn ab nb on ; do
	for j in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_aCNEs/mouse/CS_default_${i}_results.out ; do
		for k in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_aCNEs/human/CS_default_${i}_results.out ; do
			cat $j $k | # cat human and mosue files
			awk '($8 < 0.05 )' | # multiple test correction <0.05
			sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' |
			sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
			awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9}' OFS='\t' | # merge the 5 cols, separated by colon
			sort -k4,4 -rn | uniq -f 1 |
			sed $'s/;/\t/g' |
			awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="1a",$11="0.5"}' OFS='\t' >> $i-CS_default_fimo.out ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
			# cut -f1 $i-CS_default_fimo.out | sort -u > $i-CS_default_fimo.CUT.out # cut out merged first col and create no duplicates file
			# use cut file to match only the first hits in the original file which should take the top hit q-val
		done
	done
done

# Confidence level 1b - FDR corrected interactions resulting from the FIMO scans using the extrapolated non-species specific matrices: 0.3 # {DONE}

for i in mz pn ab nb on ; do
	for j in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_aCNEs/mouse/CW_default_${i}_results.out ; do
		for k in /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_aCNEs/human/CW_default_${i}_results.out ; do
			cat $j $k | # cat human and mouse files
			awk '($8 < 0.05 )' | # multiple test correction <0.05
			sed 's/RG-cich_GTRDdata_mouse_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level.ig_1//g' |
			sed 's/RG-cich_GTRDhuman_human_sites_CichLevel_TF_ig_//g' | sed 's/_cichlid_level//g' >> $i-CW_default_fimo.out1 ;
		done
	done
done # {DONE}

for i in mz pn ab nb on ; do
	for j in $i-CW_default_fimo.out1 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $ENShsmmpval $j |
		sort -k1,1 -k2,2 -k3,3 -k4,4 | # sort based on on TF, TG, start, stop to find redundant lines
		awk '{print $1";"$2";"$3";"$4";"$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' | # merge the 5 cols, separated by colon
		sort -k4,4 -rn | uniq -f 1 |
		sed $'s/;/\t/g' |
		awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19="1b",$20="0.3"}' OFS='\t' >> $i-CW_default_fimo.out2 ; # sort qval ascending and report only unique hits based on col1, add confidence and save to new file
	done
done # {DONE}


# Confidence level 1c - FDR corrected interactions resulting from the FIMO scans using the Jaspar matrices: 0.2
# this is done with TF mapping below

# Table format to output:
# motif_pattern
# GeneA_OGID
# GeneA
# GeneA_Symbol
# GeneB
# start
# stop
# strand
# score
# p-value
# q-value
# sequence
# conf_level
# conf_score

# 3. Map the first column (TF) to gene symbol/cichlid gene ID (where relevant - for 1a, 1b to TF symbol and 1c to cichlid ID)

# CS TFs
for i in mz pn ab nb on ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/merged_mat_qual_pvals_ALL.out1.TF.$i $i-CS_default_fimo.out | awk '{print $18,$20,$1,$17,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > $i-CS_default_fimo.out1
done
rm *-CS_default_fimo.out

# CW TFs # {DONE}
# just select required columns from files created above, then remove old files
awk '{print $16,$18,$11,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' mz-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > mz-CW_default_fimo.out3
awk '{print $16,$18,$12,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' pn-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > pn-CW_default_fimo.out3
awk '{print $16,$18,$13,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' ab-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > ab-CW_default_fimo.out3
awk '{print $16,$18,$14,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' nb-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > nb-CW_default_fimo.out3
awk '{print $16,$18,$15,$17,$2,$3,$4,$5,$6,$7,$8,$9,$19,$20}' OFS='\t' on-CW_default_fimo.out2 | awk '$4=tolower($4)' OFS='\t' > on-CW_default_fimo.out3
rm *-CW_default_fimo.out1
for i in mz pn ab nb on ; do
	mv $i-CW_default_fimo.out3 $i-CW_default_fimo.out1
done
rm *-CW_default_fimo.out2 # {DONE}

# JASPAR TFs # {DONE}
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_aCNEs/human/JASPAR_default_mz_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > mz-JASPARTF_interm.out1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_aCNEs/human/JASPAR_default_pn_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > pn-JASPARTF_interm.out1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_aCNEs/human/JASPAR_default_ab_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > ab-JASPARTF_interm.out1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_aCNEs/human/JASPAR_default_nb_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > nb-JASPARTF_interm.out1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $TFhsmmpval /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput_aCNEs/human/JASPAR_default_on_results.out | grep -v 'NA' | awk '($8 < 0.05 )' > on-JASPARTF_interm.out1

awk '{print $16,$18,$11,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' mz-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > mz-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $16,$18,$12,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' pn-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > pn-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $16,$18,$13,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' ab-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > ab-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $16,$18,$14,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' nb-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > nb-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)
awk '{print $16,$18,$15,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' on-JASPARTF_interm.out1 |awk '{print $0,$13="1c",$14=0.2}' OFS='\t' | awk '$4=tolower($4)' OFS='\t' > on-JASPAR_default_fimo.out1 # this file is ready for TG mapping (as first column)

rm *-JASPARTF_interm.out1

# filter, for the presence of genes in modules only - need to do for GeneA only
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Mz-speciesspecnames_clusterassign.txt mz-speciesspecnames_clusterassign.txt
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Pn-speciesspecnames_clusterassign.txt pn-speciesspecnames_clusterassign.txt
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Ab-speciesspecnames_clusterassign.txt ab-speciesspecnames_clusterassign.txt
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Nb-speciesspecnames_clusterassign.txt nb-speciesspecnames_clusterassign.txt
cp /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/On-speciesspecnames_clusterassign.txt on-speciesspecnames_clusterassign.txt

for i in mz pn ab nb on ; do
	for j in $i-*_fimo.out1 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$3]{print $0}' $i-speciesspecnames_clusterassign.txt $j > "$(basename "$j" 1)2"
	done
done

# at this point, all TF mapped files are *.out2, then merge
cat *-*_fimo.out2 > All_aCNEs_map.txt # {DONE}

# add colheaders
printf 'motif_pattern\tGeneA_OGID\tGeneA\tGeneA_Symbol\tGeneB\tstart\tstop\tstrand\tscore\tp-value\tq-value\tsequence\tconf_level\tconf_score\n' > colheads_acnes # prepare colheaders file
cat colheads_acnes All_aCNEs_map.txt > All_aCNEs_map9.txt # add column headers - ### This is the final file for TF-aCNE interactions ##
# add all other columns of interaction, effect etc. - add all of these into one line, test

awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}' All_aCNEs_map9.txt | # add effect column
awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"PD_directed"}' |
awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}' | # add interaction column
awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"unknown"}' | # add directness column
awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"unknown"}' | # add direction column
awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Transcriptional_regulation"}' | # add layer column
awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"CNE_motif"}' > All_aCNEs_map10.txt # add source column
rm All_aCNEs_map.txt
mv All_aCNEs_map10.txt ../

# run the above
sbatch generate_aCNEs.sh


## A-B final files
# All_aCNEs_map10.txt # add column headers - ### This is the final file for TF-aCNE interactions ##
# All_528k_map10.txt # add column headers - ### This is the final file for TF-CNE interactions ##

# C-D. aCNE and CNE proximity to gene ##{DONE}##
# Then, create CNE proximity to gene asociations, GeneA (CNE) to GeneB (PROXIMAL GENE)

## HERE YOU SHOULD ADD NA TO EMPTY CELLS
## Then re-do the xargs grep
## alter the mapping below to use the new orthology file ../OGIDS.txt5
## can use the same awk mapping below - do not need the R script

cat ../7.CNEs/Dataset_7.txt | tr " " "_"  | perl -pe 's/\t(?=\t)/\tNA/g' > CNE_proximity.txt  # add underscore to all spaces and then use perl script to fill all empty spaces, tab delimited with NA
# Amend order of files for awk matching - the xargs works fine below.
#awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' CNE_proximity.txt > CNE_proximity.txt2
#awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' All_528k_map10.txt > All_528k_map10.txt2
#awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' OFS='\t' All_aCNEs_map10.txt > All_aCNEs_map10.txt2
#awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$4;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}'

# C. aCNE proximity ##{DONE}##
ln -s ../7.CNEs/list_species2
while read F; do cut -f1 $F-aCNEs.txt| xargs -i grep -wiF {} CNE_proximity.txt | sort | uniq > $F-aCNEs_proximity.txt ; done < list_species2 # map presence of ALL CNEs from CNE_proximity.txt file (column 2) in each species aCNE file e.g. brichardi-aCNEs.txt > output species specific files

# map cichlid gene IDs to the species specific files generated above using OGIDS.txt5 by mapping On Ensembl IDs (column 4 in proximity file and column 13 in OGIDS.txt5)
ln -s $OGIDS6
awk 'FNR==NR{a[$13]="\t"$5"\t"$10"\t"$12;next;} {print $0,a[$4]?a[$4]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 brichardi-aCNEs_proximity.txt > brichardi-aCNEs_proximitymap.txt # match column 4 of brichardi-aCNEs_proximity (file1) to column 13 (key) of ../OGIDS.txt5 (file2) and append field 5, 10 and 12 (value) of OGIDS.txt5 (file2) to brichardi-aCNEs_proximity (file1)
awk 'FNR==NR{a[$13]="\t"$4"\t"$10"\t"$12;next;} {print $0,a[$4]?a[$4]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 burtoni-aCNEs_proximity.txt > burtoni-aCNEs_proximitymap.txt  # match column 4 of burtoni-aCNEs_proximity (file1) to column 13 (key) of ../OGIDS.txt5 (file2) and append field 4, 10 and 12 (value) of OGIDS.txt5 (file2) to burtoni-aCNEs_proximity (file1)
awk 'FNR==NR{a[$13]="\t"$2"\t"$10"\t"$12;next;} {print $0,a[$4]?a[$4]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 mzebra-aCNEs_proximity.txt > mzebra-aCNEs_proximitymap.txt  # match column 4 of mzebra-aCNEs_proximity (file1) to column 13 (key) of ../OGIDS.txt5 (file2) and append field 2, 10 and 12 (value) of OGIDS.txt5 (file2) to mzebra-aCNEs_proximity (file1)
awk 'FNR==NR{a[$13]="\t"$3"\t"$10"\t"$12;next;} {print $0,a[$4]?a[$4]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 nyererei-aCNEs_proximity.txt > nyererei-aCNEs_proximitymap.txt  # match column 4 of nyererei-aCNEs_proximity (file1) to column 13 (key) of ../OGIDS.txt5 (file2) and append field 3, 10 and 12 (value) of OGIDS.txt5 (file2) to nyererei-aCNEs_proximity (file1)
awk 'FNR==NR{a[$13]="\t"$6"\t"$10"\t"$12;next;} {print $0,a[$4]?a[$4]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 tilapia-aCNEs_proximity.txt > tilapia-aCNEs_proximitymap.txt  # match column 4 of tilapia-aCNEs_proximity (file1) to column 13 (key) of ../OGIDS.txt5 (file2) and append field 6, 10 and 12 (value) of OGIDS.txt5 (file2) to tilapia-aCNEs_proximity (file1)

printf 'S/N\tGeneA\tLength\tOn_EnsemblID\tGeneB_Symbol\tDescription_of_associated_gene\tBri_identity-Med-Til\tBur_identity-Med-Til\tNye_identity-Med-Til\tMze_identity-Med-Til\tp-value(adj)\tGeneB\tGeneB_SymbolDr\tGeneB_SymbolGa\n' > colheads_cnes3 # prepare colheaders file
while read F ; do cat colheads_cnes3 $F-aCNEs_proximitymap.txt > $F-aCNEs_proximitymap2.txt ; done < list_species2 # add all colheaders
# add all other columns of interaction, effect etc. - add all of these into one line, test
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}' $F-aCNEs_proximitymap2.txt | awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"PD_undirected"}' > $F-aCNEs_proximitymap3.txt ; done < list_species2 # add effect column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}' $F-aCNEs_proximitymap3.txt > $F-aCNEs_proximitymap4.txt ; done < list_species2 # add interaction column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"indirect"}' $F-aCNEs_proximitymap4.txt > $F-aCNEs_proximitymap5.txt ; done < list_species2 # add directness column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"undirected"}' $F-aCNEs_proximitymap5.txt > $F-aCNEs_proximitymap6.txt ; done < list_species2 # add direction column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Transcriptional_regulation"}' $F-aCNEs_proximitymap6.txt > $F-aCNEs_proximitymap7.txt ; done < list_species2 # add layer column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"aCNE_Proximal"}' $F-aCNEs_proximitymap7.txt > $F-aCNEs_proximitymap8.txt ; done < list_species2 # add source column
head -1 brichardi-aCNEs_proximitymap8.txt > colheads_cnes4 # save the column headers from one of the five files to merge later
while read F ; do echo "$(tail -n +2 $F-aCNEs_proximitymap8.txt)" > $F-aCNEs_proximitymap9.txt ; done < list_species2 # remove first row from all files
cat *-aCNEs_proximitymap9.txt > All_aCNEs_proximitymap9.txt #Join all files
rm *-aCNEs_proximitymap*.txt # remove all intermediate files
cat colheads_cnes4 All_aCNEs_proximitymap9.txt > All_aCNEs_proximitymap10.txt # add column headers - ### This is the final file for aCNE-Gene interactions ##

# D. CNE proximity - DATASET7 IS ONLY FOR aCNE PROXIMAL, thus, for CNE proximal, do a bedtools closest to promoters ##{DONE}##

ml bedtools/2.25.0
ml GCC
ml zlib

# run bedtools closest to annotate location of peaks - but only to promoters
# -D a maps to both 5' and 3' genes, but no overlaps
# -t all maps to everything, 5' and 3' genes plus overlaps
# -d is closest distance to 5' or 3' genes

# sort the promoter files
for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed ; do sort -k1,1 -k2,2n $i > "$(basename "$i" .bed)_sorted.bed" ; done
sort -k1,1 -k2,2n /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.bed > Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed

# create symbolic link to CNE bed files - details of how created in NetworkReconstruction_v5.sh
for i in /tgac/scratch/mehtat/Cichlid_GRNs/CNEs/CNE.528k*bed ; do ln -s $i ; done
for i in CNE.528k_*-Onliftover.bed ; do sort -k1,1 -k2,2n $i > "$(basename "$i" .bed)_sorted.bed" ; done

# This will find the closest promoter (in case of ties, this will report all - tested with -t all and it provides the same number, and -t first gives the original input number)

bedtools closest -d -a CNE.528k.edit.bed -b Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed > On-CNE2prom.bed
bedtools closest -d -a CNE.528k_Aburtoni-Onliftover_sorted.bed -b Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Ab-CNE2prom.bed
bedtools closest -d -a CNE.528k_Mzebra1.1-Onliftover_sorted.bed -b Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Mz-CNE2prom.bed
bedtools closest -d -a CNE.528k_Nbrichardi-Onliftover_sorted.bed -b Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Nb-CNE2prom.bed
bedtools closest -d -a CNE.528k_Pnyererei-Onliftover_sorted.bed -b Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Pn-CNE2prom.bed

awk 'FNR==NR{a[$5]="\t"$10"\t"$12;next;} {print $0,a[$11]?a[$11]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 Nb-CNE2prom.bed | awk '{print $4,$1,$2,$3,$11,$14,$15}' OFS='\t' > brichardi-CNEs_proximitymap.txt # match column 11 of Nb-CNE2prom.bed (file1) to column 5 (key) of ../OGIDS.txt5 (file2) and append field 10 and 12 (value) of OGIDS.txt6 (file2) to (file1)
awk 'FNR==NR{a[$2]="\t"$10"\t"$12;next;} {print $0,a[$11]?a[$11]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 Mz-CNE2prom.bed | awk '{print $4,$1,$2,$3,$11,$14,$15}' OFS='\t' > mzebra-CNEs_proximitymap.txt # match column 11 of Mz-CNE2prom.bed (file1) to column 2 (key) of ../OGIDS.txt5 (file2) and append field 10 and 12 (value) of OGIDS.txt6 (file2) to (file1)
awk 'FNR==NR{a[$3]="\t"$10"\t"$12;next;} {print $0,a[$11]?a[$11]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 Pn-CNE2prom.bed | awk '{print $4,$1,$2,$3,$11,$14,$15}' OFS='\t' > nyererei-CNEs_proximitymap.txt # match column 11 of Pn-CNE2prom.bed (file1) to column 3 (key) of ../OGIDS.txt5 (file2) and append field 10 and 12 (value) of OGIDS.txt6 (file2) to (file1)
awk 'FNR==NR{a[$4]="\t"$10"\t"$12;next;} {print $0,a[$11]?a[$11]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 Ab-CNE2prom.bed | awk '{print $4,$1,$2,$3,$11,$14,$15}' OFS='\t' > burtoni-CNEs_proximitymap.txt # match column 11 of Ab-CNE2prom.bed (file1) to column 4 (key) of ../OGIDS.txt5 (file2) and append field 10 and 12 (value) of OGIDS.txt6 (file2) to (file1)
awk 'FNR==NR{a[$6]="\t"$10"\t"$12;next;} {print $0,a[$11]?a[$11]:"\t""NA""\t""NA""\t""NA";}' $OGIDS6 On-CNE2prom.bed | awk '{print $4,$1,$2,$3,$11,$14,$15}' OFS='\t' > tilapia-CNEs_proximitymap.txt # match column 11 of On-CNE2prom.bed (file1) to column 6 (key) of ../OGIDS.txt5 (file2) and append field 10 and 12 (value) of OGIDS.txt6 (file2) to (file1)

printf 'GeneA\tChr\tStart\tEnd\tGeneB\tGeneB_SymbolDr\tGeneB_SymbolGa\n' > colheads_cnes5 # prepare colheaders file
while read F ; do cat colheads_cnes5 $F-CNEs_proximitymap.txt > $F-CNEs_proximitymap2.txt ; done < list_species2 # add all colheaders
# add all other columns of interaction, effect etc. - add all of these into one line, test
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}' $F-CNEs_proximitymap2.txt | awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"PD_undirected"}' > $F-CNEs_proximitymap3.txt ; done < list_species2 # add effect column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}' $F-CNEs_proximitymap3.txt > $F-CNEs_proximitymap4.txt ; done < list_species2 # add interaction column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"indirect"}' $F-CNEs_proximitymap4.txt > $F-CNEs_proximitymap5.txt ; done < list_species2 # add directness column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"undirected"}' $F-CNEs_proximitymap5.txt > $F-CNEs_proximitymap6.txt ; done < list_species2 # add direction column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Transcriptional_regulation"}' $F-CNEs_proximitymap6.txt > $F-CNEs_proximitymap7.txt ; done < list_species2 # add layer column
while read F ; do awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"CNE_Proximal"}' $F-CNEs_proximitymap7.txt > $F-CNEs_proximitymap8.txt ; done < list_species2 # add source column
head -1 brichardi-CNEs_proximitymap8.txt > colheads_cnes6 # save the column headers from one of the five files to merge later
while read F ; do echo "$(tail -n +2 $F-CNEs_proximitymap8.txt)" > $F-CNEs_proximitymap9.txt ; done < list_species2 # remove first row from all files
cat *-CNEs_proximitymap9.txt > All_CNEs_proximitymap9.txt #Join all files
rm *-CNEs_proximitymap*.txt # remove all intermediate files
cat colheads_cnes6 All_CNEs_proximitymap9.txt > All_CNEs_proximitymap10.txt # add column headers - ### This is the final file for CNE-Gene interactions ##


###### FINAL FILES FOR CNEs ARE

All_aCNEs_map10.txt # add column headers - ### This is the final file for TF-aCNE interactions ##
All_528k_map10.txt # add column headers - ### This is the final file for TF-CNE interactions ##

All_aCNEs_proximitymap10.txt ### This is the final file for aCNE-Gene interactions ##
All_CNEs_proximitymap10.txt ### This is the final file for CNE-Gene interactions ##

##########################################################################################################################

# 3. miRNAs - using targetscan 7.0 target prediction >--{DONE 24/07/17}--<

mkdir $homeWD/6.miRNAs
WDmiRNA=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/6.miRNAs)
cd $WDmiRNA

# targetscan results here
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/miRNA_targetpred/targetscan/filtered/combined_families_targets.txt

# cut out only the columns that you require here and add 'gene' to GeneB column

# 1) orthogroup ID
# 3) Mirbase ID
# 4) Site Type
# 28) context++ score - the sum of the contribution of the above features
# 29) context++ score percentile - percentage of sites for this miRNA with a less favorable context++ score
# 30) AIR - Affected isoform ratio; fraction of transcripts with this stop site containing this site
# 31) weighted context++ score - the sum of the contribution of the above features, taking the AIR into account
# 32) weighted context++ score percentile - percentage of sites for this miRNA with a less favorable weighted context++ score
# 36) GeneA (miRNA family - name/ID of miRNA family)
# 38) species
# 39) GeneB

awk '{print $36,$39,$3,$4,$1,$28,$29,$30,$31,$32,$38}' OFS="\t" combined_families_targets.txt | awk '{gsub("mz.","mz.gene.",$2)}1' OFS="\t" | awk '{gsub("pn.","pn.gene.",$2)}1' OFS="\t" | awk '{gsub("ab.","ab.gene.",$2)}1' OFS="\t" | awk '{gsub("nb.","nb.gene.",$2)}1' OFS="\t" | awk '{gsub("on.","on.gene.",$2)}1' OFS="\t" > combined_families_targets.txt2

# Only keep interactions where genes are present in modules
cat ../Module_genesandexpr/*-modulegenes.txt > ../Module_genesandexpr/allmodulegenes.txt
awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" combined_families_targets.txt2 > combined_families_targets.txt3 # rearrange order for awk matching - total rows: 19613903
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"REMOVE";}}' ../Module_genesandexpr/allmodulegenes.txt combined_families_targets.txt3 | grep -v REMOVE | cut -f1-11 | awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" > miRNAmodulemap_merged.txt # 65074/69989 module genes retained; 15390993/19613903 interaction rows retained
rm combined_families_targets.txt2 combined_families_targets.txt3 # remove intermediate files

# add colheaders
printf 'GeneA\tGeneB\tmiRbaseID[ts7]\tSiteType[ts7]\tOrthogroupID\tcontext++score[ts7]\tcontext++scorepercentile[ts7]\tAIR[ts7]\tweightedcontext++score[ts7]\tweightedcontext++scorepercentile[ts7]\tspecies\n' > colheads_miRNA
cat colheads_miRNA miRNAmodulemap_merged.txt > miRNAmodulemap_merged1.txt
rm miRNAmodulemap_merged.txt

# add interaction_type column, effect column, interaction column, directness column, direction column, layer column, and source column
awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"Post-transcriptional_directed"}' miRNAmodulemap_merged1.txt |
awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"inhibition"}'|
awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"inhibition"}'|
awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"direct"}'|
awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"directed"}'|
awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Post-transcription_regulation"}'|
awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"targetscan7"}' > miRNAmodulemap_merged2.txt #THIS IS YOUR FINAL MIRNA FILE - 15390994 lines

##########################################################################################################################

# 4. TF-TG co-expression >--{DONE 22/01/2018}--<

mkdir $homeWD/8.TFTGco
WDTFTG=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/8.TFTGco)
cd $WDTFTG

# files here: http://pages.discovery.wisc.edu/~sroy/cichlids/networkinference/projected_speciesnets.tgz
# col1 - TF geneID, col2 - target geneID, col3 - confidence (1 = high; 0 = low)

# merge files
cat *.txt > tftgco_merged.txt

# add colheaders
printf 'GeneA\tGeneB\tconfidence[coTFTG]\n' > colheads_tftgco
cat colheads_tftgco tftgco_merged.txt > tftgco_merged1.txt
rm tftgco_merged.txt
mv tftgco_merged1.txt tftgco_merged.txt

# add interaction_type column, effect column, interaction column, directness column, direction column, layer column, and source column
awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"coexpression"}' tftgco_merged.txt |
awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}'|
awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}'|
awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"direct"}'|
awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"directed"}'|
awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Co-expression"}'|
awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"co-expression_TF-TG"}' > tftgco_merged1.txt
rm tftgco_merged.txt
mv tftgco_merged1.txt tftgco_mergedA.txt # 6675036 lines (excl header)

## NOTE: Sushmita recommended using edges >=0.5
awk '$3>=0.5' tftgco_mergedA.txt > tftgco_merged.txt # THIS IS YOUR FINAL TF-TG_coexpression file - leaves 23113 lines

##########################################################################################################################
#
# Network Reconstruction {DONE}
#
##########################################################################################################################

# 1. CREATE SPECIES-SPECIFIC FILES OF ALL EDGES - these will also be useful when not constraining (looking in the forest [hairball network] and not the trees [mini candidate networks]). These are found in /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/[SpeciesID]-Edge_Attributes_Collated4c.txt
# 2a. THEN CREATE MODULE-SPECIFIC FILES FROM THE SPECIES-SPECIFIC FILES
# 2b. DO THIS FOR EVERYTHING INCLUDING THE TF-PROMOTER, TF-CNE AND MIRNA-GENE INTERACTIONS - THE TF-PROMOTER INTERACTIONS ESPECIALLY SHOULD BE MODULE-SPECIFIC AS THE MODULES WILL CLUSTER ACCORDING TO CO-EXPRESSION IN TISSUES. THAT WAY THE TF IS LIKELY CO-EXPRESSED WITH THE GENES IN THE SAME MODULE WITH SIMILAR EXPRESSION IN THE TISSUE OF INTEREST.

##########################################################################################################################

# 1. Edge Attributes - >--{TO DO ONCE NEW TFBSs (in promoters and CNEs) HAVE BEEN GENERATED}--<

uv2k2
mkdir /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes
EDGE=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes)
cd $EDGE

# create symbolic links of relevant files to this folder
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/motifenr_merged-TFBSs_map2d.txt .
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/6.miRNAs/miRNAmodulemap_merged2.txt .
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs_v2/All_aCNEs_map10.txt .
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs_v2/All_528k_map10.txt .
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs_v2/All_aCNEs_proximitymap10.txt .
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs_v2/All_CNEs_proximitymap10.txt .
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/8.TFTGco/tftgco_merged.txt .

#### Create scripts to collate the edge tables ####

#### 1a. Build initial Edge Attribute table with R
nano 1a_Edge_Attributes_Build.R

miRNAmodulemap_merged2 <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/miRNAmodulemap_merged2.txt", header=T, na.strings=c(""," ","NA"))
motifenr_merged_TFBSs_map2d <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/motifenr_merged-TFBSs_map2d.txt", header=T, na.strings=c(""," ","NA"))
All_aCNEs_map10 <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/All_aCNEs_map10.txt", header=T, na.strings=c(""," ","NA"))
All_aCNEs_proximitymap10 <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/All_aCNEs_proximitymap10.txt", header=T, na.strings=c(""," ","NA"))
All_528k_map10 <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/All_528k_map10.txt", header=T, na.strings=c(""," ","NA"))
All_CNEs_proximitymap10 <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/All_CNEs_proximitymap10.txt", header=T, na.strings=c(""," ","NA"))
TFTGco <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/tftgco_merged.txt", header=T, na.strings=c(""," ","NA"))

#Merge
EAB1a=plyr::rbind.fill(miRNAmodulemap_merged2, motifenr_merged_TFBSs_map2d)
EAB2a=plyr::rbind.fill(EAB1a, All_aCNEs_map10)
EAB3a=plyr::rbind.fill(EAB2a, All_aCNEs_proximitymap10)
EAB4a=plyr::rbind.fill(EAB3a, All_528k_map10)
EAB5a=plyr::rbind.fill(EAB4a, All_CNEs_proximitymap10)
EAB6a=plyr::rbind.fill(EAB5a, TFTGco)

#Write out table - after this, apply shell script '2_Gene_symbol_mapping.sh'
write.table(EAB6a, "Edge_Attributes_Collated1.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#### 1b. Map gene IDs to gene symbols

# copy gene names file to here
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/geneNamesMapping.txt .

nano 1b_Gene_symbol_mapping.sh

#!/bin/sh

# This is after step one (plyr::rbind) in R script 'Edge_Attributes_Build.R'

# create appropriate gene names mapping file
awk '{if($2 == "NONE")print $1, $2, "NONE";else print $1, $2, $3;}' OFS='\t' geneNamesMapping.txt > geneNamesMapping2.txt

# Cut two columns, one for GeneA and another for GeneB columns to map 'Gene Symbols' to IDs
awk '{print $1}' Edge_Attributes_Collated1.txt > Edge_Attributes_GeneA.txt
awk '{print $2}' Edge_Attributes_Collated1.txt > Edge_Attributes_GeneB.txt

#Compare first column from both files. If there is a match, print the corresponding value of the first column and the matched third column. If no match is found, fill with "NA" in second column, retaining ProteinID in first column and then edit column headers to "GeneA", "GeneA_Symbol", and "GeneB", "GeneB_Symbol".
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' geneNamesMapping2.txt Edge_Attributes_GeneA.txt | sed -e '1s/GeneA/GeneA/' -e '1s/NA/GeneA_Symbol/' > Edge_Attributes_GeneA_mapped.txt
rm Edge_Attributes_GeneA.txt # remove the intermediate files
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' geneNamesMapping2.txt Edge_Attributes_GeneB.txt | sed -e '1s/GeneB/GeneB/' -e '1s/NA/GeneB_Symbol/' > Edge_Attributes_GeneB_mapped.txt
rm Edge_Attributes_GeneB.txt # remove the intermediate files

#Remove three 'Gene_Symbol' columns from Edge_Atrributes_Collated1.txt file
cut -f1-20,22-24,26-37,39-48 Edge_Attributes_Collated1.txt > Edge_Attributes_Collated2.txt


#### 1c. Finally build your Edge Attributes tables with several commands

nano 1c_Edge_Attributes_Build.sh

#!\bin\sh

# column bind main edge files and the GeneA GeneB mapping
paste -d'\t' Edge_Attributes_Collated2.txt Edge_Attributes_GeneA_mapped.txt Edge_Attributes_GeneB_mapped.txt > Edge_Attributes_Collated2a.txt

# remove colheaders
sed -i '1d' Edge_Attributes_Collated2a.txt

# awk remove certain columns and rearange at the same time
awk '{print $1,$47,$20,$2,$49,$22,$23,$21,$11,$12,$13,$14,$15,$16,$17,$18,$3,$4,$5,$6,$7,$8,$9,$10,$19,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45}' OFS='\t' Edge_Attributes_Collated2a.txt > Edge_Attributes_Collated2b.txt #done

# remove intermediate files
rm Edge_Attributes_Collated2a.txt

# print new colheaders
printf 'GeneA\tGeneA_Symbol\tGeneA_OGID\tGeneB\tGeneB_SymbolSp\tGeneB_SymbolDr\tGeneB_SymbolGa\tGeneB_OGID\tspecies\tinteraction_type\teffect\tinteraction\tdirectness\tdirection\tlayer\tsource\tmiRbaseID[ts7]\tSiteType[ts7]\tOrthogroupID\tcontext++score[ts7]\tcontext++scorepercentile[ts7]\tAIR[ts7]\tweightedcontext++score[ts7]\tweightedcontext++scorepercentile[ts7]\tmotif_pattern[pm_cm]\tstart[pm_cm]\tstop[pm_cm]\tstrand[pm_cm]\tscore[pm_cm]\tp.value[pm_cm]\tq.value[pm_cm]\tsequence[pm_cm]\tconf_level[pm_cm]\tconf_score[pm_cm]\tS.N[cp]\tLength[cp]\tOn_EnsemblID[cp]\tDescription_of_associated_gene[cp]\tBri_identity.Med.Til[cp]\tBur_identity.Med.Til[cp]\tNye_identity.Med.Til[cp]\tMze_identity.Med.Til[cp]\tp.value.adj[cp]\tChr[cp]\tStart[cp]\tEnd[cp]\tconfidence[coTFTG]\n' > colheads_edge

# add new colheaders to file above
cat colheads_edge Edge_Attributes_Collated2b.txt > Edge_Attributes_Collated3.txt

# remove intermediate file
rm Edge_Attributes_Collated2b.txt

# replace all spaces with an underscore
cat Edge_Attributes_Collated3.txt | tr -s ' ' | tr ' ' '_' > Edge_Attributes_Collated4.txt

# Create species-specific files
head -1 Edge_Attributes_Collated4.txt > colheaders_edgeattr # save the column headers
grep -wiF ab.gene Edge_Attributes_Collated4.txt > Ab-Edge_Attributes_Collated4.txt
grep -wiF nb.gene Edge_Attributes_Collated4.txt > Nb-Edge_Attributes_Collated4.txt
grep -wiF pn.gene Edge_Attributes_Collated4.txt > Pn-Edge_Attributes_Collated4.txt
grep -wiF mz.gene Edge_Attributes_Collated4.txt > Mz-Edge_Attributes_Collated4.txt
grep -wiF on.gene Edge_Attributes_Collated4.txt > On-Edge_Attributes_Collated4.txt

# create module-specific and then species-specific edge files - this requires that both GeneA and GeneB are in the module
# 1. Map and create two extra columns in each file - one for GeneA_Module and the other for GeneB_module

#Cut the required columns for mapping from species-specific files
cut -f1 Ab-Edge_Attributes_Collated4.txt > Ab-Edge_Attributes_Collated4_GeneA.txt
cut -f4 Ab-Edge_Attributes_Collated4.txt > Ab-Edge_Attributes_Collated4_GeneB.txt
cut -f1 Nb-Edge_Attributes_Collated4.txt > Nb-Edge_Attributes_Collated4_GeneA.txt
cut -f4 Nb-Edge_Attributes_Collated4.txt > Nb-Edge_Attributes_Collated4_GeneB.txt
cut -f1 Mz-Edge_Attributes_Collated4.txt > Mz-Edge_Attributes_Collated4_GeneA.txt
cut -f4 Mz-Edge_Attributes_Collated4.txt > Mz-Edge_Attributes_Collated4_GeneB.txt
cut -f1 Pn-Edge_Attributes_Collated4.txt > Pn-Edge_Attributes_Collated4_GeneA.txt
cut -f4 Pn-Edge_Attributes_Collated4.txt > Pn-Edge_Attributes_Collated4_GeneB.txt
cut -f1 On-Edge_Attributes_Collated4.txt > On-Edge_Attributes_Collated4_GeneA.txt
cut -f4 On-Edge_Attributes_Collated4.txt > On-Edge_Attributes_Collated4_GeneB.txt

#Compare first column from both gene files with the module gene mapping files. If there is a match, print the corresponding value of the first column and the matched second column. If no match is found, fill with "NA" in second column, retaining ID in first column
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/Ab-speciesspecnames_clusterassign.txt Ab-Edge_Attributes_Collated4_GeneA.txt > Ab-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/Ab-speciesspecnames_clusterassign.txt Ab-Edge_Attributes_Collated4_GeneB.txt > Ab-Edge_Attributes_Collated4_GeneBmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/Nb-speciesspecnames_clusterassign.txt Nb-Edge_Attributes_Collated4_GeneA.txt > Nb-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/Nb-speciesspecnames_clusterassign.txt Nb-Edge_Attributes_Collated4_GeneB.txt > Nb-Edge_Attributes_Collated4_GeneBmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/Pn-speciesspecnames_clusterassign.txt Pn-Edge_Attributes_Collated4_GeneA.txt > Pn-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/Pn-speciesspecnames_clusterassign.txt Pn-Edge_Attributes_Collated4_GeneB.txt > Pn-Edge_Attributes_Collated4_GeneBmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/Mz-speciesspecnames_clusterassign.txt Mz-Edge_Attributes_Collated4_GeneA.txt > Mz-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/Mz-speciesspecnames_clusterassign.txt Mz-Edge_Attributes_Collated4_GeneB.txt > Mz-Edge_Attributes_Collated4_GeneBmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/On-speciesspecnames_clusterassign.txt On-Edge_Attributes_Collated4_GeneA.txt > On-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/On-speciesspecnames_clusterassign.txt On-Edge_Attributes_Collated4_GeneB.txt > On-Edge_Attributes_Collated4_GeneBmap.txt

#Remove the first column from each file
cut -f2 Ab-Edge_Attributes_Collated4_GeneAmap.txt > Ab-Edge_Attributes_Collated4_GeneAmap1.txt
cut -f2 Ab-Edge_Attributes_Collated4_GeneBmap.txt > Ab-Edge_Attributes_Collated4_GeneBmap1.txt
cut -f2 Nb-Edge_Attributes_Collated4_GeneAmap.txt > Nb-Edge_Attributes_Collated4_GeneAmap1.txt
cut -f2 Nb-Edge_Attributes_Collated4_GeneBmap.txt > Nb-Edge_Attributes_Collated4_GeneBmap1.txt
cut -f2 Pn-Edge_Attributes_Collated4_GeneAmap.txt > Pn-Edge_Attributes_Collated4_GeneAmap1.txt
cut -f2 Pn-Edge_Attributes_Collated4_GeneBmap.txt > Pn-Edge_Attributes_Collated4_GeneBmap1.txt
cut -f2 Mz-Edge_Attributes_Collated4_GeneAmap.txt > Mz-Edge_Attributes_Collated4_GeneAmap1.txt
cut -f2 Mz-Edge_Attributes_Collated4_GeneBmap.txt > Mz-Edge_Attributes_Collated4_GeneBmap1.txt
cut -f2 On-Edge_Attributes_Collated4_GeneAmap.txt > On-Edge_Attributes_Collated4_GeneAmap1.txt
cut -f2 On-Edge_Attributes_Collated4_GeneBmap.txt > On-Edge_Attributes_Collated4_GeneBmap1.txt

# column bind main edge files and the mapped modules
paste -d'\t' Ab-Edge_Attributes_Collated4.txt Ab-Edge_Attributes_Collated4_GeneAmap1.txt Ab-Edge_Attributes_Collated4_GeneBmap1.txt > Ab-Edge_Attributes_Collated4a.txt
paste -d'\t' Nb-Edge_Attributes_Collated4.txt Nb-Edge_Attributes_Collated4_GeneAmap1.txt Nb-Edge_Attributes_Collated4_GeneBmap1.txt > Nb-Edge_Attributes_Collated4a.txt
paste -d'\t' Pn-Edge_Attributes_Collated4.txt Pn-Edge_Attributes_Collated4_GeneAmap1.txt Pn-Edge_Attributes_Collated4_GeneBmap1.txt > Pn-Edge_Attributes_Collated4a.txt
paste -d'\t' Mz-Edge_Attributes_Collated4.txt Mz-Edge_Attributes_Collated4_GeneAmap1.txt Mz-Edge_Attributes_Collated4_GeneBmap1.txt > Mz-Edge_Attributes_Collated4a.txt
paste -d'\t' On-Edge_Attributes_Collated4.txt On-Edge_Attributes_Collated4_GeneAmap1.txt On-Edge_Attributes_Collated4_GeneBmap1.txt > On-Edge_Attributes_Collated4a.txt

# Add OGIDs to GeneA and GeneB
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-ab Ab-Edge_Attributes_Collated4_GeneA.txt > Ab-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-ab Ab-Edge_Attributes_Collated4_GeneB.txt > Ab-Edge_Attributes_Collated4_GeneBmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-nb Nb-Edge_Attributes_Collated4_GeneA.txt > Nb-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-nb Nb-Edge_Attributes_Collated4_GeneB.txt > Nb-Edge_Attributes_Collated4_GeneBmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-pn Pn-Edge_Attributes_Collated4_GeneA.txt > Pn-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-pn Pn-Edge_Attributes_Collated4_GeneB.txt > Pn-Edge_Attributes_Collated4_GeneBmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-mz Mz-Edge_Attributes_Collated4_GeneA.txt > Mz-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-mz Mz-Edge_Attributes_Collated4_GeneB.txt > Mz-Edge_Attributes_Collated4_GeneBmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-on On-Edge_Attributes_Collated4_GeneA.txt > On-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/OGIDS.txt5-on On-Edge_Attributes_Collated4_GeneB.txt > On-Edge_Attributes_Collated4_GeneBmap.txt2

#Remove the first column from each file
cut -f2 Ab-Edge_Attributes_Collated4_GeneAmap.txt2 > Ab-Edge_Attributes_Collated4_GeneAmap1.txt2
cut -f2 Ab-Edge_Attributes_Collated4_GeneBmap.txt2 > Ab-Edge_Attributes_Collated4_GeneBmap1.txt2
cut -f2 Nb-Edge_Attributes_Collated4_GeneAmap.txt2 > Nb-Edge_Attributes_Collated4_GeneAmap1.txt2
cut -f2 Nb-Edge_Attributes_Collated4_GeneBmap.txt2 > Nb-Edge_Attributes_Collated4_GeneBmap1.txt2
cut -f2 Pn-Edge_Attributes_Collated4_GeneAmap.txt2 > Pn-Edge_Attributes_Collated4_GeneAmap1.txt2
cut -f2 Pn-Edge_Attributes_Collated4_GeneBmap.txt2 > Pn-Edge_Attributes_Collated4_GeneBmap1.txt2
cut -f2 Mz-Edge_Attributes_Collated4_GeneAmap.txt2 > Mz-Edge_Attributes_Collated4_GeneAmap1.txt2
cut -f2 Mz-Edge_Attributes_Collated4_GeneBmap.txt2 > Mz-Edge_Attributes_Collated4_GeneBmap1.txt2
cut -f2 On-Edge_Attributes_Collated4_GeneAmap.txt2 > On-Edge_Attributes_Collated4_GeneAmap1.txt2
cut -f2 On-Edge_Attributes_Collated4_GeneBmap.txt2 > On-Edge_Attributes_Collated4_GeneBmap1.txt2

# column bind main edge files and the mapped modules
paste -d'\t' Ab-Edge_Attributes_Collated4a.txt Ab-Edge_Attributes_Collated4_GeneAmap1.txt2 Ab-Edge_Attributes_Collated4_GeneBmap1.txt2 > Ab-Edge_Attributes_Collated4a.txt2
paste -d'\t' Nb-Edge_Attributes_Collated4a.txt Nb-Edge_Attributes_Collated4_GeneAmap1.txt2 Nb-Edge_Attributes_Collated4_GeneBmap1.txt2 > Nb-Edge_Attributes_Collated4a.txt2
paste -d'\t' Pn-Edge_Attributes_Collated4a.txt Pn-Edge_Attributes_Collated4_GeneAmap1.txt2 Pn-Edge_Attributes_Collated4_GeneBmap1.txt2 > Pn-Edge_Attributes_Collated4a.txt2
paste -d'\t' Mz-Edge_Attributes_Collated4a.txt Mz-Edge_Attributes_Collated4_GeneAmap1.txt2 Mz-Edge_Attributes_Collated4_GeneBmap1.txt2 > Mz-Edge_Attributes_Collated4a.txt2
paste -d'\t' On-Edge_Attributes_Collated4a.txt On-Edge_Attributes_Collated4_GeneAmap1.txt2 On-Edge_Attributes_Collated4_GeneBmap1.txt2 > On-Edge_Attributes_Collated4a.txt2


# add the colheaders to the created files
# add four more extra colheaders - GeneA_module and GeneB_module
printf 'GeneA\tGeneA_Symbol\tGeneA_OGID\tGeneB\tGeneB_SymbolSp\tGeneB_SymbolDr\tGeneB_SymbolGa\tGeneB_OGID\tspecies\tinteraction_type\teffect\tinteraction\tdirectness\tdirection\tlayer\tsource\tmiRbaseID[ts7]\tSiteType[ts7]\tOrthogroupID\tcontext++score[ts7]\tcontext++scorepercentile[ts7]\tAIR[ts7]\tweightedcontext++score[ts7]\tweightedcontext++scorepercentile[ts7]\tmotif_pattern[pm_cm]\tstart[pm_cm]\tstop[pm_cm]\tstrand[pm_cm]\tscore[pm_cm]\tp.value[pm_cm]\tq.value[pm_cm]\tsequence[pm_cm]\tconf_level[pm_cm]\tconf_score[pm_cm]\tS.N[cp]\tLength[cp]\tOn_EnsemblID[cp]\tDescription_of_associated_gene[cp]\tBri_identity.Med.Til[cp]\tBur_identity.Med.Til[cp]\tNye_identity.Med.Til[cp]\tMze_identity.Med.Til[cp]\tp.value.adj[cp]\tChr[cp]\tStart[cp]\tEnd[cp]\tconfidence[coTFTG]\tGeneA_module\tGeneB_module\tGeneA_OGID\tGeneB_OGID\n' > colheaders_edgeattr2

cat colheaders_edgeattr2 Ab-Edge_Attributes_Collated4a.txt2 > Ab-Edge_Attributes_Collated4b.txt #done
cat colheaders_edgeattr2 Nb-Edge_Attributes_Collated4a.txt2 > Nb-Edge_Attributes_Collated4b.txt #done
cat colheaders_edgeattr2 Pn-Edge_Attributes_Collated4a.txt2 > Pn-Edge_Attributes_Collated4b.txt #done
cat colheaders_edgeattr2 Mz-Edge_Attributes_Collated4a.txt2 > Mz-Edge_Attributes_Collated4b.txt #done
cat colheaders_edgeattr2 On-Edge_Attributes_Collated4a.txt2 > On-Edge_Attributes_Collated4b.txt #done

# awk rearrange columns
awk '{print $1,$2,$50,$48,$4,$5,$6,$7,$51,$49,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47}' OFS='\t' Ab-Edge_Attributes_Collated4b.txt > Ab-Edge_Attributes_Collated4c.txt # final species-specific files at this point
awk '{print $1,$2,$50,$48,$4,$5,$6,$7,$51,$49,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47}' OFS='\t' Nb-Edge_Attributes_Collated4b.txt > Nb-Edge_Attributes_Collated4c.txt # final species-specific files at this point
awk '{print $1,$2,$50,$48,$4,$5,$6,$7,$51,$49,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47}' OFS='\t' Pn-Edge_Attributes_Collated4b.txt > Pn-Edge_Attributes_Collated4c.txt # final species-specific files at this point
awk '{print $1,$2,$50,$48,$4,$5,$6,$7,$51,$49,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47}' OFS='\t' Mz-Edge_Attributes_Collated4b.txt > Mz-Edge_Attributes_Collated4c.txt # final species-specific files at this point
awk '{print $1,$2,$50,$48,$4,$5,$6,$7,$51,$49,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47}' OFS='\t' On-Edge_Attributes_Collated4b.txt > On-Edge_Attributes_Collated4c.txt # final species-specific files at this point

# remove all intermediate files created above
rm *-Edge_Attributes_Collated4a.txt
rm *-Edge_Attributes_Collated4a.txt2
rm *-Edge_Attributes_Collated4b.txt
rm *-Edge_Attributes_Collated4_Gene*.txt
rm *-Edge_Attributes_Collated4_Gene*.txt*

## Add the following to a script to run the full Edge Attributes Build processes on UV

nano 0.edgeattributesbuild-Arboretum_GT_v3.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml R
ml GCC

R CMD BATCH 1a_Edge_Attributes_Build.R # This will merge all the edge tables - DONE SCRIPT
sh 1b_Gene_symbol_mapping.sh # Run gene symbol mapping
sh 1c_Edge_Attributes_Build.sh #build your Edge Attributes tables with several commands	- split	by species and module too

# Run all on UV - requires at least 100Gb of memory
qsub -q Test -l select=1:mem=150GB:ncpus=1 0.edgeattributesbuild-Arboretum_GT_v3.sh

## Check lowest and highest number of columns in rows of a file to ensure no blank cells or spaces
for i in *-Edge_Attributes_Collated4c.txt ; do
	awk '{print NF}' $i | sort -nu | head -n 1 # lowest
done # 49 columns in each

for i in *-Edge_Attributes_Collated4c.txt ; do
	awk '{print NF}' $i | sort -nu | tail -n 1 # highest
done # 49 columns in each

### Sept 2018 - in the the following script: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/Rewiring_analysis.sh
### We filtered the TF-TG lines for Present NULL OGIDs in GeneB column
### In total there are 4209 PRESENT NULL OGIDS - Filter the geneB column for the main edge tables in the same way:

nano filteredges_presentNULLOGID.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-4
#SBATCH --mem-per-cpu 48000
#SBATCH -t 0-23:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 *-Edge_Attributes_Collated4c.txt > sp_edgetable # create a list of all edge tables
mapfile -t sp_edgetable < sp_edgetable # assign as elements to $sp_edgetable variable

PresentNOG=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS-noCand-noTF.txt2)

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$9]){print $0,a[$9];}else{print $0,"KEEPME";}}' $PresentNOG ${sp_edgetable[${SLURM_ARRAY_TASK_ID}]} | grep 'KEEPME' | cut -f1-49 > "$(basename "${sp_edgetable[${SLURM_ARRAY_TASK_ID}]}" 4c.txt)4d.txt"

# run the above
sbatch filteredges_presentNULLOGID.sh

## create module-specific files

nano moduleEDGES.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 0-23:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

# create a list file of species
echo Ab > list2
echo Nb >> list2
echo Pn >> list2
echo Mz >> list2
echo On >> list2

# Pull out according to module based on number presence in both module columns (4 and 10) and the other interactions that have just one 'Gene' present in the module based on:
# a. GeneA for TF-CNE interactions - use the source column ($18) entry 'CNE_motif' to pull out only those interactions and then the modules
# b. GeneB for CNE proximal genes - use the source column ($18) entry 'CNE_Proximal' to pull out only those interactions and then the modules
# c. GeneB for miRNA-mRNA interactions - use the source column ($18) entry 'targetscan7' to pull out only those interactions and then the modules

while read F ; do awk '{if (($4 == "0") && ($10 == "0")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module0.txt ; awk '{if (($4 == "0") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module0.txt ; awk '{if (($10 == "0") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module0.txt ; awk '{if (($4 == "0") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module0.txt ; awk '{if (($10 == "0") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module0.txt ; awk '{if (($4 == "0") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module0.txt ; awk '{if (($10 == "0") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module0.txt ; awk '{if (($4 == "0") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module0.txt ; awk '{if (($10 == "0") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module0.txt ; awk '{if (($10 == "0") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module0.txt ; done < list2

while read F ; do awk '{if (($4 == "1") && ($10 == "1")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module1.txt ; awk '{if (($4 == "1") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module1.txt ; awk '{if (($10 == "1") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module1.txt ; awk '{if (($4 == "1") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module1.txt ; awk '{if (($10 == "1") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module1.txt ; awk '{if (($4 == "1") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module1.txt ; awk '{if (($10 == "1") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module1.txt ; awk '{if (($4 == "1") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module1.txt ; awk '{if (($10 == "1") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module1.txt ; awk '{if (($10 == "1") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module1.txt ; done < list2

while read F ; do awk '{if (($4 == "2") && ($10 == "2")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module2.txt ; awk '{if (($4 == "2") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module2.txt ; awk '{if (($10 == "2") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module2.txt ; awk '{if (($4 == "2") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module2.txt ; awk '{if (($10 == "2") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module2.txt ; awk '{if (($4 == "2") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module2.txt ; awk '{if (($10 == "2") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module2.txt ; awk '{if (($4 == "2") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module2.txt ; awk '{if (($10 == "2") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module2.txt ; awk '{if (($10 == "2") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module2.txt ; done < list2

while read F ; do awk '{if (($4 == "3") && ($10 == "3")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module3.txt ; awk '{if (($4 == "3") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module3.txt ; awk '{if (($10 == "3") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module3.txt ; awk '{if (($4 == "3") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module3.txt ; awk '{if (($10 == "3") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module3.txt ; awk '{if (($4 == "3") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module3.txt ; awk '{if (($10 == "3") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module3.txt ; awk '{if (($4 == "3") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module3.txt ; awk '{if (($10 == "3") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module3.txt ; awk '{if (($10 == "3") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module3.txt ; done < list2

while read F ; do awk '{if (($4 == "4") && ($10 == "4")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module4.txt ; awk '{if (($4 == "4") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module4.txt ; awk '{if (($10 == "4") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module4.txt ; awk '{if (($4 == "4") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module4.txt ; awk '{if (($10 == "4") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module4.txt ; awk '{if (($4 == "4") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module4.txt ; awk '{if (($10 == "4") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module4.txt ; awk '{if (($4 == "4") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module4.txt ; awk '{if (($10 == "4") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module4.txt ; awk '{if (($10 == "4") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module4.txt ; done < list2

while read F ; do awk '{if (($4 == "5") && ($10 == "5")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module5.txt ; awk '{if (($4 == "5") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module5.txt ; awk '{if (($10 == "5") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module5.txt ; awk '{if (($4 == "5") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module5.txt ; awk '{if (($10 == "5") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module5.txt ; awk '{if (($4 == "5") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module5.txt ; awk '{if (($10 == "5") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module5.txt ; awk '{if (($4 == "5") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module5.txt ; awk '{if (($10 == "5") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module5.txt ; awk '{if (($10 == "5") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module5.txt ; done < list2

while read F ; do awk '{if (($4 == "6") && ($10 == "6")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module6.txt ; awk '{if (($4 == "6") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module6.txt ; awk '{if (($10 == "6") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module6.txt ; awk '{if (($4 == "6") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module6.txt ; awk '{if (($10 == "6") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module6.txt ; awk '{if (($4 == "6") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module6.txt ; awk '{if (($10 == "6") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module6.txt ; awk '{if (($4 == "6") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module6.txt ; awk '{if (($10 == "6") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module6.txt ; awk '{if (($10 == "6") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module6.txt ; done < list2

while read F ; do awk '{if (($4 == "7") && ($10 == "7")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module7.txt ; awk '{if (($4 == "7") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module7.txt ; awk '{if (($10 == "7") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module7.txt ; awk '{if (($4 == "7") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module7.txt ; awk '{if (($10 == "7") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module7.txt ; awk '{if (($4 == "7") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module7.txt ; awk '{if (($10 == "7") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module7.txt ; awk '{if (($4 == "7") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module7.txt ; awk '{if (($10 == "7") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module7.txt ; awk '{if (($10 == "7") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module7.txt ; done < list2

while read F ; do awk '{if (($4 == "8") && ($10 == "8")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module8.txt ; awk '{if (($4 == "8") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module8.txt ; awk '{if (($10 == "8") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module8.txt ; awk '{if (($4 == "8") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module8.txt ; awk '{if (($10 == "8") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module8.txt ; awk '{if (($4 == "8") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module8.txt ; awk '{if (($10 == "8") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module8.txt ; awk '{if (($4 == "8") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module8.txt ; awk '{if (($10 == "8") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module8.txt ; awk '{if (($10 == "8") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module8.txt ; done < list2

while read F ; do awk '{if (($4 == "9") && ($10 == "9")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module9.txt ; awk '{if (($4 == "9") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module9.txt ; awk '{if (($10 == "9") && ($18 == "CNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module9.txt ; awk '{if (($4 == "9") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module9.txt ; awk '{if (($10 == "9") && ($18 == "CNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module9.txt ; awk '{if (($4 == "9") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module9.txt ; awk '{if (($10 == "9") && ($18 == "aCNE_motif")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module9.txt ; awk '{if (($4 == "9") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module9.txt ; awk '{if (($10 == "9") && ($18 == "aCNE_Proximal")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module9.txt ; awk '{if (($10 == "9") && ($18 == "targetscan7")) print $0}' $F-Edge_Attributes_Collated4d.txt >> $F-Edge_Attributes_Collated4d-Module9.txt ; done < list2


# Add colheaders to all module-specific files
printf 'GeneA\tGeneA_Symbol\tGeneA_OGID\tGeneA_module\tGeneB\tGeneB_SymbolSp\tGeneB_SymbolDr\tGeneB_SymbolGa\tGeneB_OGID\tGeneB_module\tspecies\tinteraction_type\teffect\tinteraction\tdirectness\tdirection\tlayer\tsource\tmiRbaseID[ts7]\tSiteType[ts7]\tOrthogroupID\tcontext++score[ts7]\tcontext++scorepercentile[ts7]\tAIR[ts7]\tweightedcontext++score[ts7]\tweightedcontext++scorepercentile[ts7]\tmotif_pattern[pm_cm]\tstart[pm_cm]\tstop[pm_cm]\tstrand[pm_cm]\tscore[pm_cm]\tp.value[pm_cm]\tq.value[pm_cm]\tsequence[pm_cm]\tconf_level[pm_cm]\tconf_score[pm_cm]\tS.N[cp]\tLength[cp]\tOn_EnsemblID[cp]\tDescription_of_associated_gene[cp]\tBri_identity.Med.Til[cp]\tBur_identity.Med.Til[cp]\tNye_identity.Med.Til[cp]\tMze_identity.Med.Til[cp]\tp.value.adj[cp]\tChr[cp]\tStart[cp]\tEnd[cp]\tconfidence[coTFTG]\n' > colheaders_edgeattr3
INFO=$(cat colheaders_edgeattr3)  # read the contents of colheaders_edgeattr3 into var INFO
for i in *-Edge_Attributes_Collated4d-Module*.txt ; do sed -i "1i$INFO" $i ; done  # insert $INFO, use -i to change the file inplace
# Copy the new files to the workarea to copy local
for i in *-Edge_Attributes_Collated4d.txt ; do cp $i /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/ ; done
for i in *-Edge_Attributes_Collated4d-Module*.txt ; do cp $i /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/ ; done


### Final species edge files are [SpeciesID]-Edge_Attributes_Collated4d.txt
### Final species edge files are [SpeciesID]-Edge_Attributes_Collated4d-Module*.txt


##########################################################################################################################
#
# Create a large node table of all possible nodes
#
##########################################################################################################################

# STRUCTURE
# 1. use OGIDs file as base
# 2. add all possible miRNAs along all lines
# 3. add all possible CNEs (proximal and scanned for motifs) along all lines
# 4. add expression values (awk match cichlidgeneID to OGID first, then create a master table by merging OGID lines from each species, then sort -u)

mkdir /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Node_Attributes
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Node_Attributes

# 1. use OGIDs file as base
cat ../Module_genesandexpr/OGIDS.txt5 > 5Cichlid_Node_Attributes.txt

# 2. add all possible miRNAs along all lines
ln -s ../6.miRNAs/miRNAmodulemap_merged2.txt
cut -f1 miRNAmodulemap_merged2.txt | sort -u | grep -v 'GeneA' | perl -lne 'print "$_\t" x 16' >> 5Cichlid_Node_Attributes.txt

# 3. add all possible CNEs (proximal and scanned for motifs) along all lines
cut -f5 /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/All_aCNEs_map10.txt | sort -u | grep -v 'GeneB' | perl -lne 'print "$_\t" x 16' >> 5Cichlid_Node_Attributes.txt
cut -f2 /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/All_aCNEs_proximitymap10.txt | sort -u | grep -v 'GeneA' | perl -lne 'print "$_\t" x 16' >> 5Cichlid_Node_Attributes.txt
cut -f5 /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/All_528k_map10.txt | sort -u | grep -v 'GeneB' | perl -lne 'print "$_\t" x 16' >> 5Cichlid_Node_Attributes.txt
cut -f1 /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/All_CNEs_proximitymap10.txt | sort -u | grep -v 'GeneA' | perl -lne 'print "$_\t" x 16' >> 5Cichlid_Node_Attributes.txt

# 4. add expression values (make the adding species specific)
for i in ../Module_genesandexpr/*-exprtab_a.txt ; do ln -s $i ; done
ln -s ../Module_genesandexpr/OGIDS.txt5

# awk 'BEGIN{OFS="\t"}NR==FNR{a[$2]=$1;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' OGIDS.txt5 Mz-exprtab_a.txt | awk '{print $8,$2,$3,$4,$5,$6,$7}' OFS='\t' | sed '1d' | grep -v 'NULL' > Mz-exprtab_b.txt
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$3]=$1;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' OGIDS.txt5 Pn-exprtab_a.txt | awk '{print $8,$2,$3,$4,$5,$6,$7}' OFS='\t' | sed '1d' | grep -v 'NULL' > Pn-exprtab_b.txt
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$4]=$1;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' OGIDS.txt5 Ab-exprtab_a.txt | awk '{print $8,$2,$3,$4,$5,$6,$7}' OFS='\t' | sed '1d' | grep -v 'NULL' > Ab-exprtab_b.txt
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$5]=$1;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' OGIDS.txt5 Nb-exprtab_a.txt | awk '{print $8,$2,$3,$4,$5,$6,$7}' OFS='\t' | sed '1d' | grep -v 'NULL' > Nb-exprtab_b.txt
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$6]=$1;next}{if(a[$1]){print $0,a[$1];}else{print $0,"NULL";}}' OGIDS.txt5 On-exprtab_a.txt | awk '{print $8,$2,$3,$4,$5,$6,$7}' OFS='\t' | sed '1d' | grep -v 'NULL' > On-exprtab_b.txt

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$2]){print $0,a[$2];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' Mz-exprtab_a.txt 5Cichlid_Node_Attributes.txt > 5Cichlid_Node_Attributes-Mz.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$3]){print $0,a[$3];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' Pn-exprtab_a.txt 5Cichlid_Node_Attributes-Mz.txt > 5Cichlid_Node_Attributes-MzPn.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$4]){print $0,a[$4];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' Ab-exprtab_a.txt 5Cichlid_Node_Attributes-MzPn.txt > 5Cichlid_Node_Attributes-MzPnAb.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$5]){print $0,a[$5];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' Nb-exprtab_a.txt 5Cichlid_Node_Attributes-MzPnAb.txt > 5Cichlid_Node_Attributes-MzPnAbNb.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$6]){print $0,a[$6];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' On-exprtab_a.txt 5Cichlid_Node_Attributes-MzPnAbNb.txt > 5Cichlid_Node_Attributes-MzPnAbNbOn.txt

# check column numbers are uniform
awk '{print NF}' 5Cichlid_Node_Attributes-MzPnAbNbOn.txt | sort -nu | head -n 1 #51
awk '{print NF}' 5Cichlid_Node_Attributes-MzPnAbNbOn.txt | sort -nu | tail -n 1 #51

# FINAL NODE ATTRIBUTES FILE:
# 5Cichlid_Node_Attributes-MzPnAbNbOn.txt - copied local too
cd ../
cp -r Node_Attributes/ /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/

##########################################################################################################################
#
# Create separate PPI networks to study functional processes linked to regulatory networks #{DONE - March 2018}
#
##########################################################################################################################

# COLUMN DETAILS TO ADD TO EACH INTERACTION SET

#interaction_type
# PPI_directed #String, GeneMania, Signafish
# PPI_undirected #String, GeneMania, Signafish

#effect
# unknown
# stimulation
# inhibition

#interaction - make sure this is not mixed up with other 'interaction' columns
# unknown
# stimulation
# inhibition

#directness
# direct
# indirect

#direction
# directed
# undirected

#layer
# Directed_protein-protein_interaction #String, GeneMania, Signafish
# Undirected_protein-protein_interaction #String, GeneMania, Signafish

#source
# genemania[gm]
# signafish[sf]
# string[st]

##########################################################################################################################

# 2. GeneMania >--{SCRIPT DONE - RAN ON CLUSTER}--<

# cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/2.GeneMania
mkdir /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/2.GeneMania
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/2.GeneMania

# Can just d/l the GeneMania database
# wget http://genemania.org/data/current/Danio_rerio/Physical_interactions.PPI-mapped.txt
# This is a two column file with GeneA, GeneB and Weight - all that is required
cp /usr/users/ga004/mehtat/Cichlid_GRNs/Arboretum_GT_v3/2.GeneMania/Physical_interactions.PPI-mapped.txt . # copied genemania ppi file to home and then copied here

# remove column headers from Genemania file - 1st is GeneA, 2nd GeneB and 3rd Weight
head -1 Physical_interactions.PPI-mapped.txt > Genemania_colheads.txt # save colheads
tail -n +2 Physical_interactions.PPI-mapped.txt > Physical_interactions.PPI-mapped2.txt # remove headers
cut -f1 Physical_interactions.PPI-mapped2.txt > Physical_interactions.PPI-mapped2a.txt # geneA column only
cut -f2 Physical_interactions.PPI-mapped2.txt > Physical_interactions.PPI-mapped2b.txt # geneB column only

for i in Ab Mz Pn Nb On ; do
	cut -f1,3 ../Module_genesandexpr/$i-speciesspecnames_clusterassign_Drmapped.txt > $i-n1.txt # extract cichlidID and DrEnsemblID columns
	awk '{print $2, $1}' OFS='\t' $i-n1.txt > $i-n2.txt # rearrange columns
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $i-n2.txt Physical_interactions.PPI-mapped2a.txt > $i-gm3ba_mapped.txt #Compare first column from both files
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $i-n2.txt Physical_interactions.PPI-mapped2b.txt > $i-gm3bb_mapped.txt
	paste -d'\t' $i-gm3ba_mapped.txt $i-gm3bb_mapped.txt > $i-gm3b_mapped.txt # merge the geneA and geneB files, rearrange columns and add colheaders
	awk '{print $2, $1, $4, $3}' OFS='\t' $i-gm3b_mapped.txt > $i-gm3c_mapped.txt
	printf 'GeneA\tDr_ID_A\tGeneB\tDr_ID_B\n' > Genemania_colheads2.txt
	cat Genemania_colheads2.txt $i-gm3c_mapped.txt > $i-gm3c.txt
	cut -f3 Physical_interactions.PPI-mapped.txt > Physical_interactions.PPI-mapped-weight.txt # merge both geneA and geneB files with weight
	paste -d'\t' $i-gm3c.txt Physical_interactions.PPI-mapped-weight.txt > $i-gm3d.txt 	# add relevant columns below
	awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"PPI_undirected"}' $i-gm3d.txt |
	awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}' |
	awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}' |
	awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"indirect"}' |
	awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"undirected"}' |
	awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Undirected_protein-protein_interaction"}' |
	awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"genemania"}' > $i-gm3e.txt
done

# This section of scripts is for all species - run last
# Add new column for 'interaction2' as 'pp' (protein>protein) - did a for loop to run code on all files and output new file with 1a at end
for file in *-gm3e.txt; do awk -v OFS="\t" 'NR==1{print $0, "interaction2";next}{print $0,"pp"}' "$file" > "$(basename "$file" e.txt)f.txt"; done

#Save column headers (just use one file)
head -1 Pn-gm3f.txt > Genemania_colheads3.txt
#Remove column headers and bind
for file in *-gm3f.txt; do tail -n +2 "$file" > "$(basename "$file" f.txt)g.txt"; done
cat *-gm3g.txt > genemania_edge_merged.txt
#Add all column headers
cat Genemania_colheads3.txt genemania_edge_merged.txt > genemania_edge_merged1.txt
#Remove rows with NA in them as they will be useless
grep -wiFv 'NA' genemania_edge_merged1.txt > genemania_edge_merged2.txt ### THIS IS THE FINAL FILE FOR GENEMANIA ###

##########################################################################################################################

# 3. SignaFish >--{SCRIPT DONE - RAN ON CLUSTER}--<

# A total of 338 interactions in Signafish db
mkdir /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/3.SignaFish
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/3.SignaFish

### Ran below previously, just copy SignaFish.1.0.0.TNpHCu_edge_Ensembl.txt to this folder - copied to home and then copied here
cp /usr/users/ga004/mehtat/Cichlid_GRNs/Arboretum_GT_v3/3.SignaFish/SignaFish.1.0.0.TNpHCu_edge_Ensembl.txt .
cp /usr/users/ga004/mehtat/Cichlid_GRNs/Arboretum_GT_v3/3.SignaFish/list_species .

# d/l all interactions (except post transcriptional regulation) in cytoscape format > 02-09-2016-signafish-1.0.0-TNpHCu_edge.csv
cat 02-09-2016-signafish-1.0.0-TNpHCu_edge.csv | tr -d '"' | tr " " "_" | tr "," "\t" > 02-09-2016-signafish-1.0.0-TNpHCu_edge1.txt # remove quotations, add underscores and csv to tab-delimited
head -1 02-09-2016-signafish-1.0.0-TNpHCu_edge1.txt > signafish_colheaders.txt ; printf '\tGeneA_symbol\tDr_ID_A\tGeneB_symbol\tDr_ID_B\n' >> signafish_colheaders.txt ; tr -d "\n\r" < signafish_colheaders.txt > signafish_colheaders2.txt # Take first row of colheaders and add extra headers, the tr command removes newlines
echo "$(tail -n +2 02-09-2016-signafish-1.0.0-TNpHCu_edge1.txt)" > 02-09-2016-signafish-1.0.0-TNpHCu_edge2.txt # remove first row to process the column contents
cut -f2 02-09-2016-signafish-1.0.0-TNpHCu_edge2.txt | awk -F "_" '$1=$1' OFS="\t" > geneA_effect_geneB_Signafish.txt # split the second 'canonicalName' into three columns 'GeneA' 'effect' 'GeneB'
cut -f1 geneA_effect_geneB_Signafish.txt > GeneA_signafish.txt # Take GeneA column
cut -f3 geneA_effect_geneB_Signafish.txt > GeneB_signafish.txt # Take GeneB column
# map the gene names to Dr ENSEMBL IDs to map back to cichlid module genes - R script mapping use Biomart
R CMD BATCH 1.Signafish-EnsemblMap.R #output file is 'SignaFish.1.0.0.TNpHCu_edge_Ensembl.txt'
### Ran the above previously, no need to re-run

# Use all module genes in /Arboretum_GT_v3/Module_genesandexpr/ for each species and map them to Signafish interactions - do for each species (not required for each module), so five files
	#File 1 - two columns; cichlid gene ID (value), drerio ensembl ID (key) ; match to column 2
	#File 2 - Signafish interactions; match to column 16 'Dr_ID_A' and 18 'Dr_ID_B'
# Use above module files to grep all Signafish interactions ($16 and $18) and output according to all species module genes that correspond
while read F ; do cut -f3 ../Module_genesandexpr/$F-speciesspecnames_clusterassign_Drmapped.txt |xargs -i grep -wiF {} SignaFish.1.0.0.TNpHCu_edge_Ensembl.txt > $F-SignaFish.1.0.0.TNpHCu_edge_Ensembl.txt ; done < list_species

# then map column 16 and column 18 to cichlid gene IDs for each species - R script
ml R
ml GCC
R CMD BATCH 2.Signafish-CichlidMap.R

# add colheaders here
printf 'SUID\tcanonicalName\tdirection\tdirectness\teffect\tinteraction\tinteraction_type\tlayer\tname\treference\tselected\tshared_interaction\tshared_name\tsfsource\tGeneA\tGeneA_Symbol\tDr_ID_A\tGeneB\tGeneB_Symbol\tDr_ID_B\n' > signafish_colheaders3.txt
cat signafish_colheaders3.txt Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl.txt > Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl1.txt
awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"signafish"}' Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl1.txt > Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl2.txt # adds source column
# when filtering out by NA rows, make sure it is only in columns 15 or 18 (as some other columns have NA, like col 14)
gawk -F"\t" '$15 != "NA"{ print}' Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl2.txt | gawk -F"\t" '$18 != "NA"{ print}' > Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl3.txt #### THIS IS THE FINAL SIGNAFISH FILE
rm Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl.txt Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl1.txt Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl2.txt
rm *-SignaFish.1.0.0.TNpHCu_edge_Ensembl.txt

##########################################################################################################################

# 4. StringInteractions >--{SCRIPT DONE - RAN ON CLUSTER}--<

# ONLY COPIED THE RELEVANT FILES TO THE CLUSTER, COPY ALL OTHER FILES.

# Include - Of the 575203 interactions in previous data, 238338 were from String.

mkdir /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/4.String
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/4.String

### Below previously ran, just copy string_merged2.txt file over and run scripts below that

# might have to re-run STRING for new modules
# best to have all STRING interactions in zebrafish > Dr map specifically for each species > Separate genes and Dr mapping for each species > pull out string interactions for each species and each module at the modularisation step
# You may look at the one with combined score > 0.6 or be more stringent if you like.
# Focus on those with literature (textmining) and experimental evidences and ignore those that have only colocalization and coexpression

# wget http://string-db.org/download/protein.links.detailed.v10/7955.protein.links.detailed.v10.txt.gz > unzip > 7955.protein.links.detailed.v10.txt

# Prepare the string interactions file
tr ' ' '\t' < 7955.protein.links.detailed.v10.txt > 7955.protein.links.detailed.v10_delimit.txt ; rm 7955.protein.links.detailed.v10.txt # convert to tab delimited file and rm old file
head -1 7955.protein.links.detailed.v10_delimit.txt > colheaders_string # save colheaders if you need them
awk '(NR==1) || ($10 > 599 )' 7955.protein.links.detailed.v10_delimit.txt > 7955.protein.links.detailed.v10_delimit_600.txt # Only interactions with combined_score more than 600; 18291254 > 867736
awk '(NR==1) || ($7 > 1 )' 7955.protein.links.detailed.v10_delimit_600.txt | awk '(NR==1) || ($9 > 1 )' > 7955.protein.links.detailed.v10_delimit_600_exp-lit-1.txt # Only interactions with combined score >600 + experimental and literature text-mining score more than 1; 867736 > 466254
sed -e 's/7955.//g' 7955.protein.links.detailed.v10_delimit_600_exp-lit-1.txt > 7955.protein.links.detailed.v10_delimit_600_exp-lit-1_a.txt #remove all instances of '7955.' in front of Dr protein ID
#cut out the protein1 and protein2 columns for conversions
cut -f1 7955.protein.links.detailed.v10_delimit_600_exp-lit-1_a.txt > s1proteinA.txt
cut -f2 7955.protein.links.detailed.v10_delimit_600_exp-lit-1_a.txt > s1proteinB.txt


#convert protein ID to Dr ENSDARG gene ID - a total of 70115 rows being inputted (after the combined score filter). Remember though that all entries are duplicated so mapping value will be a lot lower
#Input the above files into BioMart to map Dr protein ID to Dr gene ID
ml R
ml GCC
R CMD BATCH 1.string_mapping_pt1.R
# Each file wrote out from above script as 's1geneA.txt' and 's1geneB.txt', now reorder columns so that Protein ID is first column:
awk '{print $2, $1}' OFS='\t' s1geneA.txt > s1geneA_.txt
awk '{print $2, $1}' OFS='\t' s1geneB.txt > s1geneB_.txt
#Compare first column from both files. If there is a match, print the corresponding value of the first column and the matched second column. If no match is found, fill with "NA" in second column, retaining ProteinID in first column
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' s1geneA_.txt s1proteinA.txt > s1proteinA_mapped.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' s1geneB_.txt s1proteinB.txt > s1proteinB_mapped.txt
# column bind merged file and gene IDs from mapping
paste -d'\t' 7955.protein.links.detailed.v10_delimit_600_exp-lit-1_a.txt s1proteinA_mapped.txt s1proteinB_mapped.txt > string_merged2.txt

### Above previously ran, just copy string_merged2.txt and ../Module_genesandexpr/OGIDS.txt5 file over and run below

# map to cichlid gene IDs with R script > string_merged3.txt
R CMD BATCH 2.string_mapping_pt2.R # Amended this according to the use of new orthology file 'OGIDS.txt5'

#Select the columns you require for each species, remove columns where mapping is NA or NULL
cut -f1-10,12,14,17,22 string_merged3.txt | awk '$13 !="NA" && $14 !="NA"' OFS='\t' | awk '$13 !="NULL" && $14 !="NULL"' OFS='\t' > Ab_string_merged3.txt
cut -f1-10,12,14,15,20 string_merged3.txt | awk '$13 !="NA" && $14 !="NA"' OFS='\t' | awk '$13 !="NULL" && $14 !="NULL"' OFS='\t' > Mz_string_merged3.txt
cut -f1-10,12,14,16,21 string_merged3.txt | awk '$13 !="NA" && $14 !="NA"' OFS='\t' | awk '$13 !="NULL" && $14 !="NULL"' OFS='\t' > Pn_string_merged3.txt
cut -f1-10,12,14,18,23 string_merged3.txt | awk '$13 !="NA" && $14 !="NA"' OFS='\t' | awk '$13 !="NULL" && $14 !="NULL"' OFS='\t' > Nb_string_merged3.txt
cut -f1-10,12,14,19,24 string_merged3.txt | awk '$13 !="NA" && $14 !="NA"' OFS='\t' | awk '$13 !="NULL" && $14 !="NULL"' OFS='\t' > On_string_merged3.txt

# delete all column headers to row bind
for file in *_string_merged3.txt ; do tail -n +2 "$file" > "$(basename "$file" .txt)a.txt"; done

# row bind all files
cat *_string_merged3a.txt > string_merged4.txt

# rearrange columns
awk 'BEGIN {OFS="\t"} {print $1,$11,$13,$2,$12,$14,$3,$4,$5,$6,$7,$8,$9,$10}' string_merged4.txt > string_merged4a.txt

# create colheaders file and add column headersto main file
printf 'ProteinA\tDr_ID_A\tGeneA\tProteinB\tDr_ID_B\tGeneB\tneighborhood\tfusion\tcooccurence\tcoexpression\texperimental\tdatabase\ttextmining\tcombined_score\n' > colheaders_string2 ; cat colheaders_string2 string_merged4a.txt > string_merged5.txt

# add interaction_type column, effect column, interaction column, directness column, direction column, layer column, and source column
awk -v OFS="\t" 'NR==1{print $0, "interaction_type";next}{print $0,"PPI_undirected"}' string_merged5.txt | awk -v OFS="\t" 'NR==1{print $0, "effect";next}{print $0,"unknown"}' | awk -v OFS="\t" 'NR==1{print $0, "interaction";next}{print $0,"unknown"}' | awk -v OFS="\t" 'NR==1{print $0, "directness";next}{print $0,"indirect"}' | awk -v OFS="\t" 'NR==1{print $0, "direction";next}{print $0,"undirected"}' | awk -v OFS="\t" 'NR==1{print $0, "layer";next}{print $0,"Undirected_protein-protein_interaction"}' | awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"string"}' > string_merged6.txt

# Map to Module Genes for all species
head -1 string_merged6.txt > colheads_string3 # save the colheads for merging later

# grep out all species genes from string_merged7.txt first and then run the xargs grep on those file for each species
grep -wiF ab.gene string_merged6.txt > Ab-string_merged6.txt
grep -wiF nb.gene string_merged6.txt > Nb-string_merged6.txt
grep -wiF mz.gene string_merged6.txt > Mz-string_merged6.txt
grep -wiF pn.gene string_merged6.txt > Pn-string_merged6.txt
grep -wiF on.gene string_merged6.txt > On-string_merged6.txt

while read F ; do cut -f3 ../Module_genesandexpr/$F-speciesspecnames_clusterassign_Drmapped.txt |xargs -i grep -wiF {} $F-string_merged6.txt > $F-string_merged7.txt ; done < list_species

# Remove all duplicate entries from the String mapping files
awk '!a[$0]++' Ab-string_merged7.txt > Ab-string_merged7a.txt
awk '!a[$0]++' Pn-string_merged7.txt > Pn-string_merged7a.txt
awk '!a[$0]++' Mz-string_merged7.txt > Mz-string_merged7a.txt
awk '!a[$0]++' Nb-string_merged7.txt > Nb-string_merged7a.txt
awk '!a[$0]++' On-string_merged7.txt > On-string_merged7a.txt

# Join all the files in all the folders
cat *-string_merged7a.txt > string_merged8.txt # row bind all files
cat colheads_string3 string_merged8.txt > string_merged9.txt # THESE IS THE FINAL STRING FILE
rm string_merged8.txt *_string_merged3* *-string_merged6* *-string_merged7.txt # remove intermediate files

# This creates the FINAL STRING FILES > string_merged9.txt

##########################################################################################################################

# 1. PPI Edge Attributes - >--{DONE}--<

uv
mkdir /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/PPI_Edge_Attributes
PPIEDGE=(/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/PPI_Edge_Attributes)
cd $PPIEDGE

# create symbolic links of relevant files to this folder
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/2.GeneMania/genemania_edge_merged2.txt .
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/3.SignaFish/Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl3.txt .
ln -s /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/4.String/string_merged9.txt .

#### Create scripts to collate the edge tables ####

#### 1a. Build initial PPI Edge Attribute table with R
nano 1a_PPIEdge_Attributes_Build.R #>--{DONE}--<

genemania_edge_merged2 <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/2.GeneMania/genemania_edge_merged2.txt", header=T, na.strings=c(""," ","NA"))
signafish <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/3.SignaFish/Cichlid_all_SignaFish.1.0.0.TNpHCu_edge_Ensembl3.txt", header=T, na.strings=c(""," ","NA"))
string_merged9 <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/4.String/string_merged9.txt", header=T, na.strings=c(""," ","NA"))

#Merge
EAB1a=plyr::rbind.fill(genemania_edge_merged2, signafish)
EAB2a=plyr::rbind.fill(EAB1a, string_merged9)

#Write out table - after this, apply shell script '2_Gene_symbol_mapping.sh'
write.table(EAB2a, "PPIEdge_Attributes_Collated1.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#### 1b. Map PPI gene IDs to gene symbols

# copy gene names file to here
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/geneNamesMapping.txt .

nano 1b_PPIGene_symbol_mapping.sh #>--{DONE}--<

#!/bin/sh

# This is after step one (plyr::rbind) in R script 'Edge_Attributes_Build.R'

# create appropriate gene names mapping file
awk '{if($2 == "NONE")print $1, $2, "NONE";else print $1, $2, $3;}' OFS='\t' geneNamesMapping.txt > geneNamesMapping2.txt

# Cut two columns, one for GeneA and another for GeneB columns to map 'Gene Symbols' to IDs
cut -f1 PPIEdge_Attributes_Collated1.txt > PPIEdge_Attributes_GeneA.txt
cut -f3 PPIEdge_Attributes_Collated1.txt > PPIEdge_Attributes_GeneB.txt

#Compare first column from both files. If there is a match, print the corresponding value of the first column and the matched third column. If no match is found, fill with "NA" in second column, retaining ProteinID in first column and then edit column headers to "GeneA", "GeneA_Symbol", and "GeneB", "GeneB_Symbol".
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' geneNamesMapping2.txt PPIEdge_Attributes_GeneA.txt | sed -e '1s/GeneA/GeneA/' -e '1s/NA/GeneA_Symbol/' > PPIEdge_Attributes_GeneA_mapped.txt
rm PPIEdge_Attributes_GeneA.txt # remove the intermediate files
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' geneNamesMapping2.txt PPIEdge_Attributes_GeneB.txt | sed -e '1s/GeneB/GeneB/' -e '1s/NA/GeneB_Symbol/' > PPIEdge_Attributes_GeneB_mapped.txt
rm PPIEdge_Attributes_GeneB.txt # remove the intermediate files

#Remove two 'Gene_Symbol' columns from PPIEdge_Atrributes_Collated1.txt file
cut -f1-21,24-33 PPIEdge_Attributes_Collated1.txt > PPIEdge_Attributes_Collated2.txt


#### 1c. Finally build your PPI Edge Attributes tables with several commands

nano 1c_PPIEdge_Attributes_Build.sh #>--{DONE}--<

#!\bin\sh

# column bind main edge files and the GeneA GeneB mapping
paste -d'\t' PPIEdge_Attributes_Collated2.txt PPIEdge_Attributes_GeneA_mapped.txt PPIEdge_Attributes_GeneB_mapped.txt > PPIEdge_Attributes_Collated2a.txt

# remove colheaders
sed -i '1d' PPIEdge_Attributes_Collated2a.txt

# awk remove certain columns and rearange at the same time
awk '{print $1,$33,$2,$3,$35,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31}' OFS='\t' PPIEdge_Attributes_Collated2a.txt > PPIEdge_Attributes_Collated2b.txt #done

# remove intermediate files
rm PPIEdge_Attributes_Collated2a.txt

# print new colheaders
printf 'GeneA\tGeneA_Symbol\tDr_ID_A\tGeneB\tGeneB_Symbol\tDr_ID_B\tWeight[gm]\tinteraction_type\teffect\tinteraction\tdirectness\tdirection\tlayer\tsource\tinteraction2[gm]\tSUID[sf]\tcanonicalName[sf]\tname[sf]\treference[sf]\tselected[sf]\tshared_interaction[sf]\tshared_name[sf]\tsfsource[sf]\tProteinA[st]\tProteinB[st]\tneighborhood[st]\tfusion[st]\tcooccurence[st]\tcoexpression[st]\texperimental[st]\tdatabase[st]\ttextmining[st]\tcombined_score[st]\n'  > colheads_PPIedge

# add new colheaders to file above
cat colheads_PPIedge PPIEdge_Attributes_Collated2b.txt > PPIEdge_Attributes_Collated3.txt

# remove intermediate file
rm PPIEdge_Attributes_Collated2b.txt

# replace all spaces with an underscore
cat PPIEdge_Attributes_Collated3.txt | tr -s ' ' | tr ' ' '_' > PPIEdge_Attributes_Collated4.txt

# Create species-specific files
head -1 PPIEdge_Attributes_Collated4.txt > colheaders_PPIedgeattr # save the column headers
for i in ab nb pn mz on ; do
	for j in Ab Nb Pn Mz On ; do
		grep -wiF $i.gene PPIEdge_Attributes_Collated4.txt > $j-PPIEdge_Attributes_Collated4.txt
	done
done

# create module-specific and then species-specific edge files - this requires that both GeneA and GeneB are in the module
# 1. Map and create two extra columns in each file - one for GeneA_Module and the other for GeneB_module

#Cut the required columns for mapping from species-specific files
for i in Ab Nb Pn Mz On ; do
	cut -f1 $i-PPIEdge_Attributes_Collated4.txt > $i-PPIEdge_Attributes_Collated4_GeneA.txt
	cut -f4 $i-PPIEdge_Attributes_Collated4.txt > $i-PPIEdge_Attributes_Collated4_GeneB.txt
done

#Compare first column from both gene files with the module gene mapping files. If there is a match, print the corresponding value of the first column and the matched second column. If no match is found, fill with "NA" in second column, retaining ID in first column
for i in Ab Nb Pn Mz On ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/$i-speciesspecnames_clusterassign.txt $i-PPIEdge_Attributes_Collated4_GeneA.txt > $i-PPIEdge_Attributes_Collated4_GeneAmap.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../Module_genesandexpr/$i-speciesspecnames_clusterassign.txt $i-PPIEdge_Attributes_Collated4_GeneB.txt > $i-PPIEdge_Attributes_Collated4_GeneBmap.txt
done

#Remove the first column from each file
for i in Ab Nb Pn Mz On ; do
	cut -f2 $i-PPIEdge_Attributes_Collated4_GeneAmap.txt > $i-PPIEdge_Attributes_Collated4_GeneAmap1.txt
	cut -f2 $i-PPIEdge_Attributes_Collated4_GeneBmap.txt > $i-PPIEdge_Attributes_Collated4_GeneBmap1.txt
done


# column bind main edge files and the mapped modules
for i in Ab Nb Pn Mz On ; do
	paste -d'\t' $i-PPIEdge_Attributes_Collated4.txt $i-PPIEdge_Attributes_Collated4_GeneAmap1.txt $i-PPIEdge_Attributes_Collated4_GeneBmap1.txt > $i-PPIEdge_Attributes_Collated4a.txt
done

# add the colheaders to the created files
# add two more extra colheaders - GeneA_module and GeneB_module

printf 'GeneA\tGeneA_Symbol\tDr_ID_A\tGeneB\tGeneB_Symbol\tDr_ID_B\tWeight[gm]\tinteraction_type\teffect\tinteraction\tdirectness\tdirection\tlayer\tsource\tinteraction2[gm]\tSUID[sf]\tcanonicalName[sf]\tname[sf]\treference[sf]\tselected[sf]\tshared_interaction[sf]\tshared_name[sf]\tsfsource[sf]\tProteinA[st]\tProteinB[st]\tneighborhood[st]\tfusion[st]\tcooccurence[st]\tcoexpression[st]\texperimental[st]\tdatabase[st]\ttextmining[st]\tcombined_score[st]\tGeneA_module\tGeneB_module\n' > colheaders_edgeattr2

for i in Ab Nb Pn Mz On ; do
	cat colheaders_edgeattr2 $i-PPIEdge_Attributes_Collated4a.txt > $i-PPIEdge_Attributes_Collated4b.txt #done
done

# awk rearrange columns

for i in Ab Nb Pn Mz On ; do
	awk '{print $1,$2,$34,$3,$4,$5,$35,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33}' OFS='\t' $i-PPIEdge_Attributes_Collated4b.txt > $i-PPIEdge_Attributes_Collated4c.txt # final species-specific files at this point
done

# remove all intermediate files created above
rm *-PPIEdge_Attributes_Collated4a.txt
rm *-PPIEdge_Attributes_Collated4b.txt
rm *-PPIEdge_Attributes_Collated4_Gene*.txt

# create a list file of species
echo Ab > list2
echo Nb >> list2
echo Pn >> list2
echo Mz >> list2
echo On >> list2

# create module-specific files
# Pull out according to module based on number presence in both module columns (3 and 7):
while read F ; do awk '{if (($3 == "0") && ($7 == "0")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module0.txt ; done < list2
while read F ; do awk '{if (($3 == "1") && ($7 == "1")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module1.txt ; done < list2
while read F ; do awk '{if (($3 == "2") && ($7 == "2")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module2.txt ; done < list2
while read F ; do awk '{if (($3 == "3") && ($7 == "3")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module3.txt ; done < list2
while read F ; do awk '{if (($3 == "4") && ($7 == "4")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module4.txt ; done < list2
while read F ; do awk '{if (($3 == "5") && ($7 == "5")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module5.txt ; done < list2
while read F ; do awk '{if (($3 == "6") && ($7 == "6")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module6.txt ; done < list2
while read F ; do awk '{if (($3 == "7") && ($7 == "7")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module7.txt ; done < list2
while read F ; do awk '{if (($3 == "8") && ($7 == "8")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module8.txt ; done < list2
while read F ; do awk '{if (($3 == "9") && ($7 == "9")) print $0}' $F-PPIEdge_Attributes_Collated4c.txt > $F-PPIEdge_Attributes_Collated4c-Module9.txt ; done < list2

# Add colheaders to all module-specific files
printf 'GeneA\tGeneA_Symbol\tGeneA_module\tDr_ID_A\tGeneB\tGeneB_Symbol\tGeneB_module\tDr_ID_B\tWeight[gm]\tinteraction_type\teffect\tinteraction\tdirectness\tdirection\tlayer\tsource\tinteraction2[gm]\tSUID[sf]\tcanonicalName[sf]\tname[sf]\treference[sf]\tselected[sf]\tshared_interaction[sf]\tshared_name[sf]\tsfsource[sf]\tProteinA[st]\tProteinB[st]\tneighborhood[st]\tfusion[st]\tcooccurence[st]\tcoexpression[st]\texperimental[st]\tdatabase[st]\ttextmining[st]\tcombined_score[st]\n' > colheaders_edgeattr3
INFO=$(cat colheaders_edgeattr3)  # read the contents of colheaders_edgeattr3 into var INFO
for i in *-PPIEdge_Attributes_Collated4c-Module*.txt ; do sed -i "1i$INFO" $i ; done  # insert $INFO, use -i to change the file inplace
# Copy the new files to the workarea to copy local
for i in *-PPIEdge_Attributes_Collated4c.txt ; do cp $i /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/PPIEdge_Attributes/ ; done
for i in *-PPIEdge_Attributes_Collated4c-Module*.txt ; do cp $i /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/PPIEdge_Attributes/ ; done


## Add the following to a script to run the full PPI Edge Attributes Build processes on UV

nano 0.PPIedgeattributesbuild-Arboretum_GT_v3.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml R
ml GCC

R CMD BATCH 1a_PPIEdge_Attributes_Build.R # This will merge all the edge tables - DONE SCRIPT, make sure miRNAmodulemap_merged2 is same file
sh 1b_PPIGene_symbol_mapping.sh # Run gene symbol mapping - DONE MOST OF SCRIPT BUT YOU NEED TO CHECK COL NUMBERS ETC. ALSO AMEND THE WEIRD NOT ADDING 'NA', MAYBE NEED TO ADD ANOTHER 'NA'
sh 1c_PPIEdge_Attributes_Build.sh

# Run all on UV - requires at least 50Gb of memory
qsub -q Test -l select=1:mem=50GB:ncpus=2 0.PPIedgeattributesbuild-Arboretum_GT_v3.sh #>--{DONE}--<

### Final species PPI edge files are /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/PPI_Edge_Attributes/[SpeciesID]-PPIEdge_Attributes_Collated4c.txt
### Final species, per module PPI edge files are /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/PPI_Edge_Attributes/[SpeciesID]-PPIEdge_Attributes_Collated4c-Module*.txt

## Check lowest and highest number of columns in rows of a file to ensure no blank cells or spaces
for i in *-PPIEdge_Attributes_Collated4c.txt ; do
	awk '{print NF}' $i | sort -nu | head -n 1 # lowest
done # 35 columns in each

for i in *-PPIEdge_Attributes_Collated4c.txt ; do
	awk '{print NF}' $i | sort -nu | tail -n 1 # highest
done # 35 columns in each

##########################################################################################################################
#
# Cytoscape Steps
#
##########################################################################################################################

# 1. Load edge attributes as network; GeneA (col1) is Source, Interaction is 'InteractionType' (col11) and Target is GeneB (col4)
# 2. Load cichlids-catch2 file  as node attributes and map to gene ID to get gene symbols [NOT NEEDED FOR NEW EDGE FILES]
# 3. Then load in appropriate expression (see below)
# 4. Edit->Remove Duplicated Edges tool to remove the duplicate edges and check the box to "Create an edge table column with number of duplicated edges".
# 5. Load in the expression values from the files Module_genesandexpr/*_exprtab.txt files
	# a. In cytoscape > Tools > NetworkAnalyzer > NetworkAnalysis > Analyse Network > Visualise Parameters >
		# i. Map node size to: Degree >
		# ii. Mape edge size: No of underlying edges
	# b. Make the nodes have 'expression' heat-strip barcharts to visualise expression in all six tissues:
		# i. http://opentutorials.cgl.ucsf.edu/index.php/Tutorial:Basic_Expression_Analysis_in_Cytoscape_3 > go to fun with charts section and follow
# 6. Font is 'avenirNext-demibold', font size 12, black text
# 7. PLACE ALL NODES WHERE YOU WANT THEM BEFORE DOING THIS: Merging edges > Layout > Bundle edges > all nodes and edges (will do all of them)
# 8. Play around with all node shapes and labelling as see fit - circle (normal genes); square (TFs); triangle (miRNAs); diamond (CNEs)
# 9. For the key: Edge thickness - number of underlying edges
# 10. For the key: Edge darkness - dark = more edges from various sources, light = fewer sources
# 10a. Style > Network > Background Paint > White

# 11. Save as pdf (others like png and svg do not come out properly)

# ~~ Functional Landscape Analysis ~~

# 12. Compare interactions of candidate genes between species and analyse the tissue-specific state of interactions - how does the gene interact in one tissue to another
	# a. Study partners of 'candidate-genes' by right clicking on the gene (in visual network) > Select > Select first neighbours (undirected)
	# b. find ways to compare networks - especially upon highlighting primary interactions of candidate genes (Functional landscape)
	# c. compare GO/pathway enrichments of 'regulators' and 'targets' between species
