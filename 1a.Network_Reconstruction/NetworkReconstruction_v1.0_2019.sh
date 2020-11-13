#!/bin/sh

##########################################################################################################################
#
# Network Reconstruction constrained by Arboretum output - 5 cichlids (O. niloticus, P. nyererei, A. burtoni, M .zebra, N. brichardi)
#
# Interaction inputs - TFBSs and TF-TG co-expression
#
# By Tarang Mehta, Earlham, UK
# Version 1.0 2019
#
##########################################################################################################################

# Run script: ./NetworkReconstruction_v1.0_2019.sh [workingdir] Module_genesandexpr TFBSs TFTGco Edge_Attributes
# NOTE: This will not run end to end as there are instances of copying code into scripts using nano
# Files in [workingdir]:
  # - Hs-Mm_EnsemblOrthology.txt # Ensembl orthology of Human and Mouse gene IDs
  # - full_orthologs_map3c.txt # contains the published orthologs gene IDs for mapping used previously
  # - list_k10 # list of k10 module numbers
  # - list_species # list of species
  # - OGIDS.txt # orthogroups
  # - OGIDS.txt5 # processed orthogroups file
  # - *_07.extrap.annotations.blast_out # extrapolation of motifs in cichlid gene promoters using human and mouse orthology
  # - mz-speciesspecnames_clusterassign.txt # Arboretum module assignment of M. zebra genes
  # - pn-speciesspecnames_clusterassign.txt # Arboretum module assignment of P. nyererei genes
  # - ab-speciesspecnames_clusterassign.txt # Arboretum module assignment of A. burtoni genes
  # - nb-speciesspecnames_clusterassign.txt # Arboretum module assignment of N. brichardi genes
  # - on-speciesspecnames_clusterassign.txt # Arboretum module assignment of O. niloticus genes
  # - geneNamesMapping.txt # gene orthology of cichlid IDs, Ensembl IDs and gene descriptions
  # - MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS-noCand-noTF.txt2 # 4209 orthogroups with NULL orthogroup IDs for one or species are actually present in the genome (mis-annotations etc.) and not lost/absent

##########################################################################################################################

# COLUMN DETAILS TO ADD TO EACH INTERACTION SET

#interaction_type
# coexpression #TF-TG coexpression
# PD_directed #TFBSs
# PD_undirected #TFBSs

#effect
# unknown
# stimulation
# inhibition

#interaction
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
# Co-expression #TF-TG: PGG algorithm
# Transcriptional_regulation #TFBSs

#source
# co-expression_TF-TG[coTFTG]
# promoter_motif[pm]


homeWD=($1)
cd $homeWD

##########################################################################################################################

# 1. Genes and expression in each module >--{DONE - RAN ON SCRATCH > COPIED LOCAL}--<
# cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr
mkdir $homeWD/$2
cd $homeWD/$2

# Extra files required:
FO=(full_orthologs_map3c.txt) # contains the published orthologs gene IDs for mapping used previously
k10=(list_k10) # list of k10 module numbers
list_species=(list_species) . # list of species

# get extra files from Arboretum output
wget  -r -nH --cut-dirs=4 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~ckoch/cichlids/prediction_mode_results/results/
wget  -r -nH --cut-dirs=3 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~ckoch/cichlids/prediction_mode_results/relevant_ogs_new.txt
wget  -r -nH --cut-dirs=3 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~ckoch/cichlids/prediction_mode_results/relevant_ogs_old.txt
wget  -r -nH --cut-dirs=3 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~ckoch/cichlids/prediction_mode_results/OGIDS.txt

OGIDS=(OGIDS.txt) # the most recent OGids file for the new prediction mode modules
WD=($homeWD/Module_genesandexpr/)

# Map the OGIDS file to Dr symbols using full_orthology file
# already have below files where species gene placed first for awk mapping
awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' $FO > full_orthologs_map3c_mz.txt
fullorthologsmz=full_orthologs_map3c_mz.txt
awk '{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' $FO > full_orthologs_map3c_pn.txt
fullorthologspn=full_orthologs_map3c_pn.txt
awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' $FO > full_orthologs_map3c_ab.txt
fullorthologsab=full_orthologs_map3c_ab.txt
awk '{print $5,$1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' $FO > full_orthologs_map3c_nb.txt
fullorthologsnb=full_orthologs_map3c_nb.txt
awk '{print $6,$1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' $FO > full_orthologs_map3c_on.txt
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
HsMmENS=(Hs-Mm_EnsemblOrthology.txt)
geneNamesTree=(geneNamesTree)
perl -pe 's/\t/\tNULL/g' $HsMmENS | sed 's/NULLENSMU/ENSMU/g' > Hs-Mm_EnsemblOrthology.txt2 # replace all empty spaces with NULL
# copied to scratch dir
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL";}}' Hs-Mm_EnsemblOrthology.txt2 EnHsOGID | cut -f2 > EnMmOGID
paste -d'\t' OGIDS.txt2 EnOlOGID EnTnOGID EnGaOGID GaOGID EnDrOGID DrOGID EnOnOGID EnHsOGID HsOGID EnMmOGID | awk '{print $1,$4,$5,$3,$6,$2,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' > OGIDS.txt5
OGIDS5=$WD$(echo $(basename $OGIDS) | sed -e 's/.txt/.txt5/')
paste -d'\t' $OGIDS2 EnOlOGID EnTnOGID EnGaOGID GaOGID EnDrOGID DrOGID EnOnOGID EnHsOGID HsOGID EnMmOGID | awk '{print $1,$4,$5,$3,$6,$2,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' > $OGIDS5 ##### THIS IS YOUR FINAL ORTHOLOGY FILE - OGIDS.txt5 #####

## Finally, leverage the gene names from the Nature publication for the ab gene (as it has most genes annotated) to have a comprehensive set
awk '{OFS="\t"} {if ($2=="NONE") $3="NONE"; print $0}' $geneNamesTree | sed 's/ /_/g' | sed 's/Q4SR56_TETNG/sws2/g' | sed 's/OPSG_ORYLA/rh2/g' | sed 's/Q2L6A1_ORYLA/sws1/g' > geneNamesTree2
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

##########################################################################################################################

# 2. TFBSs - Motif Discovery in gene promoter regions

mkdir $homeWD/$3
WDpromtfbs=($homeWD/$3)

# Motif Scanning scripts are shared in the '2a.TF_motif_scanning' folder
# Summary of motif scanning
# Species-specific promoter sequences used as 0-order markov backgrounds
# FIMO ran with default p-value (1e-4)
# Files will be of Mouse and Human Extrapolations, Cichlid-specific (CS) and Cichlid-wide (CW) and Vertebrate JASPAR motifs - NOTE: with the exception of extrapolation files (*_07.extrap.annotations.blast_out), all other raw output files (CS and CW) are available upon requests
# Otherwise, the final file 'motifenr_merged-TFBSs_map2d.txt' is found with all the other files in the same folder

# 1. Create variables for relevant folders and mapping files
cd $WDpromtfbs
hspval=($homeWD/mat_qual_pvals_ALL.out2)
mmpval=($homeWD/mm10_mat_qual_pvals_ALL.out1)
OGIDab=($homeWD/$2/OGIDS.txt5-ab)
OGIDmz=($homeWD/$2/OGIDS.txt5-mz)
OGIDpn=($homeWD/$2/OGIDS.txt5-pn)
OGIDnb=($homeWD/$2/OGIDS.txt5-nb)
OGIDon=($homeWD/$2/OGIDS.txt5-on)

# 1a. amend the p-vals files so that they are species-specific for awk matching #
for i in $hspval $mmpval ; do
	awk '{print $4,$1,$2,$3,$5,$6,$7,$8}' OFS='\t' $i | sed 's/;/\t/' | sed 's/;/\t/' | sed 's/;/\t/' | sed 's/;/\t/' > $i.mz
	awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' $i.mz > $i.pn
	awk '{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' $i.mz > $i.ab
	awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' $i.mz > $i.nb
	awk '{print $5,$1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12}' OFS='\t' $i.mz > $i.on
done
for i in $mmpval/*out1.* ; do mv $i . ; done
for i in $hspval/*out2.* ; do mv $i . ; done

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

# Species-specific extrapolation files generated by the motif prediction pipeline are:
# ab-mm10_07.extrap.annotations.blast_out
# mz-mm10_07.extrap.annotations.blast_out
# nb-mm10_07.extrap.annotations.blast_out
# on-mm10_07.extrap.annotations.blast_out
# pn-mm10_07.extrap.annotations.blast_out
# ab-hg38_07.extrap.annotations.blast_out
# mz-hg38_07.extrap.annotations.blast_out
# nb-hg38_07.extrap.annotations.blast_out
# on-hg38_07.extrap.annotations.blast_out
# pn-hg38_07.extrap.annotations.blast_out

# find the intersect by catting, sorting and retaining duplicates
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

# Confidence level 1b - interactions in cichlids extrapolated from mouse ONLY: 0.2

for i in mz pn ab nb on ; do
	for j in ${i}-mm10hg38_07.extrap.annotations.blast_out.intersect ; do
		for k in ${i}-mm10_07.extrap.annotations.blast_out ; do
			cat $j $k | sort | uniq -u | # remove the intersected lines from the original mouse file
			awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="1b",$11="0.2"}' OFS='\t' > $i-mm10_07.extrap.annotations.blast_out.unique
		done
	done
done

# Confidence level 1c - interactions in cichlids extrapolated from human ONLY: 0.15

for i in mz pn ab nb on ; do
	for j in ${i}-mm10hg38_07.extrap.annotations.blast_out.intersect ; do
		for k in ${i}-hg38_07.extrap.annotations.blast_out ; do
			cat $j $k | sort | uniq -u | # remove the intersected lines from the original human file
			awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10="1c",$11="0.15"}' OFS='\t' > $i-hg38_07.extrap.annotations.blast_out.unique
		done
	done
done

# merge the extrapolated files
for i in mz pn ab nb on ; do
	cat $i-mm10hg38_07.extrap.annotations.blast_out.intersect $i-mm10_07.extrap.annotations.blast_out.unique $i-hg38_07.extrap.annotations.blast_out.unique >> $i-mm10hg38_07.extrap.annotations.blast_out.merged
done

# Confidence level 2 #
# Confidence level 2a - FDR corrected interactions resulting from the FIMO scans using the extrapolated cichlid species specific matrices: 0.125 #

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

# Confidence level 2b - FDR corrected interactions resulting from the FIMO scans using the extrapolated non-species specific matrices: 0.110 #

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


# Confidence level 2c - FDR corrected interactions resulting from the FIMO scans using the Jaspar matrices: 0.115 #
# this is done with TF mapping below

# 3. Map the first column (TF) to gene symbol/cichlid gene ID (where relevant - for 1a, 1b, 1c, 2a, 2b)

# Extrapolated TFs #
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

# CS TFs #
for i in mz pn ab nb on ; do
	awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' merged_mat_qual_pvals_ALL.out1.TF.$i $i-CS_default_fimo.out | awk '{print $2,$18,$20,$1,$17,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS='\t' | awk '$5=tolower($5)' OFS='\t' > $i-CS_default_fimo.out1
done
rm *-CS_default_fimo.out

# CW TFs #
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

# JASPAR TFs
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

echo $homeWD/$2/OGIDS.txt5-mz > OGIDs_list
echo $homeWD/$2/OGIDS.txt5-pn >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-ab >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-nb >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-on >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-mz >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-pn >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-ab >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-nb >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-on >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-mz >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-pn >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-ab >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-nb >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-on >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-mz >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-pn >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-ab >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-nb >> OGIDs_list
echo $homeWD/$2/OGIDS.txt5-on >> OGIDs_list

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

for i in mz pn ab nb on ; do
	for j in $i-*out3 ; do
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$3]{print $0}' $homeWD/$i-speciesspecnames_clusterassign.txt $j > $j.2
		awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;next}a[$6]{print $0}' $homeWD/$i-speciesspecnames_clusterassign.txt $j.2 > $j.3
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
awk -v OFS="\t" 'NR==1{print $0, "source";next}{print $0,"promoter_motif"}' > motifenr_merged-TFBSs_map2d.txt ### This is the final TF-Promoter interaction file

# check all rows have same number of columns
awk '{print NF}' motifenr_merged-TFBSs_map2d.txt | sort -nu | head -n 1 #25
awk '{print NF}' motifenr_merged-TFBSs_map2d.txt | sort -nu | tail -n 1 #25

# remove intermediate files
rm motifenr_merged-TFBSs_map2.txt
rm motifenr_merged-TFBSs_map2a.txt

##########################################################################################################################

# 3. TF-TG co-expression - PGG algorithm method

mkdir $homeWD/$4
WDTFTG=($homeWD/$4)
cd $WDTFTG

# download the TF-TG co-expression files generated by the PGG algorithm
wget  -r -nH --cut-dirs=4 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~sroy/cichlids/networkinference/projected_speciesnets.tgz
tar -xvzf projected_speciesnets.tgz
# structure: col1 - TF geneID, col2 - target geneID, col3 - confidence (1 = high; 0 = low)

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
awk '$3>=0.5' tftgco_mergedA.txt > tftgco_merged.txt # THIS IS THE FINAL TF-TG_coexpression file - leaves 23113 lines

##########################################################################################################################
#
# Network Reconstruction
#
##########################################################################################################################

# 1. CREATE SPECIES-SPECIFIC FILES OF ALL EDGES
# 2. THEN CREATE MODULE-SPECIFIC FILES FROM THE SPECIES-SPECIFIC FILES

##########################################################################################################################

# 1. Create Edge Attributes

mkdir $homeWD/$5
EDGE=($homeWD/$4)
cd $EDGE

# create symbolic links of relevant files to this folder
ln -s $homeWD/$2/motifenr_merged-TFBSs_map2d.txt .
ln -s $homeWD/$4/tftgco_merged.txt .

#### Create scripts to collate the edge tables ####

#### 1a. Build initial Edge Attribute table with R

echo 'motifenr_merged_TFBSs_map2d <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/motifenr_merged-TFBSs_map2d.txt", header=T, na.strings=c(""," ","NA"))' > 1a_Edge_Attributes_Build.R
echo 'TFTGco <- read.delim("/tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/tftgco_merged.txt", header=T, na.strings=c(""," ","NA"))' >> 1a_Edge_Attributes_Build.R
printf '\n' >> 1a_Edge_Attributes_Build.R
echo 'EAB1a=plyr::rbind.fill(miRNAmodulemap_merged2, motifenr_merged_TFBSs_map2d)' >> 1a_Edge_Attributes_Build.R
echo 'EAB6a=plyr::rbind.fill(EAB1a, TFTGco)' >> 1a_Edge_Attributes_Build.R
printf '\n' >> 1a_Edge_Attributes_Build.R
echo 'write.table(EAB6a, "Edge_Attributes_Collated1.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)' >> 1a_Edge_Attributes_Build.R


#### 1b. Map gene IDs to gene symbols
geneNamesMapping=($homeWD/geneNamesMapping.txt)

echo '#!/bin/sh' > 1b_Gene_symbol_mapping.sh
echo '# create appropriate gene names mapping file' >> 1b_Gene_symbol_mapping.sh
echo "awk '{if("'$2 == "NONE")print $1, $2, "NONE";else print $1, $2, $3;}'"' OFS='\t' $geneNamesMapping > geneNamesMapping2.txt" >> 1b_Gene_symbol_mapping.sh
echo "# Cut two columns, one for GeneA and another for GeneB columns to map 'Gene Symbols' to IDs" >> 1b_Gene_symbol_mapping.sh
echo "awk '{print "'$1}'"' Edge_Attributes_Collated1.txt > Edge_Attributes_GeneA.txt" >> 1b_Gene_symbol_mapping.sh
echo "awk '{print "'$2}'"' Edge_Attributes_Collated1.txt > Edge_Attributes_GeneB.txt" >> 1b_Gene_symbol_mapping.sh
echo '#Compare first column from both files. If there is a match, print the corresponding value of the first column and the matched third column. If no match is found, fill with "NA" in second column, retaining ProteinID in first column and then edit column headers to "GeneA", "GeneA_Symbol", and "GeneB", "GeneB_Symbol".' >> 1b_Gene_symbol_mapping.sh
echo "awk 'BEGIN{OFS="'"\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}'"' geneNamesMapping2.txt Edge_Attributes_GeneA.txt | sed -e '1s/GeneA/GeneA/' -e '1s/NA/GeneA_Symbol/' > Edge_Attributes_GeneA_mapped.txt" >> 1b_Gene_symbol_mapping.sh
echo 'rm Edge_Attributes_GeneA.txt # remove the intermediate files' >> 1b_Gene_symbol_mapping.sh
echo "awk 'BEGIN{OFS="'"\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}'"' geneNamesMapping2.txt Edge_Attributes_GeneB.txt | sed -e '1s/GeneB/GeneB/' -e '1s/NA/GeneB_Symbol/' > Edge_Attributes_GeneB_mapped.txt" >> 1b_Gene_symbol_mapping.sh
echo 'rm Edge_Attributes_GeneB.txt # remove the intermediate files' >> 1b_Gene_symbol_mapping.sh
echo "#Remove three 'Gene_Symbol' columns from Edge_Atrributes_Collated1.txt file" >> 1b_Gene_symbol_mapping.sh
echo "cut -f1-20,22-24,26-37,39-48 Edge_Attributes_Collated1.txt > Edge_Attributes_Collated2.txt" >> 1b_Gene_symbol_mapping.sh


#### 1c. Finally build your Edge Attributes table

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
printf 'GeneA\tGeneA_Symbol\tGeneA_OGID\tGeneB\tGeneB_SymbolSp\tGeneB_SymbolDr\tGeneB_SymbolGa\tGeneB_OGID\tspecies\tinteraction_type\teffect\tinteraction\tdirectness\tdirection\tlayer\tsource\tmotif_pattern[pm]\tstart[pm]\tstop[pm]\tstrand[pm]\tscore[pm]\tp.value[pm]\tq.value[pm]\tsequence[pm]\tconf_level[pm]\tconf_score[pm]\tconfidence[coTFTG]\n' > colheads_edge

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
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/ab-speciesspecnames_clusterassign.txt Ab-Edge_Attributes_Collated4_GeneA.txt > Ab-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/ab-speciesspecnames_clusterassign.txt Ab-Edge_Attributes_Collated4_GeneB.txt > Ab-Edge_Attributes_Collated4_GeneBmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/nb-speciesspecnames_clusterassign.txt Nb-Edge_Attributes_Collated4_GeneA.txt > Nb-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/nb-speciesspecnames_clusterassign.txt Nb-Edge_Attributes_Collated4_GeneB.txt > Nb-Edge_Attributes_Collated4_GeneBmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/pn-speciesspecnames_clusterassign.txt Pn-Edge_Attributes_Collated4_GeneA.txt > Pn-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/pn-speciesspecnames_clusterassign.txt Pn-Edge_Attributes_Collated4_GeneB.txt > Pn-Edge_Attributes_Collated4_GeneBmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/mz-speciesspecnames_clusterassign.txt Mz-Edge_Attributes_Collated4_GeneA.txt > Mz-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/mz-speciesspecnames_clusterassign.txt Mz-Edge_Attributes_Collated4_GeneB.txt > Mz-Edge_Attributes_Collated4_GeneBmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/on-speciesspecnames_clusterassign.txt On-Edge_Attributes_Collated4_GeneA.txt > On-Edge_Attributes_Collated4_GeneAmap.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/on-speciesspecnames_clusterassign.txt On-Edge_Attributes_Collated4_GeneB.txt > On-Edge_Attributes_Collated4_GeneBmap.txt

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
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-ab Ab-Edge_Attributes_Collated4_GeneA.txt > Ab-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-ab Ab-Edge_Attributes_Collated4_GeneB.txt > Ab-Edge_Attributes_Collated4_GeneBmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-nb Nb-Edge_Attributes_Collated4_GeneA.txt > Nb-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-nb Nb-Edge_Attributes_Collated4_GeneB.txt > Nb-Edge_Attributes_Collated4_GeneBmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-pn Pn-Edge_Attributes_Collated4_GeneA.txt > Pn-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-pn Pn-Edge_Attributes_Collated4_GeneB.txt > Pn-Edge_Attributes_Collated4_GeneBmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-mz Mz-Edge_Attributes_Collated4_GeneA.txt > Mz-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-mz Mz-Edge_Attributes_Collated4_GeneB.txt > Mz-Edge_Attributes_Collated4_GeneBmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-on On-Edge_Attributes_Collated4_GeneA.txt > On-Edge_Attributes_Collated4_GeneAmap.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $homeWD/Module_genesandexpr/OGIDS.txt5-on On-Edge_Attributes_Collated4_GeneB.txt > On-Edge_Attributes_Collated4_GeneBmap.txt2

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
printf 'GeneA\tGeneA_Symbol\tGeneA_OGID\tGeneB\tGeneB_SymbolSp\tGeneB_SymbolDr\tGeneB_SymbolGa\tGeneB_OGID\tspecies\tinteraction_type\teffect\tinteraction\tdirectness\tdirection\tlayer\tsource\tmotif_pattern[pm]\tstart[pm]\tstop[pm]\tstrand[pm]\tscore[pm]\tp.value[pm]\tq.value[pm]\tsequence[pm]\tconf_level[pm]\tconf_score[pm]\tconfidence[coTFTG]\tGeneA_module\tGeneB_module\tGeneA_OGID\tGeneB_OGID\n' > colheaders_edgeattr2

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

## Add the following to a script to run the full Edge Attributes Build processes
echo '#!/bin/bash' > 0.edgeattributesbuild-Arboretum_GT_v3.sh
printf '\n' >> 0.edgeattributesbuild-Arboretum_GT_v3.sh
echo 'R CMD BATCH 1a_Edge_Attributes_Build.R # This will merge all the edge tables' >> 0.edgeattributesbuild-Arboretum_GT_v3.sh
echo 'sh 1b_Gene_symbol_mapping.sh # Run gene symbol mapping' >> 0.edgeattributesbuild-Arboretum_GT_v3.sh
echo 'sh 1c_Edge_Attributes_Build.sh #build your Edge Attributes tables with several commands	- split	by species and module too' >> 0.edgeattributesbuild-Arboretum_GT_v3.sh

# run the above
sh 0.edgeattributesbuild-Arboretum_GT_v3.sh

## Check lowest and highest number of columns in rows of a file to ensure no blank cells or spaces
for i in *-Edge_Attributes_Collated4c.txt ; do
	awk '{print NF}' $i | sort -nu | head -n 1 # lowest
done # 49 columns in each

for i in *-Edge_Attributes_Collated4c.txt ; do
	awk '{print NF}' $i | sort -nu | tail -n 1 # highest
done # 49 columns in each

### Based on tblastx, we determined that 4209 orthogroups with NULL orthogroup IDs for one or species are actually present in the genome (mis-annotations etc.) and not lost/absent
### Therefore, filter the geneB column of the main edge tables for the 4209 Present NULL OGIDs:

nano filteredges_presentNULLOGID.sh

#!/bin/bash

ls -1 *-Edge_Attributes_Collated4c.txt > sp_edgetable # create a list of all edge tables
mapfile -t sp_edgetable < sp_edgetable # assign as elements to $sp_edgetable variable

PresentNOG=($homeWD/MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS-noCand-noTF.txt2)

for ele in "${sp_edgetable[@]}" ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$9]){print $0,a[$9];}else{print $0,"KEEPME";}}' $PresentNOG "$ele" | grep 'KEEPME' | cut -f1-49 > "$(basename "$ele" 4c.txt)4d.txt"
done

# run the above
sh filteredges_presentNULLOGID.sh

## create module-specific files

nano moduleEDGES.sh

#!/bin/bash

# create a list file of species
echo Ab > list2
echo Nb >> list2
echo Pn >> list2
echo Mz >> list2
echo On >> list2

# Pull out according to module based on number presence in both module columns (4 and 10)

while read F ; do awk '{if (($4 == "0") && ($10 == "0")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module0.txt ; done < list2
while read F ; do awk '{if (($4 == "1") && ($10 == "1")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module1.txt ; done < list2
while read F ; do awk '{if (($4 == "2") && ($10 == "2")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module2.txt ; done < list2
while read F ; do awk '{if (($4 == "3") && ($10 == "3")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module3.txt ; done < list2
while read F ; do awk '{if (($4 == "4") && ($10 == "4")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module4.txt ; done < list2
while read F ; do awk '{if (($4 == "5") && ($10 == "5")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module5.txt ; done < list2
while read F ; do awk '{if (($4 == "6") && ($10 == "6")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module6.txt ; done < list2
while read F ; do awk '{if (($4 == "7") && ($10 == "7")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module7.txt ; done < list2
while read F ; do awk '{if (($4 == "8") && ($10 == "8")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module8.txt ; done < list2
while read F ; do awk '{if (($4 == "9") && ($10 == "9")) print $0}' $F-Edge_Attributes_Collated4d.txt > $F-Edge_Attributes_Collated4d-Module9.txt ; done < list2

# Add colheaders to all module-specific files
printf 'GeneA\tGeneA_Symbol\tGeneA_OGID\tGeneA_module\tGeneB\tGeneB_SymbolSp\tGeneB_SymbolDr\tGeneB_SymbolGa\tGeneB_OGID\tGeneB_module\tspecies\tinteraction_type\teffect\tinteraction\tdirectness\tdirection\tlayer\tsource\tmotif_pattern[pm]\tstart[pm]\tstop[pm]\tstrand[pm]\tscore[pm]\tp.value[pm]\tq.value[pm]\tsequence[pm]\tconf_level[pm]\tconf_score[pm]\tconfidence[coTFTG]\n' > colheaders_edgeattr3
INFO=$(cat colheaders_edgeattr3)  # read the contents of colheaders_edgeattr3 into var INFO
for i in *-Edge_Attributes_Collated4d-Module*.txt ; do sed -i "1i$INFO" $i ; done  # insert $INFO, use -i to change the file inplace
# Copy the new files to the workarea to copy local
for i in *-Edge_Attributes_Collated4d.txt ; do cp $i /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/ ; done
for i in *-Edge_Attributes_Collated4d-Module*.txt ; do cp $i /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/ ; done

### Final species edge files are [SpeciesID]-Edge_Attributes_Collated4d.txt
### Final species edge files are [SpeciesID]-Edge_Attributes_Collated4d-Module*.txt
