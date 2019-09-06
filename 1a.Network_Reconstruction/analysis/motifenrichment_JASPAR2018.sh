########################## Oct/Nov 2017 ###############################################################################################
#
# New motif enrichment analysis of new promoters of Arboretum module genes using JASPAR 2018 verterbate motifs release
#
##################################################################################################################################

# run this part local
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs
curl http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt > JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt
curl http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.txt > JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.txt

# check that all matrix rows sum to 1
awk '{for(i=t=0;i<NF;) t+=$++i; $0=t}1' JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt > JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme_SUM.txt ; paste -d'\t' JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme_SUM.txt | awk -F'\t' '$2!=1' | less # they all look fine now
# convert the TF names in the meme file all to uppercase to match the mapping file - for this you need to grab the second word after MOTIF in each line and convert to uppercase
awk '{for(i=1;i<=NF;i++) if ($i=="MOTIF") print $(i+2)=toupper;print}' OFS='\t' JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt | sed '/MOTIF/{N;s/\n.*//;}' > JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt2 # the long way - awk adds the line above the original so delete the original changed line retaining the uppercase with sed
mv JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt2 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt # rename file then copy to /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters_JASPAR2018

slurm
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/
mkdir promoters_JASPAR2018
cd promoters_JASPAR2018

# get the scripts from Sushmita group to run fimo, motif enrichment and plotting
software
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters_JASPAR2018
wget -r -nH --cut-dirs=3 --no-parent --reject="index.html*" -e robots=off http://pages.discovery.wisc.edu/~sknaack/cichlids/update_motif_enrichemnts_10-9-2017/cichlids_motif_enrichments.tgz # this contains all the scripts for running FIMO and generating plots from the output
exit

# unzip scripts - also includes all the data so delete that for redundancy
tar zxvf cichlids_motif_enrichments.tgz

# WD=(/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters_JASPAR2018)

# ============================================
# Step 1: generate FIMO motif finding results
# ============================================

uv2 # run on uv2
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters_JASPAR2018

# use FIMO in the MEME suite to scan for motif models using the JASPAR file across the genomes at a level of statistical significance (use stringent default pval <1e-4)

# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-Mz.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# M. zebra
fimo --o M_zebra_1e-4_JASPAR2018 --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt cichlids_motif_enrichments/data/promoter_sequences/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-Pn.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# P. nyererei
fimo --o P_nyererei_1e-4_JASPAR2018 --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt cichlids_motif_enrichments/data/promoter_sequences/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-Ab.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# A. burtoni
fimo --o A_burtoni_1e-4_JASPAR2018 --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt cichlids_motif_enrichments/data/promoter_sequences/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-Nb.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# N. brichardi
fimo --o N_brichardi_1e-4_JASPAR2018 --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt cichlids_motif_enrichments/data/promoter_sequences/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-On.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# O. niloticus
fimo --o O_niloticus_1e-4_JASPAR2018 --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt cichlids_motif_enrichments/data/promoter_sequences/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta


# Run all on uv2 - RERAN ALL 20/10/17 WITH NEW PROMOTER SEQUENCES!
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-Mz.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-Pn.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-Ab.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-Nb.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-On.sh #{DONE}

mkdir cichlids_motif_enrichments/results/FIMO/
mv *_1e-4_JASPAR2018 cichlids_motif_enrichments/results/FIMO/

# rename output files to run with motif enrichment script
nano list
O_niloticus_1e-4_JASPAR2018
N_brichardi_1e-4_JASPAR2018
A_burtoni_1e-4_JASPAR2018
P_nyererei_1e-4_JASPAR2018
M_zebra_1e-4_JASPAR2018

while read F ; do mv cichlids_motif_enrichments/results/FIMO/$F/fimo.txt cichlids_motif_enrichments/results/FIMO/$F/allfimo.txt ; done < list


######################## the next steps are in the cichlids_motif_enrichments/scripts/README file > you need to rename the outputs to "allfimo.txt"

# ============================================
# Step 2: refine the motif results across species.
# ============================================

# Here PrepAllFiles.sh is used to take the raw FIMO results and prepare the motifs for the enrichment analysis of the expression consensus modules for the cichlids species.
#
# Next, the ../code/go_filtering/ancestralAssign code is used to filter the motif instances for evolutionary conservation in a parsimony-based algorithm.
#
# The following files are used as input to this program:
#
# score_gainloss.txt - score matrix to compute cost of motif conservation or gain/loss across species
# SpeciesTree.txt - species tree of the cichlids species
# OGIDS_new_cichlid_only_latest_restored.txt - as I understand the latest OGIDS file from Christopher that I had available.
# SpeciesOrder.txt - list of species order of the genes given in the above OGIDS file.
# conf.txt - configuration input file for inputting the motif instances into the ancestralAssign program
#
# The input for ancestralAssign is saved in the ../results/prepRefinementFiles directory.
#
# The filtered motif instances are written out in the ../results/refined_motif_Ab directory.

cd cichlids_motif_enrichments/scripts/

nano PrepAllFiles.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

ml gcc
ml glib
ml glibc

D[0]=A_burtoni_1e-4_JASPAR2018
D[1]=G_aculeatus_1e-4_JASPAR2018
D[2]=M_zebra_1e-4_JASPAR2018
D[3]=N_brichardi_1e-4_JASPAR2018
D[4]=O_niloticus_1e-4_JASPAR2018
D[5]=P_nyererei_1e-4_JASPAR2018

S[0]=Ab
S[1]=Ga
S[2]=Mz
S[3]=Nb
S[4]=On
S[5]=Pn

mkdir ../results/prepRefinementFiles
rm conf.txt
for i in 0 2 3 4 5
do
	cat ../results/FIMO/${D[$i]}/allfimo.txt | awk '{printf("%s\t%s\n",$1,$2)}' | sort -u > ../results/prepRefinementFiles/${S[$i]}.txt
	printf "../results/prepRefinementFiles/${S[$i]}.txt\t${S[$i]}\n" >> conf.txt
done
mkdir ../results/refined_motifs_Ab

# Run all on uv2
qsub -q Test -l select=1:mem=20GB:ncpus=1 PrepAllFiles.sh #{DONE}

## copied relevant files to local to finish last step as glibc not updated on cluster (requires GLIBCXX_3.4.21)
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/scripts

# ancestralAssign code compiled with
# cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/code/ancestralAssign
# g++ *.C -g -o ancestralAssign

nano PrepAllFiles_local.sh

#!/bin/bash
mkdir ../results/refined_motifs_Ab/
export exe=../code/go_filtering/ancestralAssign
export SCORE=score_gainloss.txt
export OGIDS=OGIDS_new_cichlid_only_latest_restored.txt #/mnt/ws/sysbio/roygroup/shared/projects/cichlid_reg_net/results/Arboretum/OG_analysis/OGIDS_new_6Species.txt #/mnt/ws/sysbio/roygroup/shared/projects/cichlid_reg_net/results/Arboretum/OG_analysis/OGIDS_new_cichlid_only_new_annotations.txt

# eval ${exe} SpeciesOrder.txt <(sed 's/NONE/NULL/g' ${OGIDS}) SpeciesTree.txt conf.txt ${SCORE} ../results/refined_motifs_Ab/out_ Ab
${exe} SpeciesOrder.txt ${OGIDS} SpeciesTree.txt conf.txt ${SCORE} ../results/refined_motifs_Ab/out_ Ab

# run the above
sh PrepAllFiles_local.sh

# copy all files from /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/results/refined_motifs_Ab > to >
# /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/results/refined_motifs_Ab

# ============================================
# Step 3: rename motifs for the TF's they represent.
# ============================================

# The motif ID and TF name list is in ../data/JASPAR_Vertebrates_2018_motifname.txt. Here a program called selMerge is used to convert the motif ID's in the refined motifs sets to TF names. This is done in batchMapMotifNames.sh. The files with the renamed motifs are also in the ../results/refined_motif_Ab directory.

# continue local as next program requires it
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/data
# prepare the motifname file
grep MOTIF JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt | sed 's/MOTIF //g' | awk -F' ' '{print $1,$2}' OFS='\t' > JASPAR_Vertebrates_2018_motifname.txt

# selMerge code compiled with
# cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/code/selmerge
# g++ *.C -g -o selMerge

cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/scripts

nano batchMapMotifNames_local.sh

#!/bin/bash
MOTIFNAME=../data/JASPAR_Vertebrates_2018_motifname.txt
for SPECIES in Ab Mz Nb On Pn Anc1 Anc2 Anc3 Anc4
do
	export MNAME=../results/refined_motifs_Ab/out__${SPECIES}_regnet.txt
	cut -f2 $MNAME > temp.txt
	echo "selMerge temp.txt $MOTIFNAME temp_inst.txt"
	../code/selmerge/selMerge temp.txt $MOTIFNAME temp_inst.txt
	head temp_inst.txt
	echo "paste temp_inst.txt $MNAME | awk '{printf("%s\t%s\n",$3,$2)}' > ../results/refined_motifs_Ab/${SPECIES}_motifnames_regnet.txt"
	paste temp_inst.txt $MNAME | awk '{printf("%s\t%s\n",$3,$2)}' > ../results/refined_motifs_Ab/${SPECIES}_motifnames_regnet.txt
done

rm temp*

# run the above
sh batchMapMotifNames_local.sh

# ============================================
# Step 4: Prepare the filtered/refined and renamed motif instances as inputs for the enrichment program.
# ============================================

# This is done with the PrepEAFiles.sh script and the final motif annotation files are in ../results/refined_motif_inputs.

cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/scripts

nano PrepEAfiles_local.sh

#!/bin/bash
mkdir ../results/refined_motif_inputs

for SPC in Ab Pn On Mz Nb Anc1 Anc2 Anc3 Anc4
do
	export NET=../results/refined_motifs_Ab/${SPC}_motifnames_regnet.txt
	#head ${NET}
	cat ${NET} | grep -v NULL | awk 'NF==2{printf("%s\t%s\t1\n",$1,$2)}' > ../results/refined_motif_inputs/${SPC}_gotermap.txt
	cat ${NET} | grep -v NULL  | cut -f2 | sort | uniq -c | head
	cat ${NET} | grep -v NULL  | cut -f2 | sort | uniq -c | awk '{printf("%s\t%s\t1\n",$2,$1)}' > ../results/refined_motif_inputs/${SPC}_genecnt.txt
done

# for local unix
rename 's/Ab_/Asbu_/g' ../results/refined_motif_inputs/Ab_*
rename 's/Mz_/Meze_/g' ../results/refined_motif_inputs/Mz_*
rename 's/Pn_/Puny_/g' ../results/refined_motif_inputs/Pn_*
rename 's/On_/Orni_/g' ../results/refined_motif_inputs/On_*
rename 's/Nb_/Nebr_/g' ../results/refined_motif_inputs/Nb_*

# for the cluster
# rename Ab_ Asbu_ ../results/refined_motif_inputs/Ab_*
# rename Mz_ Meze_ ../results/refined_motif_inputs/Mz_*
# rename Pn_ Puny_ ../results/refined_motif_inputs/Pn_*
# rename On_ Orni_ ../results/refined_motif_inputs/On_*
# rename Nb_ Nebr_ ../results/refined_motif_inputs/Nb_*

# run the above
sh PrepEAfiles_local.sh

# ============================================
# Step 5: Run the enrichment analysis
# ============================================

# The tool for applying the enrichment analysis is ../code/enrichanalyzer_Nongraph_Qval/enrichAnalyzer, and this is implemented with the steps in the EnrichmentAnalysis.sh

# The clustering results are in ../results/k10_arb_output/prediction, and the inputs and outputs for the enrichment program are in ../results/Motif_enrichments/Output.

# compile relevant codes
# cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/code/mergegoattrib
# g++ *.C -g -o mergeEnrichment

# cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/code/enrichanalyzer_Nongraph_Qval
# make
# the above enrichAnalyzer did not 'make' so using another version here that is already executable /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer


nano EnrichmentAnalysis_local.sh

#!/bin/bash
export GOEXE=/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer
export GOMERGOEXE=../code/mergegoattrib/mergeEnrichment
export PYGOINPUT=ReworkClusterAssign.py

#changed the following lines manually
export Species_List_All="Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4"
export CSPECIES=Ab
export OUTDIR=../results/k10_arb_output/prediction
export EAOUTDIR=../results/Motif_enrichments/Output

echo -e "Starting GO analysis"
mkdir -p ${EAOUTDIR}

for SN in Meze Puny Asbu Nebr Orni
do
        CAFILE=${OUTDIR}/${SN}_speciesspecnames_clusterassign.txt
        EAINFILE=${EAOUTDIR}/${SN}_MOTIFS_INPUT.txt
        python ${PYGOINPUT} ${CAFILE} > ${EAINFILE}
done

for SN in Anc1 Anc2 Anc3 Anc4
do
        CAFILE=${OUTDIR}/${SN}_clusterassign.txt
        EAINFILE=${EAOUTDIR}/${SN}_MOTIFS_INPUT.txt
        python ${PYGOINPUT} ${CAFILE} > ${EAINFILE}
done

for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
do
	export EAINFILE=${EAOUTDIR}/${SN}_MOTIFS_INPUT.txt
	export EAOUTPUT=${EAOUTDIR}/${SN}_MOTIFS_OUTPUT
	export CAFILE=${OUTDIR}/${SN}_clusterassign.txt
	if [[ -e ${OUTDIR}/${SN}_speciesspecnames_clusterassign.txt ]]
	then
		CAFILE=$(echo ${OUTDIR}/${SN}_speciesspecnames_clusterassign.txt | cut -f1 | sort -u)
	fi
	EADAT=../results/refined_motif_inputs/${SN}_
        ${GOEXE} ${EAINFILE} ${CAFILE} ${EADAT} 1 ${EAOUTPUT} persg
done

# run the above
sh EnrichmentAnalysis_local.sh


## Run enrichment analysis of single tissue modules

# Create relevant files by renaming based on longer species name
# cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results
# for i in br ey ht kd ms ts; do
#   pax -wrs'/Mz-speciesspecnames_clusterassign.txt/Meze_speciesspecnames_clusterassign.txt/' $i/Mz-speciesspecnames_clusterassign.txt .
#   pax -wrs'/Pn-speciesspecnames_clusterassign.txt/Puny_speciesspecnames_clusterassign.txt/' $i/Pn-speciesspecnames_clusterassign.txt .
#   pax -wrs'/Ab-speciesspecnames_clusterassign.txt/Asbu_speciesspecnames_clusterassign.txt/' $i/Ab-speciesspecnames_clusterassign.txt .
#   pax -wrs'/Nb-speciesspecnames_clusterassign.txt/Nebr_speciesspecnames_clusterassign.txt/' $i/Nb-speciesspecnames_clusterassign.txt .
#   pax -wrs'/On-speciesspecnames_clusterassign.txt/Orni_speciesspecnames_clusterassign.txt/' $i/On-speciesspecnames_clusterassign.txt .
#   pax -wrs'/Mz-clusterassign.txt/Meze_clusterassign.txt/' $i/Mz-clusterassign.txt .
#   pax -wrs'/Pn-clusterassign.txt/Puny_clusterassign.txt/' $i/Pn-clusterassign.txt .
#   pax -wrs'/Ab-clusterassign.txt/Asbu_clusterassign.txt/' $i/Ab-clusterassign.txt .
#   pax -wrs'/Nb-clusterassign.txt/Nebr_clusterassign.txt/' $i/Nb-clusterassign.txt .
#   pax -wrs'/On-clusterassign.txt/Orni_clusterassign.txt/' $i/On-clusterassign.txt .
# done

cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/scripts

nano EnrichmentAnalysis_local_ST.sh

#!/bin/bash
for i in br ey ht kd ms ts; do
  export GOEXE=/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer
  export GOMERGOEXE=../code/mergegoattrib/mergeEnrichment
  export PYGOINPUT=ReworkClusterAssign.py
  #changed the following lines manually
  export Species_List_All="Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4"
  export CSPECIES=Ab
  export OUTDIR=/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/$i
  export EAOUTDIR=../results/Motif_enrichments/Output/$i
  echo -e "Starting GO analysis"
  mkdir -p ${EAOUTDIR}
  for SN in Meze Puny Asbu Nebr Orni; do
    CAFILE=${OUTDIR}/${SN}_speciesspecnames_clusterassign.txt
    EAINFILE=${EAOUTDIR}/${SN}_MOTIFS_INPUT.txt
    python ${PYGOINPUT} ${CAFILE} > ${EAINFILE}
    for SN in Anc1 Anc2 Anc3 Anc4; do
      CAFILE=${OUTDIR}/${SN}_clusterassign.txt
      EAINFILE=${EAOUTDIR}/${SN}_MOTIFS_INPUT.txt
      python ${PYGOINPUT} ${CAFILE} > ${EAINFILE}
      for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4; do
        export EAINFILE=${EAOUTDIR}/${SN}_MOTIFS_INPUT.txt
        export EAOUTPUT=${EAOUTDIR}/${SN}_MOTIFS_OUTPUT
        export CAFILE=${OUTDIR}/${SN}_clusterassign.txt
        if [[ -e ${OUTDIR}/${SN}_speciesspecnames_clusterassign.txt ]]
        then
          CAFILE=$(echo ${OUTDIR}/${SN}_speciesspecnames_clusterassign.txt | cut -f1 | sort -u)
        fi
        EADAT=../results/refined_motif_inputs/${SN}_
              ${GOEXE} ${EAINFILE} ${CAFILE} ${EADAT} 1 ${EAOUTPUT} persg
            done
          done
        done
      done

# run the above
bash EnrichmentAnalysis_local_ST.sh


# ============================================
# Step 6a: Prepare enrichment plots
# ============================================

# The enrichment results for motifs in the promoter sequences are in ../results/Motif_enrichments/Output, and the MakeEnrichmentPlots_k10.sh script generate the enrichment plots prepared in ../results/Motif_enrichments/Figs_out

# needed to change the shebang at the top of Heatmap.awk and TT.awk to #!/usr/bin/env gawk -f so that certain functions, unique to gawk (gensub) can be used

# Several of the plots are made with a matrix that is generated from merging the enrichment results with the ../code/mergegoattrib/mergeEnrichment tool. The conf2.txt file is an input to this program and out.txt is the output matrix summarizing the enrichments across species.

# These plots are then the end of the analysis pipeline. There is a separate README in ../results/Motif_enrichments/Figs_out to describe those plots.

# this is for the multi-tissue k10 modules
nano MakeEnrichmentPlots_k10_local.sh

#!/bin/bash
GOOUTDIR=../results/Motif_enrichments/Output
outdir=../results/Motif_enrichments/Figs_out
mkdir -p ${outdir}
LBL=k10_r10

rm out.txt list_0.05.txt list_0.001.txt list_Top5.txt list_Extant_0.05.txt list_Extant_0.001.txt

for c in {0..9}
do
	export C=$c
	export i=$c
	printf "$c\t$i\n"
	export D=$C
	for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
	do
		export SPC=${SN}
		echo ${SN}
		export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | sort -k2g | awk '$4<=0.05{printf("%s:%s\t%s\t%.1f|%d|%.1f\n",$2,ENVIRON["D"],ENVIRON["SPC"],-1*log($4)/log(10),ENVIRON["D"]+1,-1*log($4)/log(10))}' >> out.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.05.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.001.txt
	done
	for SN in Meze Puny Asbu Nebr Orni
        do
                export SPC=${SN}
		echo ${SN}
                export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.05.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.001.txt
	done

	if [[ ${c} -lt 9 ]]
	then
		echo | awk '{print "|Spacer||mu" ENVIRON["C"], "|-"}' >> out.txt
	fi
done

export GOMERGOEXE=../code/mergegoattrib/mergeEnrichment
export HMAWK=../code/Heatmap.awk
export TTAWK=../code/TT.awk

# prepare color map for 10 clusters
R[0]=0
R[1]=56
R[2]=113
R[3]=170
R[4]=226
R[5]=246
R[6]=228
R[7]=210
R[8]=192
R[9]=175

G[0]=175
G[1]=192
G[2]=210
G[3]=228
G[4]=246
G[5]=226
G[6]=170
G[7]=113
G[8]=56
G[9]=0

COLORMAPS=`printf "1.3:(255,255,255);5:(${R[0]},${G[0]},0)"`
for x in 1 2 3 4 5 6 7 8 9
do
	COLORMAPS=`printf "${COLORMAPS} 1.3:(255,255,255);5:(${R[$x]},${G[$x]},0)"`
done

echo ${COLORMAPS}

cat out.txt | grep -F -w -f <(cat list_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4 Cluster5 Cluster6 Cluster7 Cluster8 Cluster9" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.05.svg

cat out.txt | grep -F -w -f <(cat list_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4 Cluster5 Cluster6 Cluster7 Cluster8 Cluster9" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.001.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4 Cluster5 Cluster6 Cluster7 Cluster8 Cluster9" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.05.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4 Cluster5 Cluster6 Cluster7 Cluster8 Cluster9" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.001.svg

cat list_0.05.txt | sort | uniq -c | awk '$1==9{print $2}' > list_Conserved_0.05.txt

cat out.txt | grep -F -w -f <(cat list_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4 Cluster5 Cluster6 Cluster7 Cluster8 Cluster9" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Conserved_q0.05.svg


cat list_Extant_0.05.txt | sort | uniq -c | awk '$1==5{print $2 }' >  list_Extant_Conserved_0.05.txt

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4 Cluster5 Cluster6 Cluster7 Cluster8 Cluster9" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_Conserved_q0.05.svg

rm list* out.txt

rm conf2.txt

for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
do
	awk '$4<0.05{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt > ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt
	printf "${SN}\t${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt\n" >> conf2.txt
done

export Species_List="Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4"
export CSPECIES="Meze"
GOMERGEOUTPUT=${GOOUTDIR}/merged_MOTIFS_OUTPUT_
eval ../code/mergegoattrib/mergeEnrichment conf2.txt ${GOMERGEOUTPUT}
TWISE=${GOMERGEOUTPUT}_termwise.txt
CWISE=${GOMERGEOUTPUT}_clusterwise.txt
TPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_termwise.svg
CPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_clusterwise.svg
CCPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_conserved_clusterwise.svg
COLORMAP="-1:(200,200,200);0:(0,175,0);4.5:(255,255,0);9:(175,0,0)"
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${CWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="11" -v L="Cluster" > ${CPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${TWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="11" -v L="Cluster" > ${TPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 <(grep -v "\-1" ${CWISE}) | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="11" -v L="Cluster" > ${CCPLOT}

# run the above (run with 'bash' so that it calls bash directly and does not invoke posix which does not allow '<(XX)' commands - throws a syntax error )
bash MakeEnrichmentPlots_k10_local.sh

# ============================================
# Step 6b: Prepare enrichment plots for Single Tissue modules
# ============================================

# this is for the single tissue modules
nano MakeEnrichmentPlots_ST_local.sh


#!/bin/bash
GOOUTDIR=../results/Motif_enrichments/Output/br
outdir=../results/Motif_enrichments/Figs_out/br
mkdir -p ${outdir}
LBL=k5_br

rm out.txt list_0.05.txt list_0.001.txt list_Top5.txt list_Extant_0.05.txt list_Extant_0.001.txt

for c in {0..4}
do
	export C=$c
	export i=$c
	printf "$c\t$i\n"
	export D=$C
	for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
	do
		export SPC=${SN}
		echo ${SN}
		export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | sort -k2g | awk '$4<=0.05{printf("%s:%s\t%s\t%.1f|%d|%.1f\n",$2,ENVIRON["D"],ENVIRON["SPC"],-1*log($4)/log(10),ENVIRON["D"]+1,-1*log($4)/log(10))}' >> out.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.05.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.001.txt
	done
	for SN in Meze Puny Asbu Nebr Orni
        do
                export SPC=${SN}
		echo ${SN}
                export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.05.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.001.txt
	done

	if [[ ${c} -lt 9 ]]
	then
		echo | awk '{print "|Spacer||mu" ENVIRON["C"], "|-"}' >> out.txt
	fi
done

export GOMERGOEXE=../code/mergegoattrib/mergeEnrichment
export HMAWK=../code/Heatmap.awk
export TTAWK=../code/TT.awk

# prepare color map for 5 clusters
R[0]=0
R[1]=56
R[2]=113
R[3]=170
R[4]=226

G[0]=175
G[1]=192
G[2]=210
G[3]=228
G[4]=246

COLORMAPS=`printf "1.3:(255,255,255);5:(${R[0]},${G[0]},0)"`
for x in 1 2 3 4
do
	COLORMAPS=`printf "${COLORMAPS} 1.3:(255,255,255);5:(${R[$x]},${G[$x]},0)"`
done

echo ${COLORMAPS}

cat out.txt | grep -F -w -f <(cat list_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.05.svg

cat out.txt | grep -F -w -f <(cat list_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.001.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.05.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.001.svg

cat list_0.05.txt | sort | uniq -c | awk '$1==9{print $2}' > list_Conserved_0.05.txt

cat out.txt | grep -F -w -f <(cat list_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Conserved_q0.05.svg


cat list_Extant_0.05.txt | sort | uniq -c | awk '$1==5{print $2 }' >  list_Extant_Conserved_0.05.txt

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_Conserved_q0.05.svg

rm list* out.txt

rm conf2.txt

for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
do
	awk '$4<0.05{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt > ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt
	printf "${SN}\t${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt\n" >> conf2.txt
done

export Species_List="Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4"
export CSPECIES="Meze"
GOMERGEOUTPUT=${GOOUTDIR}/merged_MOTIFS_OUTPUT_
eval ../code/mergegoattrib/mergeEnrichment conf2.txt ${GOMERGEOUTPUT}
TWISE=${GOMERGEOUTPUT}_termwise.txt
CWISE=${GOMERGEOUTPUT}_clusterwise.txt
TPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_termwise.svg
CPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_clusterwise.svg
CCPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_conserved_clusterwise.svg
COLORMAP="-1:(200,200,200);0:(0,175,0);2:(255,255,0);4:(175,0,0)"
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${CWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${TWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${TPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 <(grep -v "\-1" ${CWISE}) | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CCPLOT}


GOOUTDIR=../results/Motif_enrichments/Output/ey
outdir=../results/Motif_enrichments/Figs_out/ey
mkdir -p ${outdir}
LBL=k5_ey

rm out.txt list_0.05.txt list_0.001.txt list_Top5.txt list_Extant_0.05.txt list_Extant_0.001.txt

for c in {0..4}
do
	export C=$c
	export i=$c
	printf "$c\t$i\n"
	export D=$C
	for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
	do
		export SPC=${SN}
		echo ${SN}
		export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | sort -k2g | awk '$4<=0.05{printf("%s:%s\t%s\t%.1f|%d|%.1f\n",$2,ENVIRON["D"],ENVIRON["SPC"],-1*log($4)/log(10),ENVIRON["D"]+1,-1*log($4)/log(10))}' >> out.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.05.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.001.txt
	done
	for SN in Meze Puny Asbu Nebr Orni
        do
                export SPC=${SN}
		echo ${SN}
                export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.05.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.001.txt
	done

	if [[ ${c} -lt 9 ]]
	then
		echo | awk '{print "|Spacer||mu" ENVIRON["C"], "|-"}' >> out.txt
	fi
done

export GOMERGOEXE=../code/mergegoattrib/mergeEnrichment
export HMAWK=../code/Heatmap.awk
export TTAWK=../code/TT.awk

# prepare color map for 5 clusters
R[0]=0
R[1]=56
R[2]=113
R[3]=170
R[4]=226

G[0]=175
G[1]=192
G[2]=210
G[3]=228
G[4]=246

COLORMAPS=`printf "1.3:(255,255,255);5:(${R[0]},${G[0]},0)"`
for x in 1 2 3 4
do
	COLORMAPS=`printf "${COLORMAPS} 1.3:(255,255,255);5:(${R[$x]},${G[$x]},0)"`
done

echo ${COLORMAPS}

cat out.txt | grep -F -w -f <(cat list_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.05.svg

cat out.txt | grep -F -w -f <(cat list_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.001.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.05.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.001.svg

cat list_0.05.txt | sort | uniq -c | awk '$1==9{print $2}' > list_Conserved_0.05.txt

cat out.txt | grep -F -w -f <(cat list_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Conserved_q0.05.svg


cat list_Extant_0.05.txt | sort | uniq -c | awk '$1==5{print $2 }' >  list_Extant_Conserved_0.05.txt

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_Conserved_q0.05.svg

rm list* out.txt

rm conf2.txt

for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
do
	awk '$4<0.05{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt > ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt
	printf "${SN}\t${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt\n" >> conf2.txt
done

export Species_List="Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4"
export CSPECIES="Meze"
GOMERGEOUTPUT=${GOOUTDIR}/merged_MOTIFS_OUTPUT_
eval ../code/mergegoattrib/mergeEnrichment conf2.txt ${GOMERGEOUTPUT}
TWISE=${GOMERGEOUTPUT}_termwise.txt
CWISE=${GOMERGEOUTPUT}_clusterwise.txt
TPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_termwise.svg
CPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_clusterwise.svg
CCPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_conserved_clusterwise.svg
COLORMAP="-1:(200,200,200);0:(0,175,0);2:(255,255,0);4:(175,0,0)"
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${CWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${TWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${TPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 <(grep -v "\-1" ${CWISE}) | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CCPLOT}


GOOUTDIR=../results/Motif_enrichments/Output/ht
outdir=../results/Motif_enrichments/Figs_out/ht
mkdir -p ${outdir}
LBL=k5_ht

rm out.txt list_0.05.txt list_0.001.txt list_Top5.txt list_Extant_0.05.txt list_Extant_0.001.txt

for c in {0..4}
do
	export C=$c
	export i=$c
	printf "$c\t$i\n"
	export D=$C
	for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
	do
		export SPC=${SN}
		echo ${SN}
		export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | sort -k2g | awk '$4<=0.05{printf("%s:%s\t%s\t%.1f|%d|%.1f\n",$2,ENVIRON["D"],ENVIRON["SPC"],-1*log($4)/log(10),ENVIRON["D"]+1,-1*log($4)/log(10))}' >> out.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.05.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.001.txt
	done
	for SN in Meze Puny Asbu Nebr Orni
        do
                export SPC=${SN}
		echo ${SN}
                export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.05.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.001.txt
	done

	if [[ ${c} -lt 9 ]]
	then
		echo | awk '{print "|Spacer||mu" ENVIRON["C"], "|-"}' >> out.txt
	fi
done

export GOMERGOEXE=../code/mergegoattrib/mergeEnrichment
export HMAWK=../code/Heatmap.awk
export TTAWK=../code/TT.awk

# prepare color map for 5 clusters
R[0]=0
R[1]=56
R[2]=113
R[3]=170
R[4]=226

G[0]=175
G[1]=192
G[2]=210
G[3]=228
G[4]=246

COLORMAPS=`printf "1.3:(255,255,255);5:(${R[0]},${G[0]},0)"`
for x in 1 2 3 4
do
	COLORMAPS=`printf "${COLORMAPS} 1.3:(255,255,255);5:(${R[$x]},${G[$x]},0)"`
done

echo ${COLORMAPS}

cat out.txt | grep -F -w -f <(cat list_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.05.svg

cat out.txt | grep -F -w -f <(cat list_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.001.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.05.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.001.svg

cat list_0.05.txt | sort | uniq -c | awk '$1==9{print $2}' > list_Conserved_0.05.txt

cat out.txt | grep -F -w -f <(cat list_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Conserved_q0.05.svg


cat list_Extant_0.05.txt | sort | uniq -c | awk '$1==5{print $2 }' >  list_Extant_Conserved_0.05.txt

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_Conserved_q0.05.svg

rm list* out.txt

rm conf2.txt

for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
do
	awk '$4<0.05{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt > ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt
	printf "${SN}\t${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt\n" >> conf2.txt
done

export Species_List="Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4"
export CSPECIES="Meze"
GOMERGEOUTPUT=${GOOUTDIR}/merged_MOTIFS_OUTPUT_
eval ../code/mergegoattrib/mergeEnrichment conf2.txt ${GOMERGEOUTPUT}
TWISE=${GOMERGEOUTPUT}_termwise.txt
CWISE=${GOMERGEOUTPUT}_clusterwise.txt
TPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_termwise.svg
CPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_clusterwise.svg
CCPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_conserved_clusterwise.svg
COLORMAP="-1:(200,200,200);0:(0,175,0);2:(255,255,0);4:(175,0,0)"
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${CWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${TWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${TPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 <(grep -v "\-1" ${CWISE}) | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CCPLOT}


GOOUTDIR=../results/Motif_enrichments/Output/kd
outdir=../results/Motif_enrichments/Figs_out/kd
mkdir -p ${outdir}
LBL=k5_kd

rm out.txt list_0.05.txt list_0.001.txt list_Top5.txt list_Extant_0.05.txt list_Extant_0.001.txt

for c in {0..4}
do
	export C=$c
	export i=$c
	printf "$c\t$i\n"
	export D=$C
	for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
	do
		export SPC=${SN}
		echo ${SN}
		export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | sort -k2g | awk '$4<=0.05{printf("%s:%s\t%s\t%.1f|%d|%.1f\n",$2,ENVIRON["D"],ENVIRON["SPC"],-1*log($4)/log(10),ENVIRON["D"]+1,-1*log($4)/log(10))}' >> out.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.05.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.001.txt
	done
	for SN in Meze Puny Asbu Nebr Orni
        do
                export SPC=${SN}
		echo ${SN}
                export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.05.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.001.txt
	done

	if [[ ${c} -lt 9 ]]
	then
		echo | awk '{print "|Spacer||mu" ENVIRON["C"], "|-"}' >> out.txt
	fi
done

export GOMERGOEXE=../code/mergegoattrib/mergeEnrichment
export HMAWK=../code/Heatmap.awk
export TTAWK=../code/TT.awk

# prepare color map for 5 clusters
R[0]=0
R[1]=56
R[2]=113
R[3]=170
R[4]=226

G[0]=175
G[1]=192
G[2]=210
G[3]=228
G[4]=246

COLORMAPS=`printf "1.3:(255,255,255);5:(${R[0]},${G[0]},0)"`
for x in 1 2 3 4
do
	COLORMAPS=`printf "${COLORMAPS} 1.3:(255,255,255);5:(${R[$x]},${G[$x]},0)"`
done

echo ${COLORMAPS}

cat out.txt | grep -F -w -f <(cat list_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.05.svg

cat out.txt | grep -F -w -f <(cat list_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.001.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.05.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.001.svg

cat list_0.05.txt | sort | uniq -c | awk '$1==9{print $2}' > list_Conserved_0.05.txt

cat out.txt | grep -F -w -f <(cat list_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Conserved_q0.05.svg


cat list_Extant_0.05.txt | sort | uniq -c | awk '$1==5{print $2 }' >  list_Extant_Conserved_0.05.txt

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_Conserved_q0.05.svg

rm list* out.txt

rm conf2.txt

for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
do
	awk '$4<0.05{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt > ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt
	printf "${SN}\t${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt\n" >> conf2.txt
done

export Species_List="Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4"
export CSPECIES="Meze"
GOMERGEOUTPUT=${GOOUTDIR}/merged_MOTIFS_OUTPUT_
eval ../code/mergegoattrib/mergeEnrichment conf2.txt ${GOMERGEOUTPUT}
TWISE=${GOMERGEOUTPUT}_termwise.txt
CWISE=${GOMERGEOUTPUT}_clusterwise.txt
TPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_termwise.svg
CPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_clusterwise.svg
CCPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_conserved_clusterwise.svg
COLORMAP="-1:(200,200,200);0:(0,175,0);2:(255,255,0);4:(175,0,0)"
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${CWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${TWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${TPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 <(grep -v "\-1" ${CWISE}) | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CCPLOT}


GOOUTDIR=../results/Motif_enrichments/Output/ms
outdir=../results/Motif_enrichments/Figs_out/ms
mkdir -p ${outdir}
LBL=k5_ms

rm out.txt list_0.05.txt list_0.001.txt list_Top5.txt list_Extant_0.05.txt list_Extant_0.001.txt

for c in {0..4}
do
	export C=$c
	export i=$c
	printf "$c\t$i\n"
	export D=$C
	for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
	do
		export SPC=${SN}
		echo ${SN}
		export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | sort -k2g | awk '$4<=0.05{printf("%s:%s\t%s\t%.1f|%d|%.1f\n",$2,ENVIRON["D"],ENVIRON["SPC"],-1*log($4)/log(10),ENVIRON["D"]+1,-1*log($4)/log(10))}' >> out.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.05.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.001.txt
	done
	for SN in Meze Puny Asbu Nebr Orni
        do
                export SPC=${SN}
		echo ${SN}
                export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.05.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.001.txt
	done

	if [[ ${c} -lt 9 ]]
	then
		echo | awk '{print "|Spacer||mu" ENVIRON["C"], "|-"}' >> out.txt
	fi
done

export GOMERGOEXE=../code/mergegoattrib/mergeEnrichment
export HMAWK=../code/Heatmap.awk
export TTAWK=../code/TT.awk

# prepare color map for 5 clusters
R[0]=0
R[1]=56
R[2]=113
R[3]=170
R[4]=226

G[0]=175
G[1]=192
G[2]=210
G[3]=228
G[4]=246

COLORMAPS=`printf "1.3:(255,255,255);5:(${R[0]},${G[0]},0)"`
for x in 1 2 3 4
do
	COLORMAPS=`printf "${COLORMAPS} 1.3:(255,255,255);5:(${R[$x]},${G[$x]},0)"`
done

echo ${COLORMAPS}

cat out.txt | grep -F -w -f <(cat list_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.05.svg

cat out.txt | grep -F -w -f <(cat list_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.001.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.05.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.001.svg

cat list_0.05.txt | sort | uniq -c | awk '$1==9{print $2}' > list_Conserved_0.05.txt

cat out.txt | grep -F -w -f <(cat list_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Conserved_q0.05.svg


cat list_Extant_0.05.txt | sort | uniq -c | awk '$1==5{print $2 }' >  list_Extant_Conserved_0.05.txt

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_Conserved_q0.05.svg

rm list* out.txt

rm conf2.txt

for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
do
	awk '$4<0.05{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt > ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt
	printf "${SN}\t${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt\n" >> conf2.txt
done

export Species_List="Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4"
export CSPECIES="Meze"
GOMERGEOUTPUT=${GOOUTDIR}/merged_MOTIFS_OUTPUT_
eval ../code/mergegoattrib/mergeEnrichment conf2.txt ${GOMERGEOUTPUT}
TWISE=${GOMERGEOUTPUT}_termwise.txt
CWISE=${GOMERGEOUTPUT}_clusterwise.txt
TPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_termwise.svg
CPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_clusterwise.svg
CCPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_conserved_clusterwise.svg
COLORMAP="-1:(200,200,200);0:(0,175,0);2:(255,255,0);4:(175,0,0)"
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${CWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${TWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${TPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 <(grep -v "\-1" ${CWISE}) | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CCPLOT}


GOOUTDIR=../results/Motif_enrichments/Output/ts
outdir=../results/Motif_enrichments/Figs_out/ts
mkdir -p ${outdir}
LBL=k5_ts

rm out.txt list_0.05.txt list_0.001.txt list_Top5.txt list_Extant_0.05.txt list_Extant_0.001.txt

for c in {0..4}
do
	export C=$c
	export i=$c
	printf "$c\t$i\n"
	export D=$C
	for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
	do
		export SPC=${SN}
		echo ${SN}
		export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | sort -k2g | awk '$4<=0.05{printf("%s:%s\t%s\t%.1f|%d|%.1f\n",$2,ENVIRON["D"],ENVIRON["SPC"],-1*log($4)/log(10),ENVIRON["D"]+1,-1*log($4)/log(10))}' >> out.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.05.txt
		grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<=0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_0.001.txt
	done
	for SN in Meze Puny Asbu Nebr Orni
        do
                export SPC=${SN}
		echo ${SN}
                export GOOUTPUT=${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.05{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.05.txt
                grep -F -w Cluster${i} ${GOOUTPUT} | awk '$4<0.001{printf("%s:%s\n",$2,ENVIRON["D"])}' >> list_Extant_0.001.txt
	done

	if [[ ${c} -lt 9 ]]
	then
		echo | awk '{print "|Spacer||mu" ENVIRON["C"], "|-"}' >> out.txt
	fi
done

export GOMERGOEXE=../code/mergegoattrib/mergeEnrichment
export HMAWK=../code/Heatmap.awk
export TTAWK=../code/TT.awk

# prepare color map for 5 clusters
R[0]=0
R[1]=56
R[2]=113
R[3]=170
R[4]=226

G[0]=175
G[1]=192
G[2]=210
G[3]=228
G[4]=246

COLORMAPS=`printf "1.3:(255,255,255);5:(${R[0]},${G[0]},0)"`
for x in 1 2 3 4
do
	COLORMAPS=`printf "${COLORMAPS} 1.3:(255,255,255);5:(${R[$x]},${G[$x]},0)"`
done

echo ${COLORMAPS}

cat out.txt | grep -F -w -f <(cat list_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.05.svg

cat out.txt | grep -F -w -f <(cat list_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_q0.001.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.05.svg

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_0.001.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_q0.001.svg

cat list_0.05.txt | sort | uniq -c | awk '$1==9{print $2}' > list_Conserved_0.05.txt

cat out.txt | grep -F -w -f <(cat list_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Conserved_q0.05.svg


cat list_Extant_0.05.txt | sort | uniq -c | awk '$1==5{print $2 }' >  list_Extant_Conserved_0.05.txt

cat out.txt | grep -v Anc | grep -F -w -f <(cat list_Extant_Conserved_0.05.txt <(printf "Spacer\n") ) | ${HMAWK} -v L="Cluster0 Cluster1 Cluster2 Cluster3 Cluster4" -v D="space" -v LegendNum="5" -v C="${COLORMAPS}" > ${outdir}/${LBL}_MOTIFS_Extant_Conserved_q0.05.svg

rm list* out.txt

rm conf2.txt

for SN in Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4
do
	awk '$4<0.05{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_details.txt > ${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt
	printf "${SN}\t${GOOUTDIR}/${SN}_MOTIFS_OUTPUT_q0.05_details.txt\n" >> conf2.txt
done

export Species_List="Meze Puny Anc1 Asbu Anc2 Nebr Anc3 Orni Anc4"
export CSPECIES="Meze"
GOMERGEOUTPUT=${GOOUTDIR}/merged_MOTIFS_OUTPUT_
eval ../code/mergegoattrib/mergeEnrichment conf2.txt ${GOMERGEOUTPUT}
TWISE=${GOMERGEOUTPUT}_termwise.txt
CWISE=${GOMERGEOUTPUT}_clusterwise.txt
TPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_termwise.svg
CPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_clusterwise.svg
CCPLOT=${outdir}/merged_${LBL}_MOTIFS_PLOT_conserved_clusterwise.svg
COLORMAP="-1:(200,200,200);0:(0,175,0);2:(255,255,0);4:(175,0,0)"
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${CWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 ${TWISE} | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${TPLOT}
${TTAWK} -v IF=m -v IY="${Species_List}" -v OFS="\t" -v HasYH=0 <(grep -v "\-1" ${CWISE}) | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | awk -F"\t" -v OFS="\t" '$1 == "Dummy" && $2 == ENVIRON["CSPECIES"] {print "|Spacer||" NR, "|-"}; $1 == "Dummy"{next}; 1' | ${HMAWK} -v C="$COLORMAP" -v D="space" -v LegendNum="6" -v L="Cluster" > ${CCPLOT}

# run the above (run with 'bash' so that it calls bash directly and does not invoke posix which does not allow '<(XX)' commands - throws a syntax error )
bash MakeEnrichmentPlots_ST_local.sh

# ============================================
# Step 7: Convert all SVG files to pdf
# ============================================

cd ../results/Motif_enrichments/Figs_out/

# for the multi-tissue module files
nano svgtopdf.sh

#!/bin/bash
mkdir "$PWD"/pdf
for file in $PWD/*.svg
do
  filename=$(basename "$file")
  /Applications/Inkscape.app/Contents/Resources/bin/inkscape "$file" -d 1200 -A "$PWD"/pdf/"${filename%.svg}.pdf"
done

# run the above
bash svgtopdf.sh

# for the single-tissue module files
nano svgtopdf_ST.sh

#!/bin/bash
for i in br ey ht kd ms ts; do
  mkdir "$PWD"/"$i"/pdf
done
for i in br ey ht kd ms ts; do
  for file in $PWD/$i/*.svg; do
    filename=$(basename "$file")
    /Applications/Inkscape.app/Contents/Resources/bin/inkscape ${file} -d 1200 -A "$PWD"/"$i"/pdf/"${filename%.svg}.pdf"
  done
done

# run the above
bash svgtopdf_ST.sh
