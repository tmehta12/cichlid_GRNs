#Feb2017 - Cichlid GRNs functional variant analysis on cis-reg regions (Promoters and CNEs)
# Wilfried, Tomasz and Will

mkdir /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis

# Provide them with the 5-way alignment (as published from Brawand paper) - Luca is still running the new ones with improved assemblies
# found here:
/tgac/workarea/group-vh/cichlids/data/Broad/all/Data/5_way_MAFs

# use
cichild_5way.maf

# or if you require by scaffold
split_by_scaffold.tgz

# They will be initially looking at substitutions in the promoter region of EGR3 targets - only in the set with 1:1 orthologs in each of the five species
# Then looking at the list of all targets of all TFs
	# look at chromatin domains
	# use phastcons score based on teloest alignments
	# run FunSeq on tilapia genome to characterise motif break and gain

# Use old 5-way alignment to extract every subsitution vs Tilapia - put this into a vcf file
/tgac/workarea/group-vh/cichlids/data/Broad/all/Data/5_way_MAFs/cichild_5way.maf
# Provide a lsit of all target genes plus do the same analysis on a background set which is either the 5' or 3' gene adjacent to the TF target gene
	# Then plot distribution of foreground vs background


### Files that you need to generate to pass to Tomasz to run FunSeq
# Other files are present - see email sent on 13/02/17

###### 1. BED files of promoter regions in each genome - Will prepared these, found:
#The raw 5kb promoter region GTF files are here:
/tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/5kb_promoter/

#Sequences here:
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/5kb_promoter_seqs/

#The VCF files for the substitutions for each pairwise alignment are here:
/tgac/workarea/group-eg/project_cichlids/substitutions_VCFs/


###### 2. Gene prioritisation list (one per line, genes of interest)

# use cichlid gene IDs here
# prepare a master file with cichlid gene IDs, common gene name, Dr ensembl ID and 'class' of accelerated_morphogenesis, opsins, TFs and duplicates
# create one for each species, trying to retain any ortholog ordering

## NOTE: that immediately below is for the old modules (v1.1), below that is for the new modules (v3) - these are the files to be used

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_networks_v1.1/forFunSeq/
mkdir gene_prioritisation_lists

# accelerated_morphogenesis
cut -f1-6,10 Candidates_IDs.txt | awk -v OFS="\t" 'NR==1{print $0, "accelerated_morphogenesis";next}{print $0,"accelerated_morphogenesis"}' > Candidates_IDs.txt2
cut -f2 Candidates_IDs.txt2 | grep -v NA > Mz-accelerated_morphogenesis.txt
cut -f3 Candidates_IDs.txt2 | grep -v NA > Pn-accelerated_morphogenesis.txt
cut -f4 Candidates_IDs.txt2 | grep -v NA > Ab-accelerated_morphogenesis.txt
cut -f5 Candidates_IDs.txt2 | grep -v NA > Nb-accelerated_morphogenesis.txt
cut -f6 Candidates_IDs.txt2 | grep -v NA > On-accelerated_morphogenesis.txt

# opsins
cut -f2-7,11 opsin_gene_ids.txt | awk -v OFS="\t" 'NR==1{print $0, "opsins";next}{print $0,"opsins"}' > opsin_gene_ids.txt2
cut -f2 opsin_gene_ids.txt2 | grep -v NA > Mz-opsin.txt
cut -f3 opsin_gene_ids.txt2 | grep -v NA > Pn-opsin.txt
cut -f4 opsin_gene_ids.txt2 | grep -v NA > Ab-opsin.txt
cut -f5 opsin_gene_ids.txt2 | grep -v NA > Nb-opsin.txt
cut -f6 opsin_gene_ids.txt2 | grep -v NA > On-opsin.txt

# TFs
cut -f1-7 TFBS_519_full_orthologs_map7.txt | awk -v OFS="\t" 'NR==1{print $0, "TF";next}{print $0,"TF"}' > TFBS_519_full_orthologs_map7.txt2
cut -f3 TFBS_519_full_orthologs_map7.txt2 | grep -v NA > Mz-TF.txt
cut -f4 TFBS_519_full_orthologs_map7.txt2 | grep -v NA > Pn-TF.txt
cut -f5 TFBS_519_full_orthologs_map7.txt2 | grep -v NA > Ab-TF.txt
cut -f6 TFBS_519_full_orthologs_map7.txt2 | grep -v NA > Nb-TF.txt
cut -f7 TFBS_519_full_orthologs_map7.txt2 | grep -v NA > On-TF.txt

#For dup pairs put dup gene after anc gene and mark as either
	# Dup_pairs_common
	# Dup_pairs_novel

awk '{print $1 "\n"  $2}' Dup_pairs_common.txt | awk -v OFS="\t" 'NR==1{print $0, "Dup_pairs_common";next}{print $0,"Dup_pairs_common"}' > Dup_pairs_common.txt2
grep mz.gene Dup_pairs_common.txt2 > Mz-Dup_pairs_common.txt
grep pn.gene Dup_pairs_common.txt2 > Pn-Dup_pairs_common.txt
grep ab.gene Dup_pairs_common.txt2 > Ab-Dup_pairs_common.txt
grep nb.gene Dup_pairs_common.txt2 > Nb-Dup_pairs_common.txt
grep on.gene Dup_pairs_common.txt2 > On-Dup_pairs_common.txt
awk '{print $1 "\n"  $2}' Dup_pairs_novel_all.txt | awk -v OFS="\t" 'NR==1{print $0, "Dup_pairs_novel";next}{print $0,"Dup_pairs_novel"}'> Dup_pairs_novel_all.txt2
grep mz.gene Dup_pairs_novel_all.txt2 > Mz-Dup_pairs_novel.txt
grep pn.gene Dup_pairs_novel_all.txt2 > Pn-Dup_pairs_novel.txt
grep ab.gene Dup_pairs_novel_all.txt2 > Ab-Dup_pairs_novel.txt
grep nb.gene Dup_pairs_novel_all.txt2 > Nb-Dup_pairs_novel.txt
grep on.gene Dup_pairs_novel_all.txt2 > On-Dup_pairs_novel.txt

# {DONE} - can be shared with Tomasz, folder is (DO NOT SHARE THESE - SHARE BELOW)
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_networks_v1.1/forFunSeq/gene_prioritisation_lists
# share the following, genes in each file can be prioritised in this order for each species
*-accelerated_morphogenesis.txt
*-opsin.txt
*-TF.txt
*-Dup_pairs_common.txt # Note: Nb has no entries
*-Dup_pairs_novel.txt

## This is for the new modules (v3)
mkdir /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/forFunSeq
mkdir /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/forFunSeq/gene_prioritisation_lists
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/forFunSeq/gene_prioritisation_lists

# In each file use the gene with the most representation e.g. ab.gene

# accelerated_morphogenesis
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_networks_v1.1/forFunSeq/gene_prioritisation_lists/Candidates_IDs.txt .
awk -v OFS="\t" 'NR==1{print $0, "accelerated_morphogenesis";next}{print $0,"accelerated_morphogenesis"}' Candidates_IDs.txt > Candidates_IDs.txt2
awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" Candidates_IDs.txt2 > Candidates_IDs-ab.txt2
awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS="\t" OGIDS.txt5 > OGIDS-ab.txt5
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' OGIDS-ab.txt5 Candidates_IDs-ab.txt2 | awk '{print $2,$13,$14,$15,$12,$16,$17,$18,$19,$20,$22,$24,$25,$26,$11}' OFS='\t' > Candidates_IDs.txt3
cut -f3 Candidates_IDs.txt3 | grep -v NULL > Mz-accelerated_morphogenesis.txt
cut -f4 Candidates_IDs.txt3 | grep -v NULL > Pn-accelerated_morphogenesis.txt
cut -f5 Candidates_IDs.txt3 | grep -v NULL > Ab-accelerated_morphogenesis.txt
cut -f6 Candidates_IDs.txt3 | grep -v NULL > Nb-accelerated_morphogenesis.txt
cut -f7 Candidates_IDs.txt3 | grep -v NULL > On-accelerated_morphogenesis.txt

# opsins
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_networks_v1.1/forFunSeq/gene_prioritisation_lists/opsin_gene_ids.txt .
cut -f2-7,11 opsin_gene_ids.txt | awk -v OFS="\t" 'NR==1{print $0, "opsins";next}{print $0,"opsins"}' > opsin_gene_ids.txt2
awk '{print $4,$1,$2,$3,$5,$6,$7,$8}' OFS="\t" opsin_gene_ids.txt2 > opsin_gene_ids-ab.txt2
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' OGIDS-ab.txt5 opsin_gene_ids-ab.txt2 | awk '{print $2,$10,$11,$12,$9,$13,$14,$15,$16,$17,$19,$21,$22,$23,$8}' OFS='\t' > opsin_gene_ids.txt3
cut -f3 opsin_gene_ids.txt3 | grep -v NULL > Mz-opsin.txt
cut -f4 opsin_gene_ids.txt3 | grep -v NULL > Pn-opsin.txt
cut -f5 opsin_gene_ids.txt3 | grep -v NULL > Ab-opsin.txt
cut -f6 opsin_gene_ids.txt3 | grep -v NULL > Nb-opsin.txt
cut -f7 opsin_gene_ids.txt3 | grep -v NULL > On-opsin.txt

# TFs
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_networks_v1.1/forFunSeq/gene_prioritisation_lists/TFBS_OGIDs3.txt .
cut -f1-8 TFBS_OGIDs3.txt | awk -v OFS="\t" 'NR==1{print $0, "TF";next}{print $0,"TF"}' > TFBS_OGIDs3.txt2
cut -f4 TFBS_OGIDs3.txt2 | grep -v NULL > Mz-TF.txt
cut -f5 TFBS_OGIDs3.txt2 | grep -v NULL > Pn-TF.txt
cut -f6 TFBS_OGIDs3.txt2 | grep -v NULL > Ab-TF.txt
cut -f7 TFBS_OGIDs3.txt2 | grep -v NULL > Nb-TF.txt
cut -f8 TFBS_OGIDs3.txt2 | grep -v NULL > On-TF.txt

#For dup pairs put dup gene after anc gene and mark as either
	# Dup_pairs_common
	# Dup_pairs_novel
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_networks_v1.1/forFunSeq/gene_prioritisation_lists/Dup_pairs_common.txt .
awk '{print $1 "\n"  $2}' Dup_pairs_common.txt | awk -v OFS="\t" 'NR==1{print $0, "Dup_pairs_common";next}{print $0,"Dup_pairs_common"}' > Dup_pairs_common.txt2
grep mz.gene Dup_pairs_common.txt2 > Mz-Dup_pairs_common.txt
grep pn.gene Dup_pairs_common.txt2 > Pn-Dup_pairs_common.txt
grep ab.gene Dup_pairs_common.txt2 > Ab-Dup_pairs_common.txt
grep nb.gene Dup_pairs_common.txt2 > Nb-Dup_pairs_common.txt
grep on.gene Dup_pairs_common.txt2 > On-Dup_pairs_common.txt

cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_networks_v1.1/forFunSeq/gene_prioritisation_lists/Dup_pairs_novel_all.txt .
awk '{print $1 "\n"  $2}' Dup_pairs_novel_all.txt | awk -v OFS="\t" 'NR==1{print $0, "Dup_pairs_novel";next}{print $0,"Dup_pairs_novel"}'> Dup_pairs_novel_all.txt2
grep mz.gene Dup_pairs_novel_all.txt2 > Mz-Dup_pairs_novel.txt
grep pn.gene Dup_pairs_novel_all.txt2 > Pn-Dup_pairs_novel.txt
grep ab.gene Dup_pairs_novel_all.txt2 > Ab-Dup_pairs_novel.txt
grep nb.gene Dup_pairs_novel_all.txt2 > Nb-Dup_pairs_novel.txt
grep on.gene Dup_pairs_novel_all.txt2 > On-Dup_pairs_novel.txt

# {DONE} - can be shared with Tomasz, folder is
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/forFunSeq/gene_prioritisation_lists
# share the following, genes in each file can be prioritised in this order for each species
*-accelerated_morphogenesis.txt
*-opsin.txt
*-TF.txt
*-Dup_pairs_common.txt # Note: Nb has no entries
*-Dup_pairs_novel.txt

###### 3. Network information (two columns, GeneA GeneB)
# create per species and per interaction type
# use the old Arboretum v1.1 network edges

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_networks_v1.1
mkdir forFunSeq
cd forFunSeq

# change the source in column 19 from CNE_JasparTFBS to CNE_Proximal if there is a value (other than NA) in column 47 (S.N[cp]) - this is becasue I did not change the source column when generating the new files
# Use the unique S/N column (#47) to identify the rows of CNE_proximal and then change the source column (#19) for those rows
for i in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt3
do awk 'BEGIN {
     IFS = OFS = "\t"
  }
  {
     for (column = 47; column <= NF; ++column) {
        if ($column !="NA") {
            $19 = "CNE_Proximal"
        }
     }
     {print $0}
  }
' $i > "$(basename "$i" .txt3).txt4" ; done

# The above, for some reason changes all string interactions to CNE_Proximal too
# Use the unique combined_score column (#67) to identify the rows that do not contain 'NA' and then change the source column (#19) for those rows back to string
# Run this on above files which will be in this forFunSeq folder but output new files in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/ as SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4
for i in SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4
do awk 'BEGIN {
     IFS = OFS = "\t"
  }
  {
     for (column = 67; column <= NF; ++column) {
        if ($column !="NA") {
            $19 = "string"
        }
     }
     {print $0}
  }
' $i > ../Edge_Attributes/0.FEB2017-FINALEdgeTables/${i} ; done


# remove the intermediate files from this folder
rm SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4

# Cut the geneA and geneB columns, per species and per interaction

# paccmit[pm_]
# genemania[gm_]
# signafish[sf_]
# string[st_]
# mrnetb[mr_]
# promoter_JasparTFBS[pj_]
# CNE_JasparTFBS[cj_]
# CNE_Proximal[cp]

for i in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4 ; do grep -wiF paccmit $i | cut -f1,5 > "$(basename "$i" .txt4)-funseq_paccmit.txt" ; done
for i in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4 ; do grep -wiF genemania $i | cut -f1,5 > "$(basename "$i" .txt4)-funseq_genemania.txt" ; done
for i in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4 ; do grep -wiF signafish $i | cut -f1,5 > "$(basename "$i" .txt4)-funseq_signafish.txt" ; done
for i in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4 ; do grep -wiF string $i | cut -f1,5 > "$(basename "$i" .txt4)-funseq_string.txt" ; done
for i in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4 ; do grep -wiF mrnetb $i | cut -f1,5 > "$(basename "$i" .txt4)-funseq_mrnetb.txt" ; done
for i in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4 ; do grep -wiF promoter_JasparTFBS $i | cut -f1,5 > "$(basename "$i" .txt4)-funseq_promoter_JasparTFBS.txt" ; done
for i in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4 ; do grep -wiF CNE_JasparTFBS $i | cut -f1,5 > "$(basename "$i" .txt4)-funseq_CNE_JasparTFBS.txt" ; done
for i in ../Edge_Attributes/0.FEB2017-FINALEdgeTables/SpeciesFinal_*-Edge_Attributes_Collated3NEW3.txt4 ; do grep -wiF CNE_Proximal $i | cut -f1,5 > "$(basename "$i" .txt4)-funseq_CNE_Proximal.txt" ; done

# remove geneA and geneB from each string file
for i in SpeciesFinal_*-Edge_Attributes_Collated3NEW3-funseq_string.txt ; do echo "$(tail -n +2 $i)" > "$(basename "$i" .txt)2.txt" ; done
# delete old string files
rm SpeciesFinal_*-Edge_Attributes_Collated3NEW3-funseq_string.txt
# rename new string files to old name
rename string2 string *string2.txt

# move all to specific subdir
mkdir network_edges
mv *.txt network_edges
# {DONE} - can be shared with Tomasz, folder is
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_networks_v1.1/forFunSeq/network_edges


###### 4. Transcription factor (TF) binding motifs under TF peaks (BED file) - we do not have TF peaks so just have to go with best hit, basically coordinates of motifs located along the genome

##### NOTE: THE BELOW IS FOR 1E-6, BEST TO USE 1E-4 (AT END OF THIS SCRIPT)


# run on uv2
uv2
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/
mkdir motif_discovery
cd motif_discovery
mkdir whole_genome
cd whole_genome

# use FIMO in the MEME suite to scan for motif models using the JASPAR file
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/7.CNEs/JASPAR_CORE_2016_vertebrates.meme
# across the genomes at a level of statistical significance (use stringent pval of 1e-6, default is <1e-4)

#The main trick to scanning whole genomes with FIMO is to set a stringent p-value threshold. In most cases you really don't want 10 million matches almost all of which are statistically insignificant after correcting for multiple testing. When we were generating the custom tracks mentioned above we used the following settings:
#--thresh 1e-6
#--max-stored-scores 500000

# Then filter FIMO output afterwards so that we only included matches with a q-value less than 0.01.

# Take the search result from FIMO and convert it to a UCSC BED-formatted file. These files now contain the putative binding sites of all motifs from the motif models across the whole genome.

## Unzip all the fasta files first
gunzip /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta.gz
gunzip /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa.gz
#gunzip /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v2.assembly.fasta.gz
# use /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta
gunzip /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta.gz
gunzip /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta.gz
gunzip /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta.gz

# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-OnIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# O. niloticus – this is the old Illumina assembly [USE THIS]
fimo --o oniloticus_illumina --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-OnPacBio.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)
# O. niloticus - this is the improved Illumina+PacBio assembly [DO NOT USE]
fimo --o oniloticus_illumina-pacbio --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-MzPacBio.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# M. zebra – this is the improved Illumina+PacBio assembly [DO NOT USE]
fimo --o mzebra_illumina-pacbio --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-MzIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# M. zebra – this is the old Illumina assembly [USE THIS]
fimo --o mzebra_illumina --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-AbIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# A. burtoni
fimo --o aburtoni_illumina --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-PnIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# P. nyererei
fimo --o pnyererei_illumina --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-NbIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# N. brichardi
fimo --o nbrichardi_illumina --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta


# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=3 wholegenome_motifdiscovery-OnIllumina.sh #{DONE}
qsub -q Prod -l select=1:mem=100GB:ncpus=3 wholegenome_motifdiscovery-MzIllumina.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=3 wholegenome_motifdiscovery-AbIllumina.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=3 wholegenome_motifdiscovery-PnIllumina.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=3 wholegenome_motifdiscovery-NbIllumina.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=3 wholegenome_motifdiscovery-OnPacBio.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=3 wholegenome_motifdiscovery-MzPacBio.sh #{DONE}

# create a list of the outputted folders
nano list1
#add the following to above
oniloticus_illumina
oniloticus_illumina-pacbio
mzebra_illumina-pacbio
mzebra_illumina
aburtoni_illumina
pnyererei_illumina
nbrichardi_illumina

# Then filter FIMO output afterwards so that we only include matches with a q-value less than 0.05 (this is your FDR-based multiple testing correction). Note: FIMO output is 1-based to convert to 0-based by taking -1 from start e.g. 5-20 turns into 4-20
# Then convert this filtered file to a UCSC BED-formatted file (0-based). These files now contain the putative binding sites of all motifs from the motif models across the whole genome.
# BED format: chrom, chromStart, chromEnd, name, score, strand
while read F ; do awk '(NR==1) || ($8 < 0.05 )' $F/fimo.txt | echo "$(tail -n +2)" | awk '{print $2, $3-1, $4, $1, $6, $5}' OFS='\t' | sort -k1,1 -k2,2n | awk -v OFS="\t" 'NR==1{print $0, ".";next}{print $0,"."}' | awk '{print $1, $2, $3, $4, $7, $6, $4}' OFS='\t' > $F/$F-fimo.bed ; done < list1 # this filters, gets rid of header, rearranges column, converts start coordinate to 0-based, adds a row of '.' and then sorts based on chrom and then chromStart column

## Did the oniloticus where we filtered with q-val of 0.1 to maintain the biological diversity
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/whole_genome/oniloticus_illumina
awk '(NR==1) || ($8 < 0.1 )' fimo.txt | echo "$(tail -n +2)" | awk '{print $2, $3-1, $4, $1, $6, $5}' OFS='\t' | sort -k1,1 -k2,2n | awk -v OFS="\t" 'NR==1{print $0, ".";next}{print $0,"."}' | awk '{print $1, $2, $3, $4, $7, $6, $4}' OFS='\t' > oniloticus_illumina_q0.1-fimo.bed
cp oniloticus_illumina_q0.1-fimo.bed ../all_BED_files/ ## share this with Tomasz

#copy all bed files to one folder to share with Tomasz
mkdir all_BED_files
while read F ; do cp $F/$F-fimo.bed all_BED_files ; done < list1


###### NOTE: THESE HAVE BEEN AMENDED BELOW TO ADD 50 TO START AND END OF ALL COORDINATES - USE THE FILES BELOW TO SHARE WITH TOMASZ

# {DONE} - can be shared with Tomasz, folder is
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/whole_genome/all_BED_files #only need to use those that are *_illumina-fimo.bed (do not use any of the *_illumina-pacbio-fimo.bed)


###### 5. TF-target (Regulon files) - here we need regulons (TF-target gene) interactions of all genes in the genome of the five species
# Will created a gtf for the promoter regions - he will use this gtf file to extract all 5kb sequences upstream of the TSS of each gene
# GTF files which have apparent isoforms removed, and only the longest gene versions remaining, as well as the GTFs for the 5kb promoters are here: /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/5kb_promoter
# promoter sequences found here:/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/5kb_promoter_seqs/
# Also, GTF files describing genes that partially or completely overlap on different strands, or genes that are intronic, are here: /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/protein_coding/cichlid_overlap/
# Then, scan for motifs in these promoter sequences (remember, previously we used motif enichment results because not all genes in the genome were present in modules, not required here BUT if you use these for the module genes, you will need to run an enrichment)
# After that, prepare a regulon file for each species, all interactions

# run on uv2
uv2
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery
mkdir promoters
cd promoters

# convert the TF names in the JASPAR_CORE_2016_vertebrates.meme file all to uppercase to match the mapping file - for this you need to grab the second word after MOTIF in each line and convert to uppercase
# RUNNING LOCAL - ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/
awk '{for(i=1;i<=NF;i++) if ($i=="MOTIF") print $(i+2)=toupper;print}' OFS='\t' JASPAR_CORE_2016_vertebrates.meme | sed '/MOTIF/{N;s/\n.*//;}' > JASPAR_CORE_2016_vertebrates.meme2 # the long way - awk adds the line above the original so delete the original changed line retaining the uppercase with sed
mv JASPAR_CORE_2016_vertebrates.meme2 JASPAR_CORE_2016_vertebrates.meme # rename file
cp JASPAR_CORE_2016_vertebrates.meme ../6.TFBSs # copy to other folders
# also copy to workarea folders:
#/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/7.CNEs/
#/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/6.TFBSs/
#/usr/users/ga004/mehtat/Cichlid_GRNs/Arboretum_GT_v2/7.CNEs/

# use FIMO in the MEME suite to scan for motif models using the JASPAR file
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/7.CNEs/JASPAR_CORE_2016_vertebrates.meme
# across the genomes at a level of statistical significance (use stringent default pval <1e-4)

#The main trick to scanning with FIMO is to set a stringent p-value threshold. In most cases you really don't want 10 million matches almost all of which are statistically insignificant after correcting for multiple testing. When scanning the whole genome I used the following settings:
#--thresh 1e-4
#--max-stored-scores 500000
# Then filter FIMO output afterwards so that we only included matches with a q-value less than 0.1.

# create symbolic links to the promoter sequences:
for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/new_v1_5kb_promoters/* ; do ln -s $i ; done

# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-Mz.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# M. zebra
fimo --o mzebra_promoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Metriaclima_zebra.BROADMZ1.5kb_promoters.0817.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-Pn.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# P. nyererei
fimo --o pnyererei_promoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Pundamilia_nyererei.BROADPN1.5kb_promoters.0817.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-Ab.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# A. burtoni
fimo --o aburtoni_promoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Astatotilapia_burtoni.BROADAB1.5kb_promoters.0817.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-Nb.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# N. brichardi
fimo --o nbrichardi_promoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Neolamprologus_brichardi.BROADNB1.5kb_promoters.0817.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_motifdiscovery-On.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# O. niloticus
fimo --o oniloticus_promoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Oreochromis_niloticus.BROADON1.5kb_promoters.0817.fasta


# Run all on uv2 - RERAN ALL 26/07/17 WITH NEW PROMOTER SEQUENCES!
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-Mz.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-Pn.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-Ab.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-Nb.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_motifdiscovery-On.sh #{DONE}

### Running some on slurm as UV2 are overloaded!!!
nano promoter_motifdiscovery-Ab_slurm.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# A. burtoni
fimo --o aburtoni_promoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Astatotilapia_burtoni.BROADAB1.5kb_promoters.0817.fasta

nano sbatch-promoter_motifdiscovery-Ab_slurm.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 32000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh promoter_motifdiscovery-Ab_slurm.sh

# Run as batch script
sbatch sbatch-promoter_motifdiscovery-Ab_slurm.sh


nano promoter_motifdiscovery-Nb_slurm.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# A. burtoni
fimo --o nbrichardi_promoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Neolamprologus_brichardi.BROADNB1.5kb_promoters.0817.fasta

nano sbatch-promoter_motifdiscovery-Nb_slurm.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 32000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh promoter_motifdiscovery-Nb_slurm.sh

# Run as batch script
sbatch sbatch-promoter_motifdiscovery-Nb_slurm.sh


nano promoter_motifdiscovery-On_slurm.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# A. burtoni
fimo --o oniloticus_promoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Oreochromis_niloticus.BROADON1.5kb_promoters.0817.fasta

nano sbatch-promoter_motifdiscovery-On_slurm.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 32000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh promoter_motifdiscovery-On_slurm.sh

# Run as batch script
sbatch sbatch-promoter_motifdiscovery-On_slurm.sh


# create a list file for the filtering
nano list1
mzebra_promoters
pnyererei_promoters
aburtoni_promoters
nbrichardi_promoters
oniloticus_promoters


# Then filter FIMO output afterwards so that we only included matches with a q-value less than 0.05
while read F ; do awk '(NR==1) || ($8 < 0.1 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-fimo.txt ; done < list1 # this filters and gets rid of headers - you have to accept here that your usage of a relatively stringent q-value will filter out statistically insignificant hits, under this are CRX and OTX5 that have q-val>0.7. q-val of 0.1 implies that 10% of significant tests will result in false positives. I am using 0.1 as we have used vertebrate motifs to scan fish promoters, if we were scanning human/mouse then 0.01 would be ideal.

### NOTE: immediately below is for old (v1.1) modules, below that is for new modules (v3)

# Now, convert the outputs into per-species regulon tables (TF-target gene)
# create a mapping file with pattern name, cichlid ID and gene_symbol - create a new file where the motif mapping is better #{DONE}
# running this local here:
cd ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/6.TFBSs/
R CMD BATCH 0.fullorthologmapping_new.R # map orthologs to human ensembl ID and gene symbol - use this gene symbol for mapping TFs > full_orthologs_map3.txt

awk '{for(i=1;i<=NF;i++) if ($i=="MOTIF") print $(i+2),$(i+1)}' OFS='\t' JASPAR_CORE_2016_vertebrates.meme | awk '$1 = toupper($1)' OFS='\t' > TFBS_519_pattern-names.txt # extract all TFBS motif names, loop over 'MOTIF' word of each line and print second (gene symbol) and first (pattern) word after that, then convert TF to uppercase to match later
sed 's/(VAR.2)//g' TFBS_519_pattern-names.txt | sed 's/(VAR.3)//g' | awk '{split($1,a,/::/);$1=a[1]}1' OFS="\t" > TFBS_519_pattern-names.txt2 # remove var.2 and var.3 for now to match and then reintroduce later # in cases of co-bound TFs e.g. AHR::ARNT, just take the first TF, here that is AHR - reintroduce co-bound nomenclature later

awk -F"\t" '{print $14, $13, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' full_orthologs_map3.txt > full_orthologs_map3a.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' full_orthologs_map3a.txt TFBS_519_pattern-names.txt2 > TFBS_519_full_orthologs_map3-Hs.txt #403/519 match


# do this same type of map to the two other columns, 'dr_ensembl_gene_id' and 'ga_ensembl_gene_id' in the full_orthologs_map3.txt file
## when you paste into original script, retain the scripts to prepare TFBS_519_full_orthologs_map2.txt here
# just re-ruse the file TFBS_519_full_orthologs_map3.txt (as this will include paralogs too) and then take the intersection of all three outputs as the TFBS_519_ortholog_map
# amend both TFBS_519_full_orthologs_map3-Hs.txt and TFBS_519_full_orthologs_map3.txt files so that they follow the same structure
awk '$3!="NA"' OFS="\t" TFBS_519_full_orthologs_map3-Hs.txt | cut -f5-13 > TFBS_519_full_orthologs_map3-Hsa.txt
cut -f3-11 TFBS_519_full_orthologs_map3.txt > TFBS_519_full_orthologs_map3a.txt
cat TFBS_519_full_orthologs_map3-Hsa.txt TFBS_519_full_orthologs_map3a.txt > TFBS_519_full_orthologs_map4.txt # join the files
sort -u TFBS_519_full_orthologs_map4.txt > TFBS_519_full_orthologs_map5.txt
# awk map back to TF names
# for the Hs files use ENSGAC ids ($12) first, then ab.gene ($7), then on.gene ($9), hopefully get all - plus remove all rows where NA present in that column
awk '$3!="NA"' OFS='\t' TFBS_519_full_orthologs_map3-Hs.txt | awk -F"\t" '{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $13, $14, $15, $16}' OFS='\t' | awk '$1!="NA"' OFS='\t' > TFBS_519_full_orthologs_map3-Hsb.txt
awk '$3!="NA"' OFS='\t' TFBS_519_full_orthologs_map3-Hs.txt | awk -F"\t" '{print $7, $1, $2, $3, $4, $5, $6, $8, $9, $10, $11, $12, $13, $14, $15, $16}' OFS='\t' | awk '$1!="NA"' OFS='\t' > TFBS_519_full_orthologs_map3-Hsc.txt
awk '$3!="NA"' OFS='\t' TFBS_519_full_orthologs_map3-Hs.txt | awk -F"\t" '{print $9, $1, $2, $3, $4, $5, $6, $7, $8, $10, $11, $12, $13, $14, $15, $16}' OFS='\t' | awk '$1!="NA"' OFS='\t' > TFBS_519_full_orthologs_map3-Hsd.txt
# for the merge files use ENSGAC ids ($8) first as 501/538 are present, then ab.gene ($3), then on.gene ($5), hopefully get all
awk -F"\t" '{print $8, $1, $2, $3, $4, $5, $6, $7, $9}' OFS='\t' TFBS_519_full_orthologs_map5.txt | awk '$1!="NA"' OFS='\t' > TFBS_519_full_orthologs_map5a.txt
awk -F"\t" '{print $3, $1, $2, $4, $5, $6, $7, $8, $9}' OFS='\t' TFBS_519_full_orthologs_map5.txt | awk '$1!="NA"' OFS='\t' > TFBS_519_full_orthologs_map5b.txt
awk -F"\t" '{print $5, $1, $2, $3, $4, $6, $7, $8, $9}' OFS='\t' TFBS_519_full_orthologs_map5.txt | awk '$1!="NA"' OFS='\t' > TFBS_519_full_orthologs_map5c.txt
# awk match
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_519_full_orthologs_map3-Hsb.txt TFBS_519_full_orthologs_map5a.txt > TFBS_519_full_orthologs_map3Hsb-5a.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_519_full_orthologs_map3-Hsc.txt TFBS_519_full_orthologs_map5b.txt > TFBS_519_full_orthologs_map3Hsc-5b.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_519_full_orthologs_map3-Hsd.txt TFBS_519_full_orthologs_map5c.txt > TFBS_519_full_orthologs_map3Hsd-5c.txt

# rearrange columns of all awk outputs to a unified column ordering (motif pattern, TF, mz, pn, ab, nb, on, Ga_ensembl, Ga_genesymbol, Ol_ensembl, Tn_ensembl, Dr_ensembl, Dr_genesymbol, On_ensembl)
awk -F"\t" '{print $12, $13, $2, $3, $4, $5, $6, $1, $23, $7, $8, $9, $24, $25}' OFS='\t' TFBS_519_full_orthologs_map3Hsb-5a.txt > TFBS_519_full_orthologs_map3Hsb-5a.txt2
awk -F"\t" '{print $12, $13, $2, $3, $1, $4, $5, $8, $23, $6, $7, $9, $24, $25}' OFS='\t' TFBS_519_full_orthologs_map3Hsc-5b.txt > TFBS_519_full_orthologs_map3Hsc-5b.txt2
awk -F"\t" '{print $12, $13, $2, $3, $4, $5, $1, $8, $23, $6, $7, $9, $24, $25}' OFS='\t' TFBS_519_full_orthologs_map3Hsd-5c.txt > TFBS_519_full_orthologs_map3Hsd-5c.txt2
cat TFBS_519_full_orthologs_map3Hsb-5a.txt2 TFBS_519_full_orthologs_map3Hsc-5b.txt2 TFBS_519_full_orthologs_map3Hsd-5c.txt2 | sort -u > TFBS_519_full_orthologs_map6.txt
# created the above file and now manually editing in excel file
TFBS_519_full_orthologs_map6.xlsx > TFBS_519_full_orthologs_map7.txt

# for some manual entries, created a file to grep > greporthologs and then to grep on the full_orthologs file, maintaining grep order, ran:
for name in `cat greporthologs ` ; do grep -wiF $name full_orthologs_map3a.txt ; done > greporthologs2.txt

# then ran cut to extract the columns you need
cut -f1,4,5,12-14 greporthologs2.txt > greporthologs3.txt

# for manual entries of motif pattern and name, created a file to grep > greppatterns and then to grep on the TFBS_519_pattern-names.txt2 file, maintaining grep order, ran:
for name in `cat greppatterns `; do grep -wiF $name TFBS_519_pattern-names.txt2 | awk '{print $2, $1}' OFS='\t'; done > greppatterns2.txt
# To bring back to final file - manually edited in excel file
# var.2 and var.3 - added new rows for each variant (same TF)
# co-bound TFs - added extra rows for second co-bound TF (hopefully no problems with same motif name twice)
# final file > TFBS_519_full_orthologs_map7.txt

# Created the following mapping file for Tomasz
cut -f1,2,7 TFBS_519_full_orthologs_map7.txt | sort -k1,1 > TFBS_519_full_orthologs_map7-Onorthologs.txt

# Use the mapping file to create the per-species regulon tables (TF-target gene)
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/6.TFBSs/TFBS_519_full_orthologs_map7.txt

#pattern name   sequence name   start   stop    strand  score   p-value q-value matched sequence
MA0528.1        nb.gene.s17.72  4898    4918    +       29.8163 6.95e-14        9.35e-07        GGAGGAGGAGGGGGAGGAGGG
MA0138.2        nb.gene.s48.73  494     514     +       34.6122 1.25e-13        2.13e-05        TTCAGCACCATGGACAGCGCC

#### Instead, to match the Arboretum_GT_v3 module, created a new TFBS ortholog file:
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/6.TFBSs/TFBS_519_full_orthologs_map7.txt .
cut -f1,2,5,9,13 TFBS_519_full_orthologs_map7.txt > TFBS_519_full_orthologs_map7.txt2 # cut the columns that you will use to map and in the final table
cut -f1,2 TFBS_519_full_orthologs_map7.txt2 > col1-2 # to only add to final table
cut -f2 TFBS_519_full_orthologs_map7.txt2 > TFBShs
#cut -f3 TFBS_519_full_orthologs_map7.txt2 > TFBSab # do not use as it is bias to old orthology
cut -f4 TFBS_519_full_orthologs_map7.txt2 > TFBSga
cut -f5 TFBS_519_full_orthologs_map7.txt2 > TFBSdr
#awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' OGIDS.txt5 > OGIDSab # do not use as this match would be bias to old orthology
awk '{print $10,$1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$12,$13,$14,$15}' OFS='\t' OGIDS.txt5 > OGIDSga
awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$13,$14,$15}' OFS='\t' OGIDS.txt5 > OGIDSdr
awk '{print $15,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' OFS='\t' OGIDS.txt5 > OGIDShs
# match and then cut the orthogroup column
#awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' OGIDSab TFBSab | cut -f3 > OGab # bias to old orthology, do not use
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' OGIDSga TFBSga | cut -f3 > OGga
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' OGIDSdr TFBSdr | cut -f3 > OGdr
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' OGIDShs TFBShs | cut -f3 > OGhs

paste -d'\t' OGhs OGga OGdr | sed 's/NULL//g' | awk -v 'OFS=\t' '{print $1,$2,$3}' | awk '!NF{$0="NULL"}1' | cut -f1 > TFBS_OGIDs
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL";}}' OGIDS.txt5 TFBS_OGIDs | cut -f2-16 > TFBS_OGIDs2
paste -d'\t' col1-2 TFBS_OGIDs2 > TFBS_OGIDs3 # open in excel and manually edit some of the missing entries like hoxd3 e.g. grep hoxd3 OGIDS.txt5 > FINAL FILE IS TFBS_OGIDs3.txt
## Just copy file from /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/6.TFBSS/ to /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters {DONE}

# Created the following mapping file for Tomasz
cut -f1,2,8 TFBS_OGIDs3.txt | sort -k1,1 > TFBS_OGIDs3-Onorthologs.txt #{DONE}

# create a per-species mapping file
cut -f1,4 TFBS_OGIDs3.txt | grep mz.gene > TFBS_OGIDs3-Mz.txt #{DONE}
cut -f1,5 TFBS_OGIDs3.txt | grep pn.gene > TFBS_OGIDs3-Pn.txt #{DONE}
cut -f1,6 TFBS_OGIDs3.txt | grep ab.gene > TFBS_OGIDs3-Ab.txt #{DONE}
cut -f1,7 TFBS_OGIDs3.txt | grep nb.gene > TFBS_OGIDs3-Nb.txt #{DONE}
cut -f1,8 TFBS_OGIDs3.txt | grep on.gene > TFBS_OGIDs3-On.txt #{DONE}





## To run once motif scan done ## {DONE}
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters

# Then filter FIMO output afterwards so that we only included matches with a q-value less than 0.05
while read F ; do awk '(NR==1) || ($8 < 0.1 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-fimo.txt ; done < list1 # this filters and gets rid of headers - you have to accept here that your usage of a relatively stringent q-value will filter out statistically insignificant hits, under this are CRX and OTX5 that have q-val>0.7. q-val of 0.1 implies that 10% of significant tests will result in false positives. I am using 0.05 as we have used vertebrate motifs to scan fish promoters, if we were scanning human/mouse then 0.01 would be ideal.

# awk match creating the whole-genome regulon files
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-Mz.txt mzebra_promoters/mzebra_promoters-fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > mzebra_promoters/Mz-wholegenome_q0.1Regulon.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-Pn.txt pnyererei_promoters/pnyererei_promoters-fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > pnyererei_promoters/Pn-wholegenome_q0.1Regulon.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-Ab.txt aburtoni_promoters/aburtoni_promoters-fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > aburtoni_promoters/Ab-wholegenome_q0.1Regulon.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-Nb.txt nbrichardi_promoters/nbrichardi_promoters-fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > nbrichardi_promoters/Nb-wholegenome_q0.1Regulon.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-On.txt oniloticus_promoters/oniloticus_promoters-fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > oniloticus_promoters/On-wholegenome_q0.1Regulon.txt

# Final Regulon files to share with Tomasz are:
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/mzebra_promoters/Mz-wholegenome_q0.1Regulon.txt
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/pnyererei_promoters/Pn-wholegenome_q0.1Regulon.txt
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/aburtoni_promoters/Ab-wholegenome_q0.1Regulon.txt
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/nbrichardi_promoters/Nb-wholegenome_q0.1Regulon.txt
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/oniloticus_promoters/On-wholegenome_q0.1Regulon.txt





## Once done this:
# 1. copy full_orthologs_map3.txt to all relevant folders - use this from now on {DONE}
## Created a new fullorthologs file that has OGid assigned
cp /usr/users/ga004/mehtat/Cichlid_GRNs/Arboretum_networks_v1.1/march2017_regulonanalysis/full_orthologs_map3c.txt /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/

# 2. instead of amending main NetworkReconstruction_v3.sh script to reflect this new method (in 6.TFBSs and 7.CNEs - retaining R CMD BATCH 0.fullorthologmapping.R and add above below that, create a NEW file, NetworkReconstruction_v4.sh {DONE}
# 3. Amend usage of 'TFBS_519_full_orthologs_map2.txt' to now use new file in every position where it is used (TFBS_519_full_orthologs_map7.txt) - be careful though, the columns are different now, need to amend in each script {DONE}
	# TFBS_519_full_orthologs_map2.txt: mz, pn, ab, nb, on, Ol_ensembl, Tn_ensembl, Ga_ensembl, Dr_ensembl, Ga_symbol, Dr_symbol
	# TFBS_519_full_orthologs_map7.txt: pattern, TF, mz, pn, ab, nb, on, Ga_ensembl, Ga_symbol, Ol_ensembl, Tn_ensembl, Dr_ensembl, Dr_symbol, On_ensembl
	# mz$1>$3, pn$2>$4, ab$3>$5, nb$4>$6, on$5>$7, Ol_ensembl$6>$10, Tn_ensembl$7>$11, Ga_ensembl$8>$8, Dr_ensembl$9>$12, Ga_symbol$10>$9, Dr_symbol$11>$13, On_ensembl>$14
	# this also includes in the prioritisation lists that you shared with Tomasz, amend the processing to now process this new file {DONE}
	# Need to amend the TFBS file used in # Candidate genes >--{DONE - RAN ON CLUSTER + COPIED LOCAL}--< section {DONE}


###### 6. Sort out Position Frequency Matrix (PFM) file for Tomasz - needs to be fasta format

# Tomasz will prepare a new PFM file according to the ENCODE format that has IUPAC codes (http://www.bioinformatics.org/sms/iupac.html). The threshold is 0.16 for each.



###### 7. Conservation scores e.g. PhastCons
# generate this based on the new multiple alignments




#### 12/04/17 We need to amend all annotation files so that -50 and +50 to start and stop coordinates as funSeq scans both directions (29bp) on the start and end of a scaffold where a SNP is present. Thus, it spits out an error if it goes into a negative or non-existant value value

## Will is working on amending annotations for CDS BED, Intron BED, CDS Interval BED, Promoter BED and UTR BED - adding 50 to start and end coordinates plus amending any headers

## 5. Transcription factor (TF) binding motifs under TF peaks (BED file)
## Tarang needs to add 50 to col2 and col3 of each file in /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/whole_genome/all_BED_files/
for i in /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/whole_genome/all_BED_files/*fimo.bed ; do awk '{print $1,$2+50,$3+50,$4,$5,$6,$7}' OFS="\t" $i > all_BED_files/"$(basename "$i" .bed)-NEW50bp.bed" ; done

mkdir /tgac/workarea/group-vh/cichlids/funSeq/cichlids/illumina-assemblies/Tarang50bpamend/
mkdir /tgac/workarea/group-vh/cichlids/funSeq/cichlids/illumina-assemblies/Tarang50bpamend/5.TFbindingmotifs-WG # dump files here
while read F ; do cp all_BED_files/$F-fimo-NEW50bp.bed /tgac/workarea/group-vh/cichlids/funSeq/cichlids/illumina-assemblies/Tarang50bpamend/5.TFbindingmotifs-WG ; done < list1
## share folder with Tomasz /tgac/workarea/group-vh/cichlids/funSeq/cichlids/illumina-assemblies/Tarang50bpamend/5.TFbindingmotifs-WG/*-fimo-NEW50bp.bed
# do not use *-pacbio-fimo-NEW50bp.bed


## 14. Genome sequences - need to add 50Ns to start and end of each scaffold

uv

mkdir /tgac/workarea/group-vh/cichlids/funSeq/cichlids/illumina-assemblies/Tarang50bpamend/14.Genome_sequences

# for i in {1..50}; do printf "N" ; done # print 50 Ns and copy below

# awk '1;/>/{ printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"}' /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta # adds to the line below a '>' and not a new line

# awk '/>/{print "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"}1' /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta | awk '{if (NR!=1) {print}}' # adds to line above, then remove the very first line in file which will just have 50 Ns

# awk '1;/>/{ printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"}' /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta | awk '/>/{print "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"}1' | awk '{if (NR!=1) {print}}'

# awk '{if (NR!=1) {print}}' | awk '!/^>/ { printf "%s", $0; n = "\n" } ; /^>/ { print n $0; n = "" } ; END { printf "%s", n }' # this will convert everything except '>' to one single line
# fold -w 60 | grep '>LG2' -A 10 -B 10 # the fold will make each new line 60 nucleotides long (requirement for fasta), the grep will grep around to check - this is not great as any header longer than 60 will be a problem

nano list1
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta

#nano python-singlelinetomultilinefasta.py - do not use below python script, far too slow!!!!

#!/usr/bin/python

import sys

inputfile = sys.argv[1]
length = int(sys.argv[2])

outfile = open(inputfile.split(".fasta")[0] + '_multi-line.fasta', 'w') #open outfile for writing

with open(inputfile, 'r') as f:
     for line in f:
        if line.startswith(">"):
                print >> outfile, line.strip()
        else:
                sequence = line.strip()
                while len(sequence) > 0:
                        print >>outfile, sequence[:length]
                        sequence = sequence[length:]

# usage python python-singlelinetomultilinefasta.py [fasta] [nt# - normal is 60]
#python python-singlelinetomultilinefasta.py 20120125_MapAssembly.anchored.assembly_50bpN.fasta 60
#python python-singlelinetomultilinefasta.py M_zebra_v1.1_unscreened_3750.assembly_50bpN.fasta 60
#python python-singlelinetomultilinefasta.py H_burtoni_v1.assembly_50bpN.fasta 60
#python python-singlelinetomultilinefasta.py P_nyererei_v1.assembly_50bpN.fasta 60
#python python-singlelinetomultilinefasta.py N_brichardi_v1.assembly_50bpN.fasta 60


# create a shell script to run this
nano amendfasta.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

ml python
ml fastx_toolkit

while read F ; do awk '1;/>/{ printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"}' $F | awk '/>/{print "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"}1' | awk '{if (NR!=1) {print}}' | awk '!/^>/ { printf "%s", $0; n = "\n" } ; /^>/ { print n $0; n = "" } ; END { printf "%s", n }' > "$(basename "$F" .fasta)_50bpN.fasta" ; done < list1

fasta_formatter -i 20120125_MapAssembly.anchored.assembly_50bpN.fasta -o 20120125_MapAssembly.anchored.assembly_50bpN_multi-line.fasta -w 60 -e
fasta_formatter -i M_zebra_v1.1_unscreened_3750.assembly_50bpN.fasta -o M_zebra_v1.1_unscreened_3750.assembly_50bpN_multi-line.fasta -w 60 -e
fasta_formatter -i H_burtoni_v1.assembly_50bpN.fasta -o H_burtoni_v1.assembly_50bpN_multi-line.fasta -w 60 -e
fasta_formatter -i P_nyererei_v1.assembly_50bpN.fasta -o P_nyererei_v1.assembly_50bpN_multi-line.fasta -w 60 -e
fasta_formatter -i N_brichardi_v1.assembly_50bpN.fasta -o N_brichardi_v1.assembly_50bpN_multi-line.fasta -w 60 -e

# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=1 amendfasta.sh #{DONE}


## share folder with Tomasz /tgac/workarea/group-vh/cichlids/funSeq/cichlids/illumina-assemblies/Tarang50bpamend/14.Genome_sequences/*_50bpN_multi-line.fasta


###### 8. in silico motif validation
# 1. Randomize promoter sequences - Paddy running
# 2. Run motif scanning on these randomized promoter sequences
# 3. Test for enrichment in randomized vs experimental motif scanning

# run on uv2k2
uv2k2
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters
mkdir motifvalidation
cd motifvalidation

# copied randomised promoter sequences to this folder - Paddy generated these randomised sequences
ln -s ../JASPAR_CORE_2016_vertebrates.meme # create symbolic link to JASPAR file

# create a shell script to run motif scan on the promoters of one species
nano randomisedpromoter_motifdiscovery-Mz.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# M. zebra
fimo --o mzebra_shuffledpromoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Metriaclima_zebra.BROADMZ1.5kb_promoters_shuffled.fasta


# create a shell script to run motif scan on the promoters of one species
nano randomisedpromoter_motifdiscovery-Pn.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# P. nyererei
fimo --o pnyererei_shuffledpromoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Pundamilia_nyererei.BROADPN1.5kb_promoters_shuffled.fasta


# create a shell script to run motif scan on the promoters of one species
nano randomisedpromoter_motifdiscovery-Ab.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# A. burtoni
fimo --o aburtoni_shuffledpromoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Astatotilapia_burtoni.BROADAB1.5kb_promoters_shuffled.fasta


# create a shell script to run motif scan on the promoters of one species
nano randomisedpromoter_motifdiscovery-Nb.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# N. brichardi
fimo --o nbrichardi_shuffledpromoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Neolamprologus_brichardi.BROADNB1.5kb_promoters_shuffled.fasta


# create a shell script to run motif scan on the promoters of one species
nano randomisedpromoter_motifdiscovery-On.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-6)

# O. niloticus
fimo --o oniloticus_shuffledpromoters --thresh 1e-6 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Oreochromis_niloticus.BROADON1.5kb_promoters_shuffled.fasta


# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=1 randomisedpromoter_motifdiscovery-Mz.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 randomisedpromoter_motifdiscovery-Pn.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 randomisedpromoter_motifdiscovery-Ab.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 randomisedpromoter_motifdiscovery-Nb.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 randomisedpromoter_motifdiscovery-On.sh #{DONE}



# create a list file for the filtering
nano list1
mzebra_shuffledpromoters
pnyererei_shuffledpromoters
aburtoni_shuffledpromoters
nbrichardi_shuffledpromoters
oniloticus_shuffledpromoters


# Then filter FIMO output afterwards so that we only included matches with a q-value less than 0.1 {DONE}
while read F ; do awk '(NR==1) || ($8 < 0.1 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-fimo.txt ; done < list1 # this filters and gets rid of headers - you have to accept here that your usage of a relatively stringent q-value will filter out statistically insignificant hits, under this are CRX and OTX5 that have q-val>0.7. q-val of 0.1 implies that 10% of significant tests will result in false positives. I am using 0.05 as we have used vertebrate motifs to scan fish promoters, if we were scanning human/mouse then 0.01 would be ideal.

# Tar up each set of files and share with Paddy
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/motifvalidation
while read F ; do tar -zcvf randomisedpromotermotifs.tar.gz $F/$F-fimo.txt ; done < list1

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/
while read F ; do tar -zcvf originalpromotermotifs.tar.gz $F/$F-fimo.txt ; done < list1

/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/motifvalidation/randomisedpromotermotifs.tar.gz
/tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/originalpromotermotifs.tar.gz


### DO A COMPARISON OF p-val <1e-4 and <1e-6 and then different q-value filtering to observe if there are many differences

cd /tgac/scratch/mehtat/Cichlid_GRNs/
mkdir motif_discovery
mkdir motif_discovery/promoters
mkdir motif_discovery/promoters/1e-4
cd /tgac/scratch/mehtat/Cichlid_GRNs/motif_discovery/promoters/1e-4

# create symbolic links to the promoTer sequences:
for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/new_v1_5kb_promoters/* ; do ln -s $i ; done

# use FIMO in the MEME suite to scan for motif models using the JASPAR file
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/7.CNEs/JASPAR_CORE_2016_vertebrates.meme

# create a shell script to run motif scan on the promoters of one species
nano promoter_1e4motifdiscovery-Mz.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# M. zebra
fimo --o mzebra_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Metriaclima_zebra.BROADMZ1.5kb_promoters.0817.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_1e4motifdiscovery-Pn.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# P. nyererei
fimo --o pnyererei_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Pundamilia_nyererei.BROADPN1.5kb_promoters.0817.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_1e4motifdiscovery-Ab.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# A. burtoni
fimo --o aburtoni_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Astatotilapia_burtoni.BROADAB1.5kb_promoters.0817.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_1e4motifdiscovery-Nb.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# N. brichardi
fimo --o nbrichardi_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Neolamprologus_brichardi.BROADNB1.5kb_promoters.0817.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_1e4motifdiscovery-On.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# O. niloticus
fimo --o oniloticus_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Oreochromis_niloticus.BROADON1.5kb_promoters.0817.fasta


# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_1e4motifdiscovery-Mz.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_1e4motifdiscovery-Pn.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_1e4motifdiscovery-Ab.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_1e4motifdiscovery-Nb.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_1e4motifdiscovery-On.sh #{DONE}

## or run on slurm
nano promoter_1e4motifdiscovery-all.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# M. zebra
fimo --o mzebra_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Metriaclima_zebra.BROADMZ1.5kb_promoters.0817.fasta
# P. nyererei
fimo --o pnyererei_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Pundamilia_nyererei.BROADPN1.5kb_promoters.0817.fasta
# A. burtoni
fimo --o aburtoni_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Astatotilapia_burtoni.BROADAB1.5kb_promoters.0817.fasta
# N. brichardi
fimo --o nbrichardi_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Neolamprologus_brichardi.BROADNB1.5kb_promoters.0817.fasta
# O. niloticus
fimo --o oniloticus_1e4promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Oreochromis_niloticus.BROADON1.5kb_promoters.0817.fasta


nano sbatch-promoter_1e4motifdiscovery-all.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 2 # number of nodes
#SBATCH -n 2 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh promoter_1e4motifdiscovery-all.sh

# Run as batch script
sbatch sbatch-promoter_1e4motifdiscovery-all.sh

# Then filter FIMO output afterwards so that we only included matches with a q-value less than 0.05
while read F ; do awk '(NR==1) || ($8 < 0.1 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-q0.1fimo.txt ; done < list1 # this filters and gets rid of headers - you have to accept here that your usage of a relatively stringent q-value will filter out statistically insignificant hits, under this are CRX and OTX5 that have q-val>0.7. q-val of 0.1 implies that 10% of significant tests will result in false positives. I am using 0.05 as we have used vertebrate motifs to scan fish promoters, if we were scanning human/mouse then 0.01 would be ideal.

ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/TFBS_OGIDs3-Ab.txt
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/TFBS_OGIDs3-Pn.txt
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/TFBS_OGIDs3-Mz.txt
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/TFBS_OGIDs3-Nb.txt
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/TFBS_OGIDs3-On.txt

# awk match creating the whole-genome regulon files
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-Mz.txt mzebra_1e4promoters/mzebra_1e4promoters-q0.1fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > mzebra_1e4promoters/Mz-wholegenome_p1e4q0.1Regulon.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-Pn.txt pnyererei_1e4promoters/pnyererei_1e4promoters-q0.1fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > pnyererei_1e4promoters/Pn-wholegenome_p1e4q0.1Regulon.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-Ab.txt aburtoni_1e4promoters/aburtoni_1e4promoters-q0.1fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > aburtoni_1e4promoters/Ab-wholegenome_p1e4q0.1Regulon.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-Nb.txt nbrichardi_1e4promoters/nbrichardi_1e4promoters-q0.1fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > nbrichardi_1e4promoters/Nb-wholegenome_p1e4q0.1Regulon.txt
awk 'BEGIN{IGNORECASE=1}{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' TFBS_OGIDs3-On.txt oniloticus_1e4promoters/oniloticus_1e4promoters-q0.1fimo.txt | awk '{print $10, $2}' OFS='\t' | grep -v NA > oniloticus_1e4promoters/On-wholegenome_p1e4q0.1Regulon.txt


## Assess the number of unique TFs and interactions - 1e-4 vs 1e-6
cd /tgac/scratch/mehtat/Cichlid_GRNs/motif_discovery/promoters/1e-4


# create a list file for the filtering
nano list1
mzebra_1e4promoters
pnyererei_1e4promoters
aburtoni_1e4promoters
nbrichardi_1e4promoters
oniloticus_1e4promoters


# Then filter FIMO output - here we want to compare different q-vals
while read F ; do awk '(NR==1) || ($8 < 0.01 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-q0.01fimo.txt ; done < list1
while read F ; do awk '(NR==1) || ($8 < 0.05 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-q0.05fimo.txt ; done < list1
while read F ; do awk '(NR==1) || ($8 < 0.1 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-q0.1fimo.txt ; done < list1
while read F ; do awk '(NR==1) || ($8 < 0.2 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-q0.2fimo.txt ; done < list1

while read F ; do cut -f1 $F/fimo.txt | sort -u | wc -l ; done < list1
while read F ; do cut -f1 $F/$F-q0.01fimo.txt | sort -u | wc -l ; done < list1
while read F ; do cut -f1 $F/$F-q0.05fimo.txt | sort -u | wc -l ; done < list1
while read F ; do cut -f1 $F/$F-q0.1fimo.txt | sort -u | wc -l ; done < list1
while read F ; do cut -f1 $F/$F-q0.2fimo.txt | sort -u | wc -l ; done < list1

while read F ; do wc -l $F/fimo.txt ; done < list1
while read F ; do wc -l $F/$F-q0.01fimo.txt ; done < list1
while read F ; do wc -l $F/$F-q0.05fimo.txt ; done < list1
while read F ; do wc -l $F/$F-q0.1fimo.txt ; done < list1
while read F ; do wc -l $F/$F-q0.2fimo.txt ; done < list1

# compare against the 1e-6 outputs
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/
# Then filter FIMO output - here we want to compare different q-vals
while read F ; do awk '(NR==1) || ($8 < 0.01 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-q0.01fimo.txt ; done < list1
while read F ; do awk '(NR==1) || ($8 < 0.05 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-q0.05fimo.txt ; done < list1
while read F ; do awk '(NR==1) || ($8 < 0.1 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-q0.1fimo.txt ; done < list1
while read F ; do awk '(NR==1) || ($8 < 0.2 )' $F/fimo.txt | echo "$(tail -n +2)" > $F/$F-q0.2fimo.txt ; done < list1

while read F ; do cut -f1 $F/fimo.txt | sort -u | wc -l ; done < list1
while read F ; do cut -f1 $F/$F-q0.01fimo.txt | sort -u | wc -l ; done < list1
while read F ; do cut -f1 $F/$F-q0.05fimo.txt | sort -u | wc -l ; done < list1
while read F ; do cut -f1 $F/$F-q0.1fimo.txt | sort -u | wc -l ; done < list1
while read F ; do cut -f1 $F/$F-q0.2fimo.txt | sort -u | wc -l ; done < list1

while read F ; do wc -l $F/fimo.txt ; done < list1
while read F ; do wc -l $F/$F-q0.01fimo.txt ; done < list1
while read F ; do wc -l $F/$F-q0.05fimo.txt ; done < list1
while read F ; do wc -l $F/$F-q0.1fimo.txt ; done < list1
while read F ; do wc -l $F/$F-q0.2fimo.txt ; done < list1


## OWING TO ABOVE, RE-RUN WHOLE GENOME MOTIF SCAN USING 1E-4 FOR P-VAL

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/whole_genome
mkdir 1e-4
cd 1e-4

# use FIMO in the MEME suite to scan for motif models using the JASPAR file
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/7.CNEs/JASPAR_CORE_2016_vertebrates.meme

# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-OnIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use default 1e-4)

# O. niloticus – this is the old Illumina assembly [USE THIS]
fimo --o oniloticus_illumina --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-OnPacBio.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use default 1e-4)
# O. niloticus - this is the improved Illumina+PacBio assembly [DO NOT USE]
fimo --o oniloticus_illumina-pacbio --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-MzPacBio.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use default 1e-4)

# M. zebra – this is the improved Illumina+PacBio assembly [DO NOT USE]
fimo --o mzebra_illumina-pacbio --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-MzIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use default 1e-4)

# M. zebra – this is the old Illumina assembly [USE THIS]
fimo --o mzebra_illumina --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-AbIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use default 1e-4)

# A. burtoni
fimo --o aburtoni_illumina --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-PnIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use default 1e-4)

# P. nyererei
fimo --o pnyererei_illumina --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta


# create a shell script to run motif scan on the single genome
nano wholegenome_motifdiscovery-NbIllumina.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use default 1e-4)

# N. brichardi
fimo --o nbrichardi_illumina --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta


# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=2 wholegenome_motifdiscovery-OnIllumina.sh #{DONE}
qsub -q Prod -l select=1:mem=100GB:ncpus=2 wholegenome_motifdiscovery-MzIllumina.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=2 wholegenome_motifdiscovery-AbIllumina.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=2 wholegenome_motifdiscovery-PnIllumina.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=2 wholegenome_motifdiscovery-NbIllumina.sh #{DONE}
#qsub -q Test -l select=1:mem=100GB:ncpus=3 wholegenome_motifdiscovery-OnPacBio.sh # NOT RAN
#qsub -q Test -l select=1:mem=100GB:ncpus=3 wholegenome_motifdiscovery-MzPacBio.sh # NOT RAN


# Then filter FIMO output afterwards so that we only include matches with a q-value less than 0.1 (this is your FDR-based multiple testing correction). Note: FIMO output is 1-based to convert to 0-based by taking -1 from start e.g. 5-20 turns into 4-20
# Then convert this filtered file to a UCSC BED-formatted file (0-based). These files now contain the putative binding sites of all motifs from the motif models across the whole genome.
# BED format: chrom, chromStart, chromEnd, name, score, strand
while read F ; do awk '(NR==1) || ($8 < 0.1 )' $F/fimo.txt | echo "$(tail -n +2)" | awk '{print $2, $3-1, $4, $1, $6, $5}' OFS='\t' | sort -k1,1 -k2,2n | awk -v OFS="\t" 'NR==1{print $0, ".";next}{print $0,"."}' | awk '{print $1, $2, $3, $4, $7, $6, $4}' OFS='\t' > $F/$F-1e4fimo.bed ; done < list1 # this filters, gets rid of header, rearranges column, converts start coordinate to 0-based, adds a row of '.' and then sorts based on chrom and then chromStart column

#copy all bed files to one folder to share with Tomasz
mkdir all_1e4BED_files
while read F ; do cp $F/$F-1e4fimo.bed all_1e4BED_files ; done < list1 ## This was shared with Will for a new funSeq run (week of 24/04/17)


### Do a test of the M. zebra promoter to see whether the amendment of a T to a C (3029) in the sws1 promoter (mz.gene.s102.69) can give it a NR2F6, RXRB and RXRA binding site
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/
mkdir mzebra_sws1_test
cd mzebra_sws1_test


grep -A 999999 mz.gene.s102.69 ../Metriaclima_zebra.BROADMZ1.5kb_promoters.0817.fasta | awk 'NR>1 && /^>/{exit} 1' > mz.gene.s102.69_sws1.fasta
# sed 's/CTTGGGTTCATTTCTCGCATTCATCTCCCACTTCCTGTGACCTCTGACTCTGA/CTTGGGTTCATTTCTCGCATTCATCTCCCACTTCCTGTGACCTCTGACCCTGA/g' mz.gene.s102.69_sws1.fasta | grep --color CTTGGGTTCATTTCTCGCATTCATCTCCCACTTCCTGTGACCTCTGACCCTGA # do this to check the change - ok
sed 's/CTTGGGTTCATTTCTCGCATTCATCTCCCACTTCCTGTGACCTCTGACTCTGA/CTTGGGTTCATTTCTCGCATTCATCTCCCACTTCCTGTGACCTCTGACCCTGA/g' mz.gene.s102.69_sws1.fasta > mz.gene.s102.69_sws1_AMENDEDnt.fasta

# create a shell script to run motif scan on the single genome
nano Mz-sws1AMENDEDntmotifdiscovery.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use default 1e-4)

# M. zebra sws1
fimo --o mzebra_sws1AMENDEDnt --thresh 1e-4 --max-stored-scores 500000 ../JASPAR_CORE_2016_vertebrates.meme mz.gene.s102.69_sws1_AMENDEDnt.fasta

# run on uv2
qsub -q Test -l select=1:mem=50GB:ncpus=1 Mz-sws1AMENDEDntmotifdiscovery.sh

# filter based on q-val
awk '(NR==1) || ($8 < 0.1 )' mzebra_sws1AMENDEDnt/fimo.txt > mzebra_sws1AMENDEDnt/Mz-sws1_fimo_qval0.1.txt
# check if nr2f6 and rxr motifs are in there
grep MA0677.1 mzebra_sws1AMENDEDnt/Mz-sws1_fimo_qval0.1.txt # nr2f6
MA0677.1	mz.gene.s102.69	3018	3031	-	20.4394	6.13e-08	0.00049	AGGGTCAGAGGTCA
grep MA0728.1 mzebra_sws1AMENDEDnt/Mz-sws1_fimo_qval0.1.txt # nr2f6 (var2)
MA0728.1	mz.gene.s102.69	3018	3032	-	7.67143	3.05e-06	0.0243	CAGGGTCAGAGGTCA
grep MA0512.2 mzebra_sws1AMENDEDnt/Mz-sws1_fimo_qval0.1.txt # rxra
MA0512.2	mz.gene.s102.69	3018	3031	-	19.122	1.53e-07	0.00122	AGGGTCAGAGGTCA
grep MA0855.1 mzebra_sws1AMENDEDnt/Mz-sws1_fimo_qval0.1.txt # rxrb
MA0855.1	mz.gene.s102.69	3018	3031	-	19.7273	1.08e-07	0.000853	AGGGTCAGAGGTCA
MA0855.1	mz.gene.s102.69	2591	2604	-	4.94545	2.19e-05	0.0863	AGGGTTAAATGTGA

## all looks good, by changing the nt from T > C, it picks up the binding site
## in the original promoter file, there are several binding sites for NR2F6 and RXRa/b but none in the sws1 gene


############ 23/10/2017 - Paddy wants us to run fimo scan on new promoters using low threshold (5e-3) to use as a comparison for downstream analyses
uv2k2

cd /tgac/scratch/mehtat/Cichlid_GRNs/motif_discovery/promoters

### THIS IS USING THE NEW JASPAR 2018 motifs
mkdir 5e-3_Oct2017_Jaspar2018
cd 5e-3_Oct2017_Jaspar2018

# create symbolic links to the promoTer sequences:
for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta ; do ln -s $i ; done
ln -s /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta

# use FIMO in the MEME suite to scan for motif models using the JASPAR file
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters_JASPAR2018/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt

# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-Mz.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# M. zebra
fimo --o mzebra_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-Pn.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# P. nyererei
fimo --o pnyererei_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-Ab.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# A. burtoni
fimo --o aburtoni_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-Nb.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# N. brichardi
fimo --o nbrichardi_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-On.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# O. niloticus
fimo --o oniloticus_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR_CORE_2016_vertebrates.meme Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta


# Run all on uv2k2
qsub -q Test -l select=1:mem=50GB:ncpus=1 promoter_5e3motifdiscovery-Mz.sh #{DONE}
qsub -q Test -l select=1:mem=50GB:ncpus=1 promoter_5e3motifdiscovery-Pn.sh #{DONE}
qsub -q Test -l select=1:mem=50GB:ncpus=1 promoter_5e3motifdiscovery-Ab.sh #{DONE}
qsub -q Test -l select=1:mem=50GB:ncpus=1 promoter_5e3motifdiscovery-Nb.sh #{DONE}
qsub -q Test -l select=1:mem=50GB:ncpus=1 promoter_5e3motifdiscovery-On.sh #{DONE}

# filter fimo outputs, qval < 0.05
nano list
aburtoni_5e3promoters
mzebra_5e3promoters
nbrichardi_5e3promoters
oniloticus_5e3promoters
pnyererei_5e3promoters

while read F ; do awk '(NR==1) || ($8 < 0.05 )' $F/fimo.txt > $F/$F-fimo_qval0.05.txt ; done < list

# copy files to shared folder
cp -r /tgac/scratch/mehtat/Cichlid_GRNs/motif_discovery/ /tgac/workarea/Research-Groups/RG-cichlids


### THIS IS USING THE OLD JASPAR 2016 MOTIFS
mkdir 5e-3_Oct2017
cd 5e-3_Oct2017

# create symbolic links to the promoTer sequences:
for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta ; do ln -s $i ; done
ln -s /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta

# use FIMO in the MEME suite to scan for motif models using the JASPAR file
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v2/7.CNEs/JASPAR_CORE_2016_vertebrates.meme

# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-Mz.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# M. zebra
fimo --o mzebra_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-Pn.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# P. nyererei
fimo --o pnyererei_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-Ab.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# A. burtoni
fimo --o aburtoni_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-Nb.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# N. brichardi
fimo --o nbrichardi_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta


# create a shell script to run motif scan on the promoters of one species
nano promoter_5e3motifdiscovery-On.sh #{DONE}

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# use FIMO in the MEME suite to scan for motif models using the JASPAR meme file across the genomes, at a level of statistical significance (use stringent parameters of any motif instances of p-value <1e-4)

# O. niloticus
fimo --o oniloticus_5e3promoters --thresh 1e-4 --max-stored-scores 500000 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta


# Run all on uv2k2
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_5e3motifdiscovery-Mz.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_5e3motifdiscovery-Pn.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_5e3motifdiscovery-Ab.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_5e3motifdiscovery-Nb.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 promoter_5e3motifdiscovery-On.sh #{DONE}

# filter fimo outputs, qval < 0.05
nano list
aburtoni_5e3promoters
mzebra_5e3promoters
nbrichardi_5e3promoters
oniloticus_5e3promoters
pnyererei_5e3promoters

while read F ; do awk '(NR==1) || ($8 < 0.05 )' $F/fimo.txt > $F/$F-fimo_qval0.05.txt ; done < list

# copy files to shared folder
cp -r /tgac/scratch/mehtat/Cichlid_GRNs/motif_discovery/ /tgac/workarea/Research-Groups/RG-cichlids

#### 9. For the p-value imputed motif scanning, we need some additional files

# First are 10 orthology files, two for each cichlid to both human and mouse
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr

# Cut the required columns, remove columns where NULL gene IDs, then assign geneID to gene name if gene name is NULL
# Mz
awk '{print $14, $2, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-MzHs # human
awk '{print $16, $2, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-MzMm # mouse

# Pn
awk '{print $14, $3, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-PnHs # human
awk '{print $16, $3, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-PnMm # mouse

# Ab
awk '{print $14, $4, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-AbHs # human
awk '{print $16, $4, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-AbMm # mouse

# Nb
awk '{print $14, $5, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-NbHs # human
awk '{print $16, $5, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-NbMm # mouse

# On
awk '{print $14, $6, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-OnHs # human
awk '{print $16, $6, $15, $10}' OFS=';' OGIDS.txt5 | awk -F';' '$1!="NULL" && $2!="NULL"' OFS=';' | awk -F';' '{ $3 = ($3 == "NULL" ? $1 : $3) } 1' OFS=';' | awk -F';' '{ $4 = ($4 == "NULL" ? $2 : $4) } 1' OFS=';' > OGIDS.txt5-OnMm # mouse

## Also for p-value imputed motif scanning, we also need to carry out domain profiling
# Domain annotation - run an interpro scan for the cichlid TFs and their orthologs in mice and human in one batch

uv2k2
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/
mkdir TFdomain_scan
cd TFdomain_scan
ln -s ../../Arboretum_GT_v3/TFBS_OGIDs3.txt

# the protein sequences are not sorted by the longest transcript being first - hence your script may not take the longest (.1), but instead, take the first hit. So sort the peptide fasta files
for i in ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/*.fa ; do awk 'BEGIN{RS=">"} NR>1 {gsub("\n", "\t"); print ">"$0}' $i | sort -k1,1 | awk '{sub("\t", "\n"); gsub("\t", ""); print $0}' > ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/"$(basename "$i" .fa).sorted.fa" ; done

# extract protein sequences - cut the columns for each species, read line by line and extract then take top transcript, then remove stars as stop as they will interfere
cut -f4 TFBS_OGIDs3.txt > TFBS_OGIDs3-Mz.txt ; while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Metriaclima_zebra.BROADMZ1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < TFBS_OGIDs3-Mz.txt | sed 's/*//g' > TFBS_OGIDs3-Mz.fa # Mz
cut -f5 TFBS_OGIDs3.txt > TFBS_OGIDs3-Pn.txt ; while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Pundamilia_nyererei.BROADPN1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < TFBS_OGIDs3-Pn.txt | sed 's/*//g' > TFBS_OGIDs3-Pn.fa # Pn
cut -f6 TFBS_OGIDs3.txt > TFBS_OGIDs3-Ab.txt ; while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Astatotilapia_burtoni.BROADAB1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < TFBS_OGIDs3-Ab.txt | sed 's/*//g' > TFBS_OGIDs3-Ab.fa # Ab
cut -f7 TFBS_OGIDs3.txt > TFBS_OGIDs3-Nb.txt ; while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Neolamprologus_brichardi.BROADNB1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < TFBS_OGIDs3-Nb.txt | sed 's/*//g' > TFBS_OGIDs3-Nb.fa # Nb
cut -f8 TFBS_OGIDs3.txt > TFBS_OGIDs3-On.txt ; while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Oreochromis_niloticus.BROADON1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < TFBS_OGIDs3-On.txt | sed 's/*//g' > TFBS_OGIDs3-On.fa # On

# download all human and protein sequences
software
cd /tgac/workarea/group-vh/Tarang/Reference_Genomes/
mkdir Musmusculus_GRCm38
cd hg38_chr/
wget ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
cd ../Musmusculus_GRCm38
wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz
gunzip Mus_musculus.GRCm38.pep.all.fa.gz
exit

# sort the fasta files by length - that way, the first hit of any ID should be the longest
# use awk to linearize the fasta; a second awk to insert a column with the length of the sequence; sort on first column; remove the first colum; convert back to fasta
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' /tgac/workarea/group-vh/Tarang/Reference_Genomes/hg38_chr/Homo_sapiens.GRCh38.pep.all.fa | awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' | sort -rn -k1,1 | cut -f 2- | tr "\t" "\n" > /tgac/workarea/group-vh/Tarang/Reference_Genomes/hg38_chr/Homo_sapiens.GRCh38.pep.all.sorted.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' /tgac/workarea/group-vh/Tarang/Reference_Genomes/Musmusculus_GRCm38/Mus_musculus.GRCm38.pep.all.fa | awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' | sort -rn -k1,1 | cut -f 2- | tr "\t" "\n" > /tgac/workarea/group-vh/Tarang/Reference_Genomes/Musmusculus_GRCm38/Mus_musculus.GRCm38.pep.all.sorted.fa

ln -s /tgac/workarea/group-vh/Tarang/Reference_Genomes/hg38_chr/Homo_sapiens.GRCh38.pep.all.sorted.fa
ln -s /tgac/workarea/group-vh/Tarang/Reference_Genomes/Musmusculus_GRCm38/Mus_musculus.GRCm38.pep.all.sorted.fa

# extract human and mouse protein sequences - taking the firstmost hit which should be the longest
cut -f16 TFBS_OGIDs3.txt > TFBS_OGIDs3-Hs.txt ; while read F ; do grep -A 999999 -wiF $F Homo_sapiens.GRCh38.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < TFBS_OGIDs3-Hs.txt | sed 's/*//g' > TFBS_OGIDs3-Hs.fa
cut -f18 TFBS_OGIDs3.txt > TFBS_OGIDs3-Mm.txt ; while read F ; do grep -A 999999 -wiF $F Mus_musculus.GRCm38.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < TFBS_OGIDs3-Mm.txt | sed 's/*//g' > TFBS_OGIDs3-Mm.fa


nano interpro5_scan-Hs.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml java
source interproscan-5

interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i TFBS_OGIDs3-Hs.fa -b HsOUT

nano sbatch-interpro5_scan-Hs.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh interpro5_scan-Hs.sh

# Run as batch script
sbatch sbatch-interpro5_scan-Hs.sh


nano interpro5_scan-Mm.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml java
source interproscan-5

interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i TFBS_OGIDs3-Mm.fa -b MmOUT

nano sbatch-interpro5_scan-Mm.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh interpro5_scan-Mm.sh

# Run as batch script
sbatch sbatch-interpro5_scan-Mm.sh


nano interpro5_scan-Mz.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml java
source interproscan-5

interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i TFBS_OGIDs3-Mz.fa -b MzOUT

nano sbatch-interpro5_scan-Mz.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh interpro5_scan-Mz.sh

# Run as batch script
sbatch sbatch-interpro5_scan-Mz.sh


nano interpro5_scan-Pn.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml java
source interproscan-5

interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i TFBS_OGIDs3-Pn.fa -b PnOUT

nano sbatch-interpro5_scan-Pn.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh interpro5_scan-Pn.sh

# Run as batch script
sbatch sbatch-interpro5_scan-Pn.sh


nano interpro5_scan-Ab.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml java
source interproscan-5

interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i TFBS_OGIDs3-Ab.fa -b AbOUT

nano sbatch-interpro5_scan-Ab.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh interpro5_scan-Ab.sh

# Run as batch script
sbatch sbatch-interpro5_scan-Ab.sh


nano interpro5_scan-Nb.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml java
source interproscan-5

interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i TFBS_OGIDs3-Nb.fa -b NbOUT

nano sbatch-interpro5_scan-Nb.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh interpro5_scan-Nb.sh

# Run as batch script
sbatch sbatch-interpro5_scan-Nb.sh


nano interpro5_scan-On.sh

# Add the following to above script
#!/bin/bash
#load the latest module
ml java
source interproscan-5

interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i TFBS_OGIDs3-On.fa -b OnOUT

nano sbatch-interpro5_scan-On.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh interpro5_scan-On.sh

# Run as batch script
sbatch sbatch-interpro5_scan-On.sh

# insert file for colheaders
nano interpro5_columns.txt

# 1. Protein Accession (e.g. P51587)
# 2. Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
# 3. Sequence Length (e.g. 3418)
# 4. Analysis (e.g. Pfam / PRINTS / Gene3D)
# 5. Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
# 6. Signature Description (e.g. BRCA2 repeat profile)
# 7. Start location
# 8. Stop location
# 9. Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
# 10. Status - is the status of the match (T: true)
# 11. Date - is the date of the run
# 12. (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
# 13. (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
# 14. (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
# 15. (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)

# filter based on Pfam and presence of interpro annotation - this means that it has passed interpro quality checks and unlikely to be a false match
# fill all empty spaces with NA (perl), then grep for pfam then grep for 'IPR' annotation, then change gene IDs so that they match those in OGIDS files
for i in *.tsv ; do perl -pe 's/\t(?=\t)/\tNA/g' $i | sed 's/ /_/g' | grep -wiF Pfam | grep 'IPR' | sed 's/mrna/gene/g' | awk '{gsub(/\.[0-9]$/,"",$1); print $0}' OFS='\t' | awk '{gsub(/^ENSP/,"ENSG",$1); print $0}' OFS='\t' | awk '{gsub(/^ENSMUSP/,"ENSMUSG",$1); print $0}' OFS='\t' > "$(basename "$i" .tsv).pfamfilt.tsv" ; done
for i in *pfamfilt.tsv ; do cut -f1-14 $i > "$(basename "$i" .tsv)2.tsv" ; done # some rows do not have pathway annotation so remove
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz MzOUT.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > MzOUT.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn PnOUT.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > PnOUT.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab AbOUT.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > AbOUT.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb NbOUT.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > NbOUT.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on OnOUT.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > OnOUT.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-hs HsOUT.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > HsOUT.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mm MmOUT.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > MmOUT.OGID.pfamfilt2.tsv
# check how many Hs and Mm genes are not orthology mapped to cichlid genes (of the 18799 orthogroups, 5888 are not mapped to human).
grep -v NA HsOUT.OGID.pfamfilt2.tsv | wc -l # 3
grep -v NA MmOUT.OGID.pfamfilt2.tsv  | wc -l # 42

## Also do an interpro scan of all genes
# extract the longest protein coding for each gene
# create a file to use as grep to take the .1 transcripts
for i in /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/*.pep.all.sorted.fa ; do sort -u -t' ' -k2,2 $i | awk -F' ' '{print $1}' > "$(basename "$i" pep.all.sorted.fa)longest" ; done
# remove asterix and dots randomly in middle of sequence
nano run1.sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 12000 # memory pool for all cores
#SBATCH -t 0-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Metriaclima_zebra.BROADMZ1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < Metriaclima_zebra.BROADMZ1.longest | sed 's/*//g' | sed '/^>/! s/\.//g' > Metriaclima_zebra.BROADMZ1.longest.fa # Mz
while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Pundamilia_nyererei.BROADPN1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < Pundamilia_nyererei.BROADPN1.longest | sed 's/*//g' | sed '/^>/! s/\.//g' > Pundamilia_nyererei.BROADPN1.longest.fa # Pn
while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Astatotilapia_burtoni.BROADAB1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < Astatotilapia_burtoni.BROADAB1.longest | sed 's/*//g' | sed '/^>/! s/\.//g' > Astatotilapia_burtoni.BROADAB1.longest.fa # Ab
while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Neolamprologus_brichardi.BROADNB1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < Neolamprologus_brichardi.BROADNB1.longest | sed 's/*//g' | sed '/^>/! s/\.//g' > Neolamprologus_brichardi.BROADNB1.longest.fa # Nb
while read F ; do grep -A 999999 -wiF $F ../../../../cichlids/data/Broad/all/Data/Annotation/protein_coding/v1/pep/Oreochromis_niloticus.BROADON1.pep.all.sorted.fa | awk 'NR>1 && /^>/{exit} 1' ; done < Oreochromis_niloticus.BROADON1.longest | sed 's/*//g' | sed '/^>/! s/\.//g' > Oreochromis_niloticus.BROADON1.longest.fa # On

sbatch run1.sh

# remove asterix from Hs and Mm
sed -i 's/*//g' Homo_sapiens.GRCh38.pep.all.sorted.fa
sed -i 's/*//g' Mus_musculus.GRCm38.pep.all.sorted.fa

# Filter the above for the longest transcript for each gene
awk -F' ' '{print $1,$4}' Homo_sapiens.GRCh38.pep.all.sorted.fa | # pull out the protein ID, gene ID and sequence
awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | # linearize the fasta
awk -F' ' '{print ">"$2,$3}' OFS='\t'| # get rid of the protein ID from each line
sed 's/gene://g' | # get rid of 'gene' from each header
tr "." "\t" | #extract transcript no. from header - creates 2 columns, 3rd is the sequence
awk -F '\t' '{printf("%s\t%d\n",$0,length($3));}' | #extract length from sequence (3rd column), and add as seq length as 4th column
sort -k1,1 -k4,4nr | #sort on name, inverse sequence length
sort -k1,1 -u -s | #sort on name, unique, stable sort keeping previous order
sed 's/\t/./' | #restore name
cut -f 1,2 | #cut name, sequence
tr "\t" "\n" | #back to fasta
fold -w 60 > Homo_sapiens.GRCh38.pep.all.sorted.longest.fa #pretty fasta - 23,043 proteins

awk -F' ' '{print $1,$4}' Mus_musculus.GRCm38.pep.all.sorted.fa | # pull out the protein ID, gene ID and sequence
awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | # linearize the fasta
awk -F' ' '{print ">"$2,$3}' OFS='\t'| # get rid of the protein ID from each line
sed 's/gene://g' | # get rid of 'gene' from each header
tr "." "\t" | #extract transcript no. from header - creates 2 columns, 3rd is the sequence
awk -F '\t' '{printf("%s\t%d\n",$0,length($3));}' | #extract length from sequence (3rd column), and add as seq length as 4th column
sort -k1,1 -k4,4nr | #sort on name, inverse sequence length
sort -k1,1 -u -s | #sort on name, unique, stable sort keeping previous order
sed 's/\t/./' | #restore name
cut -f 1,2 | #cut name, sequence
tr "\t" "\n" | #back to fasta
fold -w 60 > Mus_musculus.GRCm38.pep.all.sorted.longest.fa #pretty fasta - 22,784 proteins


cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan

### run as an array, otherwise it will take many days!
# check how many proteins in each file
for i in *.longest.fa ; do grep '>' $i | wc -l ; done # 20611-24559 - all will need five folders for splitting, so a-e

# create folders for splitting each fasta line by line but in 5000 chunks as 6000 jobs is maximal for submitting to the queues
mkdir fastasplit
mkdir fastasplit/mz
mkdir fastasplit/pn
mkdir fastasplit/ab
mkdir fastasplit/nb
mkdir fastasplit/on
mkdir fastasplit/hs
mkdir fastasplit/mm

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/mz

# Mz
nano interpro5_scan-a_dMz-array.sh

#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-4999
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-23:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 a/*.fa > mzreads_a # create a list of all fastq files
mapfile -t mzreads_a < mzreads_a # assign as elements to $mzreads_a variable
ls -1 b/*.fa > mzreads_b # create a list of all fastq files
mapfile -t mzreads_b < mzreads_b # assign as elements to $mzreads_b variable
ls -1 c/*.fa > mzreads_c # create a list of all fastq files
mapfile -t mzreads_c < mzreads_c # assign as elements to $mzreads_c variable
ls -1 d/*.fa > mzreads_d # create a list of all fastq files
mapfile -t mzreads_d < mzreads_d # assign as elements to $mzreads_d variable

#load the latest module
ml java
source interproscan-5

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mzreads_a[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mzreads_b[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mzreads_c[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mzreads_d[${SLURM_ARRAY_TASK_ID}]}

# run
sbatch interpro5_scan-a_dMz-array.sh

# create another for the remaining smaller number of files
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/mz

nano interpro5_scan-eMz-array.sh

#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-1672
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-23:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > mzreads_e # create a list of all fastq files
mapfile -t mzreads_e < mzreads_e # assign as elements to $mzreads_a variable

#load the latest module
ml java
source interproscan-5

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mzreads_e[${SLURM_ARRAY_TASK_ID}]}

# run
sbatch interpro5_scan-eMz-array.sh


## Separate out the files of other species so that the job arrays can be run at any point

nano run3.sh

#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 12000 # memory pool for all cores
#SBATCH -t 0-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/mz
mkdir a b c d e
awk '/^>mz/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Metriaclima_zebra.BROADMZ1.longest.fa # split the fasta by each protein
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #1673

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/
cd fastasplit/pn
mkdir a b c d e
awk '/^>pn/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Pundamilia_nyererei.BROADPN1.longest.fa # split the fasta by each protein
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #611

cd ../ab
mkdir a b c d e
awk '/^>ab/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Astatotilapia_burtoni.BROADAB1.longest.fa # split the fasta by each protein
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #3436

cd ../nb
mkdir a b c d e
awk '/^>nb/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Neolamprologus_brichardi.BROADNB1.longest.fa # split the fasta by each protein
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #119

cd ../on
mkdir a b c d e
awk '/^>on/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Oreochromis_niloticus.BROADON1.longest.fa # split the fasta by each protein
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #4559

cd ../hs
mkdir a b c d e
awk '/^>ENSG/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Homo_sapiens.GRCh38.pep.all.sorted.longest.fa # split the fasta by each protein
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #3043

cd ../mm
mkdir a b c d e
awk '/^>ENSMUSG/ {OUT=substr($1,2) ".fa"}; {print >> OUT; close(OUT)}' ../../Mus_musculus.GRCm38.pep.all.sorted.longest.fa # split the fasta by each protein
ls *.fa | head -5000 | xargs -I{} mv {} a/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} b/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} c/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} d/ #5000
ls *.fa | head -5000 | xargs -I{} mv {} e/ #2784


# run
sbatch run3.sh


## create one job script for the fixed 5000 jobs, and separate scripts for the smaller folder e's

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/

nano interpro5_scan-a_dPnAbNbOnHsMm-array.sh

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-4999
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-23:59
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

#load the latest module
ml java
source interproscan-5

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${pnreads_a[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${pnreads_b[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${pnreads_c[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${pnreads_d[${SLURM_ARRAY_TASK_ID}]}

cd ../ab

ls -1 a/*.fa > pnreads_a # create a list of all fastq files
mapfile -t pnreads_a < pnreads_a # assign as elements to $pnreads_a variable
ls -1 b/*.fa > pnreads_b # create a list of all fastq files
mapfile -t pnreads_b < pnreads_b # assign as elements to $pnreads_b variable
ls -1 c/*.fa > pnreads_c # create a list of all fastq files
mapfile -t pnreads_c < pnreads_c # assign as elements to $pnreads_c variable
ls -1 d/*.fa > pnreads_d # create a list of all fastq files
mapfile -t pnreads_d < pnreads_d # assign as elements to $pnreads_d variable

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${abreads_a[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${abreads_b[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${abreads_c[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${abreads_d[${SLURM_ARRAY_TASK_ID}]}

cd ../nb

ls -1 a/*.fa > nbreads_a # create a list of all fastq files
mapfile -t nbreads_a < nbreads_a # assign as elements to $nbreads_a variable
ls -1 b/*.fa > nbreads_b # create a list of all fastq files
mapfile -t nbreads_b < nbreads_b # assign as elements to $nbreads_b variable
ls -1 c/*.fa > nbreads_c # create a list of all fastq files
mapfile -t nbreads_c < nbreads_c # assign as elements to $nbreads_c variable
ls -1 d/*.fa > nbreads_d # create a list of all fastq files
mapfile -t nbreads_d < nbreads_d # assign as elements to $nbreads_d variable

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${nbreads_a[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${nbreads_b[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${nbreads_c[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${nbreads_d[${SLURM_ARRAY_TASK_ID}]}

cd ../on

ls -1 a/*.fa > onreads_a # create a list of all fastq files
mapfile -t onreads_a < onreads_a # assign as elements to $onreads_a variable
ls -1 b/*.fa > onreads_b # create a list of all fastq files
mapfile -t onreads_b < onreads_b # assign as elements to $onreads_b variable
ls -1 c/*.fa > onreads_c # create a list of all fastq files
mapfile -t onreads_c < onreads_c # assign as elements to $onreads_c variable
ls -1 d/*.fa > onreads_d # create a list of all fastq files
mapfile -t onreads_d < onreads_d # assign as elements to $onreads_d variable

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${onreads_a[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${onreads_b[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${onreads_c[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${onreads_d[${SLURM_ARRAY_TASK_ID}]}

cd ../hs

ls -1 a/*.fa > hsreads_a # create a list of all fastq files
mapfile -t hsreads_a < hsreads_a # assign as elements to $hsreads_a variable
ls -1 b/*.fa > hsreads_b # create a list of all fastq files
mapfile -t hsreads_b < hsreads_b # assign as elements to $hsreads_b variable
ls -1 c/*.fa > hsreads_c # create a list of all fastq files
mapfile -t hsreads_c < hsreads_c # assign as elements to $hsreads_c variable
ls -1 d/*.fa > hsreads_d # create a list of all fastq files
mapfile -t hsreads_d < hsreads_d # assign as elements to $hsreads_d variable

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${hsreads_a[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${hsreads_b[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${hsreads_c[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${hsreads_d[${SLURM_ARRAY_TASK_ID}]}

cd ../mm

ls -1 a/*.fa > mmreads_a # create a list of all fastq files
mapfile -t mmreads_a < mmreads_a # assign as elements to $mmreads_a variable
ls -1 b/*.fa > mmreads_b # create a list of all fastq files
mapfile -t mmreads_b < mmreads_b # assign as elements to $mmreads_b variable
ls -1 c/*.fa > mmreads_c # create a list of all fastq files
mapfile -t mmreads_c < mmreads_c # assign as elements to $mmreads_c variable
ls -1 d/*.fa > mmreads_d # create a list of all fastq files
mapfile -t mmreads_d < mmreads_d # assign as elements to $mmreads_d variable

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mmreads_a[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mmreads_b[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mmreads_c[${SLURM_ARRAY_TASK_ID}]}
srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mmreads_d[${SLURM_ARRAY_TASK_ID}]}


# create other scripts for the remaining smaller number of files
# Pn -e
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/pn

nano interpro5_scan-ePn-array.sh

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-610
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-23:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > pnreads_e # create a list of all fastq files
mapfile -t pnreads_e < pnreads_e # assign as elements to $mzreads_a variable

#load the latest module
ml java
source interproscan-5

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${pnreads_e[${SLURM_ARRAY_TASK_ID}]}

# Ab-e
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/ab

nano interpro5_scan-eAb-array.sh

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-3435
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-23:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > abreads_e # create a list of all fastq files
mapfile -t abreads_e < abreads_e # assign as elements to $mzreads_a variable

#load the latest module
ml java
source interproscan-5

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${abreads_e[${SLURM_ARRAY_TASK_ID}]}

# Nb-e
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/nb

nano interpro5_scan-eNb-array.sh

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-118
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-23:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > nbreads_e # create a list of all fastq files
mapfile -t nbreads_e < nbreads_e # assign as elements to $mzreads_a variable

#load the latest module
ml java
source interproscan-5

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${nbreads_e[${SLURM_ARRAY_TASK_ID}]}


# On-e
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/on

nano interpro5_scan-eOn-array.sh

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-4558
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-23:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > onreads_e # create a list of all fastq files
mapfile -t onreads_e < onreads_e # assign as elements to $mzreads_a variable

#load the latest module
ml java
source interproscan-5

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${onreads_e[${SLURM_ARRAY_TASK_ID}]}

# Hs-e
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/hs

nano interpro5_scan-eHs-array.sh

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-3042
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-23:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > hsreads_e # create a list of all fastq files
mapfile -t hsreads_e < hsreads_e # assign as elements to $mzreads_a variable

#load the latest module
ml java
source interproscan-5

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${hsreads_e[${SLURM_ARRAY_TASK_ID}]}

# Mm-e
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/mm

nano interpro5_scan-eMm-array.sh

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --array=0-2783
#SBATCH --mem-per-cpu 8240
#SBATCH -t 0-23:59
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

ls -1 e/*.fa > mmreads_e # create a list of all fastq files
mapfile -t mmreads_e < mmreads_e # assign as elements to $mzreads_a variable

#load the latest module
ml java
source interproscan-5

srun interproscan.sh -mode cluster -clusterrunid uniqueName -iprlookup --goterms --pathways -i ${mmreads_e[${SLURM_ARRAY_TASK_ID}]}


# to run
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/
sbatch interpro5_scan-a_dPnAbNbOnHsMm-array.sh # 5000 * 6 - {DONE} ran 01/02/18 10.00

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/mz/
sbatch interpro5_scan-eMz-array.sh # 1673

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/pn/
sbatch interpro5_scan-ePn-array.sh #611 - {DONE}

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/
sbatch interpro5_scan-a_dHsMm-array.sh # {DONE}

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/ab/
sbatch interpro5_scan-eAb-array.sh #3436 #

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/hs/
sbatch interpro5_scan-eHs-array.sh #1043 # {DONE} ran 29/10 17:06

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/nb/
sbatch interpro5_scan-eNb-array.sh #119 # {DONE} generated 99 tsv

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/on/
sbatch interpro5_scan-eOn-array.sh #4559 #

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/mm/
sbatch interpro5_scan-eMm-array.sh #784 # {DONE} ran 29/10 17:08


# join all individual files for each species
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/TFdomain_scan/fastasplit/
nano list
ab
mz
pn
nb
on
hs
mm


nano interpro5results_cat_tar_filt.sh

#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 12480 # memory pool for all cores
#SBATCH -t 0-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

# cat all reads for each species
while read F ; do cat $F/*.fa.tsv > $F-interpro_collated.tsv ; done < list
# tar up all the individual runs and fasta files
while read F ; do tar -zcvf $F/$F-interpro_individual.tar.gz $F/ ; done < list

# filter based on Pfam and presence of interpro annotation - this means that it has passed interpro quality checks and unlikely to be a false match
# fill all empty spaces with NA (perl), then grep for pfam then grep for 'IPR' annotation, then change gene IDs so that they match those in OGIDS files
for i in *.tsv ; do perl -pe 's/\t(?=\t)/\tNA/g' $i | sed 's/ /_/g' | grep -wiF Pfam | grep 'IPR' | sed 's/mrna/gene/g' | awk '{gsub(/\.[0-9]$/,"",$1); print $0}' OFS='\t' | awk '{gsub(/^ENSP/,"ENSG",$1); print $0}' OFS='\t' | awk '{gsub(/^ENSMUSP/,"ENSMUSG",$1); print $0}' OFS='\t' > "$(basename "$i" .tsv).pfamfilt.tsv" ; done
for i in *pfamfilt.tsv ; do cut -f1-14 $i > "$(basename "$i" .tsv)2.tsv" ; done # some rows do not have pathway annotation so remove
for i in *.pfamfilt.tsv ; do sed -i 's/prot/gene/g' $i ; done
for i in *.pfamfilt2.tsv ; do sed -i 's/prot/gene/g' $i ; done


awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz mz-interpro_collated.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > mz-interpro_collated.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-pn pn-interpro_collated.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > pn-interpro_collated.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-ab ab-interpro_collated.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > ab-interpro_collated.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb nb-interpro_collated.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > nb-interpro_collated.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-on on-interpro_collated.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > on-interpro_collated.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-hs hs-interpro_collated.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > hs-interpro_collated.OGID.pfamfilt2.tsv
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../../Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mm mm-interpro_collated.pfamfilt2.tsv | awk '{print $15,$1,$5,$9,$12,$13}' OFS='\t' > mm-interpro_collated.OGID.pfamfilt2.tsv

# run the above
sbatch interpro5results_cat_tar_filt.sh

# problem is that several human genes are not orthology mapped to cichlid genes (of the 18799 orthogroups, 5888 are not mapped to human).
grep -v NA hs-interpro_collated.OGID.pfamfilt2.tsv | wc -l #4188
grep -v NA mm-interpro_collated.OGID.pfamfilt2.tsv  | wc -l #666

# filter for NA
for i in *OGID.pfamfilt2.tsv ; do grep -v NA $i > "$(basename "$i" 2.tsv)3.tsv" ; done

# copy files to shared folder
mkdir /tgac/workarea/Research-Groups/RG-cichlids/domain_scan
cp *-interpro_collated.tsv /tgac/workarea/Research-Groups/RG-cichlids/domain_scan/ # native output
cp *-interpro_collated.pfamfilt2.tsv /tgac/workarea/Research-Groups/RG-cichlids/domain_scan/ # filtered output for presence of pfam and intepro annotation - this means that it has passed interpro quality checks and unlikely to be a false match
cp *-interpro_collated.OGID.pfamfilt3.tsv /tgac/workarea/Research-Groups/RG-cichlids/domain_scan/ # parsed the filtered outputs above for only co-expressed genes in modules (and hence they have an orthogroup) - this significantly drops the Hs and Mmm numbers due to lack of orthogroup information

# The TSV format presents the match data in columns as follows:
#     1.Protein Accession (e.g. P51587)
#     2.Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
#     3.Sequence Length (e.g. 3418)
#     4.Analysis (e.g. Pfam / PRINTS / Gene3D)
#     5.Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
#     6.Signature Description (e.g. BRCA2 repeat profile)
#     7.Start location
#     8.Stop location
#     9.Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
#     10.Status - is the status of the match (T: true)
#     11.Date - is the date of the run
#     12.(InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
#     13.(InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
#     14.(GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)


##########################################################################################################################

## Since we extracted cichlid sequences from 8-way MAF, we need to deal with (gene) strand for when we run FIMO (see NetworkReconstruction_v5.sh - Evolutionary rate section)
# Tilapia used as reference in +ve strand, so if the Tilapia gene is -ve strand (and so is the promoter), we reverse complement all sequences in the file
# BED and MAF files here:
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/5kbpromoter_extractfrom8wayMAF/oniloticus_annot

# Cichlid FASTA files here:
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/5kbpromoter_extractfrom8wayMAF/oniloticus_annot/fasta

##########################################################################################################################

# 22/08/17 - trying another approach since with the extracted MAF approach we will have problems using the correct sequence when scanning for motifs
###### This is also in your NetworkReconstruction_v5.sh script

# 1. Take the cds of the first exon of the longest transcript of each gene in O. niloticus
# 2. Use the above to BLAT against the genomes of the other species to find the corresponding start site of its actual orthologous gene
# 3. See if this corresponds to the original annotations and OGIDs
# 4. Using the match as a reference, take upto 5kb upstream of the CDS sequence as promoter and then align > PAML > evol rate
# 5. The same sequences can be used for motif scanning

cd /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids
mkdir 1stexonBLAT
cd 1stexonBLAT

# 1. Take the cds of the first exon of the longest transcript of each gene in O. niloticus

# This calculates the longest transcript of each parent gene in the O. niloticus genome, then grepping on the gff, prepares a BED file of just the first exon of the longest transcript
# - THIS WILL BE INCORRECT AS IT TAKES THE LONGEST MRNA ONLY AND NOT NECESSARILY THE LONGEST PROTEIN, DO NOT USE!
awk '$3=="mRNA" && $9~/^ID/{sub(/_T.*/,"",$9);L[substr($9,4)]+=$5-$4+1}END{for(i in L){print i,L[i]}}' OFS='\t' /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/final_annotation/Oreochromis_niloticus.BROADON1.gff3 | sort -k2,2 -rn | sed 's/\;/\t/g' | sort -k2,2 -u > On_genome_longest_transcripts.txt
cut -f1 On_genome_longest_transcripts.txt > On_genome_longest_transcripts.txt2
split -l 1000 On_genome_longest_transcripts.txt2 # split the file every 1000 lines to make the grep quicker
### MAYBE DO NOT USE ANY OF THIS - THIS TAKES THE LONGEST mRNA BUT NOT THE LONGEST PROTEIN - FOR THAT YOU WOULD HAVE TO ADD UP THE EXONS AND GO FROM THERE

nano On_1stexon-bed.sh

#!/bin/sh

# grep the longest transcript information
for i in xa* ; do grep -wiFf $i /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/final_annotation/Oreochromis_niloticus.BROADON1.gff3 | awk '{print $1, $4-1, $5, $9, "0", $7, $3}' OFS='\t' | sed 's/ID.*Parent=//g' > $i.grep.bed ; done
cat *.grep.bed > On_genome_longest_transcript.bed
# extract the cds corresponding to first exon +/- 100bp, then replace any negative values with a 0 and check for the ends of scaffolds; taking the end coordinate if your end coordinate goes over
# awk '$2<0 {$2=0} 1' OFS='\t'


# then extract from FASTA using bedtools
# ml bedtools/2.25.0
#bedtools getfasta -fo On_genome_longest_transcript_1stexon-and100bp.fasta -name -fi ../Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -bed On_genome_longest_transcript_1stexon-and100bp_2.bed
# bedtools getfasta -fo On_genome_longest_transcript_1stexon.fasta -name -fi ../Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -bed On_genome_longest_transcript_1stexon_2.bed

nano sbatch-On_1stexon-bed.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh On_1stexon-bed.sh

# Run as batch script
sbatch sbatch-On_1stexon-bed.sh


## NOTE:Will is actually doing the above, file is here:
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta

# 2. Use the above to BLAT against the genomes of the other species to find the corresponding start site of its orthologous gene

nano 1stexonBLAT.sh

#!/bin/sh

source blat-35

blat ../MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta mzebra_1stexonBLAT_new.psl
blat ../Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta pnyererei_1stexonBLAT_new.psl
blat ../Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta aburtoni_1stexonBLAT_new.psl
blat ../Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta nbrichardi_1stexonBLAT_new.psl
blat ../Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta oniloticus_1stexonBLAT_new.psl

nano sbatch-1stexonBLAT.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 32000 # memory pool for all cores
#SBATCH -t 20-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh 1stexonBLAT.sh

# Run as batch script
sbatch sbatch-1stexonBLAT.sh

# Also running a BLAST too just in case we require it

nano 1stexonBLAST.sh

#!/bin/sh

ml blast

blastn -db ../MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -outfmt 6 -evalue 1e-1 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 -out mzebra_1stexonBLAST_new.txt -query Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta
blastn -db ../Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -outfmt 6 -evalue 1e-1 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 -out pnyererei_1stexonBLAST_new.txt -query Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta
blastn -db ../Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -outfmt 6 -evalue 1e-1 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 -out aburtoni_1stexonBLAST_new.txt -query Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta
blastn -db ../Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -outfmt 6 -evalue 1e-1 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 -out nbrichardi_1stexonBLAST_new.txt -query Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta
blastn -db ../Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -outfmt 6 -evalue 1e-1 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 -out oniloticus_1stexonBLAST_new.txt -query Oreochromis_niloticus.BROADON1.firstExon100nt_new.fasta

nano sbatch-1stexonBLAST.sh

# Add the following to above script
#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 0-5:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

sh 1stexonBLAST.sh

# Run as batch script
sbatch sbatch-1stexonBLAST.sh

# sort the files so that you select the best match for each gene in each genome
# in psl output - col1+col2+col6 = query size (col11)
head -4 mzebra_1stexonBLAT_new.psl > psl_colheader # save column header of one output file if you want to cat later
for i in *_new.psl ; do sort -k1,1 -rn $i | sort -u -k10,10 | grep -v 'psLayout' | grep -v 'match' > "$(basename "$i" .psl)_2.psl" ; done

# 1. create a BED file for the match in the corresponding species
# prepare a BED (start is from position 0) file for the 1st exon match in each genome (chr, chrStart, chrEnd, name, strand)
for i in *_new_2.psl ; do awk '{print $14, $16-1, $17, $10, $9}' OFS='\t' $i | awk '$2<0 {$2=0} 1' OFS='\t' > "$(basename "$i" BLAT_new_2.psl)_new_annotation.bed" ; done


# 2. map this to a gene in the species either by
  # a. mapping to closest gene (bed intersect) in annotation file (requiring an overlap - this is why we took extra on both ends), or
  # b. simply mapping an ortholog based on OGIDs - NOT the best approach

# Only pull out the mRNA lines of annotation
for i in /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/final_annotation/*.gff3 ; do awk '$3=="mRNA"' $i > "$(basename "$i" .gff3)_mRNA.gff3" ; done


nano bedintersect-transcript.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

ml bedtools/2.25.0
ml GCC
ml zlib

Ab_bed=(aburtoni_1stexon_new_annotation.bed)
Mz_bed=(mzebra_1stexon_new_annotation.bed)
Pn_bed=(pnyererei_1stexon_new_annotation.bed)
Nb_bed=(nbrichardi_1stexon_new_annotation.bed)
On_bed=(oniloticus_1stexon_new_annotation.bed)
Ab_annot=(Astatotilapia_burtoni.BROADAB1_mRNA.gff3)
Mz_annot=(Metriaclima_zebra.BROADMZ1_mRNA.gff3)
Pn_annot=(Pundamilia_nyererei.BROADPN1_mRNA.gff3)
Nb_annot=(Neolamprologus_brichardi.BROADNB1_mRNA.gff3)
On_annot=(Oreochromis_niloticus.BROADON1_mRNA.gff3)

bedtools intersect -a $Ab_bed -b $Ab_annot -wo > Ab_bedintersect_new.txt
bedtools intersect -a $Mz_bed -b $Mz_annot -wo > Mz_bedintersect_new.txt
bedtools intersect -a $Pn_bed -b $Pn_annot -wo > Pn_bedintersect_new.txt
bedtools intersect -a $Nb_bed -b $Nb_annot -wo > Nb_bedintersect_new.txt
bedtools intersect -a $On_bed -b $On_annot -wo > On_bedintersect_new.txt

# Run on uv2k2
qsub -q Test -l select=1:mem=50GB:ncpus=1 bedintersect-transcript.sh

# 3. determine the correct start site based on both the BLAT and the closest gene overlap, preparing a BED file of output

# a. First, map the gene ID in the bedintersect to gene in original annotation and output the 'end' coordinate (will be col4)
for i in /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/final_annotation/*.gff3 ; do awk '$3=="gene"' OFS='\t' $i | awk '{print $9,$1,$4-1,$5,$7}' OFS='\t' | sed 's/ID=//g' > "$(basename "$i" .gff3)_Brawandgff3.bed" ; done
for i in *_bedintersect_new.txt ; do awk '{print $14,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$15}' OFS='\t' $i | sed 's/\;/\t/g' | sed 's/ID=//g' | sed 's/Parent=//g' | awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' > $i.2 ; done
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Astatotilapia_burtoni.BROADAB1_Brawandgff3.bed Ab_bedintersect_new.txt.2 > Ab_bedintersect_new.txt.3
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Metriaclima_zebra.BROADMZ1_Brawandgff3.bed Mz_bedintersect_new.txt.2 > Mz_bedintersect_new.txt.3
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Neolamprologus_brichardi.BROADNB1_Brawandgff3.bed Nb_bedintersect_new.txt.2 > Nb_bedintersect_new.txt.3
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pundamilia_nyererei.BROADPN1_Brawandgff3.bed Pn_bedintersect_new.txt.2 > Pn_bedintersect_new.txt.3
## IMPORTANT: for positive strand, use col4 (start) and col20 (end) as annotation for the longest transcript col 1 (that has been mapped from gff3 gene) - this is only for calling a 5kb promoter upstream and using as annotation of nearest gene
## IMPORTANT: for negative strand, use col19 (start) and col5 (end) as annotation for the longest transcript col 1 (that has been mapped from gff3 gene) - really this is col3 as a start for the gene and col9 as the end
# for i in *_bedintersect_new.txt.3 ; do awk '{if($7 == "+")print $8, $4, $20, $1, $7, $2, $6;else print $8, $19, $5, $1, $7, $2, $6;}' OFS='\t' $i | sort -u -k7,7 > "$(basename "$i" .3).4" ; done
# for i in *_bedintersect_new.txt.4 ; do cut -f7 $i | sed 's/on.cds/on.gene/g' | cut -d. -f1-4 > "$(basename "$i" .4).5" ; done
# for i in *_bedintersect_new.txt ; do paste -d'\t' $i.4 $i.5 | sed 's/bn.gene/nb.gene/g' | sed 's/bn.mrna/nb.mrna/g' > "$(basename "$i" _bedintersect_new.txt)_LongestmRNA_OnBLAT_annot_new.bed" ; done
## Some of the BLAT hits are in the opposite strand to the gene! This means that you need to use the strand from the original annotation for calling start and end coordinates as well as using that strand!
for i in *_bedintersect_new.txt.3 ; do awk '{if($21 == "+")print $8, $4, $20, $1, $21, $2, $6;else print $8, $19, $5, $1, $21, $2, $6;}' OFS='\t' $i | sort -u -k7,7 > "$(basename "$i" .3).4" ; done
for i in *_bedintersect_new.txt.4 ; do cut -f7 $i | sed 's/on.cds/on.gene/g' | cut -d. -f1-4 > "$(basename "$i" .4).5" ; done
for i in *_bedintersect_new.txt ; do paste -d'\t' $i.4 $i.5 | sed 's/bn.gene/nb.gene/g' | sed 's/bn.mrna/nb.mrna/g' > "$(basename "$i" _bedintersect_new.txt)_LongestmRNA_OnBLAT_annot_new.bed" ; done


# determine how many instances in each species where the O. niloticus 1st exon (+/- 100bp) match is overlapping a transcript but the start is at a distance >10 kb away
for i in *_bedintersect_new.txt ; do awk '$16=$2-$9' OFS='\t' $i | awk '$16>=10000' | awk '$5=="+"' OFS='\t' > $i-morethan10kbdist-plustrand ; done
for i in *_bedintersect_new.txt ; do awk '$16=$10-$3' OFS='\t' $i | awk '$16>=10000' | awk '$5=="-"' OFS='\t' > $i-morethan10kbdist-negativestrand ; done
cat Ab_bedintersect_new.txt-morethan10kbdist* > Ab_bedintersect_new.txt-morethan10kbdist-ALL
cat Mz_bedintersect_new.txt-morethan10kbdist* > Mz_bedintersect_new.txt-morethan10kbdist-ALL
cat Pn_bedintersect_new.txt-morethan10kbdist* > Pn_bedintersect_new.txt-morethan10kbdist-ALL
cat Nb_bedintersect_new.txt-morethan10kbdist* > Nb_bedintersect_new.txt-morethan10kbdist-ALL
cat On_bedintersect_new.txt-morethan10kbdist* > On_bedintersect_new.txt-morethan10kbdist-ALL
# ~~~~~~~~~~~~~ ##################### ~~~~~~~~~~~~~~~~ ##########################

## IMPORTANT: WILL NASH TO DO (this is required for funseq annotation) - As of 24/08/17

# From the following files (these are BLAT output files filtered for best top hit per unique O. niloticus transcript)
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_1stexonBLAT_2.psl

# 1. Sum up all small blocks from col19 upto the largest number block is reached > extract that from start coordinate (col19)
# 2. Sum up all small blocks from col19 after the largest number block > extract from end coordinate (col17)
# 3. Use the above annotation and intersects here (/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_bedintersect_new.txt), to
	# a. add the first exon as above coordinates
	# b. go back to original annotation and add remaining exons (not any that overlap i.e. original annotated first exon) to our annotation of the first exon
	# c. these will be the new annotation files

# ~~~~~~~~~~~~~ ##################### ~~~~~~~~~~~~~~~~ ##########################

# 4. filter the output for whether ortholog is correct based on OGIDs.txt5 file, do this by
  # a. mapping cichlid gene to OGID then comparing against the OGID of the original O. niloticus transcript gene that it was indetified from

for i in *_LongestmRNA_OnBLAT_annot_new.bed ; do cut -f4 $i > $i.cichlidgene ; done
for i in *_LongestmRNA_OnBLAT_annot_new.bed ; do cut -f8 $i > $i.Ongene ; done

for i in *_new.bed.Ongene ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_On-rearrange $i > $i.OGID ;
done

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Ab-rearrange Ab_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene > Ab_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Mz-rearrange Mz_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene > Mz_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Pn-rearrange Pn_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene > Pn_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Nb-rearrange Nb_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene > Nb_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene.OGID

paste -d'\t' Ab_LongestmRNA_OnBLAT_annot_new.bed Ab_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene.OGID Ab_LongestmRNA_OnBLAT_annot_new.bed.Ongene.OGID | grep -v NA > Ab_LongestmRNA_OnBLAT_annot_new.bed2
paste -d'\t' Mz_LongestmRNA_OnBLAT_annot_new.bed Mz_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene.OGID Mz_LongestmRNA_OnBLAT_annot_new.bed.Ongene.OGID | grep -v NA > Mz_LongestmRNA_OnBLAT_annot_new.bed2
paste -d'\t' Pn_LongestmRNA_OnBLAT_annot_new.bed Pn_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene.OGID Pn_LongestmRNA_OnBLAT_annot_new.bed.Ongene.OGID | grep -v NA > Pn_LongestmRNA_OnBLAT_annot_new.bed2
paste -d'\t' Nb_LongestmRNA_OnBLAT_annot_new.bed Nb_LongestmRNA_OnBLAT_annot_new.bed.cichlidgene.OGID Nb_LongestmRNA_OnBLAT_annot_new.bed.Ongene.OGID | grep -v NA > Nb_LongestmRNA_OnBLAT_annot_new.bed2

# then match the orthogroups, only outputting lines where they match
for i in *_LongestmRNA_OnBLAT_annot_new.bed2 ; do awk '($10==$12)' OFS='\t' $i > "$(basename "$i" .bed2)_FINAL.bed" ; done

# number of genes per species
10050 Ab_LongestmRNA_OnBLAT_annot_new_FINAL.bed # from 17887
10654 Mz_LongestmRNA_OnBLAT_annot_new_FINAL.bed # from 18116
 8464 Nb_LongestmRNA_OnBLAT_annot_new_FINAL.bed # from 14894
10030 Pn_LongestmRNA_OnBLAT_annot_new_FINAL.bed # from 16867

# how many 1:1 orthologs are present - this will only work for pairs though (match to niloticus, but need not mean that all orthologs are correct)
cat *_LongestmRNA_OnBLAT_annot_new_FINAL.bed | cut -f10 | sort -u | xargs -i grep -wiF {} /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS1to1 > OGIDS1to1_LongestmRNA_OnBLAT_annot_new_FINAL.txt #6593/6844
# those that are absent - extreme cases where not a single ortholog to On is correct
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS1to1_LongestmRNA_OnBLAT_annot_new_FINAL.txt /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS1to1 | grep NA > missingOGIDS1to1_LongestmRNA_OnBLAT_annot_new_FINAL.txt
# to reconstruct whether the orthology is ok we need to match separately look for the non-matches then sort -u
for i in *_LongestmRNA_OnBLAT_annot_new.bed2 ; do awk '($10!=$12)' OFS='\t' $i | cut -f10 | sort -u > anymissingOGIDS1to1_LongestmRNA_OnBLAT_annot_new_FINAL.txt ; done # 1030


# how many of the >10kb overlaps match 1:1 orthologs
cut -f14 Ab_bedintersect_new.txt-morethan10kbdist-ALL | sed 's/ID=.*;Parent=//g' | sort -u | xargs -i grep -wiF {} /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS1to1 | wc -l # 0
cut -f14 Mz_bedintersect_new.txt-morethan10kbdist-ALL | sed 's/ID=.*;Parent=//g' | sort -u | xargs -i grep -wiF {} /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS1to1 | wc -l # 0
cut -f14 Pn_bedintersect_new.txt-morethan10kbdist-ALL | sed 's/ID=.*;Parent=//g' | sort -u | xargs -i grep -wiF {} /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS1to1 | wc -l # 0
cut -f14 Nb_bedintersect_new.txt-morethan10kbdist-ALL | sed 's/ID=.*;Parent=//g' | sort -u | xargs -i grep -wiF {} /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS1to1 | wc -l # 0

# since none of the problematic areas overlap 1:1 OGs, we can proceed!

### HOWEVER, THE ABOVE IS NOT FINAL SINCE IT ONLY CONTAINS ORTHOLOGS in ORTHOGROUPING AND NOT ALL GENES
# Also the files before this are only for BLAT hits and therefore not ALL genes in the genome (which we require for boundaries)
# Need to use and filter for following:
  # 1. Take BLAT + OGID overlap (incl 1:1) as base - take this
  # 2. BLAT + OGID overlap (incl 1:1) map to BLAT duplicates removed only - take the non-overlap and cat with OGID overlap
  # 3. Above overlap with Brawand annotations - take the non-overlap
  # 4. cat the above (no duplicates should be here)

# 1) Pull out gene information for the longest genes and create a BED output (0-based)
for i in /tgac/workarea/group-vh/cichlids/data/Broad/all/Data/Annotation/final_annotation/*.gff3 ; do awk '$3=="gene"' OFS='\t' $i | awk '{print $1,$4-1,$5,$7,$9}' OFS='\t' | sed 's/ID=//g' > "$(basename "$i" .gff3)_Brawand.bed" ; done

# 2) Cat BLAT1to1 files with the Brawand annotation to complete for other genes (we don't want to use other BLAT hits as they could be dubious but will rectify any odd BLAT hits by scanning in 10kb windows)
# also repalce any negative values with 0
for i in *_LongestmRNA_OnBLAT_annot_new_FINAL.bed ; do awk '{print $1,$2,$3,$5,$4,$4}' OFS='\t' $i > $i.rearrange ; done # BLAT 1:1 only
for i in *_Brawand.bed ; do awk '{print $1,$2,$3,$4,$5,$5}' OFS='\t' $i | sed 's/bn.gene/nb.gene/g' > $i.2 ; done # Original annotation
cat Mz_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Metriaclima_zebra.BROADMZ1_Brawand.bed.2 | sort -u -k5,5 | sort -k1,1 -k2,2n | awk '$2<0 {$2=0} 1' OFS='\t' > Mz_GeneAnnotation_11092017.bed
cat Pn_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Pundamilia_nyererei.BROADPN1_Brawand.bed.2 | sort -u -k5,5 | sort -k1,1 -k2,2n | awk '$2<0 {$2=0} 1' OFS='\t' > Pn_GeneAnnotation_11092017.bed
cat Ab_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Astatotilapia_burtoni.BROADAB1_Brawand.bed.2 | sort -u -k5,5 | sort -k1,1 -k2,2n | awk '$2<0 {$2=0} 1' OFS='\t' > Ab_GeneAnnotation_11092017.bed
cat Nb_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Neolamprologus_brichardi.BROADNB1_Brawand.bed.2 | sort -u -k5,5 | sort -k1,1 -k2,2n | awk '$2<0 {$2=0} 1' OFS='\t' > Nb_GeneAnnotation_11092017.bed
# ABOVE ARE THE FINAL GENE ANNOTATION FILES


# Rearrange ordering of three files for awk matching, then awk matching
# for i in *_LongestmRNA_OnBLAT_annot_new_FINAL.bed ; do awk '{print $4,$1,$2,$3,$5}' OFS='\t' $i > $i.rearrange ; done # BLAT 1:1 only
# for i in *_LongestmRNA_OnBLAT_annot_new.bed ; do awk '{print $4,$1,$2,$3,$5}' OFS='\t' $i > $i.rearrange ; done # BLAT + BLAT 1:1
# for i in *_Brawand.bed ; do awk '{print $5,$1,$2,$3,$4}' OFS='\t' $i | sed 's/bn.gene/nb.gene/g' > $i.2 ; done # Original annotation

# the *_LongestmRNA_OnBLAT_annot_new.bed.rearrange files will have duplicate entries so we need to ensure that the non-duplicate entries that do not overlap those with ortholog mapping are retained
# for i in *_LongestmRNA_OnBLAT_annot_new.bed.rearrange ; do awk 'NR==FNR{s[$1]++;next} (s[$1]<2)' $i $i > $i.dupremoved ; done # remove duplicates from BLAT hits

# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Ab_LongestmRNA_OnBLAT_annot_new.bed.rearrange.dupremoved | grep NA | cut -f1-5 > Ab_BLAT_no1to1.bed
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Mz_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Mz_LongestmRNA_OnBLAT_annot_new.bed.rearrange.dupremoved | grep NA | cut -f1-5 > Mz_BLAT_no1to1.bed
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pn_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Pn_LongestmRNA_OnBLAT_annot_new.bed.rearrange.dupremoved | grep NA | cut -f1-5 > Pn_BLAT_no1to1.bed
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Nb_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Nb_LongestmRNA_OnBLAT_annot_new.bed.rearrange.dupremoved | grep NA | cut -f1-5 > Nb_BLAT_no1to1.bed

# Cat the above and get the non-overlap from the BLAT hits
# cat Ab_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Ab_BLAT_no1to1.bed > Ab_BLAT_plus1to1.bed
# cat Mz_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Mz_BLAT_no1to1.bed > Mz_BLAT_plus1to1.bed
# cat Pn_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Pn_BLAT_no1to1.bed > Pn_BLAT_plus1to1.bed
# cat Nb_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange Nb_BLAT_no1to1.bed > Nb_BLAT_plus1to1.bed
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_BLAT_plus1to1.bed Astatotilapia_burtoni.BROADAB1_Brawand.bed.2 | grep NA | awk '{print $2,$3,$4,$5,$1,$1}' OFS='\t' > Ab_Brawand.bed.3
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Mz_BLAT_plus1to1.bed Metriaclima_zebra.BROADMZ1_Brawand.bed.2 | grep NA | awk '{print $2,$3,$4,$5,$1,$1}' OFS='\t' > Mz_Brawand.bed.3
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pn_BLAT_plus1to1.bed Pundamilia_nyererei.BROADPN1_Brawand.bed.2 | grep NA | awk '{print $2,$3,$4,$5,$1,$1}' OFS='\t' > Pn_Brawand.bed.3
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Nb_BLAT_plus1to1.bed Neolamprologus_brichardi.BROADNB1_Brawand.bed.2 | grep NA | awk '{print $2,$3,$4,$5,$1,$1}' OFS='\t' > Nb_Brawand.bed.3

# We have dropout of some genes based on coherence of final gene numbers with original gff3 - this is because genes have been parsed out by Will that overlap other genes or are inside them
# awk '{print $5,$1,$2$3,$4,$6}' OFS='\t' Ab_LongestmRNA_OnBLAT_annot_new.5kb_promoters.stranded.GENEBED_all.bed > Ab_LongestmRNA_OnBLAT_annot_new.5kb_promoters.stranded.GENEBED_all.bed.DROPOUT1
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_LongestmRNA_OnBLAT_annot_new.5kb_promoters.stranded.GENEBED_all.bed.DROPOUT1 Ab_LongestmRNA_OnBLAT_annot_new.bed.rearrange | grep NA

# Good to go ahead so create a complete file of the final annotations by concatenating all non-overlap files incl longest mRNA annotation files (there should be no duplicate entries)
# for i in *_BLAT_plus1to1.bed ; do awk '{print $2,$3,$4,$5,$1,$1}' OFS='\t' $i > $i.2 ; done
# cat Ab_BLAT_plus1to1.bed.2 Ab_Brawand.bed.3 | sort -k1,1 -k2,2n > Ab_GeneAnnotation_11092017.bed
# cat Mz_BLAT_plus1to1.bed.2 Mz_Brawand.bed.3 | sort -k1,1 -k2,2n > Mz_GeneAnnotation_11092017.bed
# cat Nb_BLAT_plus1to1.bed.2 Nb_Brawand.bed.3 | sort -k1,1 -k2,2n > Nb_GeneAnnotation_11092017.bed
# cat Pn_BLAT_plus1to1.bed.2 Pn_Brawand.bed.3 | sort -k1,1 -k2,2n > Pn_GeneAnnotation_11092017.bed
# ABOVE ARE THE FINAL GENE ANNOTATION FILES

# Will uses the above annotations for below to CREATE FILES THAT CONTAIN promoter regions, removing any overlapping gene bodies for annotation, these are:
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017.5kb_promoters.stranded.GENEBED_all.bed
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017.5kb_promoters.stranded.GENEBED_all.fasta

## For Will, create a file with all species BLAT hit (no duplicates) regions so that Will (when aligning) can dump alignments that are non-BLAT hits into a separate folder
# Sort out gene name and remove transcript ID
grep -v on.gene Oreochromis_niloticus.BROADON1.longestCDS_new.bed | awk '{gsub(/\.[0-9]$/,"",$5); print $0}' OFS='\t' | awk '{gsub(/\.[0-9]$/,"",$6); print $0}' OFS='\t' | sed 's/on.mrna/on.gene/g' > Oreochromis_niloticus.BROADON1.longestCDS_new.bed.2
grep on.gene Oreochromis_niloticus.BROADON1.longestCDS_new.bed > Oreochromis_niloticus.BROADON1.longestCDS_new.bed.3
cat Oreochromis_niloticus.BROADON1.longestCDS_new.bed.2 Oreochromis_niloticus.BROADON1.longestCDS_new.bed.3 | sort -k1,1 -k2,2n > Oreochromis_niloticus.BROADON1.longestCDS_new.bed.4
rm Oreochromis_niloticus.BROADON1.longestCDS_new.bed Oreochromis_niloticus.BROADON1.longestCDS_new.bed.2 Oreochromis_niloticus.BROADON1.longestCDS_new.bed.3
mv Oreochromis_niloticus.BROADON1.longestCDS_new.bed.4 Oreochromis_niloticus.BROADON1.longestCDS_new.bed
awk '{print $5,$1,$2,$3,$4}' OFS='\t' Oreochromis_niloticus.BROADON1.longestCDS_new.bed > Oreochromis_niloticus.BROADON1.longestCDS_new.bed.rearrange # this is NOT needed below!!!
for i in *_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange ; do awk '{print $5,$1,$2,$3,$4}' OFS='\t' $i > $i.2 ; done
cat *_LongestmRNA_OnBLAT_annot_new_FINAL.bed.rearrange.2 Oreochromis_niloticus.BROADON1.longestCDS_new.bed > all_BLAT_plus1to1_plusOn.bed ### ADDED O. NILOTICUS - Will uses this to put alignments from non-BLAT hits into separate folder
# new file where Will added the additional annotations that were taken from 10kb windows in BLAT adhering to orthogrouping > all_BLAT_plus1to1_plusOn_plusSalvaged.bed
# Also create a file to check how many OGs are fully represented by BLAT 1:1 hits

awk '{print $5,$1,$2,$3,$4,$6}' OFS='\t' all_BLAT_plus1to1_plusOn_plusSalvaged.bed > all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Mz all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange > all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange.Mz
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Pn all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange > all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange.Pn
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Ab all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange > all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange.Ab
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Nb all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange > all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange.Nb
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-On all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange > all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange.On

cat all_BLAT_plus1to1_plusOn_plusSalvaged.bed.rearrange.* | grep -v NA | cut -f7 > allspecies_BLAT1to1-only-1to1OGIDs
awk 'NR==FNR{s[$1]++;next} (s[$1]==1)' allspecies_BLAT1to1-only-1to1OGIDs allspecies_BLAT1to1-only-1to1OGIDs | sort -u | wc -l #251 BLAT1to1OGIDs only appear once
awk 'NR==FNR{s[$1]++;next} (s[$1]==2)' allspecies_BLAT1to1-only-1to1OGIDs allspecies_BLAT1to1-only-1to1OGIDs | sort -u | wc -l #156 BLAT1to1OGIDs only appear twice
awk 'NR==FNR{s[$1]++;next} (s[$1]==3)' allspecies_BLAT1to1-only-1to1OGIDs allspecies_BLAT1to1-only-1to1OGIDs | sort -u | wc -l #252 BLAT1to1OGIDs only appear thrice
awk 'NR==FNR{s[$1]++;next} (s[$1]==4)' allspecies_BLAT1to1-only-1to1OGIDs allspecies_BLAT1to1-only-1to1OGIDs | sort -u | wc -l #1058 BLAT1to1OGIDs only appear four times
awk 'NR==FNR{s[$1]++;next} (s[$1]==5)' allspecies_BLAT1to1-only-1to1OGIDs allspecies_BLAT1to1-only-1to1OGIDs | sort -u | wc -l #5127 BLAT1to1OGIDs appear the required five times





##### Several genes are not providing good alignments as the O. niloticus 1st exon BLAT did not intersect with the correct gene e.g. is on.gene.LG8-24.649 intersected nb.gene.s27.143 instead of on.gene.LG8-24.650 (only ~1.2kb off intersecting!)
# Strategy to fix these:
# 1. Get the on.gene in OGIDs file that have a 1:1 ortholog in corresponding species that are not present in the bed intersect
awk '$2!="NULL" && $6!="NULL"' OFS='\t' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 | awk '{print $6,$2,$1}' OFS='\t' > OGID-Mz-On
awk '$3!="NULL" && $6!="NULL"' OFS='\t' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 | awk '{print $6,$3,$1}' OFS='\t' > OGID-Pn-On
awk '$4!="NULL" && $6!="NULL"' OFS='\t' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 | awk '{print $6,$4,$1}' OFS='\t' > OGID-Ab-On
awk '$5!="NULL" && $6!="NULL"' OFS='\t' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/OGIDS.txt5 | awk '{print $6,$5,$1}' OFS='\t' > OGID-Nb-On

for i in *_bedintersect_new.txt ; do awk '{print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' $i > $i.rearrangeforOGfix ; done
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_bedintersect_new.txt.rearrangeforOGfix OGID-Ab-On | grep -wiF NA | cut -f1-3 > OGID-Ab-On_noBEDintersect
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Mz_bedintersect_new.txt.rearrangeforOGfix OGID-Mz-On | grep -wiF NA | cut -f1-3 > OGID-Mz-On_noBEDintersect
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pn_bedintersect_new.txt.rearrangeforOGfix OGID-Pn-On | grep -wiF NA | cut -f1-3 > OGID-Pn-On_noBEDintersect
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Nb_bedintersect_new.txt.rearrangeforOGfix OGID-Nb-On | grep -wiF NA | cut -f1-3 > OGID-Nb-On_noBEDintersect

# Shared the above with Will who is now preparing a script that iterates through the BLAT to see if there are any hits within around 10kb (that correspond to the correct ortholog), then take the correct orthologous BLAT hit
# Then, feed in the gff3 file to get the rest of the gene coordinates
# Once we have the correct regions, amend all incorrect gene annotation lines with the correct region (below) then re-call promoters and align
cat OGID-Mz-On_noBEDintersect_correctBLAThits.bed Mz_GeneAnnotation_11092017.bed | sort -u -k5,5 | sort -k1,1 -k2,2n > Mz_GeneAnnotation_11092017_FINALcorrected.bed
cat OGID-Pn-On_noBEDintersect_correctBLAThits.bed Pn_GeneAnnotation_11092017.bed | sort -u -k5,5 | sort -k1,1 -k2,2n > Pn_GeneAnnotation_11092017_FINALcorrected.bed
cat OGID-Ab-On_noBEDintersect_correctBLAThits.bed Ab_GeneAnnotation_11092017.bed | sort -u -k5,5 | sort -k1,1 -k2,2n > Ab_GeneAnnotation_11092017_FINALcorrected.bed
cat OGID-Nb-On_noBEDintersect_correctBLAThits.bed Nb_GeneAnnotation_11092017.bed | sort -u -k5,5 | sort -k1,1 -k2,2n > Nb_GeneAnnotation_11092017_FINALcorrected.bed


###### AS OF 13/09/2017, THESE ARE THE FINAL PROMOTER AND GENE ANNOTATION FILES TO USE
## FINAL LONGEST PROTEIN-CODING TRANSCRIPT ANNOTATIONS ARE:
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.bed
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.bed

## Longest cds sequences
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.fasta

# Oreochromis 1st exon ANNOTATIONS
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.firstExon100nt_new.bed

## FINAL PROMOTER ANNOTATIONS ARE (HOWEVER, THESE DO NOT TAKE SCAFFOLD LENGTHS INTO CONSIDERATION, SO SOME MAY BE OUT OF BOUNDS!!)
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.bed

# Check how many 1to1 orthologous groups have promoter annotations
cut -f2 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed | cut -f5 > Mz_prom_annotation_present #6711
cut -f3 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed | cut -f5 > Pn_prom_annotation_present #6758
cut -f4 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed | cut -f5 > Ab_prom_annotation_present #6711
cut -f5 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.bed | cut -f5 > Nb_prom_annotation_present #6768
cut -f6 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.bed | cut -f5 > On_prom_annotation_present #6784
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Mz_prom_annotation_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Mz | grep -v NA | cut -f2 > Mz_prom_annotation_present.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pn_prom_annotation_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Pn | grep -v NA | cut -f2 > Pn_prom_annotation_present.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_prom_annotation_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Ab | grep -v NA | cut -f2 > Ab_prom_annotation_present.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Nb_prom_annotation_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Nb | grep -v NA | cut -f2 > Nb_prom_annotation_present.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' On_prom_annotation_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-On | grep -v NA | cut -f2 > On_prom_annotation_present.OGID
cat *_prom_annotation_present.OGID | sort -u | wc -l # A total of 6842/6844 orthogroups are (in at least one species) represented in promoter annotations (but the whole orthogroup may not be represented, as in annotations for all species in an orthogroup)
cat *_prom_annotation_present.OGID > all_prom_annotation_present.OGID
awk 'NR==FNR{s[$1]++;next} (s[$1]==1)' all_prom_annotation_present.OGID all_prom_annotation_present.OGID | sort -u | wc -l #13 OGIDs only appear once
awk 'NR==FNR{s[$1]++;next} (s[$1]==2)' all_prom_annotation_present.OGID all_prom_annotation_present.OGID | sort -u | wc -l #18 OGIDs only appear twice
awk 'NR==FNR{s[$1]++;next} (s[$1]==3)' all_prom_annotation_present.OGID all_prom_annotation_present.OGID | sort -u | wc -l #58 OGIDs only appear thrice
awk 'NR==FNR{s[$1]++;next} (s[$1]==4)' all_prom_annotation_present.OGID all_prom_annotation_present.OGID | sort -u | wc -l #287 OGIDs only appear four times
awk 'NR==FNR{s[$1]++;next} (s[$1]==5)' all_prom_annotation_present.OGID all_prom_annotation_present.OGID | sort -u | wc -l #6466 OGIDs appear the required five times
# a total of 376/6844 1to1OGIDs are not fully represented in annotations (as in annotation present for all five species)

# FINAL PROMOTER FASTA SEQUENCES ARE:
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta
# create a tar ball to share with Sushmita
tar -zcvf FINALcorrected.5kb_promoters.tar.gz *_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta

# Check how many 1to1 orthologous groups have promoter sequences
for i in *_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta ; do grep '>' $i | sed 's/>//g' > "$(basename "$i" GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta)_PromoterGenesinFasta" ; done
grep '>' Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta | sed 's/>//g' > On__PromoterGenesinFasta # just pull out the gene IDs to analyse
cut -f2 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} Mz__PromoterGenesinFasta > Mz_prom_fasta_present #6711
cut -f3 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} Pn__PromoterGenesinFasta > Pn_prom_fasta_present #6754
cut -f4 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} Ab__PromoterGenesinFasta > Ab_prom_fasta_present #6708
cut -f5 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} Nb__PromoterGenesinFasta > Nb_prom_fasta_present #6764
cut -f6 /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | xargs -i grep -wiF {} On__PromoterGenesinFasta > On_prom_fasta_present #6783
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Mz_prom_fasta_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Mz | grep -v NA | cut -f2 > Mz_prom_fasta_present.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pn_prom_fasta_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Pn | grep -v NA | cut -f2 > Pn_prom_fasta_present.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_prom_fasta_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Ab | grep -v NA | cut -f2 > Ab_prom_fasta_present.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Nb_prom_fasta_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-Nb | grep -v NA | cut -f2 > Nb_prom_fasta_present.OGID
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' On_prom_fasta_present /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1-On | grep -v NA | cut -f2 > On_prom_fasta_present.OGID
cat *_prom_fasta_present.OGID | sort -u | wc -l # A total of 6842/6844 orthogroups are (in at least one species) represented in promoter fastas (but the whole orthogroup may not be represented, as in fastas for all species in an orthogroup)
cat *_prom_fasta_present.OGID > all_prom_fasta_present.OGID
awk 'NR==FNR{s[$1]++;next} (s[$1]==1)' all_prom_fasta_present.OGID all_prom_fasta_present.OGID | sort -u | wc -l #13 OGIDs only appear once
awk 'NR==FNR{s[$1]++;next} (s[$1]==2)' all_prom_fasta_present.OGID all_prom_fasta_present.OGID | sort -u | wc -l #18 OGIDs only appear twice
awk 'NR==FNR{s[$1]++;next} (s[$1]==3)' all_prom_fasta_present.OGID all_prom_fasta_present.OGID | sort -u | wc -l #58 OGIDs only appear thrice
awk 'NR==FNR{s[$1]++;next} (s[$1]==4)' all_prom_fasta_present.OGID all_prom_fasta_present.OGID | sort -u | wc -l #299 OGIDs only appear four times
awk 'NR==FNR{s[$1]++;next} (s[$1]==5)' all_prom_fasta_present.OGID all_prom_fasta_present.OGID | sort -u | wc -l #6454 OGIDs appear the required five times
# a total of 388/6844 1to1OGIDs are not fully represented in annotations (as in annotation present for all five species)

# Promoter alignments are here
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/alignments/
# separated by those with a BLAT1to1OGID
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/alignments/*.fa
# and those without, around 20% are reliable
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/alignments/no_BLAT_hit/*.fa




# 22/09/17
# Wilfried shared the promoter only outputs - they need parsing

# Wilfried ran it on both folders, first with BLAT that are 1to1 OGIDs
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/alignments/rates/output/*.out #4926 files (from 5421) - here we used a filter of 100nt minimum alignment and 10% aligned (as opposed to 50% that results in only 2992 files)
# Then, those that are not - no BLAT hit and used the original annotation
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/alignments/no_BLAT_hit/rates/output/*.out
# parse output files - 1:On, 2:Anc3, 3:Nb, 4:Anc2, 5:Ab, 6:Anc1, 7:Mz, 8:Pn; then added whole tree onto end


### For those with a BLAT 1:1 OGID
cd /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/alignments/rates/output/
for i in *.out ; do
  awk 'BEGIN {OFS="\t"} {print $0,FILENAME}' $i | grep -A 4 lnL | tail -3 | awk 'FNR!=1{print l}{l=$0};END{ORS="";print l}' ORS='\t' | cut -f1,2,5 | sed 's/_aln_general.out//g' | sed 's/tree length*.=//g' | sed 's/^ *//' | tr ' ' \\t | awk '{print $15,$1,$2,$3,$4,$5,$6,$7,$8,$16}' OFS='\t' | grep -v 'nan' >> promoter_BLAT1to1OGIDcollated.out ;
done # from 4926, this leaves 2892 promoter values

# then remove rows where any value >1
# NOTE: owing to the short evolutionary distances, we would not expect branch lengths of >1 - therefore, we will remove any where branch lengths are >1 (this is the same for 4fold sites)
awk '$2<"1" && $3<"1" && $4<"1" && $5<"1" && $6<"1" && $7<"1" && $8<"1" && $9<"1" && $10<"1"' promoter_BLAT1to1OGIDcollated.out > promoter_BLAT1to1OGIDcollated.out2 # leaves 4875 promoter evol rate rows

printf 'OGID\tOn\tAnc3\tNb\tAnc2\tAb\tAnc1\tMz\tNb\tWholeTree\n' > promoter_colheaders ; cat promoter_colheaders promoter_BLAT1to1OGIDcollated.out2 > 220917_promoter_ratesOGIDs.txt # add colheaders

# separate by branch and node (6..7 - Anc3; 7..8 - Anc2; 8..9 - Anc1)
cut -f1,2 220917_promoter_ratesOGIDs.txt > 220917_promoter_ratesOGIDs.txt-On
cut -f1,3 220917_promoter_ratesOGIDs.txt > 220917_promoter_ratesOGIDs.txt-Anc3
cut -f1,4 220917_promoter_ratesOGIDs.txt > 220917_promoter_ratesOGIDs.txt-Nb
cut -f1,5 220917_promoter_ratesOGIDs.txt > 220917_promoter_ratesOGIDs.txt-Anc2
cut -f1,6 220917_promoter_ratesOGIDs.txt > 220917_promoter_ratesOGIDs.txt-Ab
cut -f1,7 220917_promoter_ratesOGIDs.txt > 220917_promoter_ratesOGIDs.txt-Anc1
cut -f1,8 220917_promoter_ratesOGIDs.txt > 220917_promoter_ratesOGIDs.txt-Mz
cut -f1,9 220917_promoter_ratesOGIDs.txt > 220917_promoter_ratesOGIDs.txt-Pn
cut -f1,10 220917_promoter_ratesOGIDs.txt > 220917_promoter_ratesOGIDs.txt-WholeTree

# combine with 4fold rates that can be found here: /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/
# create anc evol rates for the 4fold (col2) that Will recently provided (email Saturday, 5 August 2017 at 11:55 "EvolRate") - note these contain promoter values (col3) that we will NOT use
cat /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/no_switch_ortho_MzPnvsAbNbOn_new.out /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/switch_ortho_MzPnvsAbNbOn_new.out > /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/Anc1_MzPn_evolRate.out
cat /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/no_switch_ortho_MzPnAbvsNbOn_new.out /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/switch_ortho_MzPnAbvsNbOn_new.out > /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/Anc2_MzPnAb_evolRate.out
cat /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/no_switch_ortho_MzPnAbNbvsOn_new.out /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/switch_ortho_MzPnAbNbvsOn_new.out > /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/Anc3_MzPnAbNb_evolRate.out

mkdir /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0917_USETHIS-BLAT
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0917_USETHIS-BLAT
mkdir BLAT_promrate
for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/alignments/rates/output/220917_* ; do cp $i BLAT_promrate/ ; done # copy files that you need to here

# Wilfried shared the 4fold only outputs - they need parsing (as opposed to what Will provided which will remove OGID's according to failed promoter evol rate on that run)
# files here /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/4_fold/output/*.out #6410 files
# parse output files - 1:On, 2:Anc3, 3:Nb, 4:Anc2, 5:Ab, 6:Anc1, 7:Mz, 8:Pn; then added whole tree onto end
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/4_fold/output/
for i in *.out ; do
  awk 'BEGIN {OFS="\t"} {print $0,FILENAME}' $i | grep -A 4 lnL | tail -3 | awk 'FNR!=1{print l}{l=$0};END{ORS="";print l}' ORS='\t' | cut -f1,2,5 | sed 's/_four.*out//g' | sed 's/tree length*.=//g' | sed 's/^ *//' | tr ' ' \\t | awk '{print $15,$1,$2,$3,$4,$5,$6,$7,$8,$16}' OFS='\t' | grep -v 'nan' >> 4fold_collated.out ;
done # this leaves 6306 4fold values

printf 'OGID\tOn\tAnc3\tNb\tAnc2\tAb\tAnc1\tMz\tNb\tWholeTree\n' > 4fold_colheaders ; cat 4fold_colheaders 4fold_collated.out > 4fold_collated.out2 # add colheaders


cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0917_USETHIS-BLAT
mkdir 0817_4foldrate
cd 0817_4foldrate
cp /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/4_fold/output/4fold_collated.out . # cols are 1:OGID, 2:On, 3:Anc3, 4:Nb, 5:Anc2, 6:Ab, 7:Anc1, 8:Mz, 9:Pn, 10:WholeTree
cut -f1,2 4fold_collated.out > 4fold_collated.out-On
cut -f1,3 4fold_collated.out > 4fold_collated.out-Anc3
cut -f1,4 4fold_collated.out > 4fold_collated.out-Nb
cut -f1,5 4fold_collated.out > 4fold_collated.out-Anc2
cut -f1,6 4fold_collated.out > 4fold_collated.out-Ab
cut -f1,7 4fold_collated.out > 4fold_collated.out-Anc1
cut -f1,8 4fold_collated.out > 4fold_collated.out-Mz
cut -f1,9 4fold_collated.out > 4fold_collated.out-Pn
cut -f1,10 4fold_collated.out > 4fold_collated.out-WholeTree


# now join the 4fold with the new promoters (BLAT approach version) - cols will be ($1 -OGID; $2 -4fold; $3 -promoter from MAF extract)
cd ../
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoter_ratesOGIDs.txt-On 0817_4foldrate/4fold_collated.out-On | grep -v NA > On_4fold_Prom-BLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoter_ratesOGIDs.txt-Nb 0817_4foldrate/4fold_collated.out-Nb | grep -v NA > Nb_4fold_Prom-BLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoter_ratesOGIDs.txt-Ab 0817_4foldrate/4fold_collated.out-Ab | grep -v NA > Ab_4fold_Prom-BLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoter_ratesOGIDs.txt-Mz 0817_4foldrate/4fold_collated.out-Mz | grep -v NA > Mz_4fold_Prom-BLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoter_ratesOGIDs.txt-Pn 0817_4foldrate/4fold_collated.out-Pn | grep -v NA > Pn_4fold_Prom-BLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoter_ratesOGIDs.txt-Anc3 0817_4foldrate/4fold_collated.out-Anc3 | grep -v NA > Anc3_4fold_Prom-BLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoter_ratesOGIDs.txt-Anc2 0817_4foldrate/4fold_collated.out-Anc2 | grep -v NA > Anc2_4fold_Prom-BLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoter_ratesOGIDs.txt-Anc1 0817_4foldrate/4fold_collated.out-Anc1 | grep -v NA > Anc1_4fold_Prom-BLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoter_ratesOGIDs.txt-WholeTree 0817_4foldrate/4fold_collated.out-WholeTree | grep -v NA > WholeTree_4fold_Prom-BLAT.evolrate.txt
# After joing 4fold and promoter rates, there 4588 orthogroups however, there still could be saturated 4fold values

# check for 4fold values > 1
for i in *4fold_Prom-BLAT.evolrate.txt ; do awk '$2>"1"' $i > $i.remove ; done
cat *4fold_Prom-BLAT.evolrate.txt.remove | cut -f1 | sort -u > OGIDstoremove # 26 in total to remove from each
# remove saturated 4fold orthogroups from all branches and nodes
for i in *4fold_Prom-BLAT.evolrate.txt ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDstoremove $i | grep 'NA' | cut -f1-3 > "$(basename $i .txt).txta" ;
done
## ## ## ~ IN TOTAL, THERE ARE 4562 BLAT 1:1 orthogroups left for evolutionary rate analysis ~ ## ## ## < Use these to cat with no BLAT hit below



########## For those with no BLAT hit ############

cd /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/alignments/no_BLAT_hit/rates/output/
for i in *.out ; do
  awk 'BEGIN {OFS="\t"} {print $0,FILENAME}' $i | grep -A 4 lnL | tail -3 | awk 'FNR!=1{print l}{l=$0};END{ORS="";print l}' ORS='\t' | cut -f1,2,5 | sed 's/_aln_general.out//g' | sed 's/tree length*.=//g' | sed 's/^ *//' | tr ' ' \\t | awk '{print $15,$1,$2,$3,$4,$5,$6,$7,$8,$16}' OFS='\t' | grep -v 'nan' >> promoter_noBLATcollated.out ;
done # from 1033, this leaves 236 promoter values

# then remove rows where any value >1
# NOTE: owing to the short evolutionary distances, we would not expect branch lengths of >1 - therefore, we will remove any where branch lengths are >1 (this is the same for 4fold sites)
awk '$2<"1" && $3<"1" && $4<"1" && $5<"1" && $6<"1" && $7<"1" && $8<"1" && $9<"1" && $10<"1"' promoter_noBLATcollated.out > promoter_noBLATcollated.out2 # leaves 62 promoter evol rate rows

printf 'OGID\tOn\tAnc3\tNb\tAnc2\tAb\tAnc1\tMz\tNb\tWholeTree\n' > promoter_colheaders ; cat promoter_colheaders promoter_noBLATcollated.out2 > 220917_promoternoBLAT_ratesOGIDs.txt # add colheaders

# separate by branch and node (6..7 - Anc3; 7..8 - Anc2; 8..9 - Anc1)
cut -f1,2 220917_promoternoBLAT_ratesOGIDs.txt > 220917_promoternoBLAT_ratesOGIDs.txt-On
cut -f1,3 220917_promoternoBLAT_ratesOGIDs.txt > 220917_promoternoBLAT_ratesOGIDs.txt-Anc3
cut -f1,4 220917_promoternoBLAT_ratesOGIDs.txt > 220917_promoternoBLAT_ratesOGIDs.txt-Nb
cut -f1,5 220917_promoternoBLAT_ratesOGIDs.txt > 220917_promoternoBLAT_ratesOGIDs.txt-Anc2
cut -f1,6 220917_promoternoBLAT_ratesOGIDs.txt > 220917_promoternoBLAT_ratesOGIDs.txt-Ab
cut -f1,7 220917_promoternoBLAT_ratesOGIDs.txt > 220917_promoternoBLAT_ratesOGIDs.txt-Anc1
cut -f1,8 220917_promoternoBLAT_ratesOGIDs.txt > 220917_promoternoBLAT_ratesOGIDs.txt-Mz
cut -f1,9 220917_promoternoBLAT_ratesOGIDs.txt > 220917_promoternoBLAT_ratesOGIDs.txt-Pn
cut -f1,10 220917_promoternoBLAT_ratesOGIDs.txt > 220917_promoternoBLAT_ratesOGIDs.txt-WholeTree

# combine with 4fold rates that can be found here: /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0817/
# all files already created for merging (above, BLAT hit entries) > /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0917_USETHIS-BLAT/0817_4foldrate/*
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0917_USETHIS-BLAT
for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/alignments/no_BLAT_hit/rates/output/220917_promoternoBLAT* ; do cp $i BLAT_promrate/ ; done # copy files that you need to here


# now join the 4fold with the no BLAT promoter (BLAT approach versions) - cols will be ($1 -OGID; $2 -4fold; $3 -promoter from MAF extract)
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoternoBLAT_ratesOGIDs.txt-On 0817_4foldrate/4fold_collated.out-On | grep -v NA > On_4fold_Prom-noBLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoternoBLAT_ratesOGIDs.txt-Nb 0817_4foldrate/4fold_collated.out-Nb | grep -v NA > Nb_4fold_Prom-noBLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoternoBLAT_ratesOGIDs.txt-Ab 0817_4foldrate/4fold_collated.out-Ab | grep -v NA > Ab_4fold_Prom-noBLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoternoBLAT_ratesOGIDs.txt-Mz 0817_4foldrate/4fold_collated.out-Mz | grep -v NA > Mz_4fold_Prom-noBLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoternoBLAT_ratesOGIDs.txt-Pn 0817_4foldrate/4fold_collated.out-Pn | grep -v NA > Pn_4fold_Prom-noBLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoternoBLAT_ratesOGIDs.txt-Anc3 0817_4foldrate/4fold_collated.out-Anc3 | grep -v NA > Anc3_4fold_Prom-noBLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoternoBLAT_ratesOGIDs.txt-Anc2 0817_4foldrate/4fold_collated.out-Anc2 | grep -v NA > Anc2_4fold_Prom-noBLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoternoBLAT_ratesOGIDs.txt-Anc1 0817_4foldrate/4fold_collated.out-Anc1 | grep -v NA > Anc1_4fold_Prom-noBLAT.evolrate.txt
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' BLAT_promrate/220917_promoternoBLAT_ratesOGIDs.txt-WholeTree 0817_4foldrate/4fold_collated.out-WholeTree | grep -v NA > WholeTree_4fold_Prom-noBLAT.evolrate.txt
# After joing 4fold and no BLAT promoter rates, there 60 orthogroups however, there still could be saturated 4fold values

# check for 4fold values > 1
for i in *Prom-noBLAT.evolrate.txt ; do awk '$2>"1"' $i > $i.remove ; done
cat *Prom-noBLAT.evolrate.txt.remove | cut -f1 | sort -u > noBLATOGIDstoremove # 0 in total to remove from each
# Normally here, we remove saturated 4fold orthogroups from all branches and nodes - however, none to remove, so leave!!
# for i in *Prom-noBLAT.evolrate.txt ; do
#  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' noBLATOGIDstoremove $i | grep 'NA' | cut -f1-3 > "$(basename $i .txt).txtb" ;
# done
## ## ## ~ IN TOTAL, THERE ARE 60 NO BLAT 1:1 orthogroups left for evolutionary rate analysis ~ ## ## ## < as these are now filtered for saturated values, cat with the BLAT ones above to create a final set

cat Ab_4fold_Prom-BLAT.evolrate.txta Ab_4fold_Prom-noBLAT.evolrate.txt > Ab_4fold_Prom-BLAT.evolrate.txt2
cat Anc1_4fold_Prom-BLAT.evolrate.txta Anc1_4fold_Prom-noBLAT.evolrate.txt > Anc1_4fold_Prom-BLAT.evolrate.txt2
cat Anc2_4fold_Prom-BLAT.evolrate.txta Anc2_4fold_Prom-noBLAT.evolrate.txt > Anc2_4fold_Prom-BLAT.evolrate.txt2
cat Anc3_4fold_Prom-BLAT.evolrate.txta Anc3_4fold_Prom-noBLAT.evolrate.txt > Anc3_4fold_Prom-BLAT.evolrate.txt2
cat Mz_4fold_Prom-BLAT.evolrate.txta Mz_4fold_Prom-noBLAT.evolrate.txt > Mz_4fold_Prom-BLAT.evolrate.txt2
cat Nb_4fold_Prom-BLAT.evolrate.txta Nb_4fold_Prom-noBLAT.evolrate.txt > Nb_4fold_Prom-BLAT.evolrate.txt2
cat On_4fold_Prom-BLAT.evolrate.txta On_4fold_Prom-noBLAT.evolrate.txt > On_4fold_Prom-BLAT.evolrate.txt2
cat Pn_4fold_Prom-BLAT.evolrate.txta Pn_4fold_Prom-noBLAT.evolrate.txt > Pn_4fold_Prom-BLAT.evolrate.txt2
cat WholeTree_4fold_Prom-BLAT.evolrate.txta WholeTree_4fold_Prom-noBLAT.evolrate.txt > WholeTree_4fold_Prom-BLAT.evolrate.txt2

## ## ## ##
## ## ## ## ~ IN TOTAL, THERE ARE 4622 1to1 BLAT (4562) + no BLAT (60) orthogroups left for evolutionary rate analysis ~ ## ## ## ##
## ## ## ##


# first, use the per species orthogroup and gene IDs to pull out module assignments
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr
# already prepared the majority of the files that you need when looking at switching
cut -f1,2,17 OGIDS.txt5-clusterassign > OGIDS.txt5-clusterassign_Mz
cut -f1,3,18 OGIDS.txt5-clusterassign > OGIDS.txt5-clusterassign_Pn
cut -f1,4,19 OGIDS.txt5-clusterassign > OGIDS.txt5-clusterassign_Ab
cut -f1,5,20 OGIDS.txt5-clusterassign > OGIDS.txt5-clusterassign_Nb
cut -f1,6,21 OGIDS.txt5-clusterassign > OGIDS.txt5-clusterassign_On

# for each branch, seperate the files according to module assignment to plot the evolutionary rate per module and per species
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0917_USETHIS-BLAT

# create symbolic link to the files created above
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Mz
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Pn
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Ab
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Nb
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_On

# using the first column from each, awk match and pull out the rows
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Mz Mz_4fold_Prom-BLAT.evolrate.txt2 | cut -f1,2,3,5,6 > Mz_4fold_Prom-BLAT.evolrate.txt2_modules
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Pn Pn_4fold_Prom-BLAT.evolrate.txt2 | cut -f1,2,3,5,6 > Pn_4fold_Prom-BLAT.evolrate.txt2_modules
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab Ab_4fold_Prom-BLAT.evolrate.txt2 | cut -f1,2,3,5,6 > Ab_4fold_Prom-BLAT.evolrate.txt2_modules
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Nb Nb_4fold_Prom-BLAT.evolrate.txt2 | cut -f1,2,3,5,6 > Nb_4fold_Prom-BLAT.evolrate.txt2_modules
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_On On_4fold_Prom-BLAT.evolrate.txt2 | cut -f1,2,3,5,6 > On_4fold_Prom-BLAT.evolrate.txt2_modules

ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/list_k10

# split by module
while read F ; do awk -v module="$F" '{if ($5 == module) print $0;}' Mz_4fold_Prom-BLAT.evolrate.txt2_modules > $F-Mz_4fold_Prom-BLAT.evolrate.txt2_modules ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($5 == module) print $0;}' Pn_4fold_Prom-BLAT.evolrate.txt2_modules > $F-Pn_4fold_Prom-BLAT.evolrate.txt2_modules ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($5 == module) print $0;}' Ab_4fold_Prom-BLAT.evolrate.txt2_modules > $F-Ab_4fold_Prom-BLAT.evolrate.txt2_modules ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($5 == module) print $0;}' Nb_4fold_Prom-BLAT.evolrate.txt2_modules > $F-Nb_4fold_Prom-BLAT.evolrate.txt2_modules ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($5 == module) print $0;}' On_4fold_Prom-BLAT.evolrate.txt2_modules > $F-On_4fold_Prom-BLAT.evolrate.txt2_modules ; done < list_k10 # need to add -v to pass the variable as module to awk

# Add species name to each file
for i in *-Mz_4fold_Prom-BLAT.evolrate.txt2_modules ; do awk -v OFS="\t" 'NR==1{print $0, "M. zebra";next}{print $0,"M. zebra"}' $i > $i.2 ; done
for i in *-Pn_4fold_Prom-BLAT.evolrate.txt2_modules ; do awk -v OFS="\t" 'NR==1{print $0, "P. nyererei";next}{print $0,"P. nyererei"}' $i > $i.2 ; done
for i in *-Ab_4fold_Prom-BLAT.evolrate.txt2_modules ; do awk -v OFS="\t" 'NR==1{print $0, "A. burtoni";next}{print $0,"A. burtoni"}' $i > $i.2 ; done
for i in *-Nb_4fold_Prom-BLAT.evolrate.txt2_modules ; do awk -v OFS="\t" 'NR==1{print $0, "N. brichardi";next}{print $0,"N. brichardi"}' $i > $i.2 ; done
for i in *-On_4fold_Prom-BLAT.evolrate.txt2_modules ; do awk -v OFS="\t" 'NR==1{print $0, "O. niloticus";next}{print $0,"O. niloticus"}' $i > $i.2 ; done

rm *_4fold_Prom-BLAT.evolrate.txt2_modules # remove intermediate files

# merge same module files
while read F ; do cat $F-*_4fold_Prom-BLAT.evolrate.txt2_modules.2 > $F-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2 ; done < list_k10


## do similar as above but for ancestral assignments - have the ancestral files already but some lines are as OGXX_1 meaning not an extra clade in the phylo tree but meaning that orthogroup has no ab.gene and hence used the OGID

cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr

# first, change each OG to _0 from _1
for i in Anc*_clusterassign.txt ; do sed 's/_1/_0/g' $i > "$(basename $i .txt)_amend.txt" ; done

# then map each ab.gene to an OGID
awk '{print $2,$1,$3}' OFS='\t' OGIDS.txt5-clusterassign_Ab > OGIDS.txt5-clusterassign_Ab-rearrange # rearrange column order of Ab mapping file
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab-rearrange Anc1_clusterassign_amend.txt | awk -F'\t' '{ $3 = ($3 == "NA" ? $1 : $3) } 1' OFS='\t' | awk '{print $3,$1,$2}' OFS='\t' > OGIDS.txt5-clusterassign_Anc1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab-rearrange Anc2_clusterassign_amend.txt | awk -F'\t' '{ $3 = ($3 == "NA" ? $1 : $3) } 1' OFS='\t' | awk '{print $3,$1,$2}' OFS='\t' > OGIDS.txt5-clusterassign_Anc2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab-rearrange Anc3_clusterassign_amend.txt | awk -F'\t' '{ $3 = ($3 == "NA" ? $1 : $3) } 1' OFS='\t' | awk '{print $3,$1,$2}' OFS='\t' > OGIDS.txt5-clusterassign_Anc3
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab-rearrange Anc4_clusterassign_amend.txt | awk -F'\t' '{ $3 = ($3 == "NA" ? $1 : $3) } 1' OFS='\t' | awk '{print $3,$1,$2}' OFS='\t' > OGIDS.txt5-clusterassign_Anc4

# for each ancestral node, seperate the files according to module assignment to plot the evolutionary rate per module and per ancestral node
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Evolutionary_Rate/newEvolRate_0917_USETHIS-BLAT

# create symbolic link to the files created above
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Anc1
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Anc2
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Anc3
ln -s /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign_Anc4

# using the first column from each, awk match and pull out the rows
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Anc1 Anc1_4fold_Prom-BLAT.evolrate.txt2 | cut -f1,2,3,5,6 > Anc1_4fold_Prom-BLAT.evolrate.txt2_modules
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Anc2 Anc2_4fold_Prom-BLAT.evolrate.txt2 | cut -f1,2,3,5,6 > Anc2_4fold_Prom-BLAT.evolrate.txt2_modules
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Anc3 Anc3_4fold_Prom-BLAT.evolrate.txt2 | cut -f1,2,3,5,6 > Anc3_4fold_Prom-BLAT.evolrate.txt2_modules

# split by module
while read F ; do awk -v module="$F" '{if ($5 == module) print $0;}' Anc1_4fold_Prom-BLAT.evolrate.txt2_modules > $F-Anc1_4fold_Prom-BLAT.evolrate.txt2_modules ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($5 == module) print $0;}' Anc2_4fold_Prom-BLAT.evolrate.txt2_modules > $F-Anc2_4fold_Prom-BLAT.evolrate.txt2_modules ; done < list_k10 # need to add -v to pass the variable as module to awk
while read F ; do awk -v module="$F" '{if ($5 == module) print $0;}' Anc3_4fold_Prom-BLAT.evolrate.txt2_modules > $F-Anc3_4fold_Prom-BLAT.evolrate.txt2_modules ; done < list_k10 # need to add -v to pass the variable as module to awk

# Add node name to each file
for i in *-Anc1_4fold_Prom-BLAT.evolrate.txt2_modules ; do awk -v OFS="\t" 'NR==1{print $0, "Anc1 - M. zebra/P. nyererei";next}{print $0,"Anc1 - M. zebra/P. nyererei"}' $i > $i.2 ; done
for i in *-Anc2_4fold_Prom-BLAT.evolrate.txt2_modules ; do awk -v OFS="\t" 'NR==1{print $0, "Anc2 - M. zebra/P. nyererei/A. burtoni";next}{print $0,"Anc2 - M. zebra/P. nyererei/A. burtoni"}' $i > $i.2 ; done
for i in *-Anc3_4fold_Prom-BLAT.evolrate.txt2_modules ; do awk -v OFS="\t" 'NR==1{print $0, "Anc3 - M. zebra/P. nyererei/A. burtoni/N. brichardi";next}{print $0,"Anc3 - M. zebra/P. nyererei/A. burtoni/N. brichardi"}' $i > $i.2 ; done

rm *_4fold_Prom-BLAT.evolrate.txt2_modules # remove intermediate files

# merge same module files
while read F ; do cat $F-Anc*_4fold_Prom-BLAT.evolrate.txt2_modules.2 > $F-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2 ; done < list_k10


### Sort the evolutionary rate in promoter regions to rank them
for i in *.txt2 ; do sort -k3,3nr $i > "$(basename $i ).promsort" ; done

for i in *evolrate.txt2.promsort ; do awk -F'\t' '($3>"0.2")' $i > "$(basename $i ).promsort0.2" ; done # select those with the highest evol rate (>0.2)

11 Ab_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2
 1 Anc1_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2
25 Anc2_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2
66 Anc3_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2
21 Mz_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2
55 Nb_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2
70 On_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2
 7 Pn_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2
473 WholeTree_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2
729 total

# Run GO enrichment of the above genes against the 1:1 orthologs
# first map each orthgroup to its cichlid gene, then cut cichlid gene column and create GO input files
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Mz Mz_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2 | cut -f4 | sed -e '1i\
NA
' | awk '{ for(i=1;i<=NF;i++){ print $i}}' | tr '\n' '#' | sed $'s/NA/Mz_promsort0.2\t/g' > Mz_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Pn Pn_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2 | cut -f4 | sed -e '1i\
NA
' | awk '{ for(i=1;i<=NF;i++){ print $i}}' | tr '\n' '#' | sed $'s/NA/Pn_promsort0.2\t/g' > Pn_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab Ab_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2 | cut -f4 | sed -e '1i\
NA
' | awk '{ for(i=1;i<=NF;i++){ print $i}}' | tr '\n' '#' | sed $'s/NA/Ab_promsort0.2\t/g' > Ab_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Nb Nb_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2 | cut -f4 | sed -e '1i\
NA
' | awk '{ for(i=1;i<=NF;i++){ print $i}}' | tr '\n' '#' | sed $'s/NA/Nb_promsort0.2\t/g' > Nb_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_On On_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2 | cut -f4 | sed -e '1i\
NA
' | awk '{ for(i=1;i<=NF;i++){ print $i}}' | tr '\n' '#' | sed $'s/NA/On_promsort0.2\t/g' > On_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab Anc1_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2 | cut -f4 | sed -e '1i\
NA
' | awk '{ for(i=1;i<=NF;i++){ print $i}}' | tr '\n' '#' | sed $'s/NA/Anc1_promsort0.2\t/g' > Anc1_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab Anc2_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2 | cut -f4 | sed -e '1i\
NA
' | awk '{ for(i=1;i<=NF;i++){ print $i}}' | tr '\n' '#' | sed $'s/NA/Anc2_promsort0.2\t/g' > Anc2_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab Anc3_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2 | cut -f4 | sed -e '1i\
NA
' | awk '{ for(i=1;i<=NF;i++){ print $i}}' | tr '\n' '#' | sed $'s/NA/Anc3_promsort0.2\t/g' > Anc3_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt


# create genelist of 1to1 ortholog genes (with evol rate - so only 3876/6844 1:1) only to use as background (from biological process folder only)
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_4fold_Prom-BLAT.evolrate.txt /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1 | grep -v NA > OGIDS.txt5-clusterassign-1to1evolrate
cut -f2 OGIDS.txt5-clusterassign-1to1evolrate > OGIDS.txt5-clusterassign-1to1-Mz-genelist.txt
cut -f3 OGIDS.txt5-clusterassign-1to1evolrate > OGIDS.txt5-clusterassign-1to1-Pn-genelist.txt
cut -f4 OGIDS.txt5-clusterassign-1to1evolrate > OGIDS.txt5-clusterassign-1to1-Ab-genelist.txt
cut -f5 OGIDS.txt5-clusterassign-1to1evolrate > OGIDS.txt5-clusterassign-1to1-Nb-genelist.txt
cut -f6 OGIDS.txt5-clusterassign-1to1evolrate > OGIDS.txt5-clusterassign-1to1-On-genelist.txt

# copy to local folder: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917

cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917
# Run enrichment
/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer Mz_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt OGIDS.txt5-clusterassign-1to1-Mz-genelist.txt /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Meze_ 1 Mz_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2.GOOUTPUT persg

/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer Pn_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt OGIDS.txt5-clusterassign-1to1-Pn-genelist.txt /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Puny_ 1 Pn_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2.GOOUTPUT persg

/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer Ab_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt OGIDS.txt5-clusterassign-1to1-Ab-genelist.txt /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Asbu_ 1 Ab_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2.GOOUTPUT persg

/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer Nb_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt OGIDS.txt5-clusterassign-1to1-Nb-genelist.txt /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Nebr_ 1 Nb_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2.GOOUTPUT persg

/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer On_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt OGIDS.txt5-clusterassign-1to1-On-genelist.txt /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Orni_ 1 On_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2.GOOUTPUT persg

/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer Anc1_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt OGIDS.txt5-clusterassign-1to1-Ab-genelist.txt /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Asbu_ 1 Anc1_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2.GOOUTPUT persg

/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer Anc2_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt OGIDS.txt5-clusterassign-1to1-Ab-genelist.txt /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Asbu_ 1 Anc2_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2.GOOUTPUT persg

/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/EnrichAnalyzer/enrichAnalyzer Anc3_4fold_Prom-BLAT.evolrate.txt.promsort.promsort0.2.GOINPUT.txt OGIDS.txt5-clusterassign-1to1-Ab-genelist.txt /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v2/10.EdgeGOenrichment/GO_Analysis/GO_terms/filtered/biological_process_only/Asbu_ 1 Anc3_4fold_Prom-BLAT.evolrate.txt2.promsort.promsort0.2.GOOUTPUT persg

# NO ENRICHMENT!

# Select those with a q-value (adjusted p-value) of < 0.05
for i in *GOOUTPUT_details.txt ; do
	awk -F '\t' '$4 < 0.05' $i > "$(basename "$i" .txt)_filtered.txt" ;
done # left nothing in most files - exception is Ab!

# for plotting the GO (just doing based on -log10 p-val <0.05) prepare the OUTPUT file so that matches will show the GO output on each line for each gene
for i in *GOOUTPUT_details.txt ; do awk -F '\t' '$3 < 0.05' $i | grep -v biological_process > "$(basename "$i" .txt)_filtered2.txt" ; done # filter based on p-val < 0.05 and remove any matching 'biological_process' (NOTE: no mutliple test correction filtering yielded anything - ok to just use this?!?!)
for i in *details_filtered2.txt ; do sed 's/$/;/g' $i | awk '{ gsub(";", ";"$1";"$2";"$3";"$4";"$5";"$6";"$7";"$8"\n") } 1' | cut -f10 | awk '{ gsub(";", "\t") } 1' | sed '/^$/d' > "$(basename "$i" )2" ; done
# remove first column then only display unique entries as final column displays no. of genes assigned to the term if you need it!
for i in *details_filtered2.txt2 ; do cut -f2-9 $i | sort -u | sort -k3,3 > "$(basename "$i" 2)3" ; done

## Pull out switching and non-switching genes from branches and ancestral nodes to plot

# Switching (state changes) at branches

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Mzvsall Mz_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Mz_4fold_Prom-BLAT.evolrate.txt2.sw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Mzvsallns Mz_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Mz_4fold_Prom-BLAT.evolrate.txt2.nsw

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Pnvsall Pn_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Pn_4fold_Prom-BLAT.evolrate.txt2.sw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Pnvsallns Pn_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Pn_4fold_Prom-BLAT.evolrate.txt2.nsw

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Abvsall Ab_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Ab_4fold_Prom-BLAT.evolrate.txt2.sw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Abvsallns Ab_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Ab_4fold_Prom-BLAT.evolrate.txt2.nsw

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Nbvsall Nb_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Nb_4fold_Prom-BLAT.evolrate.txt2.sw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Nbvsallns Nb_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Nb_4fold_Prom-BLAT.evolrate.txt2.nsw

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Onvsall On_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > On_4fold_Prom-BLAT.evolrate.txt2.sw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_Onvsallns On_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > On_4fold_Prom-BLAT.evolrate.txt2.nsw


### Plotting against tau to see if there is any correlation
# Map cichlid gene IDs in tau tables to OGIDs
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Ab-rearrange Ab_tau_final.txt | grep -v NA > Ab_tau_final.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Mz-rearrange Mz_tau_final.txt | grep -v NA > Mz_tau_final.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Pn-rearrange Pn_tau_final.txt | grep -v NA > Pn_tau_final.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_Nb-rearrange Nb_tau_final.txt | grep -v NA > Nb_tau_final.txt2
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS.txt5-clusterassign_On-rearrange On_tau_final.txt | grep -v NA > On_tau_final.txt2


# There are two options for the ancestral switching (state changes):

  # 1. Ancestral module assignment switches that are state changes from the LCA (previous ancestral node) - this excludes -1 non-assigned < make best sense to use this owing to an integrated approach from Arboretum
  #Anc4 represents the ancestral state derived from On and hence 0 state changes
  #Anc3 - Anc3 vs Anc4
  /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc3vs4_clusterassign_noOG-allassigned_StateChange.txt #792
  #Anc2 - Anc2 vs Anc3
  /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc2vs3_clusterassign_noOG-allassigned_StateChange.txt #293
  #Anc1 - Anc1 vs Anc2
  /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc1vs2_clusterassign_noOG-allassigned_StateChange.txt #404

# the above files only have ab.gene so we need to map to OGID to map below

for i in /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc*vs*_clusterassign_noOG-allassigned_StateChange.txt ;
  do awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Module_genesandexpr/OGIDS.txt5-ab $i | awk '{print $6,$1,$2,$3,$4,$5}' OFS='\t' > /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/"$(basename $i _clusterassign_noOG-allassigned_StateChange.txt)_clusterassign_OGIDs-allassigned_StateChange.txt" ;
done

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc1vs2_clusterassign_OGIDs-allassigned_StateChange.txt Anc1_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Anc1_4fold_Prom-BLAT.evolrate.txt2.Ancsw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc1vs2_clusterassign_OGIDs-allassigned_StateChange.txt Anc1_4fold_Prom-BLAT.evolrate.txt2 | grep NA | cut -f1-3 > Anc1_4fold_Prom-BLAT.evolrate.txt2.Ancnsw

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc2vs3_clusterassign_OGIDs-allassigned_StateChange.txt Anc2_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Anc2_4fold_Prom-BLAT.evolrate.txt2.Ancsw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc2vs3_clusterassign_OGIDs-allassigned_StateChange.txt Anc2_4fold_Prom-BLAT.evolrate.txt2 | grep NA | cut -f1-3 > Anc2_4fold_Prom-BLAT.evolrate.txt2.Ancnsw

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc3vs4_clusterassign_OGIDs-allassigned_StateChange.txt Anc3_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Anc3_4fold_Prom-BLAT.evolrate.txt2.Ancsw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/Anc3vs4_clusterassign_OGIDs-allassigned_StateChange.txt Anc3_4fold_Prom-BLAT.evolrate.txt2 | grep NA | cut -f1-3 > Anc3_4fold_Prom-BLAT.evolrate.txt2.Ancnsw


  # 2. Species that share the same LCA, have the same switch (state change) compared to the other species e.g. PnMzvsAbNbOn; PnMzAbvsNbOn; PnMzAbNbvsOn
  # Pn & Mz same switch vs Ab, Nb and On
  /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_MzPnvsAbNbOn #89
  # Pn, Mz & Ab same switch vs Nb and On
  /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_MzPnAbvsNbOn #172
  # Pn, Mz, Ab & Nb same switch vs On
  /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_MzPnAbNbvsOn #983

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_MzPnvsAbNbOn Anc1_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Anc1_4fold_Prom-BLAT.evolrate.txt2.SpAncsw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_MzPnvsAbNbOn Anc1_4fold_Prom-BLAT.evolrate.txt2 | grep NA | cut -f1-3 > Anc1_4fold_Prom-BLAT.evolrate.txt2.SpAncnsw

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_MzPnAbvsNbOn Anc2_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Anc2_4fold_Prom-BLAT.evolrate.txt2.SpAncsw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_MzPnAbvsNbOn Anc2_4fold_Prom-BLAT.evolrate.txt2 | grep NA | cut -f1-3 > Anc2_4fold_Prom-BLAT.evolrate.txt2.SpAncnsw

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_MzPnAbNbvsOn Anc3_4fold_Prom-BLAT.evolrate.txt2 | grep -v NA | cut -f1-3 > Anc3_4fold_Prom-BLAT.evolrate.txt2.SpAncsw
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-clusterassign-1to1_MzPnAbNbvsOn Anc3_4fold_Prom-BLAT.evolrate.txt2 | grep NA | cut -f1-3 > Anc3_4fold_Prom-BLAT.evolrate.txt2.SpAncnsw

## NOTE: the script for plotting figures can be found here ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/evolRate_figs_v2.R; and evolRate_figs_v2.RData
## There are several 4fold alignments that are <100nt in length, all details accordingly:
  # all alignment lengths: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/4fold_aln_lengths.txt
  # all alignment lengths <100nt : /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/4fold_aln_lengths.txt2
  # per species all alignments <100nt : /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/*-4fold_aln_lengths.txt
  # per species all alignments <100nt and 4fold = 0.000004 (includes those that are not too with NA) : /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/*-4foldlessthan100nt_low4fold.txt

### RSAT-based motif imputed p-values
mkdir /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/RSAT_sept17
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/RSAT_sept17

# copied the tar balls that Will shared
cp rsat_nw.tar.gz .
cp rsat_sub_ON_temp_MA0002.2.transfac.trial.sh .

# untar the file
tar -xvzf rsat_nw.tar.gz

### DO NOT FOLLOW THIS ORDERING AS WE CHANGED SEVERAL THINGS
# compile rsat
make -f makefiles/init_rsat.mk compile_all

# did several modifications to the scrips etc.
mkdir logs/ # create a log folder otherwise scripts won't run

### Now use new JASPAR 2018
mkdir /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/RSAT_oct17_Jaspar2018
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/RSAT_oct17_Jaspar2018

# copied JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.txt from local
cp JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.txt /tgac/workarea/Research-Groups/RG-cichlids/
cp ../RSAT_sept17/run_mat_qual2.py .
sed -i 's/run2_/run3_/g' run_mat_qual2.py
mkdir logs/ # create a log folder otherwise scripts won't run
# Run just for Pn
ml python/3.5
python3.5 run_mat_qual2.py /tgac/workarea/Research-Groups/RG-cichlids/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac_v1.txt /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta 1000 ../RSAT_sept17/Pn_bg_use.txt PN # out > rsat_results/[motif-folder]
# runs 1000 permutations creating *.transfac_PN__scan_positive_set_1000perm_score_distrib.tab
# Will creating a script to parse all these permutation files
# this will output all files to /tgac/workarea/Research-Groups/RG-cichlids/ as run3_

# this was for moving results to the shared folders
# source rsync-3.1.1
# rsync -abviuP /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/motif_discovery/promoters/RSAT_sept17/rsat_results/MA*/ /tgac/workarea/Research-Groups/RG-cichlids/rsat_results/MA*/
# for i in /tgac/workarea/Research-Groups/RG-cichlids/rsat_results/MA* ; do rm $i/temp_MA* ; done

## Running matrix quality on all extrapolated cichlid-specific PSSMs (Jan 2018)
cd /tgac/workarea/Research-Groups/RG-cichlids/
source python-3.5.1

python3.5 /tgac/workarea/Research-Groups/RG-cichlids/run_mat_qual_nm_TM.py /tgac/workarea/Research-Groups/RG-cichlids/GTRDhuman/human_sites/TF_transfac/Mzeb/ /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta 1000 /tgac/workarea/Research-Groups/RG-cichlids/Mz_bg_use.txt mz /tgac/workarea/Research-Groups/RG-cichlids/GTRDhuman/human_sites/CichLevel_TF_transfac/
python3.5 /tgac/workarea/Research-Groups/RG-cichlids/run_mat_qual_nm_TM.py /tgac/workarea/Research-Groups/RG-cichlids/GTRDdata/mouse_sites/TF_transfac/Mzeb/ /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta 1000 /tgac/workarea/Research-Groups/RG-cichlids/Mz_bg_use.txt mz /tgac/workarea/Research-Groups/RG-cichlids/GTRDdata/mouse_sites/CichLevel_TF_transfac/
python3.5 /tgac/workarea/Research-Groups/RG-cichlids/run_mat_qual_nm_TM.py /tgac/workarea/Research-Groups/RG-cichlids/GTRDdata/mouse_sites/TF_transfac/Pnye/ /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta 1000 /tgac/workarea/Research-Groups/RG-cichlids/Pn_bg_use.txt pn /tgac/workarea/Research-Groups/RG-cichlids/GTRDdata/mouse_sites/CichLevel_TF_transfac/

### Check if the extrapolated TF-TG match your co-expressed TF-TG
cd /tgac/workarea/Research-Groups/RG-cichlids
mkdir GTRD_TFTGcoexpCheck
cd GTRDhuman/human_sites/
for i in hg38_tftgtfbs_*_07.extrap.output ; do sed 's/;/\t/g' $i | awk '{print $1":"$2,$3}' OFS='\t' > ../../GTRD_TFTGcoexpCheck/$i ; done
cd ../../GTRDdata/mouse_sites/
for i in mm10_tftgtfbs_*_07.extrap.output ; do sed 's/;/\t/g' $i | awk '{print $1":"$2,$3}' OFS='\t' > ../../GTRD_TFTGcoexpCheck/$i ; done
cd ../../GTRD_TFTGcoexpCheck/
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/8.TFTGco/
for i in *.txt ; do awk '{print $1":"$2,$3}' OFS='\t' $i > /tgac/workarea/Research-Groups/RG-cichlids/GTRD_TFTGcoexpCheck/$i ; done
cd /tgac/workarea/Research-Groups/RG-cichlids/GTRD_TFTGcoexpCheck
rm tftgco_merged.txt

for i in *_Ab_GeneAnnotation_11092017_07.extrap.output ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $i Ab.txt | grep -v NA > "$(basename "$i" _GeneAnnotation_11092017_07.extrap.output)_07_TFTGcoexpr.txt" ;
done

for i in *_Mz_GeneAnnotation_11092017_07.extrap.output ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $i Mz.txt | grep -v NA > "$(basename "$i" _GeneAnnotation_11092017_07.extrap.output)_07_TFTGcoexpr.txt" ;
done

for i in *_Nb_GeneAnnotation_11092017_07.extrap.output ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $i Nb.txt | grep -v NA > "$(basename "$i" _GeneAnnotation_11092017_07.extrap.output)_07_TFTGcoexpr.txt" ;
done

for i in *_Pn_GeneAnnotation_11092017_07.extrap.output ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $i Pn.txt | grep -v NA > "$(basename "$i" _GeneAnnotation_11092017_07.extrap.output)_07_TFTGcoexpr.txt" ;
done

for i in *_Oreochromis_niloticus_07.extrap.output ; do
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' $i On.txt | grep -v NA > "$(basename "$i" _GeneAnnotation_11092017_07.extrap.output)_07_TFTGcoexpr.txt" ;
done

# 4938/1240003 hg38_tftgtfbs_Ab_07_TFTGcoexpr.txt
# 4762/1122098 hg38_tftgtfbs_Mz_07_TFTGcoexpr.txt
# 4247/674254 hg38_tftgtfbs_Nb_07_TFTGcoexpr.txt
# 5632/1155912 hg38_tftgtfbs_Oreochromis_niloticus_07.extrap.output_07_TFTGcoexpr.txt
# 5720/941826 hg38_tftgtfbs_Pn_07_TFTGcoexpr.txt

# 11690/3201969 mm10_tftgtfbs_Ab_07_TFTGcoexpr.txt
# 11912/2638846 mm10_tftgtfbs_Mz_07_TFTGcoexpr.txt
# 11392/2009045 mm10_tftgtfbs_Nb_07_TFTGcoexpr.txt
# 15807/3073437 mm10_tftgtfbs_Oreochromis_niloticus_07.extrap.output_07_TFTGcoexpr.txt
# 15135/2655667 mm10_tftgtfbs_Pn_07_TFTGcoexpr.txt

# Check against module co-expressed genes

# merge all for Hs and Mm
cat hg38_*.output > hg38_tftgtfbs_merged_07.extrap.output # 5134093
sed -i 's/:/\t/g' hg38_tftgtfbs_merged_07.extrap.output
cat mm10_*.output > mm10_tftgtfbs_merged_07.extrap.output # 13578964
sed -i 's/:/\t/g' mm10_tftgtfbs_merged_07.extrap.output

echo 0 > list_modules
echo 1 >> list_modules
echo 2 >> list_modules
echo 3 >> list_modules
echo 4 >> list_modules
echo 5 >> list_modules
echo 6 >> list_modules
echo 7 >> list_modules
echo 8 >> list_modules
echo 9 >> list_modules

awk '{print $1}' hg38_tftgtfbs_merged_07.extrap.output > hg38_tftgtfbs_merged_07.extrap.output-cutf1
awk '{print $2,$3}' OFS='\t' hg38_tftgtfbs_merged_07.extrap.output > hg38_tftgtfbs_merged_07.extrap.output-cutf2-3
while read F ; do awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/$F-modulegenes.txt hg38_tftgtfbs_merged_07.extrap.output-cutf1 > $F-hg38_tftgtfbs_merged_07.extrap.output-cutf1 ; done < list_modules
while read F ; do awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/$F-modulegenes.txt hg38_tftgtfbs_merged_07.extrap.output-cutf2-3 > $F-hg38_tftgtfbs_merged_07.extrap.output-cutf2-3 ; done < list_modules
while read F ; do paste -d'\t' $F-hg38_tftgtfbs_merged_07.extrap.output-cutf1 $F-hg38_tftgtfbs_merged_07.extrap.output-cutf2-3 > $F-hg38_tftgtfbs_merged_07.extrap.output ; done < list_modules
cat *-hg38_tftgtfbs_merged_07.extrap.output > ModuleMerged_hg38_tftgtfbs_merged_07.extrap.output
rm *-hg38_tftgtfbs_merged_07.extrap.output*
rm hg38_tftgtfbs_merged_07.extrap.output-*
awk '$2!="NA" && $5!="NA"' ModuleMerged_hg38_tftgtfbs_merged_07.extrap.output | wc -l #636810/5134093 = 12.4%

awk '{print $1}' mm10_tftgtfbs_merged_07.extrap.output > mm10_tftgtfbs_merged_07.extrap.output-cutf1
awk '{print $2,$3}' OFS='\t' mm10_tftgtfbs_merged_07.extrap.output > mm10_tftgtfbs_merged_07.extrap.output-cutf2-3
while read F ; do awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/$F-modulegenes.txt mm10_tftgtfbs_merged_07.extrap.output-cutf1 > $F-mm10_tftgtfbs_merged_07.extrap.output-cutf1 ; done < list_modules
while read F ; do awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$1;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/$F-modulegenes.txt mm10_tftgtfbs_merged_07.extrap.output-cutf2-3 > $F-mm10_tftgtfbs_merged_07.extrap.output-cutf2-3 ; done < list_modules
while read F ; do paste -d'\t' $F-mm10_tftgtfbs_merged_07.extrap.output-cutf1 $F-mm10_tftgtfbs_merged_07.extrap.output-cutf2-3 > $F-mm10_tftgtfbs_merged_07.extrap.output ; done < list_modules
cat *-mm10_tftgtfbs_merged_07.extrap.output > ModuleMerged_mm10_tftgtfbs_merged_07.extrap.output
rm *-mm10_tftgtfbs_merged_07.extrap.output*
rm mm10_tftgtfbs_merged_07.extrap.output-*
awk '$2!="NA" && $5!="NA"' ModuleMerged_mm10_tftgtfbs_merged_07.extrap.output | wc -l #1769346/13578964 = 13%

## check the outputs of some opsins

cd /tgac/workarea/Research-Groups/RG-cichlids/GTRD_TFTGcoexpCheck

# sws1 - mz.gene.s102.69
grep -wiF mz.gene.s102.69 ModuleMerged_mm10_tftgtfbs_merged_07.extrap.output | cut -f1 | sort -u
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3 # local
for i in mz.gene.s120.33 mz.gene.s127.35 mz.gene.s141.46 mz.gene.s203.14 mz.gene.s28.101 mz.gene.s31.137 mz.gene.s50.65 mz.gene.s5.254 mz.gene.s54.33 mz.gene.s57.27 mz.gene.s5.74 mz.gene.s78.27 mz.gene.s9.124 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done

nkx2.4a	nkx2.4a	NKX2-1
NULL	jun	JUN
NULL	foxa2	FOXA2
pknox1.2	pknox1.2	PKNOX1
stat5a	stat5a	STAT5A
NULL	PBX1_(1_of_many)	PBX1
NULL	ar	AR
runx2b	runx2b	RUNX2
NULL	zbtb16b	ZBTB16
stat6	stat6	STAT6
NULL	fosab	FOS
NULL	nr5a2	NR5A2
tcf7l1b	tcf7l1b	TCF7L1

# sws1 - ab.gene.s2279.1
grep -wiF ab.gene.s2279.1 ModuleMerged_mm10_tftgtfbs_merged_07.extrap.output | cut -f1 | sort -u | paste -s -d " "
for i in ab.gene.s111.11 ab.gene.s130.19 ab.gene.s138.2 ab.gene.s167.12 ab.gene.s189.12 ab.gene.s19.91 ab.gene.s25.108 ab.gene.s31.77 ab.gene.s35.64 ab.gene.s35.65 ab.gene.s376.13 ab.gene.s38.12 ab.gene.s47.42 ab.gene.s50.54 ab.gene.s619.2 ab.gene.s66.35 ab.gene.s6.86 ab.gene.s7.113 ab.gene.s96.36 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done

NULL	foxa2	FOXA2
NULL	vdra	VDR
NULL	zbtb16b	ZBTB16
ascl1a	ascl1a	ASCL1
NULL	jun	JUN
runx2b	runx2b	RUNX2
nr4a1	nr4a1	NR4A1
NULL	PBX1_(1_of_many)	PBX1
stat5a	stat5a	STAT5A
stat5a	stat5a	STAT5A
pknox1.2	pknox1.2	PKNOX1
NULL	nr5a2	NR5A2
nkx2.4a	nkx2.4a	NKX2-1
tcf7l1b	tcf7l1b	TCF7L1
rbpjl	NULL	RBPJL
NULL	NR4A1_(1_of_many)	NR4A1
NULL	fosab	FOS
stat6	stat6	STAT6
NULL	ar	AR

# sws1 - nb.gene.s1.386
grep -wiF nb.gene.s1.386 ModuleMerged_mm10_tftgtfbs_merged_07.extrap.output | cut -f1 | sort -u | paste -s -d " "
for i in nb.gene.s148.46 nb.gene.s17.254 nb.gene.s22.90 nb.gene.s235.10 nb.gene.s26.65 nb.gene.s28.92 nb.gene.s31.109 nb.gene.s31.93 nb.gene.s48.52 nb.gene.s55.94 nb.gene.s64.60 nb.gene.s73.73 nb.gene.s80.43 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done

NULL	PBX1_(1_of_many)	PBX1
NULL	nr5a2	NR5A2
nr1d2b	nr1d2b	NR1D2
runx2b	runx2b	RUNX2
NULL	foxa2	FOXA2
cebpb	cebpb	CEBPB
nr4a1	nr4a1	NR4A1
NULL	fosab	FOS
stat5a	stat5a	STAT5A
NULL	ar	AR
nkx2.4a	nkx2.4a	NKX2-1
NULL	zbtb16b	ZBTB16

# sws2a, rh2 and rho absent from Hs and Mm extrapolation


### 05/03/2018 Check the very first fimo runs using Cichlid specific, wide and JASPAR motifs

# First check the counts of edges and TFs (before and after filtering)

cd /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/human
# same can also be ran on
cd /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse

for i in CS_default_*_results.out ; do wc -l $i | awk -F' ' '{print $1}' ; done
for i in CS_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in CS_default_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in CS_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in CS_matQualPval_*_results.out ; do wc -l $i | awk -F' ' '{print $1}' ; done
for i in CS_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in CS_matQualPval_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in CS_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in CW_default_*_results.out ; do wc -l $i | awk -F' ' '{print $1}' ; done
for i in CW_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in CW_default_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in CW_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in CW_matQualPval_*_results.out ; do wc -l $i | awk -F' ' '{print $1}' ; done
for i in CW_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in CW_matQualPval_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in CW_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in JASPAR_default_*_results.out ; do wc -l $i | awk -F' ' '{print $1}' ; done
for i in JASPAR_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in JASPAR_default_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in JASPAR_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in JASPAR_matQualPval_*_results.out ; do wc -l $i | awk -F' ' '{print $1}' ; done
for i in JASPAR_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in JASPAR_matQualPval_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in JASPAR_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

# get total of TFs (except CW)
cat CS_default_ab_results.out CS_matQualPval_ab_results.out JASPAR_default_ab_results.out JASPAR_matQualPval_ab_results.out | cut -f1 | sort -u | wc -l
cat CS_default_mz_results.out CS_matQualPval_mz_results.out JASPAR_default_mz_results.out JASPAR_matQualPval_mz_results.out | cut -f1 | sort -u | wc -l
cat CS_default_nb_results.out CS_matQualPval_nb_results.out JASPAR_default_nb_results.out JASPAR_matQualPval_nb_results.out | cut -f1 | sort -u | wc -l
cat CS_default_on_results.out CS_matQualPval_on_results.out JASPAR_default_on_results.out JASPAR_matQualPval_on_results.out | cut -f1 | sort -u | wc -l
cat CS_default_pn_results.out CS_matQualPval_pn_results.out JASPAR_default_pn_results.out JASPAR_matQualPval_pn_results.out | cut -f1 | sort -u | wc -l
# get total of TFs (filtered qval < 0.05)
cat CS_default_ab_results.out CS_matQualPval_ab_results.out JASPAR_default_ab_results.out JASPAR_matQualPval_ab_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l
cat CS_default_mz_results.out CS_matQualPval_mz_results.out JASPAR_default_mz_results.out JASPAR_matQualPval_mz_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l
cat CS_default_nb_results.out CS_matQualPval_nb_results.out JASPAR_default_nb_results.out JASPAR_matQualPval_nb_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l
cat CS_default_on_results.out CS_matQualPval_on_results.out JASPAR_default_on_results.out JASPAR_matQualPval_on_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l
cat CS_default_pn_results.out CS_matQualPval_pn_results.out JASPAR_default_pn_results.out JASPAR_matQualPval_pn_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l


# GTRD
cd /tgac/workarea/Research-Groups/RG-cichlids/GTRDhuman
wc -l hg38_TFTGTFBS_sites_nr.in
cut -f1 hg38_TFTGTFBS_sites_nr.in | sort -u | wc -l

cd /tgac/workarea/Research-Groups/RG-cichlids/GTRDdata
wc -l mm10_TFTGTFBS_sites_nr.in
cut -f1 mm10_TFTGTFBS_sites_nr.in | sort -u | wc -l


## Looking at numbers after extrapolation
cd /tgac/workarea/Research-Groups/RG-cichlids/GTRDhuman/human_sites/
for i in hg38_tftgtfbs_*07.extrap.output ; do wc -l $i ; done
for i in hg38_tftgtfbs_*07.extrap.output ; do awk -F';' '{print $1}' $i | sort -u | wc -l ; done

cd /tgac/workarea/Research-Groups/RG-cichlids/GTRD_TFTGcoexpCheck
for i in mm10_tftgtfbs_*_GeneAnnotation_11092017_07.extrap.output ; do wc -l $i ; done
wc -l mm10_tftgtfbs_Oreochromis_niloticus_07.extrap.output
for i in mm10_tftgtfbs_*_GeneAnnotation_11092017_07.extrap.output ; do awk -F':' '{print $1}' $i | sort -u | wc -l ; done
awk -F':' '{print $1}' mm10_tftgtfbs_Oreochromis_niloticus_07.extrap.output | sort -u | wc -l

### Check for the opsins in the fimo runs
cd /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/human

##### SWS1

# sws1 - mz.gene.s102.69
awk '($8 < 0.05 )' OFS='\t' CS_matQualPval_mz_results.out | grep -wiF mz.gene.s102.69 | cut -f1 | sort -u | paste -sd " "
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3 # local
for i in mz.gene.s10.82 mz.gene.s136.36 mz.gene.s171.22 mz.gene.s17.45 mz.gene.s186.6 mz.gene.s26.74 mz.gene.s31.45 mz.gene.s31.46 mz.gene.s4.9 mz.gene.s50.74 mz.gene.s6.87 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
NULL	GATA2_(1_of_many)	GATA2
klf5b	klf5b	KLF5
NULL	egr2b	EGR2
nrf1	nrf1	NRF1
tp53	tp53	TP53
NULL	MXI1	MXI1
NULL	tcf3b	TCF3
si:dkey-110c1.10	tcf3b	TCF3
egr1	egr1	EGR1
ncor1	NULL	NCOR1
NULL	tcf3a	TCF3

#non-filtered for q-val
for i in mz.gene.s0.105 mz.gene.s0.220 mz.gene.s10.105 mz.gene.s10.123 mz.gene.s10.160 mz.gene.s10.326 mz.gene.s10.82 mz.gene.s109.30 mz.gene.s11.168 mz.gene.s113.2 mz.gene.s120.33 mz.gene.s12.260 mz.gene.s123.8 mz.gene.s130.9 mz.gene.s134.8 mz.gene.s135.29 mz.gene.s136.36 mz.gene.s137.11 mz.gene.s13.9 mz.gene.s141.31 mz.gene.s141.46 mz.gene.s146.35 mz.gene.s149.5 mz.gene.s15.135 mz.gene.s15.176 mz.gene.s153.32 mz.gene.s154.17 mz.gene.s162.10 mz.gene.s171.22 mz.gene.s17.142 mz.gene.s17.45 mz.gene.s186.6 mz.gene.s19.76 mz.gene.s199.19 mz.gene.s202.22 mz.gene.s211.9 mz.gene.s221.14 mz.gene.s2.35 mz.gene.s23.79 mz.gene.s24.187 mz.gene.s249.26 mz.gene.s24.95 mz.gene.s25.137 mz.gene.s26.74 mz.gene.s27.101 mz.gene.s27.110 mz.gene.s31.45 mz.gene.s31.46 mz.gene.s31.80 mz.gene.s32.112 mz.gene.s321.14 mz.gene.s33.72 mz.gene.s3.450 mz.gene.s35.42 mz.gene.s39.42 mz.gene.s39.68 mz.gene.s4.32 mz.gene.s45.62 mz.gene.s4.9 mz.gene.s50.74 mz.gene.s51.58 mz.gene.s5.4 mz.gene.s55.53 mz.gene.s5.74 mz.gene.s57.92 mz.gene.s59.22 mz.gene.s60.4 mz.gene.s6.166 mz.gene.s6.241 mz.gene.s6.263 mz.gene.s63.72 mz.gene.s66.25 mz.gene.s6.87 mz.gene.s69.29 mz.gene.s7.32 mz.gene.s77.62 mz.gene.s79.13 mz.gene.s82.26 mz.gene.s8.39 mz.gene.s9.30 mz.gene.s97.34 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
smad4a	SMAD4_(1_of_many)	SMAD4
pou5f3	pou5f3	POU5F1
stat2	stat2	STAT2
NULL	BHLHE40_(1_of_many)	BHLHE40
kdm5ba	kdm5ba	KDM5B
gata2a	gata2a	GATA2
NULL	GATA2_(1_of_many)	GATA2
zeb1a	zeb1a	ZEB1
NULL	fosl1a	FOSL1
tbx21	tbx21	TBX21
nkx2.4a	nkx2.4a	NKX2-1
kdm5bb	kdm5bb	KDM5B
NULL	NFATC1_(1_of_many)	NFATC1
foxp1b	FOXP1	FOXP1
NULL	TCF12_(1_of_many)	TCF12
CABZ01066696.1	NULL	LYL1
klf5b	klf5b	KLF5
NULL	eomesb	EOMES
tfap2a	tfap2a	TFAP2A
atf3	atf3	ATF3
NULL	foxa2	FOXA2
prdm1c	PRDM1_(1_of_many)	PRDM1
creb1b	creb1a_(1_of_many)	CREB1
tfap2c	tfap2c	TFAP2C
sp1	sp1	SP1
nfatc1	nfatc1	NFATC1
NULL	nfic	NFIC
NULL	foxa1	FOXA1
NULL	egr2b	EGR2
foxm1	foxm1	FOXM1
nrf1	nrf1	NRF1
tp53	tp53	TP53
eomesa	eomesa	EOMES
irf4b	irf4b	IRF4
gata1a	gata1a	GATA1
zeb1b	zeb1b	ZEB1
pou2f2a	POU2F2_(1_of_many)	POU2F2
pax6b	pax6b	PAX6
esr2a	esr2a	ESR2
NULL	USF2_(1_of_many)	USF2
fli1b	fli1b	ERG
pou2f2a	POU2F2_(1_of_many)	POU2F2
gabpa	gabpa	GABPA
NULL	MXI1	MXI1
mycb	myca	MYC
sox17	sox17	SOX17
NULL	tcf3b	TCF3
si:dkey-110c1.10	tcf3b	TCF3
tp63	NULL	TP63
nr2f2	nr2f2	NR2F2
spi1a	spi1a	SPI1
epas1b	epas1b	EPAS1
NULL	MAX_(1_of_many)	MAX
rfx2	rfx2	RFX2
esr2b	NULL	ESR2
whsc1	whsc1	WHSC1
NULL	nr3c1	NR3C1
NULL	pax6a	PAX6
egr1	egr1	EGR1
ncor1	NULL	NCOR1
maff	maff	MAFF
esr1	esr1	ESR1
srebf1	srebf1	SREBF1
NULL	fosab	FOS
nfe2	nfe2	NFE2
stat4	stat4	STAT4
fosl2	fosl2	FOSL2
tal1	tal1	TAL1
irf4a	irf4a	IRF4
cebpd	cebpd	CEBPD
mybl2b	mybl2b	MYBL2
prdm1a	prdm1a	PRDM1
NULL	tcf3a	TCF3
foxp2	foxp2	FOXP2
tcf12	tcf12	TCF12
thap1	thap1	THAP1
sox2	sox2	SOX2
creb1a	creb1a_(1_of_many)	CREB1
six5	six5	SIX5
ncor2	ncor2	NCOR2
pax5	pax5	PAX5

## checking for correlation with extrapolated set
grep -wiF mz.gene.s102.69 ../../GTRD_TFTGcoexpCheck/hg38_tftgtfbs_merged_07.extrap.output | cut -f1 | sort -u | paste -sd " "
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3 # local
for i in mz.gene.s10.123 mz.gene.s11.168 mz.gene.s123.8 mz.gene.s127.35 mz.gene.s141.31 mz.gene.s152.11 mz.gene.s153.32 mz.gene.s17.142 mz.gene.s249.26 mz.gene.s25.137 mz.gene.s321.14 mz.gene.s60.4 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
NULL	BHLHE40_(1_of_many)	BHLHE40
NULL	fosl1a	FOSL1
NULL	NFATC1_(1_of_many)	NFATC1
NULL	jun	JUN
atf3	atf3	ATF3
NULL	e2f4	E2F4
nfatc1	nfatc1	NFATC1
foxm1	foxm1	FOXM1
fli1b	fli1b	ERG
gabpa	gabpa	GABPA
spi1a	spi1a	SPI1
fosl2	fosl2	FOSL2

# TF outputs per set - those with 'XXXXXXX' are extrapolated but not found in the non-filtered fimo output
# hg38_tftgtfbs_merged_07.extrap.output
NULL	BHLHE40_(1_of_many)	BHLHE40
NULL	fosl1a	FOSL1
NULL	NFATC1_(1_of_many)	NFATC1
NULL	jun	JUN 	XXXXXXX
atf3	atf3	ATF3
NULL	e2f4	E2F4 	XXXXXXX
nfatc1	nfatc1	NFATC1
foxm1	foxm1	FOXM1
fli1b	fli1b	ERG
gabpa	gabpa	GABPA
spi1a	spi1a	SPI1
fosl2	fosl2	FOSL2

# Filtered (Qval<0.05) CS_matQualPval_mz_results.out
NULL	GATA2_(1_of_many)	GATA2
klf5b	klf5b	KLF5
NULL	egr2b	EGR2
nrf1	nrf1	NRF1
tp53	tp53	TP53
NULL	MXI1	MXI1
NULL	tcf3b	TCF3
si:dkey-110c1.10	tcf3b	TCF3
egr1	egr1	EGR1
ncor1	NULL	NCOR1
NULL	tcf3a	TCF3

# Non-filtered (no Qval<0.05) CS_matQualPval_mz_results.out
smad4a	SMAD4_(1_of_many)	SMAD4
pou5f3	pou5f3	POU5F1
stat2	stat2	STAT2
NULL	BHLHE40_(1_of_many)	BHLHE40
kdm5ba	kdm5ba	KDM5B
gata2a	gata2a	GATA2
NULL	GATA2_(1_of_many)	GATA2
zeb1a	zeb1a	ZEB1
NULL	fosl1a	FOSL1
tbx21	tbx21	TBX21
nkx2.4a	nkx2.4a	NKX2-1
kdm5bb	kdm5bb	KDM5B
NULL	NFATC1_(1_of_many)	NFATC1
foxp1b	FOXP1	FOXP1
NULL	TCF12_(1_of_many)	TCF12
CABZ01066696.1	NULL	LYL1
klf5b	klf5b	KLF5
NULL	eomesb	EOMES
tfap2a	tfap2a	TFAP2A
atf3	atf3	ATF3
NULL	foxa2	FOXA2
prdm1c	PRDM1_(1_of_many)	PRDM1
creb1b	creb1a_(1_of_many)	CREB1
tfap2c	tfap2c	TFAP2C
sp1	sp1	SP1
nfatc1	nfatc1	NFATC1
NULL	nfic	NFIC
NULL	foxa1	FOXA1
NULL	egr2b	EGR2
foxm1	foxm1	FOXM1
nrf1	nrf1	NRF1
tp53	tp53	TP53
eomesa	eomesa	EOMES
irf4b	irf4b	IRF4
gata1a	gata1a	GATA1
zeb1b	zeb1b	ZEB1
pou2f2a	POU2F2_(1_of_many)	POU2F2
pax6b	pax6b	PAX6
esr2a	esr2a	ESR2
NULL	USF2_(1_of_many)	USF2
fli1b	fli1b	ERG
pou2f2a	POU2F2_(1_of_many)	POU2F2
gabpa	gabpa	GABPA
NULL	MXI1	MXI1
mycb	myca	MYC
sox17	sox17	SOX17
NULL	tcf3b	TCF3
si:dkey-110c1.10	tcf3b	TCF3
tp63	NULL	TP63
nr2f2	nr2f2	NR2F2
spi1a	spi1a	SPI1
epas1b	epas1b	EPAS1
NULL	MAX_(1_of_many)	MAX
rfx2	rfx2	RFX2
esr2b	NULL	ESR2
whsc1	whsc1	WHSC1
NULL	nr3c1	NR3C1
NULL	pax6a	PAX6
egr1	egr1	EGR1
ncor1	NULL	NCOR1
maff	maff	MAFF
esr1	esr1	ESR1
srebf1	srebf1	SREBF1
NULL	fosab	FOS
nfe2	nfe2	NFE2
stat4	stat4	STAT4
fosl2	fosl2	FOSL2
tal1	tal1	TAL1
irf4a	irf4a	IRF4
cebpd	cebpd	CEBPD
mybl2b	mybl2b	MYBL2
prdm1a	prdm1a	PRDM1
NULL	tcf3a	TCF3
foxp2	foxp2	FOXP2
tcf12	tcf12	TCF12
thap1	thap1	THAP1
sox2	sox2	SOX2
creb1a	creb1a_(1_of_many)	CREB1
six5	six5	SIX5
ncor2	ncor2	NCOR2
pax5	pax5	PAX5

awk '($8 < 0.05 )' OFS='\t' CS_default_mz_results.out | grep -wiF mz.gene.s102.69 | cut -f1 | sort -u | paste -sd " "
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3 # local
for i in mz.gene.s10.82 mz.gene.s136.36 mz.gene.s171.22 mz.gene.s17.45 mz.gene.s186.6 mz.gene.s26.74 mz.gene.s31.45 mz.gene.s31.46 mz.gene.s4.9 mz.gene.s50.74 mz.gene.s6.87 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
NULL	GATA2_(1_of_many)	GATA2
klf5b	klf5b	KLF5
NULL	egr2b	EGR2
nrf1	nrf1	NRF1
tp53	tp53	TP53
NULL	MXI1	MXI1
NULL	tcf3b	TCF3
si:dkey-110c1.10	tcf3b	TCF3
egr1	egr1	EGR1
ncor1	NULL	NCOR1
NULL	tcf3a	TCF3

# JASPAR
awk '($8 < 0.05 )' OFS='\t' JASPAR_matQualPval_mz_results.out | grep -wiF mz.gene.s102.69 | cut -f1 | sort -u | paste -sd " "
RREB1 ZNF263 ZNF384 #matqual filt
grep -wiF mz.gene.s102.69 JASPAR_matQualPval_mz_results.out | cut -f1 | sort -u | paste -sd " "
#matqual no filt
Ahr::Arnt Ar Arid3a Arnt Arntl ASCL1 Ascl2 Atf3 ATF4 Atoh1 BACH2 BARHL2 BATF3 BATF::JUN Bcl6 Bhlha15 BHLHE40 BHLHE41 CEBPA CEBPB CEBPD CEBPE CLOCK CREB1 CTCF CTCFL DMRT3 Dux DUX4 E2F6 EGR1 EGR2 EGR3 ELF5 ELK4 EOMES ERF ESR1 ESR2 Esrra ESRRB Esrrg ETS1 ETV2 ETV3 ETV6 FIGLA FOS FOSB::JUN FOS::JUN FOS::JUNB FOS::JUND FOSL1 FOSL2 FOXA1 Foxa2 FOXB1 FOXC1 FOXC2 FOXD1 FOXD2 Foxd3 FOXG1 FOXH1 FOXI1 Foxj2 Foxj3 FOXK1 FOXK2 FOXL1 Foxo1 FOXO3 FOXO4 FOXO6 FOXP1 FOXP2 FOXP3 Foxq1 Gabpa Gata4 GCM1 GCM2 Gfi1 Gfi1b GRHL1 GRHL2 Hes2 HES5 Hic1 HIC2 HLTF Hnf4a HNF4G HOXA5 HOXB13 HOXD13 Hoxd9 Id2 ID4 IRF1 IRF2 IRF7 IRF8 IRF9 JDP2 JDP2_v2 JUN JUNB JUND JUND_v2 JUN::JUNB JUN_v2 Klf1 KLF16 KLF4 KLF5 KLF9 LEF1 LIN54 LMX1A Mafb MAFF MAFG MAFK MAF::NFE2 MAX MAX::MYC MEF2A MEF2C MEIS1 MEIS2 MEIS3 MGA MITF Mlxip MLXIPL MNT MSC MTF1 MXI1 MYB MYC MYCN MYF6 Myod1 Myog MZF1 NEUROD1 Neurog1 NEUROG2 NFAT5 NFATC2 NFE2 Nfe2l2 NFIA NFIC NFIX NFYA NHLH1 NKX2-3 Nkx2-5_v2 NKX2-8 Nkx3-1 Npas2 NR2C2 Nr2e1 Nr2e3 NR2F1 NR2F2 Nr2f6 Nr2f6_v2 NR3C1 NR3C2 NR4A1 NR4A2 Nr5a2 NRF1 NRL OLIG1 OLIG2 OLIG3 Pax2 PAX5 PBX1 PBX2 PBX3 PKNOX1 PKNOX2 POU1F1 POU2F1 POU2F2 Pou2f3 POU3F1 POU3F2 POU3F3 POU3F4 POU4F1 POU4F2 POU5F1 POU5F1B PPARG PRDM1 RARA RARA_v2 Rarb Rarb_v2 Rarg Rarg_v2 RBPJ REST Rfx1 Rhox11 RORA RORA_v2 RREB1 RUNX1 Rxra RXRB RXRG SCRT1 SCRT2 SIX1 SIX2 Smad4 SNAI2 Sox1 SOX10 Sox11 SOX13 SOX15 Sox17 Sox2 SOX21 Sox3 SOX4 Sox5 Sox6 SOX8 SP1 SP3 SP8 SPI1 SPIB SREBF1 Srebf1_v2 SREBF2 SREBF2_v2 STAT1 STAT3 Stat4 Stat6 T TBR1 TBX1 TBX15 TBX19 TBX2 TBX20 TBX21 TBX4 TBX5 Tcf12 Tcf21 TCF3 TCF4 Tcf7 TCF7L1 TCF7L2 TFAP2A TFAP2B_v2 TFAP2C_v2 TFAP4 TFE3 TFEB TFEC TGIF1 TGIF2 THAP1 TP63 TP73 TWIST1 Twist2 USF1 USF2 VDR YY1 ZBTB18 ZBTB7A ZEB1 Zfx ZIC1 ZIC3 ZIC4 ZNF143 ZNF24 ZNF263 ZNF354C ZNF384 ZSCAN4
awk '($8 < 0.05 )' OFS='\t' JASPAR_default_mz_results.out | grep -wiF mz.gene.s102.69 | cut -f1 | sort -u | paste -sd " "
RREB1 ZNF263 ZNF384 #default filt
grep -wiF mz.gene.s102.69 JASPAR_default_mz_results.out | cut -f1 | sort -u | paste -sd " "
# default no filt
Alx1 ARNT::HIF1A Arx ASCL1 Ascl2 Atf1 Atf3 BARHL2 BATF3 CREB3 Crem Crx CTCF CTCFL CUX1 CUX2 DMRT3 DUX4 DUXA E2F6 E2F8 EGR1 EGR2 EGR3 EGR4 ELF1 ELF4 ELK1 ELK3 ELK4 ERF ERG ESR1 ESR2 Esrra ESRRB Esrrg ETS1 ETV1 ETV2 ETV3 ETV4 ETV6 EWSR1-FLI1 FEV FIGLA FLI1 FOS FOS::JUN_v2 FOXB1 FOXC1 FOXC2 Foxd3 Foxj3 Gabpa GLI2 GLIS2 GSC GSC2 HIF1A Hmx2 HNF1A HNF1B HNF4G Hoxd3 HSF2 HSF4 ID4 IRF3 IRF4 IRF5 IRF8 IRF9 JUND JUN_v2 KLF13 KLF4 KLF9 Lhx3 LIN54 LMX1A MAX::MYC mix-a MSC MXI1 MYBL2 MYC MYCN MYF6 Myod1 Myog NEUROD1 NEUROD2 NFATC1 NFATC2 NFATC3 NFIA NFIC::TLX1 NFIX NFKB1 NFYA Nkx2-5 Nkx2-5_v2 NKX2-8 Nkx3-1 NKX3-2 NR1H2::RXRA Nr1h3::Rxra NR1H4 NR2C2 NR2F1 NR2F2 Nr2f6 Nr2f6_v2 NR4A1 NR4A2 NRL ONECUT1 ONECUT2 ONECUT3 OTX1 OTX2 PAX3 Pax6 PAX7 PBX1 PBX2 PBX3 PHOX2A Phox2b Pitx1 PITX3 POU1F1 POU2F1 POU3F1 POU3F2 POU3F3 PPARA::RXRA PPARG Pparg::Rxra PROP1 PROX1 RARA::RXRA RARA::RXRG Rarb Rarb_v2 Rarg REL Rfx1 RHOXF1 RORA RORC RREB1 RUNX2 Rxra RXRA::VDR RXRB RXRG SCRT1 SIX1 SIX2 SNAI2 SOX10 Sox11 SOX15 SOX21 Sox3 SOX4 SP1 SP2 SPDEF Spz1 SREBF2 SRF STAT1::STAT2 TAL1::TCF3 TBP Tcf12 Tcf21 TCF3 TCF4 Tcf7 TCF7L1 TFAP2A TFAP4 TFDP1 TGIF1 TWIST1 USF1 USF2 VDR ZBTB18 ZBTB7A ZBTB7B ZBTB7C ZEB1 Zfx ZIC1 ZIC3 ZIC4 ZNF263 ZNF384 Znf423 ZNF740 ZSCAN4


#### sws1 - ab.gene.s2279.1
#filtered for q-val
awk '($8 < 0.05 )' OFS='\t' CS_matQualPval_ab_results.out | grep -wiF ab.gene.s2279.1 | cut -f1 | sort -u | paste -sd " "
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3 # local
for i in ab.gene.s179.34 ab.gene.s271.28 ab.gene.s373.3 ab.gene.s395.10 ab.gene.s397.20 ab.gene.s9.66 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
klf5b	klf5b	KLF5
tp53	tp53	TP53
NULL	egr2b	EGR2
NULL	GATA2_(1_of_many)	GATA2
NULL	sp2	SP2
gata2a	gata2a	GATA2

#non-filtered for q-val
grep -wiF ab.gene.s2279.1 CS_matQualPval_ab_results.out | cut -f1 | sort -u | paste -sd " "
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3 # local
for i in ab.gene.s106.29 ab.gene.s107.48 ab.gene.s108.30 ab.gene.s111.11 ab.gene.s1.112 ab.gene.s111.27 ab.gene.s112.33 ab.gene.s117.29 ab.gene.s120.15 ab.gene.s131.26 ab.gene.s13.28 ab.gene.s14.29 ab.gene.s144.13 ab.gene.s151.11 ab.gene.s15.20 ab.gene.s156.8 ab.gene.s157.18 ab.gene.s167.12 ab.gene.s169.15 ab.gene.s179.34 ab.gene.s186.12 ab.gene.s1.90 ab.gene.s191.12 ab.gene.s203.2 ab.gene.s2.112 ab.gene.s21.28 ab.gene.s216.13 ab.gene.s220.18 ab.gene.s224.3 ab.gene.s23.59 ab.gene.s246.39 ab.gene.s248.18 ab.gene.s271.28 ab.gene.s305.3 ab.gene.s31.19 ab.gene.s314.16 ab.gene.s32.25 ab.gene.s338.5 ab.gene.s34.62 ab.gene.s35.64 ab.gene.s35.65 ab.gene.s373.3 ab.gene.s382.5 ab.gene.s3.93 ab.gene.s393.13 ab.gene.s395.10 ab.gene.s397.20 ab.gene.s40.36 ab.gene.s40.58 ab.gene.s41.81 ab.gene.s42.64 ab.gene.s458.1 ab.gene.s458.13 ab.gene.s469.13 ab.gene.s47.42 ab.gene.s4.83 ab.gene.s486.7 ab.gene.s49.24 ab.gene.s49.33 ab.gene.s53.7 ab.gene.s56.17 ab.gene.s568.5 ab.gene.s57.2 ab.gene.s6.159 ab.gene.s646.1 ab.gene.s66.44 ab.gene.s66.79 ab.gene.s6.86 ab.gene.s68.8 ab.gene.s693.3 ab.gene.s69.50 ab.gene.s70.5 ab.gene.s7.47 ab.gene.s82.17 ab.gene.s82.47 ab.gene.s827.1 ab.gene.s8.71 ab.gene.s87.60 ab.gene.s9.18 ab.gene.s95.25 ab.gene.s9.66 ab.gene.s97.2 ab.gene.s986.1 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
prdm1a	prdm1a	PRDM1
tcf12	tcf12	TCF12
whsc1	whsc1	WHSC1
NULL	foxa2	FOXA2
irf4a	irf4a	IRF4
atf3	atf3	ATF3
six5	six5	SIX5
stat2	stat2	STAT2
kdm5bb	kdm5bb	KDM5B
esr2a	esr2a	ESR2
NULL	smad4a	SMAD4
NULL	EBF1_(1_of_many)	EBF1
NULL	tcf3b	TCF3
nr2f2	nr2f2	NR2F2
hif1ab	hif-1a	HIF1A
mybl2b	mybl2b	MYBL2
NULL	pax6a	PAX6
ascl1a	ascl1a	ASCL1
maff	maff	MAFF
klf5b	klf5b	KLF5
epas1b	epas1b	EPAS1
cebpd	cebpd	CEBPD
NULL	e2f4	E2F4
CABZ01066696.1	NULL	LYL1
NULL	NR2F2_(1_of_many)	NR2F2
zeb1a	zeb1a	ZEB1
NULL	NFKB1	NFKB1
irf4b	irf4b	IRF4
usf1	usf1	USF1
zeb1b	zeb1b	ZEB1
NULL	nfic	NFIC
rfx2	rfx2	RFX2
tp53	tp53	TP53
pou2f2a	POU2F2_(1_of_many)	POU2F2
tp63	NULL	TP63
smad4a	SMAD4_(1_of_many)	SMAD4
stat1a	stat1a	STAT1
tfap2a	tfap2a	TFAP2A
creb1a	creb1a_(1_of_many)	CREB1
stat5a	stat5a	STAT5A
stat5a	stat5a	STAT5A
NULL	egr2b	EGR2
nfatc1	nfatc1	NFATC1
NULL	fosl1a	FOSL1
foxp1b	FOXP1	FOXP1
NULL	GATA2_(1_of_many)	GATA2
NULL	sp2	SP2
egr1	egr1	EGR1
NULL	nr3c1	NR3C1
atf2	atf2	ATF2
prdm1c	PRDM1_(1_of_many)	PRDM1
NULL	TCF12_(1_of_many)	TCF12
spi1a	spi1a	SPI1
gabpa	gabpa	GABPA
nkx2.4a	nkx2.4a	NKX2-1
tbx21	tbx21	TBX21
fosl2	fosl2	FOSL2
sox17	sox17	SOX17
mycb	myca	MYC
NULL	irf1b	IRF1
foxm1	foxm1	FOXM1
creb1b	creb1a_(1_of_many)	CREB1
fli1b	fli1b	ERG
esr1	esr1	ESR1
esr2b	NULL	ESR2
kdm5ba	kdm5ba	KDM5B
NULL	BHLHE40_(1_of_many)	BHLHE40
NULL	fosab	FOS
pou2f2a	POU2F2_(1_of_many)	POU2F2
NULL	eomesb	EOMES
srebf1	srebf1	SREBF1
NULL	USF2_(1_of_many)	USF2
nfe2	nfe2	NFE2
pbx3b	pbx3a	PBX3
ncor2	ncor2	NCOR2
tfap2c	tfap2c	TFAP2C
pax6b	pax6b	PAX6
eomesa	eomesa	EOMES
gata1a	gata1a	GATA1
NULL	arntl1a	ARNTL
gata2a	gata2a	GATA2
sp4	sp4	SP4
NULL	foxa1	FOXA1

# JASPAR
awk '($8 < 0.05 )' OFS='\t' JASPAR_matQualPval_ab_results.out | grep -wiF ab.gene.s2279.1 | cut -f1 | sort -u | paste -sd " "
NR2C2 Nr2f6 RREB1 Rxra RXRB RXRG ZNF263 ZNF384 #matqual filt
grep -wiF ab.gene.s2279.1 JASPAR_matQualPval_ab_results.out | cut -f1 | sort -u | paste -sd " "
#matqual no filt
Ahr::Arnt Ar Arid3a Arnt ARNT::HIF1A Arntl ASCL1 Ascl2 Atf3 ATF4 Atoh1 Bach1::Mafk BACH2 BARHL2 BATF3 BATF::JUN Bcl6 Bhlha15 BHLHE40 BHLHE41 CDX1 CEBPA CEBPB CEBPD CEBPE CLOCK CREB1 CTCF CTCFL Ddit3::Cebpa Dux DUX4 E2F6 EGR1 EGR2 EGR3 ELK4 EOMES ERF ESR1 ESR2 Esrra ESRRB Esrrg ETS1 ETV2 ETV3 ETV6 EWSR1-FLI1 FIGLA FOS FOSB::JUN FOSB::JUNB FOSB::JUNB_v2 FOS::JUN FOS::JUNB FOS::JUND FOS::JUN_v2 FOSL1 FOSL1::JUN FOSL1::JUNB FOSL1::JUND FOSL1::JUN_v2 FOSL2 FOSL2::JUN FOSL2::JUNB FOSL2::JUNB_v2 FOSL2::JUND FOSL2::JUND_v2 FOSL2::JUN_v2 FOXA1 Foxa2 FOXB1 FOXC1 FOXC2 FOXD1 FOXD2 Foxd3 FOXG1 FOXH1 FOXI1 Foxj2 Foxj3 FOXK1 FOXK2 FOXL1 Foxo1 FOXO3 FOXO4 FOXO6 FOXP1 FOXP2 FOXP3 Foxq1 Gabpa GATA1::TAL1 Gata4 GCM1 GCM2 Gfi1 Gfi1b GRHL1 GRHL2 Hand1::Tcf3 Hes2 HES5 HEY1 Hic1 HIC2 HLTF Hnf4a HNF4G HOXA13 HOXA5 Hoxa9 HOXB13 HOXD13 Hoxd9 Id2 ID4 IRF1 IRF2 IRF5 IRF7 IRF9 JDP2 JDP2_v2 JUN JUNB JUND JUND_v2 JUN::JUNB JUN_v2 Klf1 KLF16 KLF4 KLF5 KLF9 LEF1 LIN54 LMX1A Mafb MAFF MAFG MAFG::NFE2L1 MAFK MAF::NFE2 MAX MAX::MYC MEF2A MEF2C MEF2D MEIS1 MEIS2 MEIS3 MGA MITF Mlxip MLXIPL MNT MSC MTF1 MXI1 MYB MYC MYCN MYF6 Myod1 Myog MZF1 NEUROD1 Neurog1 NEUROG2 NFAT5 NFATC2 NFE2 Nfe2l2 NFIA NFIC NFIC::TLX1 NFIX NFYA NFYB NHLH1 NKX2-3 Nkx2-5_v2 NKX2-8 Nkx3-1 Npas2 NR1H2::RXRA Nr1h3::Rxra NR2C2 Nr2e1 Nr2e3 NR2F1 NR2F2 Nr2f6 Nr2f6_v2 NR3C1 NR3C2 NR4A1 NR4A2 Nr5a2 NRF1 NRL OLIG1 OLIG2 OLIG3 Pax2 PAX5 PBX1 PBX2 PBX3 PKNOX1 PKNOX2 POU1F1 POU2F1 POU2F2 Pou2f3 POU3F1 POU3F2 POU3F3 POU3F4 POU4F1 POU4F2 POU5F1 POU5F1B Pou5f1::Sox2 PPARA::RXRA PPARG Pparg::Rxra PRDM1 RARA RARA::RXRA RARA::RXRG RARA_v2 Rarb Rarb_v2 Rarg Rarg_v2 RBPJ REST Rfx1 Rhox11 RORA RREB1 RUNX1 Rxra RXRB RXRG SCRT1 SCRT2 SIX1 SIX2 SMAD2::SMAD3::SMAD4 Smad4 SNAI2 Sox1 SOX10 Sox11 SOX13 SOX15 Sox17 Sox2 SOX21 Sox3 SOX4 Sox5 Sox6 SOX8 SOX9 SP1 SP3 SP8 SPI1 SPIB SREBF1 Srebf1_v2 SREBF2 SREBF2_v2 STAT1 STAT1::STAT2 STAT3 Stat4 Stat5a::Stat5b Stat6 T TAL1::TCF3 TBR1 TBX1 TBX15 TBX19 TBX2 TBX20 TBX21 TBX4 TBX5 Tcf12 Tcf21 TCF3 TCF4 Tcf7 TCF7L1 TCF7L2 TFAP2A TFAP2A_v2 TFAP2A_v3 TFAP2B TFAP2B_v2 TFAP2B_v3 TFAP2C TFAP2C_v2 TFAP2C_v3 TFAP4 TFE3 TFEB TFEC TGIF1 TGIF2 TP53 TP63 TP73 TWIST1 Twist2 USF1 USF2 VDR YY1 ZBTB18 ZBTB7A ZEB1 Zfx ZIC1 ZIC3 ZIC4 ZNF143 ZNF24 ZNF263 ZNF354C ZNF384 ZSCAN4
awk '($8 < 0.05 )' OFS='\t' JASPAR_default_ab_results.out | grep -wiF ab.gene.s2279.1 | cut -f1 | sort -u | paste -sd " "
NR2C2 Nr2f6 RREB1 Rxra RXRB RXRG ZNF263 ZNF384 #default filt
grep -wiF ab.gene.s2279.1 JASPAR_default_ab_results.out | cut -f1 | sort -u | paste -sd " "
# default no filt
Ahr::Arnt Ar Arid3a Arnt ARNT::HIF1A Arntl ASCL1 Ascl2 Atf3 ATF4 Atoh1 Bach1::Mafk BACH2 BARHL2 BATF3 BATF::JUN Bcl6 Bhlha15 BHLHE40 BHLHE41 CDX1 CEBPA CEBPB CEBPD CEBPE CLOCK CREB1 CTCF CTCFL Ddit3::Cebpa Dux DUX4 E2F6 EGR1 EGR2 EGR3 ELK4 EOMES ERF ESR1 ESR2 Esrra ESRRB Esrrg ETS1 ETV2 ETV3 ETV6 EWSR1-FLI1 FIGLA FOS FOSB::JUN FOSB::JUNB FOSB::JUNB_v2 FOS::JUN FOS::JUNB FOS::JUND FOS::JUN_v2 FOSL1 FOSL1::JUN FOSL1::JUNB FOSL1::JUND FOSL1::JUN_v2 FOSL2 FOSL2::JUN FOSL2::JUNB FOSL2::JUNB_v2 FOSL2::JUND FOSL2::JUND_v2 FOSL2::JUN_v2 FOXA1 Foxa2 FOXB1 FOXC1 FOXC2 FOXD1 FOXD2 Foxd3 FOXG1 FOXH1 FOXI1 Foxj2 Foxj3 FOXK1 FOXK2 FOXL1 Foxo1 FOXO3 FOXO4 FOXO6 FOXP1 FOXP2 FOXP3 Foxq1 Gabpa GATA1::TAL1 Gata4 GCM1 GCM2 Gfi1 Gfi1b GRHL1 GRHL2 Hand1::Tcf3 Hes2 HES5 HEY1 Hic1 HIC2 HLTF Hnf4a HNF4G HOXA13 HOXA5 Hoxa9 HOXB13 HOXD13 Hoxd9 Id2 ID4 IRF1 IRF2 IRF5 IRF7 IRF9 JDP2 JDP2_v2 JUN JUNB JUND JUND_v2 JUN::JUNB JUN_v2 Klf1 KLF16 KLF4 KLF5 KLF9 LEF1 LIN54 LMX1A Mafb MAFF MAFG MAFG::NFE2L1 MAFK MAF::NFE2 MAX MAX::MYC MEF2A MEF2C MEF2D MEIS1 MEIS2 MEIS3 MGA MITF Mlxip MLXIPL MNT MSC MTF1 MXI1 MYB MYC MYCN MYF6 Myod1 Myog MZF1 NEUROD1 Neurog1 NEUROG2 NFAT5 NFATC2 NFE2 Nfe2l2 NFIA NFIC NFIC::TLX1 NFIX NFYA NFYB NHLH1 NKX2-3 Nkx2-5_v2 NKX2-8 Nkx3-1 Npas2 NR1H2::RXRA Nr1h3::Rxra NR2C2 Nr2e1 Nr2e3 NR2F1 NR2F2 Nr2f6 Nr2f6_v2 NR3C1 NR3C2 NR4A1 NR4A2 Nr5a2 NRF1 NRL OLIG1 OLIG2 OLIG3 Pax2 PAX5 PBX1 PBX2 PBX3 PKNOX1 PKNOX2 POU1F1 POU2F1 POU2F2 Pou2f3 POU3F1 POU3F2 POU3F3 POU3F4 POU4F1 POU4F2 POU5F1 POU5F1B Pou5f1::Sox2 PPARA::RXRA PPARG Pparg::Rxra PRDM1 RARA RARA::RXRA RARA::RXRG RARA_v2 Rarb Rarb_v2 Rarg Rarg_v2 RBPJ REST Rfx1 Rhox11 RORA RREB1 RUNX1 Rxra RXRB RXRG SCRT1 SCRT2 SIX1 SIX2 SMAD2::SMAD3::SMAD4 Smad4 SNAI2 Sox1 SOX10 Sox11 SOX13 SOX15 Sox17 Sox2 SOX21 Sox3 SOX4 Sox5 Sox6 SOX8 SOX9 SP1 SP3 SP8 SPI1 SPIB SREBF1 Srebf1_v2 SREBF2 SREBF2_v2 STAT1 STAT1::STAT2 STAT3 Stat4 Stat5a::Stat5b Stat6 T TAL1::TCF3 TBR1 TBX1 TBX15 TBX19 TBX2 TBX20 TBX21 TBX4 TBX5 Tcf12 Tcf21 TCF3 TCF4 Tcf7 TCF7L1 TCF7L2 TFAP2A TFAP2A_v2 TFAP2A_v3 TFAP2B TFAP2B_v2 TFAP2B_v3 TFAP2C TFAP2C_v2 TFAP2C_v3 TFAP4 TFE3 TFEB TFEC TGIF1 TGIF2 TP53 TP63 TP73 TWIST1 Twist2 USF1 USF2 VDR YY1 ZBTB18 ZBTB7A ZEB1 Zfx ZIC1 ZIC3 ZIC4 ZNF143 ZNF24 ZNF263 ZNF354C ZNF384 ZSCAN4


##### sws1 - nb.gene.s1.386
#filtered for q-val
awk '($8 < 0.05 )' OFS='\t' CS_matQualPval_nb_results.out | grep -wiF nb.gene.s1.386 | cut -f1 | sort -u | paste -sd " "
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3 # local
# NONE

#non-filtered for q-val
grep -wiF nb.gene.s1.386 CS_matQualPval_nb_results.out | cut -f1 | sort -u | paste -sd " "
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3 # local
for i in nb.gene.s0.260 nb.gene.s0.327 nb.gene.s0.45 nb.gene.s10.26 nb.gene.s10.48 nb.gene.s11.136 nb.gene.s11.212 nb.gene.s1.223 nb.gene.s124.7 nb.gene.s12.71 nb.gene.s130.21 nb.gene.s13.31 nb.gene.s14.143 nb.gene.s15.129 nb.gene.s15.5 nb.gene.s16.63 nb.gene.s168.5 nb.gene.s17.179 nb.gene.s17.181 nb.gene.s194.14 nb.gene.s195.22 nb.gene.s197.8 nb.gene.s199.9 nb.gene.s21.97 nb.gene.s2.243 nb.gene.s2.269 nb.gene.s24.167 nb.gene.s247.17 nb.gene.s24.9 nb.gene.s2.71 nb.gene.s28.105 nb.gene.s28.92 nb.gene.s31.113 nb.gene.s3.156 nb.gene.s32.121 nb.gene.s33.134 nb.gene.s34.59 nb.gene.s35.47 nb.gene.s4.42 nb.gene.s48.52 nb.gene.s53.55 nb.gene.s54.38 nb.gene.s55.94 nb.gene.s58.59 nb.gene.s6.225 nb.gene.s63.51 nb.gene.s6.381 nb.gene.s6.425 nb.gene.s6.449 nb.gene.s66.39 nb.gene.s73.10 nb.gene.s73.73 nb.gene.s7.97 nb.gene.s8.350 nb.gene.s9.110 nb.gene.s9.280 nb.gene.s96.12 nb.gene.s9.83 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
pax6b	pax6b	PAX6
tcf12	tcf12	TCF12
NULL	NR2F2_(1_of_many)	NR2F2
cebpd	cebpd	CEBPD
irf4a	irf4a	IRF4
pax5	pax5	PAX5
CABZ01066696.1	NULL	LYL1
foxm1	foxm1	FOXM1
zeb1b	zeb1b	ZEB1
sp4	sp4	SP4
NULL	arntl1a	ARNTL
klf5b	klf5b	KLF5
smad4a	SMAD4_(1_of_many)	SMAD4
nr2f2	nr2f2	NR2F2
NULL	TCF12_(1_of_many)	TCF12
NULL	smad4a	SMAD4
NULL	GATA2_(1_of_many)	GATA2
NULL	tcf3b	TCF3
si:dkey-110c1.10	tcf3b	TCF3
stat2	stat2	STAT2
stat4	stat4	STAT4
tp53	tp53	TP53
zeb1a	zeb1a	ZEB1
prdm1a	prdm1a	PRDM1
egr1	egr1	EGR1
NULL	nr3c1	NR3C1
esr2a	esr2a	ESR2
NULL	BHLHE40_(1_of_many)	BHLHE40
esr1	esr1	ESR1
NULL	fosl1a	FOSL1
atf3	atf3	ATF3
NULL	foxa2	FOXA2
kdm5bb	kdm5bb	KDM5B
NULL	pax6a	PAX6
pou2f2a	POU2F2_(1_of_many)	POU2F2
tfap2a	tfap2a	TFAP2A
NULL	NFATC1_(1_of_many)	NFATC1
sox17	sox17	SOX17
tfap2c	tfap2c	TFAP2C
NULL	fosab	FOS
irf4b	irf4b	IRF4
srebf2	srebf2	SREBF2
stat5a	stat5a	STAT5A
NULL	rest	REST
kdm5ba	kdm5ba	KDM5B
fosl2	fosl2	FOSL2
gata2a	gata2a	GATA2
gata1a	gata1a	GATA1
foxp1b	FOXP1	FOXP1
rfx2	rfx2	RFX2
NULL	foxa1	FOXA1
nkx2.4a	nkx2.4a	NKX2-1
NULL	irf1b	IRF1
nfe2	nfe2	NFE2
whsc1	whsc1	WHSC1
hif1ab	hif-1a	HIF1A
pou2f2a	POU2F2_(1_of_many)	POU2F2
esr2b	NULL	ESR2

# JASPAR
awk '($8 < 0.05 )' OFS='\t' JASPAR_matQualPval_nb_results.out | grep -wiF nb.gene.s1.386 | cut -f1 | sort -u | paste -sd " "
NR2C2 Nr2f6 Rxra RXRB RXRG #matqual filt
grep -wiF nb.gene.s1.386 JASPAR_matQualPval_nb_results.out | cut -f1 | sort -u | paste -sd " "
#matqual no filt
Ahr::Arnt Ar Arid3a Arnt ARNT::HIF1A Arntl ASCL1 Ascl2 Atf3 ATF4 ATF7 Atoh1 Bach1::Mafk BACH2 BATF3 BATF::JUN Bcl6 Bhlha15 BHLHE40 BHLHE41 CEBPA CEBPB CEBPD CEBPE CLOCK CREB1 Creb5 CTCF CTCFL Ddit3::Cebpa Dux DUX4 E2F6 EGR2 ELF5 ELK4 EOMES ERF ESR1 ESR2 Esrra ESRRB Esrrg ETS1 ETV2 ETV3 ETV6 FIGLA FOS FOSB::JUN FOSB::JUNB FOSB::JUNB_v2 FOS::JUN FOS::JUNB FOS::JUND FOS::JUN_v2 FOSL1 FOSL1::JUN FOSL1::JUNB FOSL1::JUND FOSL1::JUN_v2 FOSL2 FOSL2::JUN FOSL2::JUNB FOSL2::JUNB_v2 FOSL2::JUND FOSL2::JUND_v2 FOSL2::JUN_v2 FOXA1 Foxa2 FOXB1 FOXC1 FOXC2 FOXD1 FOXD2 Foxd3 FOXG1 FOXH1 FOXI1 Foxj2 Foxj3 FOXK1 FOXK2 FOXL1 Foxo1 FOXO3 FOXO4 FOXO6 FOXP1 FOXP2 FOXP3 Foxq1 Gabpa GATA1::TAL1 Gata4 GCM1 GCM2 Gfi1 Gfi1b GRHL1 GRHL2 Hand1::Tcf3 Hes2 HES5 HIC2 HLTF Hnf4a HNF4G HOXA5 Hoxa9 HOXB13 HOXD13 Hoxd9 Id2 ID4 IRF1 IRF2 IRF5 IRF7 IRF9 JDP2 JDP2_v2 JUN JUNB JUND JUND_v2 JUN::JUNB JUN_v2 Klf1 KLF4 KLF5 KLF9 LEF1 LIN54 Mafb MAFF MAFG MAFG::NFE2L1 MAFK MAF::NFE2 MAX MAX::MYC MEF2A MEF2C MEIS1 MEIS2 MEIS3 MGA MITF Mlxip MLXIPL MNT MSC MTF1 MXI1 MYB MYC MYF6 Myod1 Myog MZF1 NEUROD1 Neurog1 NEUROG2 NFATC2 NFE2 Nfe2l2 NFIA NFIC NFIX NHLH1 NKX2-3 Nkx2-5_v2 NKX2-8 Nkx3-1 NKX3-2 Npas2 NR1H2::RXRA NR2C2 Nr2e1 Nr2e3 NR2F1 NR2F2 Nr2f6 Nr2f6_v2 NR3C1 NR3C2 NR4A1 NR4A2 Nr5a2 NRF1 NRL OLIG1 OLIG2 OLIG3 Pax2 PAX5 PBX2 PBX3 PKNOX1 PKNOX2 POU1F1 POU2F1 POU2F2 Pou2f3 POU3F1 POU3F2 POU3F3 POU3F4 POU4F1 POU4F2 POU5F1 POU5F1B Pou5f1::Sox2 PPARA::RXRA PPARG Pparg::Rxra PRDM1 RARA::RXRA RARA::RXRG RARA_v2 Rarb Rarb_v2 Rarg Rarg_v2 RBPJ REST Rfx1 Rhox11 RORA RREB1 RUNX1 Rxra RXRB RXRG SCRT1 SCRT2 SIX1 SIX2 SMAD2::SMAD3::SMAD4 Smad4 SNAI2 Sox1 SOX10 Sox11 SOX15 Sox17 Sox2 SOX21 Sox3 SOX4 Sox5 Sox6 SOX8 SPI1 SPIB SREBF1 Srebf1_v2 SREBF2 SREBF2_v2 STAT1 STAT1::STAT2 STAT3 Stat4 Stat5a::Stat5b Stat6 T TAL1::TCF3 TBR1 TBX1 TBX15 TBX19 TBX2 TBX20 TBX21 TBX4 TBX5 Tcf12 Tcf21 TCF3 TCF4 Tcf7 TCF7L1 TCF7L2 TFAP2A TFAP2A_v3 TFAP2B_v2 TFAP2B_v3 TFAP2C_v2 TFAP2C_v3 TFAP4 TFE3 TFEB TFEC TGIF1 TGIF2 TP53 TP63 TP73 TWIST1 Twist2 USF1 USF2 VDR YY1 ZBTB18 ZBTB7A ZEB1 Zfx ZIC1 ZIC3 ZIC4 ZNF143 ZNF24 ZNF263 ZNF354C ZNF384 ZSCAN4
awk '($8 < 0.05 )' OFS='\t' JASPAR_default_nb_results.out | grep -wiF nb.gene.s1.386 | cut -f1 | sort -u | paste -sd " "
NR2C2 Nr2f6 Rxra RXRB RXRG #default filt
grep -wiF nb.gene.s1.386 JASPAR_default_nb_results.out | cut -f1 | sort -u | paste -sd " "
# default no filt
Ahr::Arnt Ar Arid3a Arnt ARNT::HIF1A Arntl ASCL1 Ascl2 Atf3 ATF4 ATF7 Atoh1 Bach1::Mafk BACH2 BATF3 BATF::JUN Bcl6 Bhlha15 BHLHE40 BHLHE41 CEBPA CEBPB CEBPD CEBPE CLOCK CREB1 Creb5 CTCF CTCFL Ddit3::Cebpa Dux DUX4 E2F6 EGR2 ELF5 ELK4 EOMES ERF ESR1 ESR2 Esrra ESRRB Esrrg ETS1 ETV2 ETV3 ETV6 FIGLA FOS FOSB::JUN FOSB::JUNB FOSB::JUNB_v2 FOS::JUN FOS::JUNB FOS::JUND FOS::JUN_v2 FOSL1 FOSL1::JUN FOSL1::JUNB FOSL1::JUND FOSL1::JUN_v2 FOSL2 FOSL2::JUN FOSL2::JUNB FOSL2::JUNB_v2 FOSL2::JUND FOSL2::JUND_v2 FOSL2::JUN_v2 FOXA1 Foxa2 FOXB1 FOXC1 FOXC2 FOXD1 FOXD2 Foxd3 FOXG1 FOXH1 FOXI1 Foxj2 Foxj3 FOXK1 FOXK2 FOXL1 Foxo1 FOXO3 FOXO4 FOXO6 FOXP1 FOXP2 FOXP3 Foxq1 Gabpa GATA1::TAL1 Gata4 GCM1 GCM2 Gfi1 Gfi1b GRHL1 GRHL2 Hand1::Tcf3 Hes2 HES5 HIC2 HLTF Hnf4a HNF4G HOXA5 Hoxa9 HOXB13 HOXD13 Hoxd9 Id2 ID4 IRF1 IRF2 IRF5 IRF7 IRF9 JDP2 JDP2_v2 JUN JUNB JUND JUND_v2 JUN::JUNB JUN_v2 Klf1 KLF4 KLF5 KLF9 LEF1 LIN54 Mafb MAFF MAFG MAFG::NFE2L1 MAFK MAF::NFE2 MAX MAX::MYC MEF2A MEF2C MEIS1 MEIS2 MEIS3 MGA MITF Mlxip MLXIPL MNT MSC MTF1 MXI1 MYB MYC MYF6 Myod1 Myog MZF1 NEUROD1 Neurog1 NEUROG2 NFATC2 NFE2 Nfe2l2 NFIA NFIC NFIX NHLH1 NKX2-3 Nkx2-5_v2 NKX2-8 Nkx3-1 NKX3-2 Npas2 NR1H2::RXRA NR2C2 Nr2e1 Nr2e3 NR2F1 NR2F2 Nr2f6 Nr2f6_v2 NR3C1 NR3C2 NR4A1 NR4A2 Nr5a2 NRF1 NRL OLIG1 OLIG2 OLIG3 Pax2 PAX5 PBX2 PBX3 PKNOX1 PKNOX2 POU1F1 POU2F1 POU2F2 Pou2f3 POU3F1 POU3F2 POU3F3 POU3F4 POU4F1 POU4F2 POU5F1 POU5F1B Pou5f1::Sox2 PPARA::RXRA PPARG Pparg::Rxra PRDM1 RARA::RXRA RARA::RXRG RARA_v2 Rarb Rarb_v2 Rarg Rarg_v2 RBPJ REST Rfx1 Rhox11 RORA RREB1 RUNX1 Rxra RXRB RXRG SCRT1 SCRT2 SIX1 SIX2 SMAD2::SMAD3::SMAD4 Smad4 SNAI2 Sox1 SOX10 Sox11 SOX15 Sox17 Sox2 SOX21 Sox3 SOX4 Sox5 Sox6 SOX8 SPI1 SPIB SREBF1 Srebf1_v2 SREBF2 SREBF2_v2 STAT1 STAT1::STAT2 STAT3 Stat4 Stat5a::Stat5b Stat6 T TAL1::TCF3 TBR1 TBX1 TBX15 TBX19 TBX2 TBX20 TBX21 TBX4 TBX5 Tcf12 Tcf21 TCF3 TCF4 Tcf7 TCF7L1 TCF7L2 TFAP2A TFAP2A_v3 TFAP2B_v2 TFAP2B_v3 TFAP2C_v2 TFAP2C_v3 TFAP4 TFE3 TFEB TFEC TGIF1 TGIF2 TP53 TP63 TP73 TWIST1 Twist2 USF1 USF2 VDR YY1 ZBTB18 ZBTB7A ZEB1 Zfx ZIC1 ZIC3 ZIC4 ZNF143 ZNF24 ZNF263 ZNF354C ZNF384 ZSCAN4

## Check that the sws1 rxra/b/g TFBS in N. brichardi is the same as that identified previously

# Previous details
# MA0512.2	mz.gene.s102.69	619	632	-	10.4024	6.02e-06	0.291	AGAGTCAGAGGTCA - SNP at 630
# MA0512.2	nb.gene.s1.386	344	357	-	19.122	1.53e-07	0.0329	AGGGTCAGAGGTCA
# mz.gene.s102.69	scaffold_102	2339992	2342530	+ > position is 2339992+630= 2340622; 'motif' position is 2340611-2340624

# The score for the match of a position in a sequence to a motif is computed by summing the appropriate entries from each column of the position-dependent scoring matrix that represents the motif.
awk '$2=="nb.gene.s1.386"' JASPAR_matQualPval_nb_results.out | awk '$1=="Rxra"'
# motif_name sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
# Rxra	nb.gene.s1.386	344	357	-	19.716	1.14e-07	0.0256	AGGGTCAGAGGTCA
awk '$2=="mz.gene.s102.69"' JASPAR_matQualPval_mz_results.out | awk '$1=="Rxra"'
# Rxra	mz.gene.s102.69	619	632	-	10.6914	5.08e-06	0.251	AGAGTCAGAGGTCA

# predict DNA binding domain here: http://lcg.rit.albany.edu/dp-bind/
for i in NR2C2 Nr2f6 Rxra RXRB RXRG ; do grep $i OGIDS.txt5 | cut -f5,10,12,15 ; done
nb.gene.s44.59	nr2c2	nr2c2	NR2C2
NULL	rxrba	rxrba	RXRB
nb.gene.s79.100	NULL	rxrbb	RXRB
NULL	rxrgb	rxrgb	RXRG

cd /tgac/workarea/Research-Groups/RG-cichlids/domain_scan
for i in nb.gene.s44.59 nb.gene.s79.100 ; do grep -wiF $i nb-interpro_collated.pfamfilt2.tsv ; done
nb.gene.s79.100	59e765ff6dfa5e6f83af31416e48c582	444	Pfam	PF00104	Ligand-binding_domain_of_nuclear_hormone_receptor	235	424	7.7E-36	T	06-02-2018	IPR000536	Nuclear_hormone_receptor,_ligand-binding,_core	GO:0003700|GO:0003707|GO:0005634|GO:0006355|GO:0043401
nb.gene.s79.100	59e765ff6dfa5e6f83af31416e48c582	444	Pfam	PF00105	Zinc_finger,_C4_type_(two_domains)	89	157	1.3E-31	T	06-02-2018	IPR001628	Zinc_finger,_nuclear_hormone_receptor-type	GO:0003700|GO:0005634|GO:0006355|GO:0008270|GO:0043565
nb.gene.s79.100	59e765ff6dfa5e6f83af31416e48c582	444	Pfam	PF11825	Nuclear/hormone_receptor_activator_site_AF-1	4	82	7.2E-14	T	06-02-2018	IPR021780	Nuclear/hormone_receptor_activator_site_AF-1	NA
for i in nb.prot.s44.59 nb.prot.s79.100 ; do grep $i nb-interpro_collated.tsv ; done

####### sws2a

# P. nyererei SWS2a
awk '($8 < 0.05 )' OFS='\t' CS_matQualPval_pn_results.out | grep -wiF pn.gene.s177.2 | cut -f1 | sort -u | paste -sd " "
for i in pn.gene.s109.14 pn.gene.s177.51 pn.gene.s227.19 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
NULL	GATA2_(1_of_many)	GATA2
gata2a	gata2a	GATA2
NULL	egr2b	EGR2

grep -wiF pn.gene.s177.2 CS_matQualPval_pn_results.out | cut -f1 | sort -u | paste -sd " "
for i in pn.gene.s0.237 pn.gene.s0.265 pn.gene.s100.28 pn.gene.s101.70 pn.gene.s102.26 pn.gene.s103.41 pn.gene.s106.41 pn.gene.s109.14 pn.gene.s110.12 pn.gene.s111.33 pn.gene.s116.33 pn.gene.s122.38 pn.gene.s127.11 pn.gene.s128.3 pn.gene.s165.8 pn.gene.s16.71 pn.gene.s167.7 pn.gene.s177.51 pn.gene.s189.13 pn.gene.s20.29 pn.gene.s207.18 pn.gene.s22.38 pn.gene.s22.39 pn.gene.s22.69 pn.gene.s227.19 pn.gene.s233.1 pn.gene.s24.13 pn.gene.s241.3 pn.gene.s26.28 pn.gene.s28.59 pn.gene.s3.110 pn.gene.s31.71 pn.gene.s3.182 pn.gene.s3.251 pn.gene.s3.272 pn.gene.s342.15 pn.gene.s34.39 pn.gene.s40.29 pn.gene.s412.8 pn.gene.s4.150 pn.gene.s43.16 pn.gene.s44.10 pn.gene.s46.57 pn.gene.s470.1 pn.gene.s47.7 pn.gene.s49.47 pn.gene.s51.18 pn.gene.s5.34 pn.gene.s54.33 pn.gene.s56.69 pn.gene.s57.25 pn.gene.s57.33 pn.gene.s57.96 pn.gene.s58.60 pn.gene.s6.10 pn.gene.s62.18 pn.gene.s70.23 pn.gene.s7.101 pn.gene.s71.5 pn.gene.s73.40 pn.gene.s84.15 pn.gene.s85.42 pn.gene.s86.28 pn.gene.s86.36 pn.gene.s87.44 pn.gene.s89.23 pn.gene.s91.40 pn.gene.s9.171 pn.gene.s91.75 pn.gene.s91.95 pn.gene.s92.48 pn.gene.s95.2 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
NULL	nr3c1	NR3C1
egr1	egr1	EGR1
zeb1b	zeb1b	ZEB1
eomesa	eomesa	EOMES
zeb1a	zeb1a	ZEB1
foxm1	foxm1	FOXM1
pax5	pax5	PAX5
NULL	GATA2_(1_of_many)	GATA2
epas1b	epas1b	EPAS1
tcf12	tcf12	TCF12
six5	six5	SIX5
NULL	MAX_(1_of_many)	MAX
NULL	eomesb	EOMES
NULL	NFATC1_(1_of_many)	NFATC1
pax6b	pax6b	PAX6
tbx21	tbx21	TBX21
bcl11aa	bcl11aa	BCL11A
gata2a	gata2a	GATA2
gabpa	gabpa	GABPA
mybl2b	mybl2b	MYBL2
NULL	pax6a	PAX6
NULL	tcf3b	TCF3
si:dkey-110c1.10	tcf3b	TCF3
tp63	NULL	TP63
NULL	egr2b	EGR2
CABZ01066696.1	NULL	LYL1
NULL	NR2F2_(1_of_many)	NR2F2
pou2f2a	POU2F2_(1_of_many)	POU2F2
NULL	TCF12_(1_of_many)	TCF12
whsc1	whsc1	WHSC1
NULL	tcf3a	TCF3
pou5f3	pou5f3	POU5F1
tal1	tal1	TAL1
irf4a	irf4a	IRF4
cebpd	cebpd	CEBPD
irf4b	irf4b	IRF4
nr2f2	nr2f2	NR2F2
prdm1c	PRDM1_(1_of_many)	PRDM1
sp4	sp4	SP4
maff	maff	MAFF
kdm5bb	kdm5bb	KDM5B
pou2f2a	POU2F2_(1_of_many)	POU2F2
nfe2	nfe2	NFE2
NULL	foxa1	FOXA1
atf3	atf3	ATF3
rfx2	rfx2	RFX2
nrf1	nrf1	NRF1
NULL	smad4a	SMAD4
klf5b	klf5b	KLF5
NULL	nfic	NFIC
batf	batf	BATF
NULL	fosab	FOS
esr1	esr1	ESR1
srebf1	srebf1	SREBF1
tfap2a	tfap2a	TFAP2A
smad4a	SMAD4_(1_of_many)	SMAD4
stat4	stat4	STAT4
NULL	fosl1a	FOSL1
sox2	sox2	SOX2
esr2a	esr2a	ESR2
stat1a	stat1a	STAT1
fosl2	fosl2	FOSL2
sox17	sox17	SOX17
mycb	myca	MYC
nkx2.4a	nkx2.4a	NKX2-1
prdm1a	prdm1a	PRDM1
kdm5ba	kdm5ba	KDM5B
hif1ab	hif-1a	HIF1A
NULL	BHLHE40_(1_of_many)	BHLHE40
stat2	stat2	STAT2
nfatc1	nfatc1	NFATC1
tfap2c	tfap2c	TFAP2C

awk '($8 < 0.05 )' OFS='\t' JASPAR_matQualPval_pn_results.out | grep -wiF pn.gene.s177.2 | cut -f1 | sort -u | paste -sd " "
CTCF RREB1 ZNF263 ZNF384
 #matqual filt
grep -wiF pn.gene.s177.2 JASPAR_matQualPval_pn_results.out | cut -f1 | sort -u | paste -sd " "
#matqual no filt
awk '($8 < 0.05 )' OFS='\t' JASPAR_default_pn_results.out | grep -wiF pn.gene.s177.2 | cut -f1 | sort -u | paste -sd " "
CTCF RREB1 ZNF263 ZNF384
 #default filt
grep -wiF pn.gene.s177.2 JASPAR_default_pn_results.out | cut -f1 | sort -u | paste -sd " "
# default no filt
Alx1 ALX3 Alx4 Arid3b Arntl ASCL1 Ascl2 Barhl1 BARHL2 BATF::JUN Bcl6 BCL6B BHLHE40 BHLHE41 CDX1 CDX2 CEBPA CEBPB CEBPD CEBPE CEBPG CLOCK Creb3l2 Crx CTCF CTCFL Dlx2 Dlx3 Dmbx1 DMRT3 DUX4 DUXA E2F6 E2F7 E2F8 EGR1 EGR3 EGR4 ELK4 EMX2 EOMES ERF ERG ESR2 Esrra ESRRB Esrrg ESX1 ETS1 ETV3 FIGLA FLI1 FOS FOSB::JUNB FOS::JUN FOS::JUNB FOS::JUND FOSL1 FOSL1::JUN FOSL1::JUNB FOSL1::JUND FOSL2 FOSL2::JUN FOSL2::JUNB FOSL2::JUND FOXA1 Foxa2 FOXB1 FOXC1 FOXC2 Foxd3 FOXF2 FOXH1 Foxj3 Foxo1 FOXP2 Foxq1 Gabpa Gata1 GATA1::TAL1 GATA2 Gata4 GATA6 GBX1 GBX2 GCM1 GCM2 Gfi1 Gfi1b GLI2 GLIS2 GRHL1 GRHL2 GSC GSC2 GSX1 GSX2 Hes2 HESX1 HEY1 Hic1 HIC2 Hmx1 Hmx2 Hmx3 HNF1A Hnf4a HNF4G HOXA10 Hoxa11 HOXA13 HOXA2 HOXA5 Hoxa9 HOXB13 HOXB2 HOXB3 Hoxb5 HOXC11 HOXC12 Hoxc9 HOXD12 HOXD13 Hoxd3 Hoxd8 Hoxd9 HSF1 HSF2 Id2 ID4 INSM1 IRF3 IRF5 IRF7 IRF8 IRF9 ISL2 ISX JDP2 JUNB JUND JUN::JUNB JUN_v2 Klf1 KLF16 KLF4 KLF5 KLF9 LBX1 LBX2 LEF1 LHX2 Lhx3 Lhx4 LHX6 LHX9 LIN54 LMX1A LMX1B Mafb MAF::NFE2 MAX MAX::MYC Mecom MEF2A MEF2C MEIS2 MEIS3 MEOX1 MEOX2 MGA MITF mix-a MIXL1 MLX Mlxip MLXIPL MNT MNX1 MSX1 MSX2 Msx3 MTF1 MYBL2 MYCN MYF6 Myod1 Myog MZF1_v2 NEUROD1 NEUROG2 NFAT5 NFATC1 NFATC2 NFATC3 Nfe2l2 NFIA NFIC::TLX1 NFIL3 NFIX NFYA NHLH1 NKX2-3 Nkx2-5 NKX2-8 Nkx3-1 NKX3-2 NKX6-1 NKX6-2 Nobox NOTO Npas2 NR1H2::RXRA Nr1h3::Rxra NR2C2 Nr2e1 Nr2e3 NR2F1 NR2F2 NR4A1 NR4A2 NR4A2::RXRA NRF1 NRL ONECUT1 ONECUT2 ONECUT3 OTX1 OTX2 PAX3 PAX5 PAX7 PBX1 PBX2 PBX3 PHOX2A Phox2b Pitx1 PITX3 PKNOX1 PKNOX2 PLAG1 POU2F1 POU2F2 Pou2f3 POU3F4 POU4F1 POU4F2 POU4F3 POU5F1 POU5F1B Pou5f1::Sox2 POU6F1 POU6F2 PPARA::RXRA PPARG Pparg::Rxra PRDM1 PROP1 PROX1 PRRX1 Prrx2 RARA Rarb Rarg Rarg_v2 RAX RAX2 RBPJ REL RELA RELB REST Rhox11 RORA RORA_v2 RORB RORC RREB1 RUNX1 RUNX2 RUNX3 SCRT1 SCRT2 SHOX Shox2 SIX1 SMAD2::SMAD3::SMAD4 SMAD3 SNAI2 Sox1 SOX10 SOX13 Sox17 Sox2 SOX21 Sox3 SOX4 Sox6 SP1 SP2 SP3 SP8 SPI1 SPIB SPIC Spz1 SREBF1 Srebf1_v2 SREBF2 SREBF2_v2 SRY STAT1 STAT3 Stat4 Stat6 TAL1::TCF3 TBR1 TBX20 TBX21 Tcf12 Tcf21 TCF3 TCF4 Tcf7 TCF7L1 TCF7L2 TEAD1 TEAD2 TEAD3 TEAD4 TFAP4 TFCP2 TFE3 TFEB TFEC TGIF1 TGIF2 THAP1 TP53 TWIST1 USF1 USF2 VAX1 VAX2 YY1 ZBTB18 ZBTB7A ZBTB7B ZBTB7C ZEB1 Zfx ZIC1 ZIC3 ZIC4 ZNF143 ZNF24 ZNF263 ZNF282 ZNF384 ZNF740 ZSCAN4
# check for Crx
awk '$2=="pn.gene.s177.2"' JASPAR_matQualPval_pn_results.out | awk '$1=="Crx"'

awk '$2=="pn.gene.s177.2"' JASPAR_default_pn_results.out | awk '$1=="Crx"'
# Crx	pn.gene.s177.2	4410	4420	+	11.5439	5.81e-05	1	AGAGGGATTAA

# A. burtoni SWS2a
awk '($8 < 0.05 )' OFS='\t' CS_matQualPval_ab_results.out | grep -wiF ab.gene.s9.15 | cut -f1 | sort -u | paste -sd " "
for i in ab.gene.s271.28 ab.gene.s373.3 ab.gene.s395.10 ab.gene.s9.66 ; do grep -wiF $i OGIDS.txt5 | cut -f10,12,15 ; done
tp53	tp53	TP53
NULL	egr2b	EGR2
NULL	GATA2_(1_of_many)	GATA2
gata2a	gata2a	GATA2
grep -wiF ab.gene.s9.15 CS_matQualPval_ab_results.out | cut -f1 | sort -u | paste -sd " "
awk '($8 < 0.05 )' OFS='\t' JASPAR_matQualPval_ab_results.out | grep -wiF ab.gene.s9.15 | cut -f1 | sort -u | paste -sd " "
RREB1 ZNF263 ZNF384
 #matqual filt
grep -wiF ab.gene.s9.15 JASPAR_matQualPval_ab_results.out | cut -f1 | sort -u | paste -sd " "
#matqual no filt
awk '($8 < 0.05 )' OFS='\t' JASPAR_default_ab_results.out | grep -wiF ab.gene.s9.15 | cut -f1 | sort -u | paste -sd " "
RREB1 ZNF263 ZNF384
 #default filt
grep -wiF ab.gene.s9.15 JASPAR_default_ab_results.out | cut -f1 | sort -u | paste -sd " "
# default no filt

# O. niloticus SWS2a
awk '($8 < 0.05 )' OFS='\t' CS_matQualPval_on_results.out | grep -wiF on.gene.LG5.187 | cut -f1 | sort -u | paste -sd " " #NONE
grep -wiF on.gene.LG5.187 CS_matQualPval_on_results.out | cut -f1 | sort -u | paste -sd " "
awk '($8 < 0.05 )' OFS='\t' JASPAR_matQualPval_on_results.out | grep -wiF on.gene.LG5.187 | cut -f1 | sort -u | paste -sd " " # NONE
 #matqual filt
grep -wiF on.gene.LG5.187 JASPAR_matQualPval_on_results.out | cut -f1 | sort -u | paste -sd " "
#matqual no filt
awk '($8 < 0.05 )' OFS='\t' JASPAR_default_on_results.out | grep -wiF on.gene.LG5.187 | cut -f1 | sort -u | paste -sd " " # NONE
 #default filt
grep -wiF on.gene.LG5.187 JASPAR_default_on_results.out | cut -f1 | sort -u | paste -sd " "
# default no filt


# rh2
mz.gene.s18.149
pn.gene.s11.51
ab.gene.s86.69
nb.gene.s4.225
on.gene.LG5.742

# rho
mz.gene.s12.95
ab.gene.s134.8
on.gene.LG20.398




#### Creating alternative background files for fimo runs (07/03/2018)

# Providing an appropriate background model is the single most important step you can take to improve your FIMO results.
# The FIMO match score is the log-odds ratio of the probability of observing a sequence given the motif PWM to the probability of observing the same sequence given the background model.
# Ideally you'd create a background model from sequence data that is biologically similar to the sequences you are analyzing with FIMO, but that doesn't contain any instances of the motifs you are trying to match. If you are scanning a full genome then using the nucleotide frequencies for that genome as a background model is a reasonable approach.
# In this way, provide three background files
  # 1. 0-order background of promoters {DONE}
  # 2. 0-order background of whole genome {DONE - use whole fasta}
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta
  # 3. 0-order background of whole genome minus promoter regions
# Then re-run motif scanning using cichlid-specific motifs on a single species and compare the results

cd /tgac/workarea/Research-Groups/RG-cichlids/
mkdir fimoBG
cd /tgac/workarea/Research-Groups/RG-cichlids/fimoBG


ml bedtools/2.25.0
ml zlib
# create 0-order background of whole genome minus promoter regions
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed
# a. create files for scaffold lengths
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Ab_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Mz_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Nb_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > On_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Pn_scafflength
# b. create a new bedfile that just has 0 (start) to scaffold_length (end)
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Ab_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Mz_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Mz_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Nb_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Nb_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' On_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > On_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pn_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Pn_scafflength.bed
# c. Use bedtools subtract to subtract the overlapping promoter regions from the full scaffold length
bedtools subtract -a Ab_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Ab_noprom.bed
bedtools subtract -a Mz_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Mz_noprom.bed
bedtools subtract -a Nb_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Nb_noprom.bed
bedtools subtract -a On_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed > On_noprom.bed
bedtools subtract -a Pn_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Pn_noprom.bed
# d. extract regions from whole genome
nano getfasta.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml bedtools/2.25.0
ml zlib

bedtools getfasta -fo Ab_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -bed Ab_noprom.bed
bedtools getfasta -fo Mz_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -bed Mz_noprom.bed
bedtools getfasta -fo Nb_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -bed Nb_noprom.bed
bedtools getfasta -fo On_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -bed On_noprom.bed
bedtools getfasta -fo Pn_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -bed Pn_noprom.bed

# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=2 getfasta.sh #{DONE}

## Create 0-order backgrounds for whole genomes as input for fimo

nano Mz_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# M. zebra
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta Mz.wholegenome.0orderMarkovBackgrnd

nano Ab_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# Ab
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta Ab.wholegenome.0orderMarkovBackgrnd

nano Nb_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# Nb
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta Nb.wholegenome.0orderMarkovBackgrnd

nano On_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# On
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta On.wholegenome.0orderMarkovBackgrnd

nano Pn_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# Pn
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta Pn.wholegenome.0orderMarkovBackgrnd


# Run all on uv2 - RAN ALL 9/3/18
qsub -q Test -l select=1:mem=100GB:ncpus=1 Mz_WG_BG.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 Pn_WG_BG.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 Ab_WG_BG.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 Nb_WG_BG.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 On_WG_BG.sh #{DONE}

### run relevant fimo on whole genome backgrounds
scratch

nano fimo_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
source python-3.5.1

# run script
python3.5 matQualResultsParser1_hsap_PBS_wg.py /tgac/workarea/Research-Groups/RG-cichlids/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac_v2.txt /tgac/workarea/Research-Groups/RG-cichlids/Jaspar_PSSMs_hsap_1.txt /tgac/workarea/Research-Groups/RG-cichlids/Jaspar_PSSMs_mmus_1.txt /tgac/workarea/Research-Groups/RG-cichlids/OGIDS.txt5 /tgac/workarea/Research-Groups/RG-cichlids/rsat_results/ /tgac/workarea/Research-Groups/RG-cichlids/rsat_results_MZ/ /tgac/workarea/Research-Groups/RG-cichlids/rsat_results_PN/ /tgac/scratch/mehtat/

# run the above
qsub -q Test -l select=1:mem=10GB:ncpus=1 fimo_WG_BG.sh


nano fimo_WG_BGmm.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
source python-3.5.1

# run script
python3.5 matQualResultsParser1_mus_PBS_wg.py /tgac/workarea/Research-Groups/RG-cichlids/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac_v2.txt /tgac/workarea/Research-Groups/RG-cichlids/Jaspar_PSSMs_hsap_1.txt /tgac/workarea/Research-Groups/RG-cichlids/Jaspar_PSSMs_mmus_1.txt /tgac/workarea/Research-Groups/RG-cichlids/OGIDS.txt5 /tgac/workarea/Research-Groups/RG-cichlids/rsat_results/ /tgac/workarea/Research-Groups/RG-cichlids/rsat_results_MZ/ /tgac/workarea/Research-Groups/RG-cichlids/rsat_results_PN/ /tgac/scratch/mehtat/

# run the above
qsub -q Test -l select=1:mem=10GB:ncpus=1 fimo_WG_BGmm.sh

## then run the parsing scripts

nano Hsparse_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
source python-3.5.1

# run script
python3.5 mat_qual_fimo_parse_hsap_wg.py bigFIMOtings/

# run the above
qsub -q Test -l select=1:mem=50GB:ncpus=1 Hsparse_WG_BG.sh


nano Mmparse_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
source python-3.5.1

# run script
python3.5 mat_qual_fimo_parse_mmus_wg.py bigFIMOtings_Mmus/

# run the above
qsub -q Test -l select=1:mem=50GB:ncpus=1 Mmparse_WG_BG.sh

cd /tgac/scratch/mehtat
# filter for qval < 0.05 of all parsed WG background files
nano qvalfilt.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

# qval filt < 0.05
for i in matQualOutput_wg/human/*.out ; do
  awk '($8 < 0.05 )' OFS='\t' $i > $i.filt
  for j in matQualOutput_wg/mouse/*.out ; do
    awk '($8 < 0.05 )' OFS='\t' $j > $j.filt ;
  done
done

# run counts of each file and output to file
for i in matQualOutput_wg/human/*.filt ; do
  wc -l $i >> matQualOutput_wg/human/Hsfiltcounts.txt
  for j in matQualOutput_wg/mouse/*.filt ; do
    wc -l $j >> matQualOutput_wg/mouse/Mmfiltcounts.txt ;
  done
done

# run the above
qsub -q Test -l select=1:mem=50GB:ncpus=1 qvalfilt.sh


# Create new promtoer annotations where none are out of bounds and extract new strand specific promoter fasta files

ml bedtools/2.25.0
ml zlib

# a. create files for scaffold lengths
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Ab_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Mz_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Nb_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > On_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Pn_scafflength
# b. create a new bedfile that corrects for promoters out of bounds
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Ab_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Mz_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Mz_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Nb_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Nb_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' On_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > On_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pn_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Pn_scafflength.bed


# re-extract promoter regions with native promoter strand - previous ones are all +ve strands
nano getPromStrandedfasta.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml bedtools/2.25.0
ml zlib

bedtools getfasta -s -fo Ab_14032018_5kbpromoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -bed Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed

bedtools getfasta -s -fo Mz_14032018_5kbpromoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -bed Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed

bedtools getfasta -s -fo Nb_14032018_5kbpromoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -bed Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed

bedtools getfasta -s -fo On_14032018_5kbpromoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -bed Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed

bedtools getfasta -s -fo Pn_14032018_5kbpromoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -bed Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed

# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=1 getPromStrandedfasta.sh #{DONE}





# create 0-order background of whole genome minus promoter regions
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed
/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed
# a. create files for scaffold lengths
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Ab_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Mz_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Nb_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > On_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Pn_scafflength
# b. create a new bedfile that just has 0 (start) to scaffold_length (end)
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Ab_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Mz_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Mz_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Nb_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Nb_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' On_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > On_scafflength.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pn_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{print $1,$2=0,$7,$4,$5,$6}' OFS='\t' > Pn_scafflength.bed
# c. Use bedtools subtract to subtract the overlapping promoter regions from the full scaffold length
bedtools subtract -a Ab_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Ab_noprom.bed
bedtools subtract -a Mz_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Mz_noprom.bed
bedtools subtract -a Nb_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Nb_noprom.bed
bedtools subtract -a On_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed > On_noprom.bed
bedtools subtract -a Pn_scafflength.bed -b /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed > Pn_noprom.bed
# d. extract regions from whole genome
nano getfasta.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml bedtools/2.25.0
ml zlib

bedtools getfasta -fo Ab_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -bed Ab_noprom.bed
bedtools getfasta -fo Mz_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -bed Mz_noprom.bed
bedtools getfasta -fo Nb_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -bed Nb_noprom.bed
bedtools getfasta -fo On_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -bed On_noprom.bed
bedtools getfasta -fo Pn_noprom.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -bed Pn_noprom.bed

# Run all on uv2
qsub -q Test -l select=1:mem=100GB:ncpus=2 getfasta.sh #{DONE}

## Create 0-order backgrounds for whole genomes as input for fimo

nano Mz_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# M. zebra
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta Mz.wholegenome.0orderMarkovBackgrnd

nano Ab_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# Ab
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta Ab.wholegenome.0orderMarkovBackgrnd

nano Nb_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# Nb
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta Nb.wholegenome.0orderMarkovBackgrnd

nano On_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# On
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta On.wholegenome.0orderMarkovBackgrnd

nano Pn_WG_BG.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml meme/4.11.1

# create 0-order markov background for fimo runs

# Pn
fasta-get-markov -dna /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta Pn.wholegenome.0orderMarkovBackgrnd


# Run all on uv2 - RAN ALL 9/3/18
qsub -q Test -l select=1:mem=100GB:ncpus=1 Mz_WG_BG.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 Pn_WG_BG.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 Ab_WG_BG.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 Nb_WG_BG.sh #{DONE}
qsub -q Test -l select=1:mem=100GB:ncpus=1 On_WG_BG.sh #{DONE}


# Create new promoter annotations where none are out of bounds
cd /tgac/workarea/Research-Groups/RG-cichlids/fimoBG

nano getPromStrandedfasta.sh

#!/bin/bash
cd $PBS_O_WORKDIR;
#load the latest module
ml bedtools/2.25.0
ml zlib

# a. create files for scaffold lengths
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Ab_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Mz_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Nb_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > On_scafflength
cat /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Pn_scafflength
# b. create new promoter bedfile annotations that correct for promoters out of bounds
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Ab_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{OFS="\t"} {if ($1==$7 && $3>$8) $3=$8; print $0}' | awk '{OFS="\t"} {if ($1==$7 && $2>$8) $2=$8; print $0}' | awk '$3!=0' | sort -k1,1 -k2,2n -k3,3n | awk '{print $1,$2,$3,$5,$5="1",$4}' OFS='\t' > /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Ab_14032018.5kb_sorted.promoters.stranded.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Mz_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{OFS="\t"} {if ($1==$7 && $3>$8) $3=$8; print $0}' | awk '{OFS="\t"} {if ($1==$7 && $2>$8) $2=$8; print $0}' | awk '$3!=0' | sort -k1,1 -k2,2n -k3,3n | awk '{print $1,$2,$3,$5,$5="1",$4}' OFS='\t' > /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Mz_14032018.5kb_sorted.promoters.stranded.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Nb_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{OFS="\t"} {if ($1==$7 && $3>$8) $3=$8; print $0}' | awk '{OFS="\t"} {if ($1==$7 && $2>$8) $2=$8; print $0}' | awk '$3!=0' | sort -k1,1 -k2,2n -k3,3n | awk '{print $1,$2,$3,$5,$5="1",$4}' OFS='\t' > /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Nb_14032018.5kb_sorted.promoters.stranded.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' On_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{OFS="\t"} {if ($1==$7 && $3>$8) $3=$8; print $0}' | awk '{OFS="\t"} {if ($1==$7 && $2>$8) $2=$8; print $0}' | awk '$3!=0' | sort -k1,1 -k2,2n -k3,3n | awk '{print $1,$2,$3,$5,$5="1",$4}' OFS='\t' > /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/On_14032018.5kb_sorted.promoters.stranded.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' Pn_scafflength /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all_sorted.bed | awk '{OFS="\t"} {if ($1==$7 && $3>$8) $3=$8; print $0}' | awk '{OFS="\t"} {if ($1==$7 && $2>$8) $2=$8; print $0}' | awk '$3!=0' | sort -k1,1 -k2,2n -k3,3n | awk '{print $1,$2,$3,$5,$5="1",$4}' OFS='\t' > /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Pn_14032018.5kb_sorted.promoters.stranded.bed
# if you also want to extract as fasta to check then
# c. create symbolic links to new promoter annotations
for i in /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/*_14032018.5kb_sorted.promoters.stranded.bed ; do ln -s $i ; done

# d. re-extract promoter regions with native promoter strand - previous ones are all +ve strands
bedtools getfasta -s -name -fo Ab_14032018.5kb_sorted.promoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta -bed Ab_14032018.5kb_sorted.promoters.stranded.bed

bedtools getfasta -s -name -fo Mz_14032018.5kb_sorted.promoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/MetZeb1.1/M_zebra_v1.1_unscreened_3750.assembly.fasta -bed Mz_14032018.5kb_sorted.promoters.stranded.bed

bedtools getfasta -s -name -fo Nb_14032018.5kb_sorted.promoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta -bed Nb_14032018.5kb_sorted.promoters.stranded.bed

bedtools getfasta -s -name -fo On_14032018.5kb_sorted.promoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Orenil1.1/20120125_MapAssembly.anchored.assembly.fasta -bed On_14032018.5kb_sorted.promoters.stranded.bed

bedtools getfasta -s -name -fo Pn_14032018.5kb_sorted.promoters.stranded.fasta -fi /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta -bed Pn_14032018.5kb_sorted.promoters.stranded.bed

# e. copy to group vh
for i in *_14032018.5kb_sorted.promoters.stranded.fasta ; do cp $i /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/ ; done

#### preliminary look at the SWS1 networks

## SWS1 - small reconstruction of networks

mkdir /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Opsins_prelimMar2018
cd /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Opsins_prelimMar2018

nano 1.MzNbSWS1.sh

#!/bin/bash
cd $PBS_O_WORKDIR;

# Mz
grep -wiF mz.gene.s102.69 /tgac/workarea/Research-Groups/RG-cichlids/GTRD_TFTGcoexpCheck/ModuleMerged_mm10_tftgtfbs_merged_07.extrap.output > 1b_MzSWS1
awk '($8 < 0.05 )' OFS='\t' /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse/CS_default_mz_results.out | grep mz.gene.s102.69 > 2a_MzSWS1
awk '($8 < 0.05 )' OFS='\t' /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse/CW_default_mz_results.out | grep mz.gene.s102.69 > 2b_MzSWS1
awk '($8 < 0.05 )' OFS='\t' /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse/JASPAR_default_mz_results.out | grep mz.gene.s102.69 > 2c_MzSWS1
awk '($3>=0.5)' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/8.TFTGco/tftgco_merged.txt | grep -wiF mz.gene.s102.69 > 3_MzSWS1-tftgco
grep -wiF mz.gene.s102.69 /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/PPI_Edge_Attributes/Mz-PPIEdge_Attributes_Collated4c.txt > 4_MzSWS1-PPI

# Nb
grep -wiF nb.gene.s1.386 /tgac/workarea/Research-Groups/RG-cichlids/GTRD_TFTGcoexpCheck/ModuleMerged_mm10_tftgtfbs_merged_07.extrap.output > 1b_NbSWS1
awk '($8 < 0.05 )' OFS='\t' /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse/CS_default_nb_results.out | grep nb.gene.s1.386 > 2a_NbSWS1
awk '($8 < 0.05 )' OFS='\t' /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse/CW_default_nb_results.out | grep nb.gene.s1.386 > 2b_NbSWS1
awk '($8 < 0.05 )' OFS='\t' /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput/mouse/JASPAR_default_nb_results.out | grep nb.gene.s1.386 > 2c_NbSWS1
awk '($3>=0.5)' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/8.TFTGco/tftgco_merged.txt | grep -wiF nb.gene.s1.386 > 3_NbSWS1-tftgco
grep -wiF nb.gene.s1.386 /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/PPI_Edge_Attributes/Nb-PPIEdge_Attributes_Collated4c.txt > 4_NbSWS1-PPI

# run the above
qsub -q Test -l select=1:mem=20GB:ncpus=1 1.MzNbSWS1.sh

# Mz - prep unified format files
cut -f1,3 1b_MzSWS1 | awk '$3="1b"' OFS='\t' > 1b_MzSWS1.a ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz 1b_MzSWS1.a | cut -f1-3,13,15,18 | awk '{print $2,$1,$3,$4,$5,$6}' OFS='\t' > 1b_MzSWS1.b ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz 1b_MzSWS1.b | cut -f1-6,18 | awk '{print $2,$4,$5,$6,$1,$7,$3}' OFS='\t' > 1b_MzSWS1.c

cut -f1,2 2a_MzSWS1 | awk '$3="2a"' OFS='\t' | sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' > 2a_MzSWS1.a ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz 2a_MzSWS1.a | cut -f1-3,13,15,18 | awk '{print $2,$1,$3,$4,$5,$6}' OFS='\t' > 2a_MzSWS1.b ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz 2a_MzSWS1.b | cut -f1-6,18 | awk '{print $2,$4,$5,$6,$1,$7,$3}' OFS='\t' > 2a_MzSWS1.c

awk '{print $2,$1,$3="2c"}' OFS='\t' 2c_MzSWS1 > 2c_MzSWS1.a
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz 2c_MzSWS1.a | cut -f1-3,13 | awk '{print $2,$2,$2,$2,$1,$4,$3}' OFS='\t' > 2c_MzSWS1.c

cut -f1,2 3_MzSWS1-tftgco | awk '$3="3"' OFS='\t' > 3_MzSWS1-tftgco.a ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz 3_MzSWS1-tftgco.a | cut -f1-3,13,15,18 | awk '{print $2,$1,$3,$4,$5,$6}' OFS='\t' > 3_MzSWS1-tftgco.b ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-mz 3_MzSWS1-tftgco.b | cut -f1-6,18 | awk '{print $2,$4,$5,$6,$1,$7,$3}' OFS='\t' > 3_MzSWS1-tftgco.c

awk '{print $1,$2,$2,$2,$5,$6,$7="4"}' 4_MzSWS1-PPI > 4_MzSWS1-PPI.c
for i in *Mz*.c ; do cat $i | awk '{OFS="\t"} {if ($2=="NULL") $2=$3; print $0}' | awk '{OFS="\t"} {if ($3=="NULL") $3=$4; print $0}' | awk '{OFS="\t"} {if ($2=="NULL") $2=$3; print $0}' | awk '{OFS="\t"} {if ($4=="NULL") $4=$2; print $0}' | awk '{OFS="\t"} {if ($3=="NULL") $3=$4; print $0}' >> MzSWS1.merge ; done
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Edge_Attributes/geneNamesMapping2.txt MzSWS1.merge | awk '{OFS="\t"} {if ($8=="NONE") $8=$2; print $0}' | awk '{print $1,$8,$5,$6,$7}' OFS='\t' > MzSWS1.merge2.txt # merge to create final file


# Nb - prep unified format files
cut -f1,3 1b_NbSWS1 | awk '$3="1b"' OFS='\t' > 1b_NbSWS1.a ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb 1b_NbSWS1.a | cut -f1-3,13,15,18 | awk '{print $2,$1,$3,$4,$5,$6}' OFS='\t' > 1b_NbSWS1.b ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb 1b_NbSWS1.b | cut -f1-6,18 | awk '{print $2,$4,$5,$6,$1,$7,$3}' OFS='\t' > 1b_NbSWS1.c

cut -f1,2 2a_NbSWS1 | awk '$3="2a"' OFS='\t' | sed 's/RG-cich_GTRDdata_mouse_sites_TF_ig_//g' | sed 's/.ig_1//g' > 2a_NbSWS1.a ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb 2a_NbSWS1.a | cut -f1-3,13,15,18 | awk '{print $2,$1,$3,$4,$5,$6}' OFS='\t' > 2a_NbSWS1.b ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb 2a_NbSWS1.b | cut -f1-6,18 | awk '{print $2,$4,$5,$6,$1,$7,$3}' OFS='\t' > 2a_NbSWS1.c

awk '{print $2,$1,$3="2c"}' OFS='\t' 2c_NbSWS1 > 2c_NbSWS1.a
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb 2c_NbSWS1.a | cut -f1-3,13 | awk '{print $2,$2,$2,$2,$1,$4,$3}' OFS='\t' > 2c_NbSWS1.c

cut -f1,2 3_NbSWS1-tftgco | awk '$3="3"' OFS='\t' > 3_NbSWS1-tftgco.a ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb 3_NbSWS1-tftgco.a | cut -f1-3,13,15,18 | awk '{print $2,$1,$3,$4,$5,$6}' OFS='\t' > 3_NbSWS1-tftgco.b ;
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' /tgac/scratch/mehtat/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5-nb 3_NbSWS1-tftgco.b | cut -f1-6,18 | awk '{print $2,$4,$5,$6,$1,$7,$3}' OFS='\t' > 3_NbSWS1-tftgco.c

awk '{print $1,$2,$2,$2,$5,$6,$7="4"}' 4_NbSWS1-PPI > 4_NbSWS1-PPI.c
for i in *Nb*.c ; do cat $i | awk '{OFS="\t"} {if ($2=="NULL") $2=$3; print $0}' | awk '{OFS="\t"} {if ($3=="NULL") $3=$4; print $0}' | awk '{OFS="\t"} {if ($2=="NULL") $2=$3; print $0}' | awk '{OFS="\t"} {if ($4=="NULL") $4=$2; print $0}' | awk '{OFS="\t"} {if ($3=="NULL") $3=$4; print $0}' >> NbSWS1.merge ; done
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' ../../Edge_Attributes/geneNamesMapping2.txt NbSWS1.merge | awk '{OFS="\t"} {if ($8=="NONE") $8=$2; print $0}' | awk '{print $1,$8,$5,$6,$7}' OFS='\t' > NbSWS1.merge2.txt # merge to create final file

for i in *merge2.txt ; do cp $i /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/Opsins_prelimMar2018/ ; done # copy to workarea to copy local

### Adding in some gene names in OGIDs - based on blastp of protein > OGIDS_genes.toadd
cd /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;};NR>FNR{if($1 in a){print $0,a[$1];}else{print $0,"NA";}}' OGIDS_genes.toadd OGIDS.txt5 | awk '{OFS="\t"} {if ($10=="NULL" && $17!="NA") $10=$17; print $0}' | cut -f1-16 > OGIDS.txt6
rm OGIDS.txt5 ; mv OGIDS.txt6 OGIDS.txt5

# try to get rid of most 'si:dkey' symbols
awk '{OFS="\t"} {if ($10~"si:dkey" && $15!="NULL") $10=$15; print $0}' OGIDS.txt5 | awk '{OFS="\t"} {if ($12~"si:dkey" && $15!="NULL") $12=$15; print $0}' > OGIDS.txt6 ; rm OGIDS.txt5 ; mv OGIDS.txt6 OGIDS.txt5
