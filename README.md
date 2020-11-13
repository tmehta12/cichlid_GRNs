# Cichlid GRNs
Gene regulatory network reconstruction of five cichlid species (M. zebra, P. nyererei, A. burtoni, N. brichardi and O. niloticus)

    - Construction of cichlid gene trees - MUSCLE-3.8.31; OrthoMCL-1.4.0; TreeFix-1.1.10 
    - Inference of multi- and single- tissue transcriptional modules - Arboretum2.0 
    - Transcription factor (TF) motif scanning  - MEME-4.11.1; RSAT Update 2017/04/21; FIMO-4.11.1 
    - Variation and evolutionary rate at coding and non-coding regions - BLAT-35; MAFFT-7.271; PAML-4.9 
    - Regulatory network rewiring analysis of gene sets - DyNet-2.0; Cytoscape-3.7.1
    - Identification of segregating sites in TFBSs - bedtools-2.25.0

================================================================

## 0.Arboretum_scripts

Run as described in published literature: Roy, S. et al. Arboretum: Reconstruction and analysis of the evolutionary history of condition-specific transcriptional modules. Genome Res. 23, 1039–1050 (2013).
https://github.com/Roy-lab/Arboretum2.0

================================================================

## 1a.Network_Reconstruction_script

Script to run network reconstruction from TFBSs and TF-TG co-expression files
- NetworkReconstruction_v1.0_2019.sh 

================================================================

## 1b.Network_Reconstruction_files

In Network_Reconstruction_files.tgz are the major files required to run network reconstruction and some output files:
      - geneNamesTree # gene orthology of cichlid IDs, Ensembl IDs and gene descriptions
      - Hs-Mm_EnsemblOrthology.txt # Ensembl orthology of Human and Mouse gene IDs
      - full_orthologs_map3c.txt # contains the published orthologs gene IDs for mapping used previously
      - list_k10 # list of k10 module numbers
      - list_species # list of species
      - OGIDS.txt # orthogroups
      - OGIDS.txt5 # processed orthogroups file
      - *_07.extrap.annotations.blast_out # extrapolation of motifs in cichlid gene promoters using human and mouse orthology 
      - mz-speciesspecnames_clusterassign.txt # Arboretum module assignment of M. zebra genes
      - pn-speciesspecnames_clusterassign.txt # Arboretum module assignment of P. nyererei genes
      - ab-speciesspecnames_clusterassign.txt # Arboretum module assignment of A. burtoni genes
      - nb-speciesspecnames_clusterassign.txt # Arboretum module assignment of N. brichardi genes
      - on-speciesspecnames_clusterassign.txt # Arboretum module assignment of O. niloticus genes
      - geneNamesMapping.txt # gene orthology of cichlid IDs, Ensembl IDs and gene descriptions
      - MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS-noCand-noTF.txt2 # 4209 orthogroups with NULL orthogroup IDs for one or species are actually present in the genome (mis-annotations etc.) and not lost/absent
      - motifenr_merged-TFBSs_map2d.txt # output of processing files in TFBSs - Motif Discovery in gene promoter regions
      - tftgco_merged.txt # output of processing TF-TG co-expression files

Files are here: https://doi.org/10.6084/m9.figshare.7707437

================================================================

## 2a.TF_motif_scanning_scripts

The scripts are in the process of being collated into a pipeline tool which will be published as a companion paper. 
As such, scripts lack comments within their body. Future versions will be fully commented and released with full documentation

Author:
All scripts described below were written and tested by Dr. Will Nash
(ORCID: 0000-0002-6790-1167) during the preparation of material for this paper. 

The scripts below were written to run in the stated order, with each script
relying on the output of that prior to it. Where a and b versions exist, these are
motifs for a) human and b) mouse data and were handled separately
during the development of the approach. The scripts are identical, but deal with
species specific descriptors differently. Both versions are included here for
clarity, but scripts will be generalised in the publication version.

Scripts in 2a.TF_motif_scanning_scripts.zip:
      - 1_GTRD_parse.py
      - 2a_extrapSites_hum_v3.py
      - 2b_extrapSites_mus_v3.py
      - 3a_extrapSites_hum_v4.py
      - 3b_extrapSites_mus_v4.py
      - 4a_extrapSites_hum_v4ex.py
      - 4b_extrapSites_mus_v4ex.py
      - 5a_extrapSites_hum_afterMEME.py
      - 5b_extrapSites_mus_afterMEME.py
      - 6_TFBSextrapolation_v5.py
      - 7_TFsortSplit5.py
      - 8_run_mat_qual_nm.py
      - 9a_matQualResultsParser1_hsap.py
      - 9b_matQualResultsParser1_mus.py
      - 10a_mat_qual_fimo_parse_hsap_wg.py
      - 10b_mat_qual_fimo_parse_mmus_wg.py
      - Motifs.zip # contains MEME format position weight matrices (PWM) of cichlid-specific (CS), cichlid-wide (CW) and JASPAR vertebrate motifs for scanning cichlid gene promoters (includes a Readme file). Files are here: https://doi.org/10.6084/m9.figshare.7599293

================================================================

## 2b.TF_motif_scanning_outputs

This contains output files from the TF motif scanning scripts.
Files are here: https://doi.org/10.6084/m9.figshare.7712423

================================================================

## 2c. Co-expression (TF-TG) based regulatory networks 

These were inferred using MERLIN-P without the module and motif prior available at https://github.com/Roy-lab/merlin-p

================================================================

## 3.Edge_attribute_file_RewiringAnalysis_file

Based on the NetworkReconstruction_v1.0_2019.sh script and subsequent formation of a matrix for each edge (1 for present and 0 for absent), the file in the below zip file was used for 1.enrichment analysis and 2.Regulatory network rewiring analysis of gene sets 
- Edge_Attributes_Collated4c.coexpr_promONLY.NEWsimplified.matrix3.filteredPresentNULLOGIDS.txt.zip 
File is here: https://doi.org/10.6084/m9.figshare.7707455

================================================================

## 4. TF-TG gain and loss rates (TF binding rates)

These were estimated using code at https://github.com/Roy-lab/compfuncgen_utils

================================================================
