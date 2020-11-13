############################################################################

Network Reconstruction constrained by Arboretum output of

5 cichlids (O. niloticus, P. nyererei, A. burtoni, M .zebra, N. brichardi)

Interaction inputs - TFBSs and TF-TG co-expression

By Tarang Mehta, Earlham, UK
Version 1.0 2019

License: 
This source code is freely available to all under the Creative Commons 
Attribution-ShareAlike license (CC BY-SA) and under the standard GPL 
3.0 license from Github

############################################################################

Run script: ./NetworkReconstruction_v1.0_2019.sh [workingdir] Module_genesandexpr TFBSs TFTGco Edge_Attributes

NOTE: This will not run end to end as there are instances of copying code into scripts using nano

Files in [workingdir]:
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
- MzPnAbNbOnGenome-BLAST_PresentNULLOGIDS-noCand-noTF.txt2

############################################################################
