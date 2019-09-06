# Prepare plots of the SNPs in candidate gene TFBSs that overlap Lake species SNPs (in Malawi and Victoria)

library(ggplot2)
library(hexbin)
library(reshape2)
library(gridExtra)
library(grid)
library(magrittr)
library(dplyr)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/EMSA_candidates/")

# 1. How many times, in each species, does the genotype match Mz/Pn (Malawi/Victoria - Ref allele) or not (match the Alt or different allele)
# 2. Do the species that are closely related to Mz/Pn have higher same TFBS SNP genotypes: clade-specific plots of genotype match > is there a correlation of segregating TFBS SNPs with phylogenetic distance

mz_malawiGT2 = read.table('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3', sep='\t', header=TRUE)
mz_malawiGT3 = read.table('mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3', sep='\t', header=TRUE)
pn_victoriaGT2 = read.table('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3', sep='\t', header=TRUE)
pn_victoriaGT3 = read.table('pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3', sep='\t', header=TRUE)

### A. Lake Malawi species
# this creates a dataframe of the lake malawi species and clades
Lake_Malawi_Species <- c("Nbrichardi","Rhamphochromis_longiceps","Rhamphochromis_esox","Rhamphochromis_woodi","Diplotaxodon_limnothrissa","Diplotaxodon_macrops","Diplotaxodon_macrops_black_dorsal","Pallidochromis_tokolosh","Diplotaxodon_greenwoodi","Diplotaxodon_white_back_similis","Diplotaxodon_ngulube","Labeotropheus_trewavasae","Genyochromis_mento","Cynotilapia_axelrodi","Cynotilapia_afra","Iodotropheus_sprengerae","Tropheops_tropheops","Petrotilapia_genalutea","Metriaclima_zebra","A_calliptera_Chizumulu","A_calliptera_Luwawa","A_calliptera_Bua","A_calliptera_Salima_Father","A_calliptera_Salima_Mother","A_calliptera_Salima_Offspring","A_calliptera_Chitimba","A_calliptera_Enukweni","A_calliptera_South_Rukuru","A_calliptera_North_Rukuru","A_calliptera_Near_Kyela","A_calliptera_Songwe_River","Astatotilapia_Kingiri","A_calliptera_Mbaka_river_Female","Massoko_benthic_HC_1","Massoko_littoral_HC_1","A_calliptera_Lake_Chilwa","A_calliptera_Rovuma_river_Female","A_calliptera_Lake_Chidya","A_calliptera_Kitai_Dam","A_calliptera_Upper_Rovuma","Copadichromis_quadrimaculatus","Copadichromis_virginalis","Copadichromis_virginalis_2","Copadichromis_virginalis_3","Copadichromis_virginalis_4","Copadichromis_virginalis_5","Copadichromis_likomae","Copadichromis_trimaculatus","Copadichromis_mloto","Aulonocara_steveni","Aulonocara_stuartgranti_Father","Aulonocara_stuartgranti_Mother","Aulonocara_stuartgranti_Offspring","Alticorpus_geoffreyi","Lethrinops_longimanus_redhead","Lethrinops_gossei","Lethrinops_sp_oliveri","Aulonocara_yellow","Aulonocara_minutus","Alticorpus_macrocleithrum","Taeniolethrinops_praeorbitalis","Taeniolethrinops_macrorhynchus","Taeniolethrinops_furcicauda","Lethrinops_lethrinus","Lethrinops_lethrinus_Mazinzi_Father","Lethrinops_lethrinus_Mazinzi_Mother","Lethrinops_lethrinus_Mazinzi_Offspring","Lethrinops_auritus","Lethrinops_albus","Placidochromis_electra","Mylochromis_anaphyrmus","Mylochromis_anaphyrmus_2","Mylochromis_anaphyrmus_3","Mylochromis_anaphyrmus_4","Mylochromis_anaphyrmus_5","Placidochromis_longimanus_1","Placidochromis_longimanus_2","Placidochromis_longimanus_3","Placidochromis_longimanus_4","Placidochromis_longimanus_5","Ctenopharynx_intermedius_1","Ctenopharynx_intermedius_2","Ctenopharynx_intermedius_3","Ctenopharynx_nitidus","Ctenopharynx_nitidus_2","Otopharynx_tetrastigma","Otopharynx_speciosus_1","Otopharynx_speciosus_2","Buccochromis_rhoadesii","Buccochromis_nototaenia","Taeniochromis_holotaenia","Stigmatochromis_modestus","Stigmatochromis_guttatus","Otopharynx_brooksi_nkhata","Trematocranus_placodon_2","Trematocranus_placodon_3","Trematocranus_placodon_4","Trematocranus_placodon_5","Mylochromis_melanotaenia","Tyrannochromis_nigriventer","Champsochromis_caeruelus_1","Champsochromis_caeruelus_2","Nimbochromis_linni","Nimbochromis_polystigma","Nimbochromis_livingstoni","Dimidiochromis_kiwinge","Dimidiochromis_compressiceps","Dimidiochromis_dimidiatus","Dimidiochromis_strigatus","Dimidiochromis_strigatus_2","Mylochromis_ericotaenia","Otopharynx_lithobates","Placidochromis_milomo","Hemitilapia_oxyrhynchus","Hemitilapia_oxyrhynchus_2","Tremitochranus_placodon","Chilotilapia_rhoadesii_1","Chilotilapia_rhoadesii_2","Chilotilapia_rhoadesii_3","Chilotilapia_rhoadesii_4","Protomelas_ornatus_1","Protomelas_ornatus_2","Hemitaenichromis_spilopterus","Hemitaeniochromis_spilopterus_2","Placidochromis_johnstoni","Fossorochromis_rostratus_3","Fossorochromis_rostratus_4","Placidochromis_subocularis_1","Placidochromis_subocularis_2","Placidochromis_subocularis_3","Placidochromis_subocularis_4","Placidochromis_subocularis_5","Placidochromis_subocularis_6","Placidochromis_subocularis_7","Placidochromis_subocularis_8")
Lake_Malawi_Clades <- c("Lake_Tanganyika","Rhamphochromis","Rhamphochromis","Rhamphochromis","Diplotaxodon","Diplotaxodon","Diplotaxodon","Diplotaxodon","Diplotaxodon","Diplotaxodon","Diplotaxodon","Mbuna","Mbuna","Mbuna","Mbuna","Mbuna","Mbuna","Mbuna","Mbuna","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","Utaka","Utaka","Utaka","Utaka","Utaka","Utaka","Utaka","Utaka","Utaka","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic")
Lake_Malawi_Species_and_Clades <- data.frame(Lake_Malawi_Species, Lake_Malawi_Clades)

# these are the clades for when the species are in alphabetical order like in the dfs
Lake_Malawi_Clades2 <- c("A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","A_calliptera","Deep_benthic","Deep_benthic","A_calliptera","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Deep_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Utaka","Utaka","Utaka","Utaka","Utaka","Utaka","Utaka","Utaka","Utaka","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Mbuna","Mbuna","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Diplotaxodon","Diplotaxodon","Diplotaxodon","Diplotaxodon","Diplotaxodon","Diplotaxodon","Shallow_benthic","Shallow_benthic","Mbuna","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Mbuna","Mbuna","Shallow_benthic","Shallow_benthic","Deep_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Deep_benthic","Deep_benthic","A_calliptera","A_calliptera","Mbuna","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Diplotaxodon","Mbuna","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Rhamphochromis","Rhamphochromis","Rhamphochromis","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Shallow_benthic","Mbuna","Shallow_benthic","Lake_Tanganyika")

# this is for when the lake malawi species have the same genotype as M. zebra
mz_malawiGT2_refsame = t(as.data.frame(lapply(mz_malawiGT2, function(x){ length(which(x==1))})))
mz_malawiGT2_refsame2 = as.data.frame(mz_malawiGT2_refsame)
mz_malawiGT2_refsame2$species <- rownames(mz_malawiGT2_refsame2)
rownames(mz_malawiGT2_refsame2) <- NULL
mz_malawiGT2_refsame2$species <- gsub("_GT.PA", "", mz_malawiGT2_refsame2$species)
mz_malawiGT2_refsame3 <- mz_malawiGT2_refsame2[307:441,]
names(mz_malawiGT2_refsame3) <- c("Genotype count", "Species")
mz_malawiGT2_refsame3 <- cbind(mz_malawiGT2_refsame3, Lake_Malawi_Clades2, Source = 'Species genotype same as M. zebra')
mz_malawiGT2_refsame3$Species <- factor(mz_malawiGT2_refsame3$Species, levels=c("Nbrichardi","Rhamphochromis_longiceps","Rhamphochromis_esox","Rhamphochromis_woodi","Diplotaxodon_limnothrissa","Diplotaxodon_macrops","Diplotaxodon_macrops_black_dorsal","Pallidochromis_tokolosh","Diplotaxodon_greenwoodi","Diplotaxodon_white_back_similis","Diplotaxodon_ngulube","Labeotropheus_trewavasae","Genyochromis_mento","Cynotilapia_axelrodi","Cynotilapia_afra","Iodotropheus_sprengerae","Tropheops_tropheops","Petrotilapia_genalutea","Metriaclima_zebra","A_calliptera_Chizumulu","A_calliptera_Luwawa","A_calliptera_Bua","A_calliptera_Salima_Father","A_calliptera_Salima_Mother","A_calliptera_Salima_Offspring","A_calliptera_Chitimba","A_calliptera_Enukweni","A_calliptera_South_Rukuru","A_calliptera_North_Rukuru","A_calliptera_Near_Kyela","A_calliptera_Songwe_River","Astatotilapia_Kingiri","A_calliptera_Mbaka_river_Female","Massoko_benthic_HC_1","Massoko_littoral_HC_1","A_calliptera_Lake_Chilwa","A_calliptera_Rovuma_river_Female","A_calliptera_Lake_Chidya","A_calliptera_Kitai_Dam","A_calliptera_Upper_Rovuma","Copadichromis_quadrimaculatus","Copadichromis_virginalis","Copadichromis_virginalis_2","Copadichromis_virginalis_3","Copadichromis_virginalis_4","Copadichromis_virginalis_5","Copadichromis_likomae","Copadichromis_trimaculatus","Copadichromis_mloto","Aulonocara_steveni","Aulonocara_stuartgranti_Father","Aulonocara_stuartgranti_Mother","Aulonocara_stuartgranti_Offspring","Alticorpus_geoffreyi","Lethrinops_longimanus_redhead","Lethrinops_gossei","Lethrinops_sp_oliveri","Aulonocara_yellow","Aulonocara_minutus","Alticorpus_macrocleithrum","Taeniolethrinops_praeorbitalis","Taeniolethrinops_macrorhynchus","Taeniolethrinops_furcicauda","Lethrinops_lethrinus","Lethrinops_lethrinus_Mazinzi_Father","Lethrinops_lethrinus_Mazinzi_Mother","Lethrinops_lethrinus_Mazinzi_Offspring","Lethrinops_auritus","Lethrinops_albus","Placidochromis_electra","Mylochromis_anaphyrmus","Mylochromis_anaphyrmus_2","Mylochromis_anaphyrmus_3","Mylochromis_anaphyrmus_4","Mylochromis_anaphyrmus_5","Placidochromis_longimanus_1","Placidochromis_longimanus_2","Placidochromis_longimanus_3","Placidochromis_longimanus_4","Placidochromis_longimanus_5","Ctenopharynx_intermedius_1","Ctenopharynx_intermedius_2","Ctenopharynx_intermedius_3","Ctenopharynx_nitidus","Ctenopharynx_nitidus_2","Otopharynx_tetrastigma","Otopharynx_speciosus_1","Otopharynx_speciosus_2","Buccochromis_rhoadesii","Buccochromis_nototaenia","Taeniochromis_holotaenia","Stigmatochromis_modestus","Stigmatochromis_guttatus","Otopharynx_brooksi_nkhata","Trematocranus_placodon_2","Trematocranus_placodon_3","Trematocranus_placodon_4","Trematocranus_placodon_5","Mylochromis_melanotaenia","Tyrannochromis_nigriventer","Champsochromis_caeruelus_1","Champsochromis_caeruelus_2","Nimbochromis_linni","Nimbochromis_polystigma","Nimbochromis_livingstoni","Dimidiochromis_kiwinge","Dimidiochromis_compressiceps","Dimidiochromis_dimidiatus","Dimidiochromis_strigatus","Dimidiochromis_strigatus_2","Mylochromis_ericotaenia","Otopharynx_lithobates","Placidochromis_milomo","Hemitilapia_oxyrhynchus","Hemitilapia_oxyrhynchus_2","Tremitochranus_placodon","Chilotilapia_rhoadesii_1","Chilotilapia_rhoadesii_2","Chilotilapia_rhoadesii_3","Chilotilapia_rhoadesii_4","Protomelas_ornatus_1","Protomelas_ornatus_2","Hemitaenichromis_spilopterus","Hemitaeniochromis_spilopterus_2","Placidochromis_johnstoni","Fossorochromis_rostratus_3","Fossorochromis_rostratus_4","Placidochromis_subocularis_1","Placidochromis_subocularis_2","Placidochromis_subocularis_3","Placidochromis_subocularis_4","Placidochromis_subocularis_5","Placidochromis_subocularis_6","Placidochromis_subocularis_7","Placidochromis_subocularis_8")) # Change the lake species order according to phylogeny
# this is for when the lake malawi species have a different genotype to M. zebra
mz_malawiGT2_refdiff = t(as.data.frame(lapply(mz_malawiGT2, function(x){ length(which(x==0))})))
mz_malawiGT2_refdiff2 = as.data.frame(mz_malawiGT2_refdiff)
mz_malawiGT2_refdiff2$species <- rownames(mz_malawiGT2_refdiff2)
rownames(mz_malawiGT2_refdiff2) <- NULL
mz_malawiGT2_refdiff2$species <- gsub("_GT.PA", "", mz_malawiGT2_refdiff2$species)
mz_malawiGT2_refdiff3 <- mz_malawiGT2_refdiff2[307:441,]
names(mz_malawiGT2_refdiff3) <- c("Genotype count", "Species")
mz_malawiGT2_refdiff3 <- cbind(mz_malawiGT2_refdiff3, Lake_Malawi_Clades2, Source = 'Species genotype different to M. zebra')
mz_malawiGT2_refdiff3$Species <- factor(mz_malawiGT2_refdiff3$Species, levels=c("Nbrichardi","Rhamphochromis_longiceps","Rhamphochromis_esox","Rhamphochromis_woodi","Diplotaxodon_limnothrissa","Diplotaxodon_macrops","Diplotaxodon_macrops_black_dorsal","Pallidochromis_tokolosh","Diplotaxodon_greenwoodi","Diplotaxodon_white_back_similis","Diplotaxodon_ngulube","Labeotropheus_trewavasae","Genyochromis_mento","Cynotilapia_axelrodi","Cynotilapia_afra","Iodotropheus_sprengerae","Tropheops_tropheops","Petrotilapia_genalutea","Metriaclima_zebra","A_calliptera_Chizumulu","A_calliptera_Luwawa","A_calliptera_Bua","A_calliptera_Salima_Father","A_calliptera_Salima_Mother","A_calliptera_Salima_Offspring","A_calliptera_Chitimba","A_calliptera_Enukweni","A_calliptera_South_Rukuru","A_calliptera_North_Rukuru","A_calliptera_Near_Kyela","A_calliptera_Songwe_River","Astatotilapia_Kingiri","A_calliptera_Mbaka_river_Female","Massoko_benthic_HC_1","Massoko_littoral_HC_1","A_calliptera_Lake_Chilwa","A_calliptera_Rovuma_river_Female","A_calliptera_Lake_Chidya","A_calliptera_Kitai_Dam","A_calliptera_Upper_Rovuma","Copadichromis_quadrimaculatus","Copadichromis_virginalis","Copadichromis_virginalis_2","Copadichromis_virginalis_3","Copadichromis_virginalis_4","Copadichromis_virginalis_5","Copadichromis_likomae","Copadichromis_trimaculatus","Copadichromis_mloto","Aulonocara_steveni","Aulonocara_stuartgranti_Father","Aulonocara_stuartgranti_Mother","Aulonocara_stuartgranti_Offspring","Alticorpus_geoffreyi","Lethrinops_longimanus_redhead","Lethrinops_gossei","Lethrinops_sp_oliveri","Aulonocara_yellow","Aulonocara_minutus","Alticorpus_macrocleithrum","Taeniolethrinops_praeorbitalis","Taeniolethrinops_macrorhynchus","Taeniolethrinops_furcicauda","Lethrinops_lethrinus","Lethrinops_lethrinus_Mazinzi_Father","Lethrinops_lethrinus_Mazinzi_Mother","Lethrinops_lethrinus_Mazinzi_Offspring","Lethrinops_auritus","Lethrinops_albus","Placidochromis_electra","Mylochromis_anaphyrmus","Mylochromis_anaphyrmus_2","Mylochromis_anaphyrmus_3","Mylochromis_anaphyrmus_4","Mylochromis_anaphyrmus_5","Placidochromis_longimanus_1","Placidochromis_longimanus_2","Placidochromis_longimanus_3","Placidochromis_longimanus_4","Placidochromis_longimanus_5","Ctenopharynx_intermedius_1","Ctenopharynx_intermedius_2","Ctenopharynx_intermedius_3","Ctenopharynx_nitidus","Ctenopharynx_nitidus_2","Otopharynx_tetrastigma","Otopharynx_speciosus_1","Otopharynx_speciosus_2","Buccochromis_rhoadesii","Buccochromis_nototaenia","Taeniochromis_holotaenia","Stigmatochromis_modestus","Stigmatochromis_guttatus","Otopharynx_brooksi_nkhata","Trematocranus_placodon_2","Trematocranus_placodon_3","Trematocranus_placodon_4","Trematocranus_placodon_5","Mylochromis_melanotaenia","Tyrannochromis_nigriventer","Champsochromis_caeruelus_1","Champsochromis_caeruelus_2","Nimbochromis_linni","Nimbochromis_polystigma","Nimbochromis_livingstoni","Dimidiochromis_kiwinge","Dimidiochromis_compressiceps","Dimidiochromis_dimidiatus","Dimidiochromis_strigatus","Dimidiochromis_strigatus_2","Mylochromis_ericotaenia","Otopharynx_lithobates","Placidochromis_milomo","Hemitilapia_oxyrhynchus","Hemitilapia_oxyrhynchus_2","Tremitochranus_placodon","Chilotilapia_rhoadesii_1","Chilotilapia_rhoadesii_2","Chilotilapia_rhoadesii_3","Chilotilapia_rhoadesii_4","Protomelas_ornatus_1","Protomelas_ornatus_2","Hemitaenichromis_spilopterus","Hemitaeniochromis_spilopterus_2","Placidochromis_johnstoni","Fossorochromis_rostratus_3","Fossorochromis_rostratus_4","Placidochromis_subocularis_1","Placidochromis_subocularis_2","Placidochromis_subocularis_3","Placidochromis_subocularis_4","Placidochromis_subocularis_5","Placidochromis_subocularis_6","Placidochromis_subocularis_7","Placidochromis_subocularis_8")) # Change the lake species order according to phylogeny
# join the above
mz_malawiGT2_refsamediff <- rbind(mz_malawiGT2_refsame3,mz_malawiGT2_refdiff3)
mz_malawiGT2_refsamediff$Lake_Malawi_Clades2 <- factor(mz_malawiGT2_refsamediff$Lake_Malawi_Clades2, levels=c("Lake_Tanganyika","Rhamphochromis","Diplotaxodon","Mbuna","A_calliptera","Utaka","Deep_benthic","Shallow_benthic")) # Change the clade order according to the ASTRAL phylogeny in Malinsky 2018
# remove the M. zebra vcf entries
mz_malawiGT3_refsamediff <- mz_malawiGT2_refsamediff[-c(83,218),]


g1 <- ggplot(data=mz_malawiGT3_refsamediff, aes(x=Species, y=`Genotype count`, fill=Source)) +
  geom_bar(stat="identity") +
  labs(y="Genotype count", x="Lake Malawi species",title="Lake Malawi species genotype counts") +
  scale_fill_brewer(palette="Paired")+
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,5708,500)) + 
  theme(plot.title = element_text(hjust = 0.5, size=22,face = "bold")) +
  theme(text = element_text(size=14), axis.text.x = element_text(size=12,angle = 90, hjust = 1, vjust=0.5),axis.text.y=element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16)) +
  theme(legend.position="bottom") +
  facet_grid(.~Lake_Malawi_Clades2,scales = "free_x", space = "free_x") +
  theme(strip.background = element_rect(color="black", fill="#e8f5e9", size=0.5, linetype="solid")) +
  theme(strip.text.x = element_text(size = 14, colour = "black", angle = 90, face="bold")) +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))


### B. Lake Victoria species
Lake_Victoria_Species = read.table('~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_Victoria_Data_Ole2018/List_of_species.txt', sep='\t', header=TRUE)

# this is for when the lake victoria species have the same genotype as P. nyererei
pn_victoriaGT2_refsame = t(as.data.frame(lapply(pn_victoriaGT2, function(x){ length(which(x==1))})))
pn_victoriaGT2_refsame2 = as.data.frame(pn_victoriaGT2_refsame)
pn_victoriaGT2_refsame2$species <- rownames(pn_victoriaGT2_refsame2)
rownames(pn_victoriaGT2_refsame2) <- NULL
pn_victoriaGT2_refsame2$species <- gsub("_GT.PA", "", pn_victoriaGT2_refsame2$species)
pn_victoriaGT2_refsame2$species <- gsub("X", "", pn_victoriaGT2_refsame2$species)
pn_victoriaGT2_refsame3 <- pn_victoriaGT2_refsame2[243:345,]
names(pn_victoriaGT2_refsame3) <- c("Genotype count", "Species_ID")
# pn_victoriaGT2_refsame3a <- merge(pn_victoriaGT2_refsame3,Lake_Victoria_Species, by.x="Species_ID", by.y="fileName") # this merges and hence, changes the order
victoria_sp_match <- (Lake_Victoria_Species[match(pn_victoriaGT2_refsame3[,2],Lake_Victoria_Species[,1]),(c(3,5,6,7))]) # this matches the Species_ID to species names and habitats
pn_victoriaGT2_refsame3a <- cbind(pn_victoriaGT2_refsame3,victoria_sp_match,Source = 'Species genotype same as P. nyererei')

# this is for when the lake victoria species have a different genotype to P. nyererei
pn_victoriaGT2_refdiff = t(as.data.frame(lapply(pn_victoriaGT2, function(x){ length(which(x==0))})))
pn_victoriaGT2_refdiff2 = as.data.frame(pn_victoriaGT2_refdiff)
pn_victoriaGT2_refdiff2$species <- rownames(pn_victoriaGT2_refdiff2)
rownames(pn_victoriaGT2_refdiff2) <- NULL
pn_victoriaGT2_refdiff2$species <- gsub("_GT.PA", "", pn_victoriaGT2_refdiff2$species)
pn_victoriaGT2_refdiff2$species <- gsub("X", "", pn_victoriaGT2_refdiff2$species)
pn_victoriaGT2_refdiff3 <- pn_victoriaGT2_refdiff2[243:345,]
names(pn_victoriaGT2_refdiff3) <- c("Genotype count", "Species_ID")
# pn_victoriaGT2_refdiff3a <- merge(pn_victoriaGT2_refdiff3,Lake_Victoria_Species, by.x="Species_ID", by.y="fileName") # this merges and hence, changes the order
victoria_sp_match <- (Lake_Victoria_Species[match(pn_victoriaGT2_refdiff3[,2],Lake_Victoria_Species[,1]),(c(3,5,6,7))]) # this matches the Species_ID to species names and habitats
pn_victoriaGT2_refdiff3a <- cbind(pn_victoriaGT2_refdiff3,victoria_sp_match,Source = 'Species genotype different to P. nyererei')

# join the above
pn_victoriaGT2_refsamediff <- rbind(pn_victoriaGT2_refsame3a,pn_victoriaGT2_refdiff3a)
pn_victoriaGT2_refsamediff_a <- pn_victoriaGT2_refsamediff[- grep("Unknown", pn_victoriaGT2_refsamediff$GeneralLocation),] # remove the odd species with few counts and unknown ID

g2 <- ggplot(data=pn_victoriaGT2_refsamediff_a, aes(x=Species_name, y=`Genotype count`, fill=Source)) +
  geom_bar(stat="identity") +
  labs(y="Genotype count", x="Lake Victoria species",title="Lake Victoria species genotype counts") +
  scale_fill_brewer(palette="Paired")+
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,17009,1000)) + 
  theme(plot.title = element_text(hjust = 0.5, size=22,face = "bold")) +
  theme(text = element_text(size=14), axis.text.x = element_text(size=12,angle = 90, hjust = 1, vjust=0.5),axis.text.y=element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16)) +
  theme(legend.position="bottom") +
  facet_grid(.~Substrate,scales = "free_x", space = "free_x") +
  theme(strip.background = element_rect(color="black", fill="#e8f5e9", size=0.5, linetype="solid")) +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 90, face="bold")) +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

## Output tiff files of above ggplots

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/SNPs/FigS-R3fA_LakeMalawi_species_genotype_counts.tiff', units="in", width=20, height=10, res=200)
g1
dev.off()

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/SNPs/FigS-R3fB_LakeVictoria_species_genotype_counts.tiff', units="in", width=20, height=10, res=200)
g2
dev.off()


## over certain ecological substrates in LV species, calculating averages to determine rough distribution
pn_victoriaGT2_refsame_average_list = list()
for (i in c("Mud","Open","Plants","Rock","Sand")) {
  pn_victoriaGT2_refsame_average_list[[i]] <- with(pn_victoriaGT2_refsamediff, mean(`Genotype count`[Substrate == i & Source == "Species genotype same as P. nyererei"]) )
}
pn_victoriaGT2_refdiff_average_list = list()
for (i in c("Mud","Open","Plants","Rock","Sand")) {
  pn_victoriaGT2_refdiff_average_list[[i]] <- with(pn_victoriaGT2_refsamediff, mean(`Genotype count`[Substrate == i & Source == "Species genotype different to P. nyererei"]) )
}
# pn_victoriaGT2_refsame_average = do.call(rbind, pn_victoriaGT2_refsame_average_list)
pn_victoriaGT2_refsame_average <- as.data.frame(pn_victoriaGT2_refsame_average_list)
pn_victoriaGT2_refsame_average2 <- as.data.frame(t(pn_victoriaGT2_refsame_average))
pn_victoriaGT2_refsame_average3 <- add_rownames(pn_victoriaGT2_refsame_average2, "Ecological substrate")
names(pn_victoriaGT2_refsame_average3) <- c("Ecological substrate", "Average genotype count")
pn_victoriaGT2_refsame_average4 <- cbind(pn_victoriaGT2_refsame_average3, Source = 'Species genotype same as P. nyererei')

pn_victoriaGT2_refdiff_average <- as.data.frame(pn_victoriaGT2_refdiff_average_list)
pn_victoriaGT2_refdiff_average2 <- as.data.frame(t(pn_victoriaGT2_refdiff_average))
pn_victoriaGT2_refdiff_average3 <- add_rownames(pn_victoriaGT2_refdiff_average2, "Ecological substrate")
names(pn_victoriaGT2_refdiff_average3) <- c("Ecological substrate", "Average genotype count")
pn_victoriaGT2_refdiff_average4 <- cbind(pn_victoriaGT2_refdiff_average3, Source = 'Species genotype different to P. nyererei')

pn_victoriaGT2_refsamediff_average <- rbind(pn_victoriaGT2_refsame_average4,pn_victoriaGT2_refdiff_average4)


# 3. Prepare indvidual plots for gene-specific TFBS SNPs and genotypes that show:
# a. the phylogeny - colouring each clade (similar to Malawi sws1 nr2c2/rxrb tree)
# b. the genotype in each lake species
# c. the foraging habit
# d. the habitat
# e. in visual opsins, the wavelength palette

Mz_RefSp = read.table("mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3", header = TRUE,stringsAsFactors=FALSE)
Mz_AltSp = read.table("mz-Candidates_TFBSsSNPs_Malawi_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3", header = TRUE,stringsAsFactors=FALSE)
Pn_RefSp = read.table("pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.RefSpMatch.bed3", header = TRUE,stringsAsFactors=FALSE)
Pn_AltSp = read.table("pn-Candidates_TFBSsSNPs_Victoria_filteredSNPs_SNPsInPromAln.0718.SpGenotypes.AltSpMatch.bed3", header = TRUE,stringsAsFactors=FALSE)
                      
# created a separate file with more details including clade, foraging, habitat
LMSpecies2 = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_MalawiVCF/Lake_Malawi_species_info.txt", header = TRUE,stringsAsFactors=FALSE)

# build the LM species tree
# https://4va.github.io/biodatasci/r-ggtree.html # guide here
# newick tree manually generated
library(tidyverse)
library(ggtree)
LM_tree <- read.tree("~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_MalawiVCF/Lake_Malawi_species_phylo.nwk")

p <- ggtree(LM_tree) + 
  geom_tiplab(offset = 4.5) + 
  geom_cladelabel(node=179, label="Shallow benthic",
                  color="red", offset=20, align=TRUE) +
  geom_cladelabel(node=53, label="Utaka",
                  color="darkgreen", offset=20, align=TRUE) +
  geom_cladelabel(node=166, label="Deep benthic", 
                  color="blue", offset=20, align=TRUE) + 
  geom_cladelabel(node=160, label="Utaka", 
                  color="darkgreen", offset=20, align=TRUE) +
  geom_cladelabel(node=139, label="A. calliptera", 
                  color="lightgreen", offset=20, align=TRUE) +
  geom_cladelabel(node=131, label="Mbuna", 
                  color="purple", offset=20, align=TRUE) +
  geom_cladelabel(node=123, label="Diplotaxodon", 
                  color="orange", offset=20, align=TRUE) +
  geom_cladelabel(node=121, label="Rhamphochromis", 
                  color="brown", offset=20, align=TRUE) +
  geom_cladelabel(node=1, label="Lake Tanganyika", 
                  color="pink", offset=20, align=TRUE) +
  theme_tree2() + 
  xlim(0, 70) + 
  theme_tree()
plot(p)

# + geom_text(aes(label=node), hjust=-.3) # this will let you see the node numbers

# load in Lake Malawi species details
LM_Astral <- read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_MalawiVCF/Lake_Malawi_species_info.txt", header = TRUE,stringsAsFactors=FALSE)


## A. nr2c2>sws1
# select Mz nr2c2>sws1 and genotype columns only - other cols could be 15-TG; 21-TF; 31-RefGT; 32-AltGT
Mz_nr2c2_sws1_GT <- as.data.frame(t(Mz_RefSp[grep("nr2c2", Mz_RefSp$TFmotif_gene_symbol),c(172:306)]),stringsAsFactors=FALSE)
Mz_nr2c2_sws1_GT2 <- add_rownames(Mz_nr2c2_sws1_GT, "Species")
Mz_nr2c2_sws1_GT2$Species <- gsub("_GT", "", Mz_nr2c2_sws1_GT2$Species)
names(Mz_nr2c2_sws1_GT2) <- c("Species", "Genotype")
Mz_nr2c2_sws1_GT3 <- cbind(Mz_nr2c2_sws1_GT2, Target_Gene = 'sws1', TF = 'Nr2c2')
# remove the following species that are not in the tree:
# Labeotropheus_trewavasae
# A_calliptera_Lake_Chilwa
# A_calliptera_Rovuma_river_Female
# A_calliptera_Lake_Chidya
# A_calliptera_Kitai_Dam
# A_calliptera_Upper_Rovuma
# Copadichromis_mloto
# Tremitochranus_placodon
# Fossorochromis_rostratus_3
# Fossorochromis_rostratus_4
# Placidochromis_subocularis_1
# Placidochromis_subocularis_2
# Placidochromis_subocularis_3
# Placidochromis_subocularis_4
# Placidochromis_subocularis_5
# Placidochromis_subocularis_6
# Placidochromis_subocularis_7
# Placidochromis_subocularis_8
removesp <- c("Labeotropheus_trewavasae","A_calliptera_Lake_Chilwa","A_calliptera_Rovuma_river_Female","A_calliptera_Lake_Chidya","A_calliptera_Kitai_Dam","A_calliptera_Upper_Rovuma","Copadichromis_mloto","Tremitochranus_placodon","Fossorochromis_rostratus_3","Fossorochromis_rostratus_4","Placidochromis_subocularis_1","Placidochromis_subocularis_2","Placidochromis_subocularis_3","Placidochromis_subocularis_4","Placidochromis_subocularis_5","Placidochromis_subocularis_6","Placidochromis_subocularis_7","Placidochromis_subocularis_8")
Mz_nr2c2_sws1_GT4 <- Mz_nr2c2_sws1_GT3[!grepl(paste(removesp, collapse="|"), Mz_nr2c2_sws1_GT3$Species),] # this removes rows of the species above
# add the foraging and ecology of each species:
Mz_nr2c2_sws1_GT5 <- (LM_Astral[match(Mz_nr2c2_sws1_GT4[,1],LM_Astral[,1]),(c(2,3,4,5))]) # this matches the Species and outputs ecology
Mz_nr2c2_sws1_GT6 <- cbind(Mz_nr2c2_sws1_GT4,Mz_nr2c2_sws1_GT5)
# Create column of star markers for habitat
Mz_nr2c2_sws1_GT6$Habitat_marker[Mz_nr2c2_sws1_GT6$Habitat == "Demersal"]<-"****"
Mz_nr2c2_sws1_GT6$Habitat_marker[Mz_nr2c2_sws1_GT6$Habitat == "Benthopelagic"]<-"***"
Mz_nr2c2_sws1_GT6$Habitat_marker[Mz_nr2c2_sws1_GT6$Habitat == "Pelagic"]<-"**"
Mz_nr2c2_sws1_GT6$Habitat_marker[Mz_nr2c2_sws1_GT6$Habitat == "Rock"]<-"*"
# change all NA to N/A:
Mz_nr2c2_sws1_GT6[is.na(Mz_nr2c2_sws1_GT6)] <- "N/A"
# order the factor levels according to the phylogeny order:
Mz_nr2c2_sws1_GT6$Species <- factor(Mz_nr2c2_sws1_GT6$Species, levels=c("Nbrichardi","Rhamphochromis_longiceps","Rhamphochromis_esox","Rhamphochromis_woodi","Diplotaxodon_limnothrissa","Diplotaxodon_macrops","Diplotaxodon_macrops_black_dorsal","Pallidochromis_tokolosh","Diplotaxodon_greenwoodi","Diplotaxodon_white_back_similis","Diplotaxodon_ngulube","Genyochromis_mento","Cynotilapia_axelrodi","Cynotilapia_afra","Iodotropheus_sprengerae","Tropheops_tropheops","Petrotilapia_genalutea","Metriaclima_zebra","A_calliptera_Chizumulu","A_calliptera_Luwawa","A_calliptera_Bua","A_calliptera_Salima_Father","A_calliptera_Salima_Mother","A_calliptera_Salima_Offspring","A_calliptera_Chitimba","A_calliptera_Enukweni","A_calliptera_South_Rukuru","A_calliptera_North_Rukuru","A_calliptera_Near_Kyela","A_calliptera_Songwe_River","Astatotilapia_Kingiri","A_calliptera_Mbaka_river_Female","Massoko_benthic_HC_1","Massoko_littoral_HC_1","Copadichromis_quadrimaculatus","Copadichromis_virginalis","Copadichromis_virginalis_2","Copadichromis_virginalis_3","Copadichromis_virginalis_4","Copadichromis_virginalis_5","Copadichromis_likomae","Aulonocara_steveni","Aulonocara_stuartgranti_Father","Aulonocara_stuartgranti_Mother","Aulonocara_stuartgranti_Offspring","Alticorpus_geoffreyi","Lethrinops_longimanus_redhead","Lethrinops_gossei","Lethrinops_sp_oliveri","Aulonocara_yellow","Aulonocara_minutus","Alticorpus_macrocleithrum","Copadichromis_trimaculatus","Taeniolethrinops_praeorbitalis","Taeniolethrinops_macrorhynchus","Taeniolethrinops_furcicauda","Lethrinops_lethrinus","Lethrinops_lethrinus_Mazinzi_Father","Lethrinops_lethrinus_Mazinzi_Mother","Lethrinops_lethrinus_Mazinzi_Offspring","Lethrinops_auritus","Lethrinops_albus","Placidochromis_electra","Mylochromis_anaphyrmus","Mylochromis_anaphyrmus_2","Mylochromis_anaphyrmus_3","Mylochromis_anaphyrmus_4","Mylochromis_anaphyrmus_5","Placidochromis_longimanus_1","Placidochromis_longimanus_2","Placidochromis_longimanus_3","Placidochromis_longimanus_4","Placidochromis_longimanus_5","Ctenopharynx_intermedius_1","Ctenopharynx_intermedius_2","Ctenopharynx_intermedius_3","Ctenopharynx_nitidus","Ctenopharynx_nitidus_2","Otopharynx_tetrastigma","Otopharynx_speciosus_1","Otopharynx_speciosus_2","Buccochromis_rhoadesii","Buccochromis_nototaenia","Taeniochromis_holotaenia","Stigmatochromis_modestus","Stigmatochromis_guttatus","Otopharynx_brooksi_nkhata","Trematocranus_placodon_2","Trematocranus_placodon_3","Trematocranus_placodon_4","Trematocranus_placodon_5","Mylochromis_melanotaenia","Tyrannochromis_nigriventer","Champsochromis_caeruelus_1","Champsochromis_caeruelus_2","Nimbochromis_linni","Nimbochromis_polystigma","Nimbochromis_livingstoni","Dimidiochromis_kiwinge","Dimidiochromis_compressiceps","Dimidiochromis_dimidiatus","Dimidiochromis_strigatus","Dimidiochromis_strigatus_2","Mylochromis_ericotaenia","Otopharynx_lithobates","Placidochromis_milomo","Hemitilapia_oxyrhynchus","Hemitilapia_oxyrhynchus_2","Chilotilapia_rhoadesii_1","Chilotilapia_rhoadesii_2","Chilotilapia_rhoadesii_3","Chilotilapia_rhoadesii_4","Protomelas_ornatus_1","Protomelas_ornatus_2","Hemitaenichromis_spilopterus","Hemitaeniochromis_spilopterus_2","Placidochromis_johnstoni")) # Change the lake species order according to phylogeny
# order the genotypes and wavelength
Mz_nr2c2_sws1_GT6$Genotype <- factor(Mz_nr2c2_sws1_GT6$Genotype, levels=c("C|C","C|T","T|C","T|T","N/A"))
Mz_nr2c2_sws1_GT6$Wavelength_palette <- factor(Mz_nr2c2_sws1_GT6$Wavelength_palette, levels=c("Short","Medium","Long","N/A"))

# plot the above - legend:
# 1. phylo tree where nodes tips are:
# A. Colour = Foraging/Diet and Shape = Genotype
# B. Shape = Wavelength palette
# C. Marker = Habitat: **** = Demersal; *** = benthopelagic; ** = pelagic, * = rock
p1 <- p %<+% Mz_nr2c2_sws1_GT6 + 
  geom_text(aes(label=Habitat_marker), show.legend = FALSE, position = position_nudge(y = -0.05 ,x = 3)) +
  geom_point(aes(shape=Wavelength_palette), size=3, show.legend = TRUE, position = position_nudge(x = 1.5)) +
  geom_point(aes(color=Foraging.Diet, shape=Genotype), size=3, show.legend = TRUE, position = position_nudge(x = 0.5)) +
  scale_color_manual(breaks=c('Algae','Aufwuchs','Benthivore','Fish','Herbivore','Invertebrates','Mix','Waste','Zooplankton','N/A'),values = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', 'lavenderblush2', '#999999', 'darkolivegreen2')) +
  # scale_color_manual(values = c("Mix" = '#E69F00','Invertebrates' = '#56B4E9','Fish' = '#009E73','Zooplankton' = '#F0E442','Algae' = '#0072B2','Aufwuchs' = '#D55E00','Herbivore' = '#CC79A7','Benthivore' = 'green','Waste' = '#999999','N/A' = 'lavenderblush2')) +
  scale_shape_manual(values = c("C|C" = 16,'C|T'=18,'T|T' = 17,'T|C' = 15, 'Short' = 12, 'Medium' = 13, 'Long' = 14,'N/A' = 4)) +
  annotate("text", x=25.5, y=119, label= "A", fontface="bold") +
  annotate("text", x=26.5, y=119, label= "B", fontface="bold") +
  annotate("text", x=28, y=119, label= "C", fontface="bold") +
  theme(legend.position=c(.9, .6))
plot(p1)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/SNPs/FigS-R4fA_nr2c2_sws1_LakeMalawi-SNPs.tiff', units="in", width=14, height=16, res=300)
p1
dev.off()

## B. stat1>rho (in Ab, lost in Mz) - Not much variation so will not be using, plotted as an example
# For this, you will need to add A. burtoni to the tree and species info

# https://4va.github.io/biodatasci/r-ggtree.html # guide here
# newick tree manually generated
library(tidyverse)
library(ggtree)
LM_tree2 <- read.tree("~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_MalawiVCF/Lake_Malawi_species_phylo2.nwk")

p2 <- ggtree(LM_tree2) + 
  # geom_text(aes(label=node), hjust=-.3) + # this will let you see the node numbers
  geom_tiplab(offset = 3.5) + 
  geom_cladelabel(node=181, label="Shallow benthic",
                  color="red", offset=19, align=TRUE) +
  geom_cladelabel(node=55, label="Utaka",
                  color="darkgreen", offset=19, align=TRUE) +
  geom_cladelabel(node=168, label="Deep benthic", 
                  color="blue", offset=19, align=TRUE) + 
  geom_cladelabel(node=162, label="Utaka", 
                  color="darkgreen", offset=19, align=TRUE) +
  geom_cladelabel(node=141, label="A. calliptera", 
                  color="lightgreen", offset=19, align=TRUE) +
  geom_cladelabel(node=133, label="Mbuna", 
                  color="purple", offset=19, align=TRUE) +
  geom_cladelabel(node=125, label="Diplotaxodon", 
                  color="orange", offset=19, align=TRUE) +
  geom_cladelabel(node=123, label="Rhamphochromis", 
                  color="brown", offset=19, align=TRUE) +
  geom_cladelabel(node=3, label="Lakes and Rivers", 
                  color="grey", offset=19, align=TRUE) +
  geom_cladelabel(node=1, label="Lake Tanganyika", 
                  color="pink", offset=19, align=TRUE) +
  theme_tree2() + 
  xlim(0, 70) + 
  theme_tree()
plot(p2)

# load in Lake Malawi species details
LM_Astral2 <- read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_MalawiVCF/Lake_Malawi_species_info2.txt", header = TRUE,stringsAsFactors=FALSE)

# select Mz stat1>rho and genotype columns only - other cols could be 15-TG; 21-TF; 31-RefGT; 32-AltGT
Mz_stat1_rho_GT <- as.data.frame(Mz_RefSp[grep("stat1", Mz_RefSp$TFmotif_gene_symbol),c(1:306)])
Mz_stat1_rho_GTa <- as.data.frame(t(Mz_stat1_rho_GT[grep("rho", Mz_stat1_rho_GT$motif_genesymbolSp),c(172:306)]),stringsAsFactors=FALSE)
# add the A. burtoni genotype
Ab_stat1_rhoGT<-data.frame(row.names = "Aburtoni_GT","A/A")
names(Ab_stat1_rhoGT) <- ("2763")
Mz_stat1_rho_GTb <- rbind(Mz_stat1_rho_GTa,Ab_stat1_rhoGT)
Mz_stat1_rho_GT2 <- add_rownames(Mz_stat1_rho_GTb, "Species")
Mz_stat1_rho_GT2$Species <- gsub("_GT", "", Mz_stat1_rho_GT2$Species)
names(Mz_stat1_rho_GT2) <- c("Species", "Genotype")
Mz_stat1_rho_GT3 <- cbind(Mz_stat1_rho_GT2, Target_Gene = 'rho', TF = 'stat1')
# remove the following species that are not in the tree:
# Labeotropheus_trewavasae
# A_calliptera_Lake_Chilwa
# A_calliptera_Rovuma_river_Female
# A_calliptera_Lake_Chidya
# A_calliptera_Kitai_Dam
# A_calliptera_Upper_Rovuma
# Copadichromis_mloto
# Tremitochranus_placodon
# Fossorochromis_rostratus_3
# Fossorochromis_rostratus_4
# Placidochromis_subocularis_1
# Placidochromis_subocularis_2
# Placidochromis_subocularis_3
# Placidochromis_subocularis_4
# Placidochromis_subocularis_5
# Placidochromis_subocularis_6
# Placidochromis_subocularis_7
# Placidochromis_subocularis_8
removesp <- c("Labeotropheus_trewavasae","A_calliptera_Lake_Chilwa","A_calliptera_Rovuma_river_Female","A_calliptera_Lake_Chidya","A_calliptera_Kitai_Dam","A_calliptera_Upper_Rovuma","Copadichromis_mloto","Tremitochranus_placodon","Fossorochromis_rostratus_3","Fossorochromis_rostratus_4","Placidochromis_subocularis_1","Placidochromis_subocularis_2","Placidochromis_subocularis_3","Placidochromis_subocularis_4","Placidochromis_subocularis_5","Placidochromis_subocularis_6","Placidochromis_subocularis_7","Placidochromis_subocularis_8")
Mz_stat1_rho_GT4 <- Mz_stat1_rho_GT3[!grepl(paste(removesp, collapse="|"), Mz_stat1_rho_GT3$Species),] # this removes rows of the species above
# add the foraging and ecology of each species:
Mz_stat1_rho_GT5 <- (LM_Astral2[match(Mz_stat1_rho_GT4[,1],LM_Astral2[,1]),(c(2,3,4,5))]) # this matches the Species and outputs ecology
Mz_stat1_rho_GT6 <- cbind(Mz_stat1_rho_GT4,Mz_stat1_rho_GT5)
# Create column of star markers for habitat
Mz_stat1_rho_GT6$Habitat_marker[Mz_stat1_rho_GT6$Habitat == "Demersal"]<-"****"
Mz_stat1_rho_GT6$Habitat_marker[Mz_stat1_rho_GT6$Habitat == "Benthopelagic"]<-"***"
Mz_stat1_rho_GT6$Habitat_marker[Mz_stat1_rho_GT6$Habitat == "Pelagic"]<-"**"
Mz_stat1_rho_GT6$Habitat_marker[Mz_stat1_rho_GT6$Habitat == "Rock"]<-"*"
# change all NA to N/A:
Mz_stat1_rho_GT6[is.na(Mz_stat1_rho_GT6)] <- "N/A"
# order the factor levels according to the phylogeny order:
Mz_stat1_rho_GT6$Species <- factor(Mz_stat1_rho_GT6$Species, levels=c("Nbrichardi","Aburtoni","Rhamphochromis_longiceps","Rhamphochromis_esox","Rhamphochromis_woodi","Diplotaxodon_limnothrissa","Diplotaxodon_macrops","Diplotaxodon_macrops_black_dorsal","Pallidochromis_tokolosh","Diplotaxodon_greenwoodi","Diplotaxodon_white_back_similis","Diplotaxodon_ngulube","Genyochromis_mento","Cynotilapia_axelrodi","Cynotilapia_afra","Iodotropheus_sprengerae","Tropheops_tropheops","Petrotilapia_genalutea","Metriaclima_zebra","A_calliptera_Chizumulu","A_calliptera_Luwawa","A_calliptera_Bua","A_calliptera_Salima_Father","A_calliptera_Salima_Mother","A_calliptera_Salima_Offspring","A_calliptera_Chitimba","A_calliptera_Enukweni","A_calliptera_South_Rukuru","A_calliptera_North_Rukuru","A_calliptera_Near_Kyela","A_calliptera_Songwe_River","Astatotilapia_Kingiri","A_calliptera_Mbaka_river_Female","Massoko_benthic_HC_1","Massoko_littoral_HC_1","Copadichromis_quadrimaculatus","Copadichromis_virginalis","Copadichromis_virginalis_2","Copadichromis_virginalis_3","Copadichromis_virginalis_4","Copadichromis_virginalis_5","Copadichromis_likomae","Aulonocara_steveni","Aulonocara_stuartgranti_Father","Aulonocara_stuartgranti_Mother","Aulonocara_stuartgranti_Offspring","Alticorpus_geoffreyi","Lethrinops_longimanus_redhead","Lethrinops_gossei","Lethrinops_sp_oliveri","Aulonocara_yellow","Aulonocara_minutus","Alticorpus_macrocleithrum","Copadichromis_trimaculatus","Taeniolethrinops_praeorbitalis","Taeniolethrinops_macrorhynchus","Taeniolethrinops_furcicauda","Lethrinops_lethrinus","Lethrinops_lethrinus_Mazinzi_Father","Lethrinops_lethrinus_Mazinzi_Mother","Lethrinops_lethrinus_Mazinzi_Offspring","Lethrinops_auritus","Lethrinops_albus","Placidochromis_electra","Mylochromis_anaphyrmus","Mylochromis_anaphyrmus_2","Mylochromis_anaphyrmus_3","Mylochromis_anaphyrmus_4","Mylochromis_anaphyrmus_5","Placidochromis_longimanus_1","Placidochromis_longimanus_2","Placidochromis_longimanus_3","Placidochromis_longimanus_4","Placidochromis_longimanus_5","Ctenopharynx_intermedius_1","Ctenopharynx_intermedius_2","Ctenopharynx_intermedius_3","Ctenopharynx_nitidus","Ctenopharynx_nitidus_2","Otopharynx_tetrastigma","Otopharynx_speciosus_1","Otopharynx_speciosus_2","Buccochromis_rhoadesii","Buccochromis_nototaenia","Taeniochromis_holotaenia","Stigmatochromis_modestus","Stigmatochromis_guttatus","Otopharynx_brooksi_nkhata","Trematocranus_placodon_2","Trematocranus_placodon_3","Trematocranus_placodon_4","Trematocranus_placodon_5","Mylochromis_melanotaenia","Tyrannochromis_nigriventer","Champsochromis_caeruelus_1","Champsochromis_caeruelus_2","Nimbochromis_linni","Nimbochromis_polystigma","Nimbochromis_livingstoni","Dimidiochromis_kiwinge","Dimidiochromis_compressiceps","Dimidiochromis_dimidiatus","Dimidiochromis_strigatus","Dimidiochromis_strigatus_2","Mylochromis_ericotaenia","Otopharynx_lithobates","Placidochromis_milomo","Hemitilapia_oxyrhynchus","Hemitilapia_oxyrhynchus_2","Chilotilapia_rhoadesii_1","Chilotilapia_rhoadesii_2","Chilotilapia_rhoadesii_3","Chilotilapia_rhoadesii_4","Protomelas_ornatus_1","Protomelas_ornatus_2","Hemitaenichromis_spilopterus","Hemitaeniochromis_spilopterus_2","Placidochromis_johnstoni")) # Change the lake species order according to phylogeny
# order the genotypes and wavelength
Mz_stat1_rho_GT6$Genotype <- factor(Mz_stat1_rho_GT6$Genotype, levels=c("T|T","T|A","A/A"))
Mz_stat1_rho_GT6$Wavelength_palette <- factor(Mz_stat1_rho_GT6$Wavelength_palette, levels=c("Short","Medium","Long","N/A"))

# plot the above - legend:
# 1. phylo tree where nodes tips are:
# A. Colour = Foraging/Diet and Shape = Genotype
# B. Shape = Wavelength palette
# C. Marker = Habitat: **** = Demersal; *** = benthopelagic; ** = pelagic, * = rock
p3 <- p2 %<+% Mz_stat1_rho_GT6 + 
  geom_text(aes(label=Habitat_marker), show.legend = FALSE, position = position_nudge(y = -0.05 ,x = 2)) +
  # geom_point(aes(shape=Wavelength_palette), size=3, show.legend = TRUE, position = position_nudge(x = 1.5)) +
  geom_point(aes(color=Foraging.Diet, shape=Genotype), size=3, show.legend = TRUE, position = position_nudge(x = 0.5)) +
  scale_color_manual(breaks=c('Algae','Aufwuchs','Benthivore','Fish','Herbivore','Invertebrates','Mix','Waste','Zooplankton','N/A'),values = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', 'lavenderblush2', '#999999', 'darkolivegreen2')) +
  # scale_color_manual(values = c("Mix" = '#E69F00','Invertebrates' = '#56B4E9','Fish' = '#009E73','Zooplankton' = '#F0E442','Algae' = '#0072B2','Aufwuchs' = '#D55E00','Herbivore' = '#CC79A7','Benthivore' = 'green','Waste' = '#999999','N/A' = 'lavenderblush2')) +
  scale_shape_manual(values = c("T|T" = 16,'T|A' = 17,'A/A' = 15, 'Short' = 12, 'Medium' = 13, 'Long' = 14,'N/A' = 4)) +
  annotate("text", x=26.5, y=119, label= "A", fontface="bold") +
  annotate("text", x=28, y=119, label= "B", fontface="bold") +
  # annotate("text", x=29, y=119, label= "C", fontface="bold") +
  theme(legend.position=c(.9, .6))
plot(p3)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/SNPs/FigS-R4fC_stat1_rho_LakeMalawi-SNPs.tiff', units="in", width=14, height=16, res=300)
p3
dev.off()

## C. gata2a>rho (in On and Ab, lost in Mz)
# For this, you will need to add A. burtoni and O. niloticus to the tree and species info
# On and Ab genotype is C/C
# Mz genotype is T|T

# GATA2a	OG5052_0	rho	OG1655_0	0	0	1	0	1 - present in Ab and On (including their rho networks) < this is the one we target as there is a SNP that may have broken the site in Mz (only gata2, not gata2a targets Mz in network)
# GATA2	OG5053_0	rho	OG1655_0	1	0	1	0	1 - present in Mz, Ab and On
# scaffold_12     3797654 3797655 on11_mz11       100     -       G       A       LG20    13233110        LG20    13233110        13233132        on.gene.LG20.398        rho     rho     rho     +       motif_29        on.gene.LG5.242 gata2   2.80412 1.97e-05  0.0175  GACTGATTGTAAGTGATTGTGCT 2a      0.125   scaffold_12     3797655 .       T|T     C|C

# https://4va.github.io/biodatasci/r-ggtree.html # guide here
# newick tree manually generated
library(tidyverse)
library(ggtree)
LM_tree3 <- read.tree("~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_MalawiVCF/Lake_Malawi_species_phylo3.nwk")

# label clades according to node number
p4 <- ggtree(LM_tree3) + 
  # geom_text(aes(label=node), hjust=-.3) + # this will let you see the node numbers
  geom_tiplab(offset = 3.5) + 
  geom_cladelabel(node=183, label="Shallow benthic",
                  color="red", offset=19, align=TRUE) +
  geom_cladelabel(node=55, label="Utaka",
                  color="darkgreen", offset=19, align=TRUE) +
  geom_cladelabel(node=170, label="Deep benthic", 
                  color="blue", offset=19, align=TRUE) + 
  geom_cladelabel(node=164, label="Utaka", 
                  color="darkgreen", offset=19, align=TRUE) +
  geom_cladelabel(node=143, label="A. calliptera", 
                  color="lightgreen", offset=19, align=TRUE) +
  geom_cladelabel(node=135, label="Mbuna", 
                  color="purple", offset=19, align=TRUE) +
  geom_cladelabel(node=127, label="Diplotaxodon", 
                  color="orange", offset=19, align=TRUE) +
  geom_cladelabel(node=125, label="Rhamphochromis", 
                  color="brown", offset=19, align=TRUE) +
  geom_cladelabel(node=3, label="Lakes and Rivers", 
                  color="grey", offset=19, align=TRUE) +
  geom_cladelabel(node=2, label="Lake Tanganyika", 
                  color="pink", offset=19, align=TRUE) +
  geom_cladelabel(node=1, label="Riverine", 
                  color="gold", offset=19, align=TRUE) +
  theme_tree2() + 
  xlim(0, 70) + 
  theme_tree()
plot(p4)

# load in Lake Malawi species details
LM_Astral3 <- read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/Lake_MalawiVCF/Lake_Malawi_species_info3.txt", header = TRUE,stringsAsFactors=FALSE)

# select Mz gata2a>rho and genotype columns only - other cols could be 15-TG; 21-TF; 31-RefGT; 32-AltGT
Mz_gata2a_rho_GT <- as.data.frame(Mz_RefSp[grep("gata2", Mz_RefSp$TFmotif_gene_symbol),c(1:306)])
Mz_gata2a_rho_GTa <- as.data.frame(t(Mz_gata2a_rho_GT[grep("rho", Mz_gata2a_rho_GT$motif_genesymbolSp),c(172:306)]),stringsAsFactors=FALSE)
# add the A. burtoni and O. niloticus genotype
Ab_gata2a_rhoGT<-data.frame(row.names = "Aburtoni_GT","C/C")
names(Ab_gata2a_rhoGT) <- ("4539")
On_gata2a_rhoGT<-data.frame(row.names = "Oniloticus_GT","C/C")
names(On_gata2a_rhoGT) <- ("4539")
Mz_gata2a_rho_GTb <- rbind(Mz_gata2a_rho_GTa,Ab_gata2a_rhoGT,On_gata2a_rhoGT)
Mz_gata2a_rho_GT2 <- add_rownames(Mz_gata2a_rho_GTb, "Species")
Mz_gata2a_rho_GT2$Species <- gsub("_GT", "", Mz_gata2a_rho_GT2$Species)
names(Mz_gata2a_rho_GT2) <- c("Species", "Genotype")
Mz_gata2a_rho_GT3 <- cbind(Mz_gata2a_rho_GT2, Target_Gene = 'rho', TF = 'gata2a')
# remove the following species that are not in the tree:
# Labeotropheus_trewavasae
# A_calliptera_Lake_Chilwa
# A_calliptera_Rovuma_river_Female
# A_calliptera_Lake_Chidya
# A_calliptera_Kitai_Dam
# A_calliptera_Upper_Rovuma
# Copadichromis_mloto
# Tremitochranus_placodon
# Fossorochromis_rostratus_3
# Fossorochromis_rostratus_4
# Placidochromis_subocularis_1
# Placidochromis_subocularis_2
# Placidochromis_subocularis_3
# Placidochromis_subocularis_4
# Placidochromis_subocularis_5
# Placidochromis_subocularis_6
# Placidochromis_subocularis_7
# Placidochromis_subocularis_8
removesp <- c("Labeotropheus_trewavasae","A_calliptera_Lake_Chilwa","A_calliptera_Rovuma_river_Female","A_calliptera_Lake_Chidya","A_calliptera_Kitai_Dam","A_calliptera_Upper_Rovuma","Copadichromis_mloto","Tremitochranus_placodon","Fossorochromis_rostratus_3","Fossorochromis_rostratus_4","Placidochromis_subocularis_1","Placidochromis_subocularis_2","Placidochromis_subocularis_3","Placidochromis_subocularis_4","Placidochromis_subocularis_5","Placidochromis_subocularis_6","Placidochromis_subocularis_7","Placidochromis_subocularis_8")
Mz_gata2a_rho_GT4 <- Mz_gata2a_rho_GT3[!grepl(paste(removesp, collapse="|"), Mz_gata2a_rho_GT3$Species),] # this removes rows of the species above
# add the foraging and ecology of each species:
Mz_gata2a_rho_GT5 <- (LM_Astral3[match(Mz_gata2a_rho_GT4[,1],LM_Astral3[,1]),(c(2,3,4,5))]) # this matches the Species and outputs ecology
Mz_gata2a_rho_GT6 <- cbind(Mz_gata2a_rho_GT4,Mz_gata2a_rho_GT5)
# Create column of star markers for habitat
Mz_gata2a_rho_GT6$Habitat_marker[Mz_gata2a_rho_GT6$Habitat == "Demersal"]<-"****"
Mz_gata2a_rho_GT6$Habitat_marker[Mz_gata2a_rho_GT6$Habitat == "Benthopelagic"]<-"***"
Mz_gata2a_rho_GT6$Habitat_marker[Mz_gata2a_rho_GT6$Habitat == "Pelagic"]<-"**"
Mz_gata2a_rho_GT6$Habitat_marker[Mz_gata2a_rho_GT6$Habitat == "Rock"]<-"*"
# change all NA to N/A:
Mz_gata2a_rho_GT6[is.na(Mz_gata2a_rho_GT6)] <- "N/A"
## the M. zebra genotype is not present in VCF and hence, amend to browser genotype of T|C
Mz_gata2a_rho_GT6$Genotype[74] <- "T|C"
## update other missing genotypes based on browser


# order the factor levels according to the phylogeny order:
Mz_gata2a_rho_GT6$Species <- factor(Mz_gata2a_rho_GT6$Species, levels=c("Oniloticus","Nbrichardi","Aburtoni","Rhamphochromis_longiceps","Rhamphochromis_esox","Rhamphochromis_woodi","Diplotaxodon_limnothrissa","Diplotaxodon_macrops","Diplotaxodon_macrops_black_dorsal","Pallidochromis_tokolosh","Diplotaxodon_greenwoodi","Diplotaxodon_white_back_similis","Diplotaxodon_ngulube","Genyochromis_mento","Cynotilapia_axelrodi","Cynotilapia_afra","Iodotropheus_sprengerae","Tropheops_tropheops","Petrotilapia_genalutea","Metriaclima_zebra","A_calliptera_Chizumulu","A_calliptera_Luwawa","A_calliptera_Bua","A_calliptera_Salima_Father","A_calliptera_Salima_Mother","A_calliptera_Salima_Offspring","A_calliptera_Chitimba","A_calliptera_Enukweni","A_calliptera_South_Rukuru","A_calliptera_North_Rukuru","A_calliptera_Near_Kyela","A_calliptera_Songwe_River","Astatotilapia_Kingiri","A_calliptera_Mbaka_river_Female","Massoko_benthic_HC_1","Massoko_littoral_HC_1","Copadichromis_quadrimaculatus","Copadichromis_virginalis","Copadichromis_virginalis_2","Copadichromis_virginalis_3","Copadichromis_virginalis_4","Copadichromis_virginalis_5","Copadichromis_likomae","Aulonocara_steveni","Aulonocara_stuartgranti_Father","Aulonocara_stuartgranti_Mother","Aulonocara_stuartgranti_Offspring","Alticorpus_geoffreyi","Lethrinops_longimanus_redhead","Lethrinops_gossei","Lethrinops_sp_oliveri","Aulonocara_yellow","Aulonocara_minutus","Alticorpus_macrocleithrum","Copadichromis_trimaculatus","Taeniolethrinops_praeorbitalis","Taeniolethrinops_macrorhynchus","Taeniolethrinops_furcicauda","Lethrinops_lethrinus","Lethrinops_lethrinus_Mazinzi_Father","Lethrinops_lethrinus_Mazinzi_Mother","Lethrinops_lethrinus_Mazinzi_Offspring","Lethrinops_auritus","Lethrinops_albus","Placidochromis_electra","Mylochromis_anaphyrmus","Mylochromis_anaphyrmus_2","Mylochromis_anaphyrmus_3","Mylochromis_anaphyrmus_4","Mylochromis_anaphyrmus_5","Placidochromis_longimanus_1","Placidochromis_longimanus_2","Placidochromis_longimanus_3","Placidochromis_longimanus_4","Placidochromis_longimanus_5","Ctenopharynx_intermedius_1","Ctenopharynx_intermedius_2","Ctenopharynx_intermedius_3","Ctenopharynx_nitidus","Ctenopharynx_nitidus_2","Otopharynx_tetrastigma","Otopharynx_speciosus_1","Otopharynx_speciosus_2","Buccochromis_rhoadesii","Buccochromis_nototaenia","Taeniochromis_holotaenia","Stigmatochromis_modestus","Stigmatochromis_guttatus","Otopharynx_brooksi_nkhata","Trematocranus_placodon_2","Trematocranus_placodon_3","Trematocranus_placodon_4","Trematocranus_placodon_5","Mylochromis_melanotaenia","Tyrannochromis_nigriventer","Champsochromis_caeruelus_1","Champsochromis_caeruelus_2","Nimbochromis_linni","Nimbochromis_polystigma","Nimbochromis_livingstoni","Dimidiochromis_kiwinge","Dimidiochromis_compressiceps","Dimidiochromis_dimidiatus","Dimidiochromis_strigatus","Dimidiochromis_strigatus_2","Mylochromis_ericotaenia","Otopharynx_lithobates","Placidochromis_milomo","Hemitilapia_oxyrhynchus","Hemitilapia_oxyrhynchus_2","Chilotilapia_rhoadesii_1","Chilotilapia_rhoadesii_2","Chilotilapia_rhoadesii_3","Chilotilapia_rhoadesii_4","Protomelas_ornatus_1","Protomelas_ornatus_2","Hemitaenichromis_spilopterus","Hemitaeniochromis_spilopterus_2","Placidochromis_johnstoni")) # Change the lake species order according to phylogeny
# order the genotypes and wavelength
Mz_gata2a_rho_GT6$Genotype <- factor(Mz_gata2a_rho_GT6$Genotype, levels=c("C/C","C|C","C|T","T|C","T|T","N/A"))
Mz_gata2a_rho_GT6$Wavelength_palette <- factor(Mz_gata2a_rho_GT6$Wavelength_palette, levels=c("Short","Medium","Long","N/A"))


# plot the above - legend:
# 1. phylo tree where nodes tips are:
# A. Colour = Foraging/Diet and Shape = Genotype
# B. Shape = Wavelength palette
# C. Marker = Habitat: **** = Demersal; *** = benthopelagic; ** = pelagic, * = rock
p5 <- p4 %<+% Mz_gata2a_rho_GT6 + 
  geom_text(aes(label=Habitat_marker), show.legend = FALSE, position = position_nudge(y = -0.05 ,x = 2)) +
  # geom_point(aes(shape=Wavelength_palette), size=3, show.legend = TRUE, position = position_nudge(x = 1.5)) +
  geom_point(aes(color=Foraging.Diet, shape=Genotype), size=3, show.legend = TRUE, position = position_nudge(x = 0.5)) +
  scale_color_manual(breaks=c('Algae','Aufwuchs','Benthivore','Fish','Herbivore','Invertebrates','Mix','Waste','Zooplankton','N/A'),values = c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', 'lavenderblush2', '#999999', 'darkolivegreen2')) +
  # scale_color_manual(values = c("Mix" = '#E69F00','Invertebrates' = '#56B4E9','Fish' = '#009E73','Zooplankton' = '#F0E442','Algae' = '#0072B2','Aufwuchs' = '#D55E00','Herbivore' = '#CC79A7','Benthivore' = 'green','Waste' = '#999999','N/A' = 'lavenderblush2')) +
  scale_shape_manual(values = c('C/C' = 0,'C|C' = 15,'C|T'=20,'T|C' = 17, 'T|T' = 18,'N/A' = 4)) +
  annotate("text", x=27.5, y=121, label= "A", fontface="bold") +
  annotate("text", x=29, y=121, label= "B", fontface="bold") +
  # annotate("text", x=29, y=119, label= "C", fontface="bold") +
  theme(legend.position=c(.9, .6))
plot(p5)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/SNPs/FigS-R4fB_gata2a_rho_LakeMalawi-SNPs.tiff', units="in", width=14, height=16, res=300)
p5
dev.off()

######################################################################################################################
# ## DO NOT USE WIHOUT AMENDING
# # This way you create a separate heatmap plot for the above and put side by side but it is very offset!!!
# # legend:
# # x-axis: foraging/diet
# # y-axis: species
# # fill: genotype
# # label1: visual opsins wavelength palette
# # label2: **** = Demersal; *** = benthopelagic; ** = pelagic, * = rock
# 
# 
# p1 <- ggplot(data = Mz_nr2c2_sws1_GT6, mapping = aes(x = Genotype,y = Species,fill = Foraging.Diet)) +
#   geom_tile() +
#   xlab(label = "Genotype") +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank()) +
#   geom_text(aes(label=Wavelength_palette)) +
#   geom_text(aes(label=Habitat_marker), position = position_nudge(y = 0.2))
# 
# # use multiplot to place the phylo tree and the 'heatmap' side by side. 
# multiplot(p,p1, ncol=2)
#
# # dummy example/test of plotting the genotypes
# XX1 <- c("Nbrichardi","Rhamphochromis_longiceps","Rhamphochromis_esox","Rhamphochromis_woodi")
# XX2 <- c("A|A","T|T","T|T","T|T")
# XX3 <- c("sws1","sws1","sws1","sws1")
# XX4 <- c("Nr2c2","Nr2c2","Nr2c2","Nr2c2")
# XX5 <- c("Lake_Tanganyika","Rhamphochromis","Rhamphochromis","Rhamphochromis")
# XX6 <- c("Short","Medium","NA","Short")
# XX7 <- c("Benthivore","Fish","Zooplankton","Algae")
# XX8 <- c("Benthopelagic","Pelagic","Rock","Rock")
# XX9 <- data.frame(XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8)
# 
# # Create column of star markers for habitat
# XX9$XX10[XX9$XX8 == "Benthopelagic"]<-"***"
# XX9$XX10[XX9$XX8 == "Pelagic"]<-"**"
# XX9$XX10[XX9$XX8 == "Rock"]<-"*"
# 
# # genotype according to foraging/diet, then, label boxes according to wavelength palette and symbols for habitat
# ggplot(data = XX9, mapping = aes(x = XX7,y = XX1,fill = XX2)) +
#   geom_tile() +
#   xlab(label = "Foraging/Diet") +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank()) +
#   geom_text(aes(label=XX6), size = 6) +
#   geom_text(aes(label=XX10), position = position_nudge(y = 0.2), size = 15)


######################################################################################################################
