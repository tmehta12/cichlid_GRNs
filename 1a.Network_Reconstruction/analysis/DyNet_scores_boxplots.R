### Preparing boxplot of derived DyNet scores from TF-TG edges comparison in all five species only:
# a. edges that include only 1to1 orthologs of TFs and TGs
# b. edges that include only 1to1 orthologs of TFs in edges where the TF-TG is present and TG confirmed as absent from genome

library(ggplot2)
library(hexbin)
library(reshape2)
library(gridExtra)
library(grid)
library(magrittr)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/")

# aa. edges that include only 1to1 orthologs of TFs and TGs
A_1to15sp = read.table("/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Extended_Data_Table_S-R3A_TFTG1to1_DyNet_rewiring_allfive_species.tsv", header=TRUE)
colnames(A_1to15sp) <- c('name','SUID','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','Pn_present','Ab_present','Nb_present','On_present','Gene_Symbol')
A_1to15sp <- cbind(A_1to15sp, Source = '1-to-1 TF-TG edges - 5 species')
A_1to15sp_scores <- melt(A_1to15sp,id.vars=c('Source','Gene_Symbol'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# ab. candidate gene edges that include only 1to1 orthologs of TFs and TGs
B_1to15sp = read.table("/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Extended_Data_Table_S-R3B_TFTG1to1_DyNet_rewiring_allfive_species_candidates.tsv", header=TRUE)
colnames(B_1to15sp) <- c('name','SUID','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','Pn_present','Ab_present','Nb_present','On_present','Gene_Symbol')
B_1to15sp <- cbind(B_1to15sp, Source = 'Candidates 1-to-1 TF-TG edges - 5 species')
B_1to15sp_scores <- melt(B_1to15sp,id.vars=c('Source','Gene_Symbol'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

AB_1to15sp_scores <- rbind(A_1to15sp_scores,B_1to15sp_scores)

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/TFTGconsensusedges/A.Boxplot_of_Dn_and_rewiringscore_1to1TFTG_5sp.tiff', units="in", width=13, height=6, res=300)
ggplot(AB_1to15sp_scores, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_violin() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Source) + guides(fill=guide_legend(title="Source of score")) + 
  labs (x = "DyNet score") + labs (y = "log2 (score value)") + theme(legend.position="none") + ylim(-7, 15) +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black"))
dev.off()
# + labs(title="Boxplot of Dn-score (degree corrected) and DyNet Rewiring (Dn-score) from 1-to-1 TF-TG edges and canndidate genes") 

# ba. edges that include only 1to1 orthologs of TFs in edges where the TF-TG is present and TG confirmed as absent from genome
A_all5sp = read.table("/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Extended_Data_Table_S-R3C_TFTGalledges_DyNet_rewiring_allfive_species.tsv", header=TRUE)
colnames(A_all5sp) <- c('name','SUID','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','Pn_present','Ab_present','Nb_present','On_present','Gene_Symbol')
A_all5sp <- cbind(A_all5sp, Source = 'All TF-TG edges - 5 species')
A_all5sp_scores <- melt(A_all5sp,id.vars=c('Source','Gene_Symbol'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# bb. candidate gene # b. edges that include only all orthologs of TFs in edges where the TF-TG is present and TG confirmed as absent from genome
B_all5sp = read.table("/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Extended_Data_Table_S-R3D_TFTGalledges_DyNet_rewiring_allfive_species_candidates.tsv", header=TRUE)
colnames(B_all5sp) <- c('name','SUID','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','Pn_present','Ab_present','Nb_present','On_present','Gene_Symbol')
B_all5sp <- cbind(B_all5sp, Source = 'Candidates All TF-TG edges - 5 species')
B_all5sp_scores <- melt(B_all5sp,id.vars=c('Source','Gene_Symbol'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

AB_all5sp_scores <- rbind(A_all5sp_scores,B_all5sp_scores)

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/TFTGconsensusedges/B.Boxplot_of_Dn_and_rewiringscore_allTFTG_5sp.tiff', units="in", width=13, height=6, res=300)
ggplot(AB_all5sp_scores, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_violin() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Source) + guides(fill=guide_legend(title="Source of score")) + 
  labs (x = "DyNet score") + labs (y = "log2 (score value)") + theme(legend.position="none") + ylim(-7, 15) +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black"))
dev.off()
# + labs(title="Boxplot of Dn-score (degree corrected) and DyNet Rewiring (Dn-score) from 1-to-1 TF-TG edges and canndidate genes") 


### Preparing boxplot of derived DyNet scores from TF-TG edges using different inputs:
# a. edges that include only 1to1 orthologs of TFs and TGs
# b. edges that include 1to1 orthologs of TFs only, all TGs (no constraints)
# c. all edges (no 1to1 orthology constraints)

# above for the following sets:
# i. derived from all 5 species
# ii. derived from haplochromines (Mz, Pn, Ab)
# iii. derived from pairwise comparisons



# i. derived from all 5 species
# a. edges that include only 1to1 orthologs of TFs and TGs
ia_1to15sp = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/1.1to1_files_shared/Result_files/All_five_species_node_noExpression.csv")
colnames(ia_1to15sp) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ia_1to15sp <- cbind(ia_1to15sp, Source = 'all1to1_5species')
ia_1to15sp_scores <- melt(ia_1to15sp,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# b. edges that include 1to1 orthologs of TFs only, all TGs (no constraints)
ib_1to1TF5sp = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/4.TF1to1_TGallnodes/drive-download-20180625T092257Z-001/All_five_species_TF1to1.csv")
colnames(ib_1to1TF5sp) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ib_1to1TF5sp <- cbind(ib_1to1TF5sp, Source = 'all1to1TF_5species')
ib_1to1TF5sp_scores <- melt(ib_1to1TF5sp,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# c. all edges (no 1to1 orthology constraints)
ic_no1to15sp = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/3.alledges_no1to1/drive-download-20180625T091957Z-001/All_five_species_no1to1.csv")
colnames(ic_no1to15sp) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ic_no1to15sp <- cbind(ic_no1to15sp, Source = 'no1to1_5species')
ic_no1to15sp_scores <- melt(ic_no1to15sp,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

iabc_5sp_scores <- rbind(ia_1to15sp_scores,ib_1to1TF5sp_scores,ic_no1to15sp_scores)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/A.Boxplot_of_Dn_and_rewiringscore_of_three_sources_5sp_comp.tiff', units="in", width=13, height=6, res=300)
ggplot(iabc_5sp_scores, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Source) + labs(title="Boxplot of Dn-score (degree corrected) and DyNet Rewiring (Dn-score) from three different sources in 5-species comparisons") + guides(fill=guide_legend(title="Source of score")) + 
  labs (x = "DyNet score") + labs (y = "log2 (score value)") + theme(legend.position="none") + ylim(-7, 15)
dev.off()

# same as above, but DyNet scores derived from module constrained edges only

# a. edges that include only 1to1 orthologs of TFs and TGs - Module Constrained
ia_1to15spMC = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/5.ModuleConstrained/all1to1/allfiveSpecies_all1to1.csv")
colnames(ia_1to15spMC) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ia_1to15spMC <- cbind(ia_1to15spMC, Source = 'all1to1_5species')
ia_1to15spMC_scores <- melt(ia_1to15spMC,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# b. edges that include 1to1 orthologs of TFs only, all TGs (no constraints) - Module Constrained
ib_1to1TF5spMC = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/5.ModuleConstrained/TF1to1_allTG/allfiveSpecies_TF1to1_allTG.csv")
colnames(ib_1to1TF5spMC) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ib_1to1TF5spMC <- cbind(ib_1to1TF5spMC, Source = 'all1to1TF_5species')
ib_1to1TF5spMC_scores <- melt(ib_1to1TF5spMC,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# c. all edges (no 1to1 orthology constraints) - Module Constrained
ic_no1to15spMC = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/5.ModuleConstrained/allTFTG/allSpecies_allTFTG_moduleconstrained.csv")
colnames(ic_no1to15spMC) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ic_no1to15spMC <- cbind(ic_no1to15spMC, Source = 'no1to1_5species')
ic_no1to15spMC_scores <- melt(ic_no1to15spMC,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

iabc_5spMC_scores <- rbind(ia_1to15spMC_scores,ib_1to1TF5spMC_scores,ic_no1to15spMC_scores)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/A.Boxplot_of_Dn_and_rewiringscore_of_three_sources_5spModuleConstrained_comp.tiff', units="in", width=13, height=6, res=300)
ggplot(iabc_5spMC_scores, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Source) + labs(title="Boxplot of Dn-score (degree corrected) and DyNet Rewiring (Dn-score) from three different sources (Module constrained TF/TG) in 5-species comparisons") + guides(fill=guide_legend(title="Source of score")) + 
  labs (x = "DyNet score") + labs (y = "log2 (score value)") + theme(legend.position="none")+ ylim(-7, 15)
dev.off()


# ii. derived from haplochromines (Mz, Pn, Ab)
# a. edges that include only 1to1 orthologs of TFs and TGs
ia_1to1Haplo = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/1.1to1_files_shared/Result_files/Haplochromines_Mz_Pn_Ab_node.csv")
colnames(ia_1to1Haplo) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ia_1to1Haplo <- cbind(ia_1to1Haplo, Source = 'all1to1_Haplo')
ia_1to1Haplo_scores <- melt(ia_1to1Haplo,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# b. edges that include 1to1 orthologs of TFs only, all TGs (no constraints)
ib_1to1TFHaplo = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/4.TF1to1_TGallnodes/drive-download-20180625T092257Z-001/Haplochromines_TF1to1.csv")
colnames(ib_1to1TFHaplo) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ib_1to1TFHaplo <- cbind(ib_1to1TFHaplo, Source = 'all1to1TF_Haplo')
ib_1to1TFHaplo_scores <- melt(ib_1to1TFHaplo,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# c. all edges (no 1to1 orthology constraints)
ic_no1to1Haplo = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/3.alledges_no1to1/drive-download-20180625T091957Z-001/Haplochromines_Mz_Pn_Ab_no1to1.csv")
colnames(ic_no1to1Haplo) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ic_no1to1Haplo <- cbind(ic_no1to1Haplo, Source = 'no1to1_Haplo')
ic_no1to1Haplo_scores <- melt(ic_no1to1Haplo,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

iabc_Haplo_scores <- rbind(ia_1to1Haplo_scores,ib_1to1TFHaplo_scores,ic_no1to1Haplo_scores)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/B.Boxplot_of_Dn_and_rewiringscore_of_three_sources_Haplo_comp.tiff', units="in", width=13, height=6, res=300)
ggplot(iabc_Haplo_scores, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Source) + labs(title="Boxplot of Dn-score (degree corrected) and DyNet Rewiring (Dn-score) from three different sources in haplochromines comparisons") + guides(fill=guide_legend(title="Source of score")) + 
  labs (x = "DyNet score") + labs (y = "log2 (score value)") + theme(legend.position="none") + ylim(-7, 15)
dev.off()

# same as above, but DyNet scores derived from module constrained edges only

# a. edges that include only 1to1 orthologs of TFs and TGs - Module Constrained
ia_1to1HaploMC = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/5.ModuleConstrained/all1to1/Haplochromines_Mz_Pn_Ab_all1to1.csv")
colnames(ia_1to1HaploMC) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ia_1to1HaploMC <- cbind(ia_1to1HaploMC, Source = 'all1to1_Haplo')
ia_1to1HaploMC_scores <- melt(ia_1to1HaploMC,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# b. edges that include 1to1 orthologs of TFs only, all TGs (no constraints) - Module Constrained
ib_1to1TFHaploMC = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/5.ModuleConstrained/TF1to1_allTG/Haplochromines_Mz_Pn_Ab_TF1to1_allTG.csv")
colnames(ib_1to1TFHaploMC) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ib_1to1TFHaploMC <- cbind(ib_1to1TFHaploMC, Source = 'all1to1TF_Haplo')
ib_1to1TFHaploMC_scores <- melt(ib_1to1TFHaploMC,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

# c. all edges (no 1to1 orthology constraints) - Module Constrained
ic_no1to1HaploMC = read.csv("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/5.ModuleConstrained/allTFTG/Haplochromines_Ab_Mz_Pn_allTFTG_moduleconstrained.csv")
colnames(ic_no1to1HaploMC) <- c('SUID','Ab_present','Dn-score_degree-corrected','DyNet_Rewiring_Dn-score','EdgeCount','Mz_present','name','Nb_present','On_present','Pn_present','selected','shared_name')
ic_no1to1HaploMC <- cbind(ic_no1to1HaploMC, Source = 'no1to1_Haplo')
ic_no1to1HaploMC_scores <- melt(ic_no1to1HaploMC,id.vars=c('Source','name'), measure.vars=c('DyNet_Rewiring_Dn-score','Dn-score_degree-corrected'))

iabc_HaploMC_scores <- rbind(ia_1to1HaploMC_scores,ib_1to1TFHaploMC_scores,ic_no1to1HaploMC_scores)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/rewiring_files_DyNet_Marton/B.Boxplot_of_Dn_and_rewiringscore_of_three_sources_HaploModuleConstrained_comp.tiff', units="in", width=13, height=6, res=300)
ggplot(iabc_HaploMC_scores, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Source) + labs(title="Boxplot of Dn-score (degree corrected) and DyNet Rewiring (Dn-score) from three different sources (Module constrained TF/TG) in haplochromines comparisons") + guides(fill=guide_legend(title="Source of score")) + 
  labs (x = "DyNet score") + labs (y = "log2 (score value)") + theme(legend.position="none") + ylim(-7, 15)
dev.off()
