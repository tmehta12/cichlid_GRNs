### Plotting GO enrichment of module switched and rewired sites and edges of rewired sites 

library(ggplot2)
library(hexbin)
library(reshape2)
library(gridExtra)
library(grid)
library(magrittr)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/")

# What are the module switched genes that are also rewired in networks? - only focus on those that are switched and rewired in edges only, as well as the whole edge
### A. GO enrichment of the gene that is switched in the edge
base_dir1 <-("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/GOOUTPUT_BPonly/")
### B. GO enrichment of the edges in each file
base_dir2 <-("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/GOOUTPUT_BPonly/geneAandgeneB/")

GOOUTPUT_BPonly_rewiredsiteonly <- read.delim("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/GOOUTPUT_BPonly/GOOUTPUT_BPonly-rewiredsiteonly_cut_cat_GOOUTPUT_details_filtered.txt3",header=F)
GOOUTPUT_BPonly_allsites <- read.delim("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/1to1only/GOOUTPUT_BPonly/geneAandgeneB/GOOUTPUT_BPonly-allsites_cut_cat_GOOUTPUT_details_filtered.txt3",header=F)

# these are for individual bar plots
GOOUTPUT_BPonly_rewiredsiteonly_plot <- ggplot(`GOOUTPUT_BPonly_rewiredsiteonly`, aes(x = reorder(V2, -V9), y = -log2(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=12,face="bold")) + theme(legend.text=element_text(size=12),legend.title=element_text(size=12,face="bold")) + labs(title="Promoter TFBSs") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

GOOUTPUT_BPonly_allsites_plot <- ggplot(`GOOUTPUT_BPonly_allsites`, aes(x = reorder(V2, -V9), y = -log2(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="Promoter TFBSs") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

# produce heatmaps for the above too

library(reshape2)
library(dplyr)
melted_GOOUTPUT_BPonly_rewiredsiteonly <- melt(GOOUTPUT_BPonly_rewiredsiteonly, id.vars = 1:2, value.var = GOOUTPUT_BPonly_rewiredsiteonly$V9)
melted_GOOUTPUT_BPonly_rewiredsiteonly2 <- dplyr::filter(melted_GOOUTPUT_BPonly_rewiredsiteonly, grepl('V9', variable))

GOOUTPUT_BPonly_rewiredsiteonly_plot <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.0–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1.tiff', units="in", width=14, height=20, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot
dev.off()

# Can also split the GO enrichment into 

# A. all five species
# Ai. no cis/trans - TF (geneASW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveA <- melted_GOOUTPUT_BPonly_rewiredsiteonly2[grep("rewired", melted_GOOUTPUT_BPonly_rewiredsiteonly2$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveB <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveA[!grepl("Cis", melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveA$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveC <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveB[!grepl("Trans", melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveB$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveD <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveC[!grepl("geneBSW", melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveC$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveD$V1 <- gsub('geneAB1to1OGID-networkandmoduleSW_rewired_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveD$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveD$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveD$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveD$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveD$V1, levels=c("On-MzPnAbNb","Nb-MzPnAbOn","Ab-MzPnNbOn","Pn-MzAbNbOn","Mz-PnAbNbOn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASW <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveD, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- 5-species comparison TF module switch") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.1–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-5spTFSW.tiff', units="in", width=10, height=6, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASW
dev.off()

# Aii. cis - TF (geneASW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveE <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveA[grep("Cis", melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveA$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveF <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveE[!grepl("geneBSW", melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveE$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveF$V1 <- gsub('geneAB1to1OGID-networkandmoduleSW_Cis-rewired_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveF$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveF$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveF$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveF$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveF$V1, levels=c("On-MzPnAbNb","Nb-MzPnAbOn","Ab-MzPnNbOn","Pn-MzAbNbOn","Mz-PnAbNbOn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASWcis <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveF, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- 5-species comparison TF module switch (cis TF-TG)") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +
  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.2–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-5spTFSWcis.tiff', units="in", width=10, height=6, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASWcis
dev.off()

# Aiii. trans - TF (geneASW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveG <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveA[grep("Trans", melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveA$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveH <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveG[!grepl("geneBSW", melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveG$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveH$V1 <- gsub('geneAB1to1OGID-networkandmoduleSW_Trans-rewired_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveH$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveH$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveH$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveH$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveH$V1, levels=c("On-MzPnAbNb","Nb-MzPnAbOn","Ab-MzPnNbOn","Pn-MzAbNbOn","Mz-PnAbNbOn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASWtrans <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveH, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- 5-species comparison TF module switch (trans TF-TG)") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.3–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-5spTFSWtrans.tiff', units="in", width=10, height=6, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASWtrans
dev.off()

# Aiv. no cis/trans - TG (geneBSW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveI <- melted_GOOUTPUT_BPonly_rewiredsiteonly2[grep("geneBSW_rewired", melted_GOOUTPUT_BPonly_rewiredsiteonly2$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveJ <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveI[!grepl("Cis", melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveI$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveK <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveJ[!grepl("Trans", melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveJ$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveK$V1 <- gsub('geneAB1to1OGID-networkandmodule_geneBSW_rewired_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveK$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveK$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveK$V1)
#melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveK$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveK$V1, levels=c("On-MzPnAbNb","Nb-MzPnAbOn","Ab-MzPnNbOn","Pn-MzAbNbOn","Mz-PnAbNbOn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASW <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveK, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- 5-species comparison TG module switch") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.4–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-5spTGSW.tiff', units="in", width=7, height=4, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASW
dev.off()

# Av. cis - TG (geneBSW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveL <- melted_GOOUTPUT_BPonly_rewiredsiteonly2[grep("geneBSW_Cis-rewired", melted_GOOUTPUT_BPonly_rewiredsiteonly2$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveL$V1 <- gsub('geneAB1to1OGID-networkandmodule_geneBSW_Cis-rewired_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveL$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveL$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveL$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveL$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveL$V1, levels=c("On-MzPnAbNb","Nb-MzPnAbOn","Ab-MzPnNbOn","Pn-MzAbNbOn","Mz-PnAbNbOn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASWcis <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveL, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- 5-species comparison TG module switch (cis TF-TG)") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +
  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.5–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-5spTGSWcis.tiff', units="in", width=10, height=6, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASWcis
dev.off()

# Avi. trans - TG (geneBSW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveM <- melted_GOOUTPUT_BPonly_rewiredsiteonly2[grep("geneBSW_Trans-rewired", melted_GOOUTPUT_BPonly_rewiredsiteonly2$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveM$V1 <- gsub('geneAB1to1OGID-networkandmodule_geneBSW_Trans-rewired_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveM$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveM$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveM$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveM$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveM$V1, levels=c("On-MzPnAbNb","Nb-MzPnAbOn","Ab-MzPnNbOn","Pn-MzAbNbOn","Mz-PnAbNbOn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASWtrans <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_allfiveM, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- 5-species comparison TG module switch (trans TF-TG)") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.6–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-5spTGSWtrans.tiff', units="in", width=10, height=6, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_5spgeneASWtrans
dev.off()

# B. haplochromines
# TF, TG and cis and trans combined
melted_GOOUTPUT_BPonly_rewiredsiteonly2_haploA <- melted_GOOUTPUT_BPonly_rewiredsiteonly2[grep("haplo", melted_GOOUTPUT_BPonly_rewiredsiteonly2$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_haploA$V1 <- gsub('geneAB1to1OGID-networkandmodule_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_haploA$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_haploA$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_haploA$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_haploA$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_haploA$V1, levels=c("geneBhaploSW_Trans-AbvsMzPn","geneBhaploSW_Trans-PnvsMzAb","geneBhaploSW_Trans-MzvsPnAb","geneAhaploSW_Trans-AbvsMzPn","geneAhaploSW_Trans-PnvsMzAb","geneAhaploSW_Trans-MzvsPnAb","geneBhaploSW_Cis-AbvsMzPn","geneBhaploSW_Cis-PnvsMzAb","geneBhaploSW_Cis-MzvsPnAb","geneAhaploSW_Cis-AbvsMzPn","geneAhaploSW_Cis-PnvsMzAb","geneAhaploSW_Cis-MzvsPnAb","geneBhaploSW_AbvsMzPn","geneBhaploSW_PnvsMzAb","geneBhaploSW_MzvsPnAb","geneAhaploSW_AbvsMzPn","geneAhaploSW_PnvsMzAb","geneAhaploSW_MzvsPnAb"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_haplo <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_haploA, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across haplochromines") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.7–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-haplo.tiff', units="in", width=10, height=6, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_haplo
dev.off()


# C. pairwise comparisons
# Ci. only the vs edges with no cis and trans - TF (geneASW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseA <- melted_GOOUTPUT_BPonly_rewiredsiteonly2[grep("vs", melted_GOOUTPUT_BPonly_rewiredsiteonly2$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseB <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseA[!grepl("Cis", melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseA$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseC <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseB[!grepl("Trans", melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseB$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseD <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseC[!grepl("haplo", melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseC$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseE <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseD[!grepl("geneBSW", melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseD$V1),]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseE$V1 <- gsub('geneAB1to1OGID-networkandmoduleSW_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseE$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseE$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseE$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseE$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseE$V1, levels=c("OnvsNb","OnvsAb","OnvsPn","OnvsMz","NbvsOn","NbvsAb","NbvsPn","NbvsMz","AbvsOn","AbvsNb","AbvsPn","AbvsMz","PnvsOn","PnvsNb","PnvsAb","PnvsMz","MzvsOn","MzvsNb","MzvsAb","MzvsPn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_pairwisegeneASW <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseE, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- pairwise comparison TF module switch") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.8–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-pairwiseTFSW.tiff', units="in", width=8, height=10, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_pairwisegeneASW
dev.off()

# Cii. only cis vs edges

# a. only TF (geneASW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseF <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseA[grep("Cis", melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseA$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseG <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseF[!grepl("geneB", melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseF$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseG$V1 <- gsub('geneAB1to1OGID-networkandmoduleSW_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseG$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseG$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseG$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseG$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseG$V1, levels=c("Cis-OnvsNb","Cis-OnvsAb","Cis-OnvsPn","Cis-OnvsMz","Cis-NbvsOn","Cis-NbvsAb","Cis-NbvsPn","Cis-NbvsMz","Cis-AbvsOn","Cis-AbvsNb","Cis-AbvsPn","Cis-AbvsMz","Cis-PnvsOn","Cis-PnvsNb","Cis-PnvsAb","Cis-PnvsMz","Cis-MzvsOn","Cis-MzvsNb","Cis-MzvsAb","Cis-MzvsPn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_pairwisegeneASWcis <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseG, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- pairwise comparison TF module switch (cis TF-TG)") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.9–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-pairwiseTFSWcis.tiff', units="in", width=8, height=10, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_pairwisegeneASWcis
dev.off()

# b. only TG (geneBSW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseH <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseD[grep("geneBSW", melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseD$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseH$V1 <- gsub('geneAB1to1OGID-networkandmodule_geneBSW_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseH$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseH$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseH$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseH$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseH$V1, levels=c("OnvsNb","OnvsAb","OnvsPn","OnvsMz","NbvsOn","NbvsAb","NbvsPn","NbvsMz","AbvsOn","AbvsNb","AbvsPn","AbvsMz","PnvsOn","PnvsNb","PnvsAb","PnvsMz","MzvsOn","MzvsNb","MzvsAb","MzvsPn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_pairwisegeneBSWcis <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseH, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- pairwise comparison TG module switch") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.10–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-pairwiseTGSW.tiff', units="in", width=10, height=8, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_pairwisegeneBSWcis
dev.off()

# Civ. only trans vs edges
# only TG (geneBSW)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseJ <- melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseA[grep("geneBSW_Trans", melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseA$V1), ]
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseJ$V1 <- gsub('geneAB1to1OGID-networkandmodule_geneBSW_', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseJ$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseJ$V1 <- gsub('.txt', '', melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseJ$V1)
melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseJ$V1 <- factor(melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseJ$V1, levels=c("Trans-OnvsNb","Trans-OnvsAb","Trans-OnvsPn","Trans-OnvsMz","Trans-NbvsOn","Trans-NbvsAb","Trans-NbvsPn","Trans-NbvsMz","Trans-AbvsOn","Trans-AbvsNb","Trans-AbvsPn","Trans-AbvsMz","Trans-PnvsOn","Trans-PnvsNb","Trans-PnvsAb","Trans-PnvsMz","Trans-MzvsOn","Trans-MzvsNb","Trans-MzvsAb","Trans-MzvsPn"))

GOOUTPUT_BPonly_rewiredsiteonly_plot_pairwisegeneBSWtrans <- ggplot(data = melted_GOOUTPUT_BPonly_rewiredsiteonly2_pairwiseJ, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment of module switched and rewired in\nTF-TG collated edges of 1-to-1 orthologs across all five species\n- pairwise comparison TG module switch (trans TF-TG)") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=2.12, size=12,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.2.b.11–GO_enrichment_mod_and_networkSW_rewired_TFTG1to1-pairwiseTGSWtrans.tiff', units="in", width=8, height=8, res=300)
GOOUTPUT_BPonly_rewiredsiteonly_plot_pairwisegeneBSWtrans
dev.off()