### Plotting GO enrichment of the TF/TG (no TG's have significant enrichment) of 1to1 edges and all edges (Using 1to1 TF orthogroups and non 1to1 TGs - filtered for bonafide absence in genome) ) 

library(ggplot2)
library(hexbin)
library(reshape2)
library(gridExtra)
library(grid)
library(magrittr)
library(dplyr)


setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/GOenrichment_BPonly_alledges_and1to1/GOOUTPUT_BPonly/")

# A. load the files
collated_TF_1to1edges_GOOUTPUT <- read.delim("collated_TF_1to1edges_GOOUTPUT_details_filtered.txt3",header=F)
collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT <- read.delim("collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT_details_filtered.txt3",header=F)
# collated_TG_1to1edges_GOOUTPUT <- read.delim("collated_TG_1to1edges_GOOUTPUT_details_filtered.txt3",header=F)
# collated_TG_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT <- read.delim("collated_TG_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT_details_filtered.txt3",header=F)

# B. melt the tables to get the log2 fold enrichment values linked to the species and GO term
melted_collated_TF_1to1edges_GOOUTPUT <- melt(collated_TF_1to1edges_GOOUTPUT, id.vars = 1:2, value.var = collated_TF_1to1edges_GOOUTPUT$V9)
melted_collated_TF_1to1edges_GOOUTPUT2 <- dplyr::filter(melted_collated_TF_1to1edges_GOOUTPUT, grepl('V9', variable))

melted_collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT <- melt(collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT, id.vars = 1:2, value.var = collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT$V9)
melted_collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT2 <- dplyr::filter(melted_collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT, grepl('V9', variable))

# C. create tile heatmaps for the enrichment (will only be for the TFs in 1to1 and alledges)

# Ci. 1to1edges
melted_collated_TF_1to1edges_GOOUTPUT2$V1 <- factor(melted_collated_TF_1to1edges_GOOUTPUT2$V1, levels=c("Mz","Pn","Ab","Nb","On"))

GOOUTPUT_BPonly_collated_TF_1to1edges <- ggplot(data = melted_collated_TF_1to1edges_GOOUTPUT2, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment (FDR < 0.05) of 1-to-1 TF-TG collated edges across all five species") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 3, limit = c(0,6), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=80,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.5, vjust=20.12, size=14,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/TFTGconsensusedges/FigS-R3.1.a–GO_enrichment_of1to1_TFTGedges_allfivespecies.tiff', units="in", width=15, height=8, res=300)
GOOUTPUT_BPonly_collated_TF_1to1edges
dev.off()

# Cii. alledges
melted_collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT2$V1 <- factor(melted_collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT2$V1, levels=c("Mz","Pn","Ab","Nb","On"))

GOOUTPUT_BPonly_collated_TF_alledges <- ggplot(data = melted_collated_TF_TF1to1TGnon1to1PresentNULLOGIDedges_GOOUTPUT2, aes(x=V2, y=V1, fill=log2(value))) + 
  geom_tile(color = "black") + 
  labs(y="Edge comparison in species", x="GO term",title="GO term enrichment (FDR < 0.05) of TF-TG collated edges of 1-to-1 TF and non 1-to1 TG (absent in genome) orthologs across all five species") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 2.5, limit = c(0,5), space = "Lab", name="log2 fold enrichment\n (FDR < 0.05)") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=80,vjust=1,size=12,hjust=1),axis.text.y=element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.8, vjust=20.12, size=14,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/TFTGconsensusedges/FigS-R3.1.b–GO_enrichment_ofalledges_TF1to1TGno1to1AG-allfivespecies.tiff', units="in", width=20, height=9, res=300)
GOOUTPUT_BPonly_collated_TF_alledges
dev.off()
