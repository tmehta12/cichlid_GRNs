# Prepare plots of the topmost three node motifs in each species to observe any correlation

library(ggplot2)
library(hexbin)
library(reshape2)
library(gridExtra)
library(grid)
library(magrittr)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/")

# all edges
threenode_all = read.table("Edge_Attributes_Collated4d_miRNA-TG-TFmotifs.TF-miRNAoverlap.txt2a", header = FALSE)
threenode_all2 <- melt(threenode_all)
# threenode_all3 <- threenode_all2[1:250,] # this selects the top 250
# select the top 100 from each species and plot that
threenode_all2Pn <- threenode_all2[threenode_all2$V3=="Pn",][1:100,]
threenode_all2Mz <- threenode_all2[threenode_all2$V3=="Mz",][1:100,]
threenode_all2Ab <- threenode_all2[threenode_all2$V3=="Ab",][1:100,]
threenode_all2Nb <- threenode_all2[threenode_all2$V3=="Nb",][1:100,]
threenode_all2On <- threenode_all2[threenode_all2$V3=="On",][1:100,]
threenode_all3 <- rbind(threenode_all2Mz,threenode_all2Pn,threenode_all2Ab,threenode_all2Nb,threenode_all2On)
threenode_all3$V3 <- factor(threenode_all3$V3, levels=c("Mz","Pn","Ab","Nb","On")) # change species order according to phylogeny - put them in the reverse order

# heatmap of counts of TF-miRNA interactions in three-node motifs of 5-cichlid species 
threenode_all_heatmap <- ggplot(data = threenode_all3, aes(x=V3, y=V2, fill=value)) + 
  geom_tile(color = "black") + 
  labs(y="TF-miRNA interactions", x="Species",title="Counts of top 150 TF-miRNA overlaps in\n three-node motifs of 5-cichlid species") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 3, limit = c(0,120), space = "Lab", name="Three node motif\ncounts") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=90,vjust=0.5,size=12,hjust=1.5),axis.text.y=element_text(size=10),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.5, vjust=20.12, size=14,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.1.b.b–TF-miRNA_threenodemotifcounts_all_edges_allfivespecies.tiff', units="in", width=5, height=50, res=300)
threenode_all_heatmap
dev.off()

# 1-to-1 node edges
threenode_1to1 = read.table("Edge_Attributes_Collated4d_miRNA-TG-TFmotifs-1to1.TF-miRNAoverlap.txt2a", header = FALSE)
threenode_1to12 <- melt(threenode_1to1)
# threenode_1to13 <- threenode_1to12[1:250,] # this selects the top 250
# select the top 100 from each species and plot that
threenode_1to12Pn <- threenode_1to12[threenode_1to12$V3=="Pn",][1:100,]
threenode_1to12Mz <- threenode_1to12[threenode_1to12$V3=="Mz",][1:100,]
threenode_1to12Ab <- threenode_1to12[threenode_1to12$V3=="Ab",][1:100,]
threenode_1to12Nb <- threenode_1to12[threenode_1to12$V3=="Nb",][1:100,]
threenode_1to12On <- threenode_1to12[threenode_1to12$V3=="On",][1:100,]
threenode_1to13 <- rbind(threenode_1to12Mz,threenode_1to12Pn,threenode_1to12Ab,threenode_1to12Nb,threenode_1to12On)
threenode_1to13$V3 <- factor(threenode_1to13$V3, levels=c("Mz","Pn","Ab","Nb","On")) # change species order according to phylogeny - put them in the reverse order

# heatmap of counts of TF-miRNA interactions in three-node motifs of 5-cichlid species 
threenode_1to1_heatmap <- ggplot(data = threenode_1to13, aes(x=V3, y=V2, fill=value)) + 
  geom_tile(color = "black") + 
  labs(y="TF-miRNA interactions", x="Species",title="Counts of top 150 TF-miRNA overlaps in\n three-node motifs of 5-cichlid species") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 3, limit = c(0,120), space = "Lab", name="Three node motif\ncounts") +
  theme_minimal() + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle=90,vjust=0.5,size=12,hjust=1.5),axis.text.y=element_text(size=10),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14)) +  coord_fixed(ratio=1) + 
  theme(plot.title = element_text(hjust = 0.5, vjust=20.12, size=14,face = "bold")) +
  theme(panel.background=element_rect(fill="white", colour="black"))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3.1.b.c–TF-miRNA_threenodemotifcounts_1to1_edges_allfivespecies.tiff', units="in", width=5, height=50, res=300)
threenode_1to1_heatmap
dev.off()