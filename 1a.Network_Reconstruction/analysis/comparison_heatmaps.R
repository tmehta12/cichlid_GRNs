setwd("/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/matrix/")

### coexpr_prom - geneA (TF) ###
coexpr_prom_geneA_pairwise_percentage_matrix <- as.matrix(read.table("coexpr_prom_geneA_pairwise_percentage_matrix",header = TRUE, row.names = 1, sep = "\t"))

# simple heatmap
heatmap(as.matrix(coexpr_prom_geneA_pairwise_percentage_matrix))

# reshape matrix for ggplot
library(reshape2)
melted_coexpr_prom_geneA_pairwise_percentage_matrix <- melt(coexpr_prom_geneA_pairwise_percentage_matrix)
head(melted_coexpr_prom_geneA_pairwise_percentage_matrix)

# prepare heatmap with ggplot
library(ggplot2)
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R–Pairwise_TF_modSW_rewiredTFTG_1to1.tiff', units="in", width=10, height=10, res=300)
ggplot(data = melted_coexpr_prom_geneA_pairwise_percentage_matrix, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(color = "white") + 
  labs(y="Edge present in species", x="Edge absent in species",title="Percentage pairwise TF module switched and rewired in TF-TG collated edges of 1-to-1 orthologs across all five species") +
  scale_fill_gradient2(low = "yellow", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="% of collated\nedges") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  coord_fixed()
dev.off()

### coexpr_prom - geneB (TG) ###
coexpr_prom_geneB_pairwise_percentage_matrix <- as.matrix(read.table("coexpr_prom_geneB_pairwise_percentage_matrix",header = TRUE, row.names = 1, sep = "\t"))
melted_coexpr_prom_geneB_pairwise_percentage_matrix <- melt(coexpr_prom_geneB_pairwise_percentage_matrix)
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R–Pairwise_TG_modSW_rewiredTFTG_1to1.tiff', units="in", width=10, height=10, res=300)
ggplot(data = melted_coexpr_prom_geneB_pairwise_percentage_matrix, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(color = "white") + 
  labs(y="Edge present in species", x="Edge absent in species",title="Percentage pairwise TG module switched and rewired in TF-TG collated edges of 1-to-1 orthologs across all five species") +
  scale_fill_gradient2(low = "yellow", high = "red", mid = "orange", midpoint = 5, limit = c(0,10), space = "Lab", name="% of collated\nedges") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  coord_fixed()
dev.off()
