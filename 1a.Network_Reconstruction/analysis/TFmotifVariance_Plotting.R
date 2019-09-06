### Plotting variance over mean of TF motif enrichment and looking for correlation with state-changes 

library(ggplot2)
library(hexbin)
library(reshape2)
library(gridExtra)
library(grid)
library(magrittr)
library(cowplot)
library(egg)


### Sushmita carried out different analyses of looking at the correlation between motif enrichment and expression across tissues
# She prepared other required plots (Top5TF_scatter.pdf)
# She also prepared text files (*_perTFCC.txt) for the correlation values for each motif, across all tissues
# From this, create Fig. S-R1f: heatmap matrix for the correlation values for each motif, across all tissues
# She also shared merged_q_exp_with_header.txt that has:
  # 3rd to 8th columns are -log(qval) where qval was from the file you shared.
  # The remaining are species specific expression values. These are zero-meaned per species which was the input to Arboretum and per-gene greedy network inference.
# In unix, prepared six tissue-specific files where the correlation is mapped > XX

# For the heatmap:
# Make a heatmap showing both the original motif enrichment and expression levels for each motif-cluster pair.
# Do this for each of the six tissues (brain, eye, heart, kidney, muscle, testis) - so 6 heatmaps:
#   Matrix with 11 columns - 1-5: each species motif enrichment (-log_qval); 6-10: zero-meaned tissue expression in five species; 11: correlation. 
setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/Edits/srfigs/motif_exp_divanalysis/results")

br <- read.table("Br_perTFCC_mapped.txt", header=T)
ey <- read.table("Ey_perTFCC_mapped.txt", header=T)
ht <- read.table("Ht_perTFCC_mapped.txt", header=T)
kd <- read.table("Kd_perTFCC_mapped.txt", header=T)
ms <- read.table("Ms_perTFCC_mapped.txt", header=T)
ts <- read.table("Ts_perTFCC_mapped.txt", header=T)

### br
# create three melts for each - 1. qvalue enrichment; 2. expression; 3. pearson correlation. Then create plots (1 can have the ID) and then join, placing the legend for each on out right (as column of three rows).
br_melt_q <- melt(br, id.vars = "Cluster_TF", measure.var = c("Mz_q","Pn_q","Ab_q","Nb_q","On_q"))
br_melt_e <- melt(br, id.vars = "Cluster_TF", measure.var = c("Mz_Br","Pn_Br","Ab_Br","Nb_Br","On_Br"))
br_melt_p <- melt(br, id.vars = "Cluster_TF", measure.var = "Pearson_Correlation")

# determine max and min values
max(br_melt_q$value) # 22
min(br_melt_q$value) # 0
max(br_melt_e$value) # 4.12
min(br_melt_e$value) # -2
max(br_melt_p$value) # 1
min(br_melt_p$value) # -1

br_p1 <- ggplot(data = br_melt_q, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="Motif enrichment in\nmodule\ngene promoters") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 11, limit = c(0,22), space = "Lab", name="Motif enrichment:\n-log(q-value)\n(FDR < 0.05)") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.x=element_blank())

br_p2 <- ggplot(data = br_melt_e, aes(x=factor(variable), y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") +
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="TF\ntissue\nexpression") +
  scale_fill_gradient2(low = "#009e4a", high = "#9e0099", mid = "white", midpoint = 1.25, limit = c(-2,4.5), space = "Lab", name="TF tissue expression:\nZero-meaned\nexpression value") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

br_p3 <- ggplot(data = br_melt_p, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  labs(x="PCC of motif\nenrichment &\nTF expression") +
  scale_fill_gradient2(low = "#9e0000", high = "#03009e", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson correlation\ncoeeficient (PCC)") +
  scale_x_discrete(labels=c('PearsonCC')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

# extract legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

br_p1_legend <- g_legend(br_p1)
br_p2_legend <- g_legend(br_p2)
br_p3_legend <- g_legend(br_p3)

# store legends in one plot
br_p4 <- plot_grid(br_p1_legend,br_p2_legend,br_p3_legend,ncol=1,nrow=3, align="hv")

# Arrange the plots - first store the three heatmaps together (as ggarrange does not like working with the legend grob) and then save the file by adding the legend grob with grid.arrange
br_p5 <- ggarrange(br_p1 + theme(legend.position="none"),
          br_p2 + theme(legend.position="none"),
          br_p3 + theme(legend.position="none"),
          ncol=3)

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/ExtData_S-R1C.Brain_heatmap_pearson-correlation.tiff', units="in", width=7, height=160, res=200)
plot_grid(br_p5,br_p4,ncol=2)
dev.off()

### ey
# create three melts for each - 1. qvalue enrichment; 2. expression; 3. pearson correlation. Then create plots (1 can have the ID) and then join, placing the legend for each on out right (as column of three rows).
ey_melt_q <- melt(ey, id.vars = "Cluster_TF", measure.var = c("Mz_q","Pn_q","Ab_q","Nb_q","On_q"))
ey_melt_e <- melt(ey, id.vars = "Cluster_TF", measure.var = c("Mz_Ey","Pn_Ey","Ab_Ey","Nb_Ey","On_Ey"))
ey_melt_p <- melt(ey, id.vars = "Cluster_TF", measure.var = "Pearson_Correlation")

# determine max and min values
max(ey_melt_q$value) # 22
min(ey_melt_q$value) # 0
max(ey_melt_e$value) # 5
min(ey_melt_e$value) # -2.5
max(ey_melt_p$value) # 1
min(ey_melt_p$value) # -1

ey_p1 <- ggplot(data = ey_melt_q, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="Motif enrichment in\nmodule\ngene promoters") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 11, limit = c(0,22), space = "Lab", name="Motif enrichment:\n-log(q-value)\n(FDR < 0.05)") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.x=element_blank())

ey_p2 <- ggplot(data = ey_melt_e, aes(x=factor(variable), y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") +
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="TF\ntissue\nexpression") +
  scale_fill_gradient2(low = "#009e4a", high = "#9e0099", mid = "white", midpoint = 1.25, limit = c(-2.5,5), space = "Lab", name="TF tissue expression:\nZero-meaned\nexpression value") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

ey_p3 <- ggplot(data = ey_melt_p, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  labs(x="PCC of motif\nenrichment &\nTF expression") +
  scale_fill_gradient2(low = "#9e0000", high = "#03009e", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson correlation\ncoeeficient (PCC)") +
  scale_x_discrete(labels=c('PearsonCC')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

# extract legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

ey_p1_legend <- g_legend(ey_p1)
ey_p2_legend <- g_legend(ey_p2)
ey_p3_legend <- g_legend(ey_p3)

# store legends in one plot
ey_p4 <- plot_grid(ey_p1_legend,ey_p2_legend,ey_p3_legend,ncol=1,nrow=3, align="hv")

# Arrange the plots - first store the three heatmaps together (as ggarrange does not like working with the legend grob) and then save the file by adding the legend grob with grid.arrange
ey_p5 <- ggarrange(ey_p1 + theme(legend.position="none"),
                   ey_p2 + theme(legend.position="none"),
                   ey_p3 + theme(legend.position="none"),
                   ncol=3)

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/ExtData_S-R1D.Eye_heatmap_pearson-correlation.tiff', units="in", width=7, height=160, res=200)
plot_grid(ey_p5,ey_p4,ncol=2)
dev.off()

### ht
# create three melts for each - 1. qvalue enrichment; 2. expression; 3. pearson correlation. Then create plots (1 can have the ID) and then join, placing the legend for each on out right (as column of three rows).
ht_melt_q <- melt(ht, id.vars = "Cluster_TF", measure.var = c("Mz_q","Pn_q","Ab_q","Nb_q","On_q"))
ht_melt_e <- melt(ht, id.vars = "Cluster_TF", measure.var = c("Mz_Ht","Pn_Ht","Ab_Ht","Nb_Ht","On_Ht"))
ht_melt_p <- melt(ht, id.vars = "Cluster_TF", measure.var = "Pearson_Correlation")

# determine max and min values
max(ht_melt_q$value) # 22
min(ht_melt_q$value) # 0
max(ht_melt_e$value) # 3.5
min(ht_melt_e$value) # -2.1
max(ht_melt_p$value) # 1
min(ht_melt_p$value) # -1

ht_p1 <- ggplot(data = ht_melt_q, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="Motif enrichment in\nmodule\ngene promoters") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 11, limit = c(0,22), space = "Lab", name="Motif enrichment:\n-log(q-value)\n(FDR < 0.05)") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.x=element_blank())

ht_p2 <- ggplot(data = ht_melt_e, aes(x=factor(variable), y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") +
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="TF\ntissue\nexpression") +
  scale_fill_gradient2(low = "#009e4a", high = "#9e0099", mid = "white", midpoint = 1.25, limit = c(-2.1,3.5), space = "Lab", name="TF tissue expression:\nZero-meaned\nexpression value") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

ht_p3 <- ggplot(data = ht_melt_p, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  labs(x="PCC of motif\nenrichment &\nTF expression") +
  scale_fill_gradient2(low = "#9e0000", high = "#03009e", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson correlation\ncoeeficient (PCC)") +
  scale_x_discrete(labels=c('PearsonCC')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

# extract legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

ht_p1_legend <- g_legend(ht_p1)
ht_p2_legend <- g_legend(ht_p2)
ht_p3_legend <- g_legend(ht_p3)

# store legends in one plot
ht_p4 <- plot_grid(ht_p1_legend,ht_p2_legend,ht_p3_legend,ncol=1,nrow=3, align="hv")

# Arrange the plots - first store the three heatmaps together (as ggarrange does not like working with the legend grob) and then save the file by adding the legend grob with grid.arrange
ht_p5 <- ggarrange(ht_p1 + theme(legend.position="none"),
                   ht_p2 + theme(legend.position="none"),
                   ht_p3 + theme(legend.position="none"),
                   ncol=3)

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/ExtData_S-R1E.Heart_heatmap_pearson-correlation.tiff', units="in", width=7, height=160, res=200)
plot_grid(ht_p5,ht_p4,ncol=2)
dev.off()

### kd
# create three melts for each - 1. qvalue enrichment; 2. expression; 3. pearson correlation. Then create plots (1 can have the ID) and then join, placing the legend for each on out rigkd (as column of three rows).
kd_melt_q <- melt(kd, id.vars = "Cluster_TF", measure.var = c("Mz_q","Pn_q","Ab_q","Nb_q","On_q"))
kd_melt_e <- melt(kd, id.vars = "Cluster_TF", measure.var = c("Mz_Kd","Pn_Kd","Ab_Kd","Nb_Kd","On_Kd"))
kd_melt_p <- melt(kd, id.vars = "Cluster_TF", measure.var = "Pearson_Correlation")

# determine max and min values
max(kd_melt_q$value) # 22
min(kd_melt_q$value) # 0
max(kd_melt_e$value) # 4.5
min(kd_melt_e$value) # -2.3
max(kd_melt_p$value) # 1
min(kd_melt_p$value) # -1

kd_p1 <- ggplot(data = kd_melt_q, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="Motif enrichment in\nmodule\ngene promoters") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 11, limit = c(0,22), space = "Lab", name="Motif enrichment:\n-log(q-value)\n(FDR < 0.05)") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.x=element_blank())

kd_p2 <- ggplot(data = kd_melt_e, aes(x=factor(variable), y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") +
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="TF\ntissue\nexpression") +
  scale_fill_gradient2(low = "#009e4a", high = "#9e0099", mid = "white", midpoint = 1.25, limit = c(-2.3,4.5), space = "Lab", name="TF tissue expression:\nZero-meaned\nexpression value") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

kd_p3 <- ggplot(data = kd_melt_p, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  labs(x="PCC of motif\nenrichment &\nTF expression") +
  scale_fill_gradient2(low = "#9e0000", high = "#03009e", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson correlation\ncoeeficient (PCC)") +
  scale_x_discrete(labels=c('PearsonCC')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

# extract legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

kd_p1_legend <- g_legend(kd_p1)
kd_p2_legend <- g_legend(kd_p2)
kd_p3_legend <- g_legend(kd_p3)

# store legends in one plot
kd_p4 <- plot_grid(kd_p1_legend,kd_p2_legend,kd_p3_legend,ncol=1,nrow=3, align="hv")

# Arrange the plots - first store the three heatmaps together (as ggarrange does not like working with the legend grob) and then save the file by adding the legend grob with grid.arrange
kd_p5 <- ggarrange(kd_p1 + theme(legend.position="none"),
                   kd_p2 + theme(legend.position="none"),
                   kd_p3 + theme(legend.position="none"),
                   ncol=3)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/ExtData_S-R1F.Kidney_heatmap_pearson-correlation.tiff', units="in", width=7, height=160, res=200)
plot_grid(kd_p5,kd_p4,ncol=2)
dev.off()

### ms
# create three melts for each - 1. qvalue enrichment; 2. expression; 3. pearson correlation. Then create plots (1 can have the ID) and then join, placing the legend for each on out rigms (as column of three rows).
ms_melt_q <- melt(ms, id.vars = "Cluster_TF", measure.var = c("Mz_q","Pn_q","Ab_q","Nb_q","On_q"))
ms_melt_e <- melt(ms, id.vars = "Cluster_TF", measure.var = c("Mz_Ms","Pn_Ms","Ab_Ms","Nb_Ms","On_Ms"))
ms_melt_p <- melt(ms, id.vars = "Cluster_TF", measure.var = "Pearson_Correlation")

# determine max and min values
max(ms_melt_q$value) # 22
min(ms_melt_q$value) # 0
max(ms_melt_e$value) # 4
min(ms_melt_e$value) # -3.1
max(ms_melt_p$value) # 1
min(ms_melt_p$value) # -1

ms_p1 <- ggplot(data = ms_melt_q, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="Motif enrichment in\nmodule\ngene promoters") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 11, limit = c(0,22), space = "Lab", name="Motif enrichment:\n-log(q-value)\n(FDR < 0.05)") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.x=element_blank())

ms_p2 <- ggplot(data = ms_melt_e, aes(x=factor(variable), y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") +
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="TF\ntissue\nexpression") +
  scale_fill_gradient2(low = "#009e4a", high = "#9e0099", mid = "white", midpoint = 1.25, limit = c(-3.1,4), space = "Lab", name="TF tissue expression:\nZero-meaned\nexpression value") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

ms_p3 <- ggplot(data = ms_melt_p, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  labs(x="PCC of motif\nenrichment &\nTF expression") +
  scale_fill_gradient2(low = "#9e0000", high = "#03009e", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson correlation\ncoeeficient (PCC)") +
  scale_x_discrete(labels=c('PearsonCC')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

# extract legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

ms_p1_legend <- g_legend(ms_p1)
ms_p2_legend <- g_legend(ms_p2)
ms_p3_legend <- g_legend(ms_p3)

# store legends in one plot
ms_p4 <- plot_grid(ms_p1_legend,ms_p2_legend,ms_p3_legend,ncol=1,nrow=3, align="hv")

# Arrange the plots - first store the three heatmaps together (as ggarrange does not like working with the legend grob) and then save the file by adding the legend grob with grid.arrange
ms_p5 <- ggarrange(ms_p1 + theme(legend.position="none"),
                   ms_p2 + theme(legend.position="none"),
                   ms_p3 + theme(legend.position="none"),
                   ncol=3)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/ExtData_S-R1G.Muscle_heatmap_pearson-correlation.tiff', units="in", width=7, height=160, res=200)
plot_grid(ms_p5,ms_p4,ncol=2)
dev.off()

### ts
# create three melts for each - 1. qvalue enrichment; 2. expression; 3. pearson correlation. Then create plots (1 can have the ID) and then join, placing the legend for each on out rigts (as column of three rows).
ts_melt_q <- melt(ts, id.vars = "Cluster_TF", measure.var = c("Mz_q","Pn_q","Ab_q","Nb_q","On_q"))
ts_melt_e <- melt(ts, id.vars = "Cluster_TF", measure.var = c("Mz_Ts","Pn_Ts","Ab_Ts","Nb_Ts","On_Ts"))
ts_melt_p <- melt(ts, id.vars = "Cluster_TF", measure.var = "Pearson_Correlation")

# determine max and min values
max(ts_melt_q$value) # 22
min(ts_melt_q$value) # 0
max(ts_melt_e$value) # 3.5
min(ts_melt_e$value) # -2.1
max(ts_melt_p$value) # 1
min(ts_melt_p$value) # -1

ts_p1 <- ggplot(data = ts_melt_q, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="Motif enrichment in\nmodule\ngene promoters") +
  scale_fill_gradient2(low = "white", high = "red", mid = "orange", midpoint = 11, limit = c(0,22), space = "Lab", name="Motif enrichment:\n-log(q-value)\n(FDR < 0.05)") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.x=element_blank())

ts_p2 <- ggplot(data = ts_melt_e, aes(x=factor(variable), y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") +
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  # labs(x="TF\ntissue\nexpression") +
  scale_fill_gradient2(low = "#009e4a", high = "#9e0099", mid = "white", midpoint = 1.25, limit = c(-2.1,3.5), space = "Lab", name="TF tissue expression:\nZero-meaned\nexpression value") +
  scale_x_discrete(labels=c('M. zebra', 'P. nyererei','A. burtoni','N. brichardi','O. niloticus')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

ts_p3 <- ggplot(data = ts_melt_p, aes(x=variable, y=Cluster_TF, fill=value)) + 
  geom_tile(color = "black") + 
  coord_fixed(ratio=0.2) +
  labs(y="Module number and enriched TF motif") +
  labs(x="PCC of motif\nenrichment &\nTF expression") +
  scale_fill_gradient2(low = "#9e0000", high = "#03009e", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson correlation\ncoeeficient (PCC)") +
  scale_x_discrete(labels=c('PearsonCC')) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
        axis.text.y = element_text(size=6),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=10),
        legend.title=element_text(size=13,face="bold"),
        legend.key.size = unit(0.7,"cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank())

# extract legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

ts_p1_legend <- g_legend(ts_p1)
ts_p2_legend <- g_legend(ts_p2)
ts_p3_legend <- g_legend(ts_p3)

# store legends in one plot
ts_p4 <- plot_grid(ts_p1_legend,ts_p2_legend,ts_p3_legend,ncol=1,nrow=3, align="hv")

# Arrange the plots - first store the three heatmaps together (as ggarrange does not like working with the legend grob) and then save the file by adding the legend grob with grid.arrange
ts_p5 <- ggarrange(ts_p1 + theme(legend.position="none"),
                   ts_p2 + theme(legend.position="none"),
                   ts_p3 + theme(legend.position="none"),
                   ncol=3)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/ExtData_S-R1H.Testis_heatmap_pearson-correlation.tiff', units="in", width=7, height=160, res=200)
plot_grid(ts_p5,ts_p4,ncol=2)
dev.off()


##########################################################################################################################
# Test plottings

# These didn't work so well and either A) didn't arrange plots as required, or B) introduced a lot of space between plots?

# tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/Fig.S-R1fA.Brain_heatmap_pearson-correlation.tiff', units="in", width=5, height=160, res=200)
# ggarrange(br_p1 + theme(legend.position="none"),
#           br_p2 + theme(legend.position="none"),
#           br_p3 + theme(legend.position="none"),
#           br_p4, # consider removing, saving legends as separate file and pasting in - google ways to paste files?
#           ncol=3)
# dev.off()

# tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/Fig.S-R1fA.Brain_heatmap_pearson-correlation.tiff', units="in", width=30, height=160, res=100)
# grid.arrange(br_p1 + theme(legend.position="none"),br_p2 + theme(legend.position="none"),br_p3 + theme(legend.position="none"),br_p4,ncol=4)
# dev.off()

# tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/Fig.S-R1fA.Brain_heatmap_pearson-correlation.tiff', units="in", width=30, height=190, res=100)
# plot_grid(
#   br_p1 + theme(legend.position="none") + theme(plot.margin = unit(c(2, -1000, 2, -1000), "cm")),
#   br_p2 + theme(legend.position="none") + theme(plot.margin = unit(c(2, -1000, 2, -1000), "cm")),
#   br_p3 + theme(legend.position="none") + theme(plot.margin = unit(c(2, -1000, 2, -1000), "cm")),
#   br_p4,
#   ncol=4,
#   # nrow=1,
#   # labels = c('A', 'B', 'C'),
#   # label_size = 18,
#   align="hv"
# )
# dev.off()

#################################################################################################################################
### Sushmita carried out different analyses of looking at the correlation between motif enrichment and expression across tissues
# Did not use analyses below as a result

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/May2019_ReviewerComments/1a.TFenrich-statechange_corr")

# 1a. Load files and add colnames - NOTE: Hs and Mm will be exactly the same since they are JASPAR derived - Just generate for one
ConsHs = read.table("All_MOTIFS_OUTPUT_q0.05_details.Conserved.-log10qval-Mean-SVAR-StDev.Hs.Nsw-Sw.txt", header=FALSE)
# ConsMm = read.table("All_MOTIFS_OUTPUT_q0.05_details.Conserved.-log10qval-Mean-SVAR-StDev.Mm.Nsw-Sw.txt", header=FALSE)
NotConsHs = read.table("All_MOTIFS_OUTPUT_q0.05_details.NotConserved.-log10qval-Mean-SVAR-StDev.Hs.Nsw-Sw.txt", header=FALSE)
# NotConsMm = read.table("All_MOTIFS_OUTPUT_q0.05_details.NotConserved.-log10qval-Mean-SVAR-StDev.Mm.Nsw-Sw.txt", header=FALSE)
colnames(ConsHs) <- c('Module','TF gene ID','Mean', 'Variance', 'Standard_Deviation', 'Mz ID', 'Pn ID', 'Ab ID', 'Nb ID', 'On ID', 'motif ID', 'Hs/Mm ENSEMBL ID', 'Orthogroup','State_Change')
# colnames(ConsMm) <- c('Module','TF gene ID','Mean', 'Variance', 'Standard_Deviation', 'Mz ID', 'Pn ID', 'Ab ID', 'Nb ID', 'On ID', 'motif ID', 'Hs/Mm ENSEMBL ID', 'Orthogroup','State_Change')
colnames(NotConsHs) <- c('Module','TF gene ID','Mean', 'Variance', 'Standard_Deviation', 'Mz ID', 'Pn ID', 'Ab ID', 'Nb ID', 'On ID', 'motif ID', 'Hs/Mm ENSEMBL ID', 'Orthogroup','State_Change')
# colnames(NotConsMm) <- c('Module','TF gene ID','Mean', 'Variance', 'Standard_Deviation', 'Mz ID', 'Pn ID', 'Ab ID', 'Nb ID', 'On ID', 'motif ID', 'Hs/Mm ENSEMBL ID', 'Orthogroup','State_Change')

# function to replace nan with 0 in dataframes
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

ConsHs[is.nan(ConsHs)] <- 0 # replace nan values with 0
# ConsMm[is.nan(ConsMm)] <- 0 # replace nan values with 0
NotConsHs[is.nan(NotConsHs)] <- 0 # replace nan values with 0
# NotConsMm[is.nan(NotConsMm)] <- 0 # replace nan values with 0


# box plot for the fold enrichment of TFs that are state-changed (1) vs non state-changed TFs

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/1aA.Boxplot_of_non_and_state_change_TFmotifenrichment-ConservedTFenrich.tiff', units="in", width=8, height=6, res=100)
ggplot(ConsHs, aes(x = factor(State_Change), y = log2(Mean), fill=factor(State_Change))) + 
  geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_x_discrete(labels=c("Non-state-changed", "State-changed")) +
  labs (x = "Stability of TF") + labs (y = "log2 mean fold enrichment (-log10 q-val) of TF motifs conserved in all five species") + theme(legend.position="none")
dev.off()

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/1aB.Boxplot_of_non_and_state_change_TFmotifenrichment-NotConservedTFenrich.tiff', units="in", width=8, height=6, res=100)
ggplot(NotConsHs, aes(x = factor(State_Change), y = log2(Mean), fill=factor(State_Change))) + 
  geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_x_discrete(labels=c("Non-state-changed", "State-changed")) +
  labs (x = "Stability of TF") + labs (y = "log2 mean fold enrichment (-log10 q-val) of TF motifs not conserved in all five species") + theme(legend.position="none")
dev.off()



# Determine the difference in number of state-changed vs non-state-changed TFs at similiar FE (Variance 0-5 of -log10(q-val)), Mid-level FE (5-10 of -log10(q-val)), dissimilar FE (10+ of -log10(q-val)) 
sum(ConsHs$State_Change == "0")
sum(ConsHs$State_Change == "1")

sum(ConsHs$State_Change == "0" & ConsHs$Variance > 0 & ConsHs$Variance < 1) # 290/765 = 37.9%
sum(ConsHs$State_Change == "1" & ConsHs$Variance > 0 & ConsHs$Variance < 1) # 104/317 = 32.8%

sum(ConsHs$State_Change == "0" & ConsHs$Variance > 0 & ConsHs$Variance < 2) # 475/765 = 62.1%
sum(ConsHs$State_Change == "1" & ConsHs$Variance > 0 & ConsHs$Variance < 2) # 193/317 = 60.9%

sum(ConsHs$State_Change == "0" & ConsHs$Variance > 0 & ConsHs$Variance < 3) # 585/765 = 76.5%
sum(ConsHs$State_Change == "1" & ConsHs$Variance > 0 & ConsHs$Variance < 3) # 228/317 = 71.9%

sum(ConsHs$State_Change == "0" & ConsHs$Variance > 0 & ConsHs$Variance < 4) # 633/765 = 82.7%
sum(ConsHs$State_Change == "1" & ConsHs$Variance > 0 & ConsHs$Variance < 4) # 262/317 = 82.6%

sum(ConsHs$State_Change == "0" & ConsHs$Variance > 0 & ConsHs$Variance < 5) # 668/765 = 87.3%
sum(ConsHs$State_Change == "1" & ConsHs$Variance > 0 & ConsHs$Variance < 5) # 278/317 = 87.7%

sum(ConsHs$State_Change == "0" & ConsHs$Variance > 5 & ConsHs$Variance < 10) # 82/765 = 10.7%
sum(ConsHs$State_Change == "1" & ConsHs$Variance > 5 & ConsHs$Variance < 10) # 32/317 = 10.1%

sum(ConsHs$State_Change == "0" & ConsHs$Variance > 5) # 97/765 = 12.6%
sum(ConsHs$State_Change == "1" & ConsHs$Variance > 5) # 39/317 = 12.3%

sum(ConsHs$State_Change == "0" & ConsHs$Variance > 10) # 15/765 = 2.0%
sum(ConsHs$State_Change == "1" & ConsHs$Variance > 10) # 7/317 = 2.2%

mean(ConsHs[["Variance"]]) # 2.35 - any value above the mean considered variation on fold enrichment, any value below considered similar FE

sum(ConsHs$State_Change == "0" & ConsHs$Variance > 2.35) # 253/765 = 33%
sum(ConsHs$State_Change == "1" & ConsHs$Variance > 2.35) # 111/317 = 35%

sum(ConsHs$State_Change == "0" & ConsHs$Variance < 2.35) # 512/765 = 67%
sum(ConsHs$State_Change == "1" & ConsHs$Variance < 2.35) # 206/317 = 65%

# sum(NotConsHs$State_Change == "0") # 1441
# sum(NotConsHs$State_Change == "1") # 629
# 
# sum(NotConsHs$State_Change == "0" & NotConsHs$Variance > 0 & NotConsHs$Variance < 1)/sum(NotConsHs$State_Change == "0")*100 # 809/1441 = 56.1%
# sum(NotConsHs$State_Change == "1" & NotConsHs$Variance > 0 & NotConsHs$Variance < 1)/sum(NotConsHs$State_Change == "1")*100 # 357/629 = 56.8%
# 
# sum(NotConsHs$State_Change == "0" & NotConsHs$Variance > 0 & NotConsHs$Variance < 2)/sum(NotConsHs$State_Change == "0")*100 # 930/1441 = 64.5%
# sum(NotConsHs$State_Change == "1" & NotConsHs$Variance > 0 & NotConsHs$Variance < 2)/sum(NotConsHs$State_Change == "1")*100 # 409/629 = 65.0%
# 
# sum(NotConsHs$State_Change == "0" & NotConsHs$Variance > 0 & NotConsHs$Variance < 3)/sum(NotConsHs$State_Change == "0")*100 # 968/1441 = 67.2%
# sum(NotConsHs$State_Change == "1" & NotConsHs$Variance > 0 & NotConsHs$Variance < 3)/sum(NotConsHs$State_Change == "1")*100 # 424/629 = 67.4%
# 
# sum(NotConsHs$State_Change == "0" & NotConsHs$Variance > 0 & NotConsHs$Variance < 4)/sum(NotConsHs$State_Change == "0")*100 # 983/1441 = 68.2%
# sum(NotConsHs$State_Change == "1" & NotConsHs$Variance > 0 & NotConsHs$Variance < 4)/sum(NotConsHs$State_Change == "1")*100 # 429/629 = 68.2%
# 
# sum(NotConsHs$State_Change == "0" & NotConsHs$Variance > 0 & NotConsHs$Variance < 5)/sum(NotConsHs$State_Change == "0")*100 # 991/1441 = 68.8%
# sum(NotConsHs$State_Change == "1" & NotConsHs$Variance > 0 & NotConsHs$Variance < 5)/sum(NotConsHs$State_Change == "1")*100 # 431/629 = 68.5%
# 
# sum(NotConsHs$State_Change == "0" & NotConsHs$Variance > 5 & NotConsHs$Variance < 10)/sum(NotConsHs$State_Change == "0")*100 # 9/1441 = 0.6%
# sum(NotConsHs$State_Change == "1" & NotConsHs$Variance > 5 & NotConsHs$Variance < 10)/sum(NotConsHs$State_Change == "1")*100 # 4/629 = 0.6%
# 
# sum(NotConsHs$State_Change == "0" & NotConsHs$Variance > 10)/sum(NotConsHs$State_Change == "0")*100 # 1/1441 = 0.07%
# sum(NotConsHs$State_Change == "1" & NotConsHs$Variance > 10)/sum(NotConsHs$State_Change == "1")*100 # 1/629 = 0.16%

# You need to plot, on the same plot, the density of variance of state-changed (1) and non state-changed TFs
# ggplot(ConsHs, aes(x=Variance, fill=factor(State_Change))) + geom_density(alpha=.3)

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/1aA.Variance_distribution_of_non_and_state_change_TFmotifenrichment-ConservedTFenrich.tiff', units="in", width=13, height=6, res=300)
ggplot(ConsHs, aes(x=Variance)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, aes(fill=factor(State_Change))) + # Overlay with transparent density plot
  guides(fill=guide_legend(title="TF Module assignment\nstability")) + 
  labs (x = "Variance of TF Motif enrichment as -log10(q-value)") + labs (y = "Density") +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black")) +
  scale_fill_discrete(breaks=c("0", "1"),
                      labels=c("Non-state-changed", "State-changed")) +
  scale_x_continuous(breaks=seq(0,35,1)) +
  scale_y_continuous(breaks=seq(0,0.4,0.025)) +
  geom_vline(xintercept=2.35, linetype="dotted", colour="red", size=1)
dev.off()

# tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/1aC.Variance_distribution_of_non_and_state_change_TFmotifenrichment-ConservedTFenrich_MmDerived.tiff', units="in", width=13, height=6, res=300)
# ggplot(ConsMm, aes(x=Variance)) + 
#   geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
#                  binwidth=.5,
#                  colour="black", fill="white") +
#   geom_density(alpha=.2, aes(fill=factor(State_Change))) + # Overlay with transparent density plot
#   guides(fill=guide_legend(title="No State-Change (0)\nor State-change (1)")) + 
#   labs (x = "Variance of TF Motif enrichment as -log10(q-value))") + labs (y = "Density") +
#   theme(axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16),  
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=20),
#         strip.text.x = element_text(size = 20, colour = "black"))
# dev.off()

# tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/1aX.Variance_distribution_of_non_and_state_change_TFmotifenrichment-NotConservedTFenrich.tiff', units="in", width=13, height=6, res=300)
# ggplot(NotConsHs, aes(x=Variance)) + 
#   geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
#                  binwidth=.5,
#                  colour="black", fill="white") +
#   geom_density(alpha=.2, aes(fill=factor(State_Change))) + # Overlay with transparent density plot
#   guides(fill=guide_legend(title="TF Module assignment\nstability")) + 
#   labs (x = "Variance of TF Motif enrichment as -log10(q-value))") + labs (y = "Density") +
#   theme(axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16),  
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=20),
#         strip.text.x = element_text(size = 20, colour = "black")) +
#   scale_fill_discrete(breaks=c("0", "1"),
#                       labels=c("Non-state-changed", "State-changed")) +
#   scale_x_continuous(breaks=seq(0,20,1)) +
#   scale_y_continuous(breaks=seq(0,3,0.25))
# dev.off() # Exclude this one as it will not represent all taxa

# tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/1aD.Variance_distribution_of_non_and_state_change_TFmotifenrichment-NotConservedTFenrich_MmDerived.tiff', units="in", width=13, height=6, res=300)
# ggplot(NotConsMm, aes(x=Variance)) + 
#   geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
#                  binwidth=.5,
#                  colour="black", fill="white") +
#   geom_density(alpha=.2, aes(fill=factor(State_Change))) + # Overlay with transparent density plot
#   guides(fill=guide_legend(title="No State-Change (0)\nor State-change (1)")) + 
#   labs (x = "Variance of TF Motif enrichment as -log10(q-value))") + labs (y = "Density") +
#   theme(axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16),  
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=20),
#         strip.text.x = element_text(size = 20, colour = "black"))
# dev.off()

# Split each of the four sets above by module and then plot below
split.ConsHs <- split(ConsHs, ConsHs$Module) # split the dataframe by module and store each as elements in list
names(split.ConsHs) <- c("ConsHs_Mod1", "ConsHs_Mod2", "ConsHs_Mod3", "ConsHs_Mod6", "ConsHs_Mod7", "ConsHs_Mod8", "ConsHs_Mod9") # rename elements
# list2env(split.ConsHs,envir=.GlobalEnv) # convert all elements of list into separate dfs
# Files1 <- c("ConsHs_Mod1", "ConsHs_Mod2", "ConsHs_Mod3", "ConsHs_Mod6", "ConsHs_Mod7", "ConsHs_Mod8", "ConsHs_Mod9")
split.NotConsHs <- split(NotConsHs, NotConsHs$Module) # split the dataframe by module and store each as elements in list
names(split.NotConsHs) <- c("ConsHs_Mod1", "ConsHs_Mod2", "ConsHs_Mod3", "ConsHs_Mod4", "ConsHs_Mod5", "ConsHs_Mod6", "ConsHs_Mod7", "ConsHs_Mod8", "ConsHs_Mod9") # rename elements
# split.ConsMm <- split(ConsMm, ConsMm$Module) # split the dataframe by module and store each as elements in list
# names(split.ConsMm) <- c("ConsHs_Mod1", "ConsHs_Mod2", "ConsHs_Mod3", "ConsHs_Mod6", "ConsHs_Mod7", "ConsHs_Mod8", "ConsHs_Mod9") # rename elements
# split.NotConsMm <- split(NotConsMm, NotConsMm$Module) # split the dataframe by module and store each as elements in list
# names(split.NotConsMm) <- c("ConsHs_Mod1", "ConsHs_Mod2", "ConsHs_Mod3", "ConsHs_Mod4", "ConsHs_Mod5", "ConsHs_Mod6", "ConsHs_Mod7", "ConsHs_Mod8", "ConsHs_Mod9") # rename elements


# Make plots of the individual modules and store in a list
ConsHs_modulePlots = list()
for (i in split.ConsHs) {
  ConsHs_modulePlots1 = ggplot(i, aes(x=Variance)) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    geom_density(alpha=.2, aes(fill=factor(State_Change))) + # Overlay with transparent density plot
    guides(fill=guide_legend(title="TF Module assignment\nstability")) +
    labs (x = "Variance of TF Motif enrichment\nas -log10(q-value))") + labs (y = "Density") + labs(title=paste(i[1,1])) + # this will paste the module number as title
    theme(axis.text.x = element_text(size=4),
          axis.text.y = element_text(size=4),  
          axis.title.x = element_text(size=6),
          axis.title.y = element_text(size=6),
          strip.text.x = element_text(size = 6, colour = "black")) + 
    expand_limits(x = 0, y = 0) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.title=element_text(size=6), 
          legend.text=element_text(size=6),
          legend.background = element_rect(fill = "lightgray")) +
    scale_fill_discrete(breaks=c("0", "1"),
                        labels=c("Non-state-changed", "State-changed"))
  ConsHs_modulePlots <- c(ConsHs_modulePlots, list(ConsHs_modulePlots1))
} 

ConsHs_modulePlots_glist <- lapply(ConsHs_modulePlots, ggplotGrob) # if you don't do this and ggsave immediately, a blank first page is created
ggsave(filename = "/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/1aB.ConservedMotifs_Modules.pdf", marrangeGrob(ConsHs_modulePlots_glist, nrow=2, ncol=2, top=NULL)) # arrange them 2x2 in multiple pages

# NotConsHs_modulePlots = list()
# for (i in split.NotConsHs) {
#   NotConsHs_modulePlots1 = ggplot(i, aes(x=Variance)) + 
#     geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
#                    binwidth=.5,
#                    colour="black", fill="white") +
#     geom_density(alpha=.2, aes(fill=factor(State_Change))) + # Overlay with transparent density plot
#     guides(fill=guide_legend(title="TF Module assignment\nstability")) + 
#     labs (x = "Variance of TF Motif enrichment\nas -log10(q-value))") + labs (y = "Density") + labs(title=paste(i[1,1])) + # this will paste the module number as title
#     theme(axis.text.x = element_text(size=4),
#           axis.text.y = element_text(size=4),  
#           axis.title.x = element_text(size=6),
#           axis.title.y = element_text(size=6),
#           strip.text.x = element_text(size = 6, colour = "black")) + 
#     expand_limits(x = 0, y = 0) +
#     theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
#           legend.title=element_text(size=6), 
#           legend.text=element_text(size=6),
#           legend.background = element_rect(fill = "lightgray")) +
#     scale_fill_discrete(breaks=c("0", "1"),
#                         labels=c("Non-state-changed", "State-changed"))
#   NotConsHs_modulePlots <- c(NotConsHs_modulePlots, list(NotConsHs_modulePlots1))
# } 
# 
# NotConsHs_modulePlots_glist <- lapply(NotConsHs_modulePlots, ggplotGrob) # if you don't do this and ggsave immediately, a blank first page is created
# ggsave(filename = "/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/1aD.NotConservedMotifs_Modules.pdf", marrangeGrob(NotConsHs_modulePlots_glist, nrow=2, ncol=2, top=NULL)) # arrange them 2x2 in multiple pages


# # Save plots to tiff. Makes a separate file for each plot.
# pdf("all.pdf")
# invisible(lapply(ConsHs_modulePlots, print))
# dev.off()
# 
# invisible(mapply(ggsave, file=paste0("plot-", names(ConsHs_modulePlots), ".pdf"), plot=l)) # individual plots




