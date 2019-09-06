### 1.6. Correlation of TF-TG expression to support sws1 and rho expression variation driven by selected TFs (NR2C2 and Gata2) specifically.

library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(ggpubr)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/May2019_ReviewerComments/1.6_eye_exprcorrelation/")

# Load in the data frames - rows names are the geneID and cols are normalised expression values in brain, eye, heart, kidney, muscle, testis
mz_sws1 =read.table("Mz_sws1_network.txtD", header=FALSE, row.names=1, sep="\t")
nb_sws1 =read.table("Nb_sws1_network.txtD", header=FALSE, row.names=1, sep="\t")

mz_rho =read.table("Mz_rho_network.txtD", header=FALSE, row.names=1, sep="\t")
ab_rho =read.table("Ab_rho_network.txtD", header=FALSE, row.names=1, sep="\t")
on_rho =read.table("On_rho_network.txtD", header=FALSE, row.names=1, sep="\t")

# log2(x+1) transformation of all expression values
mz_sws1_log2 = log2(mz_sws1+1)
nb_sws1_log2 = log2(nb_sws1+1)
mz_rho_log2 = log2(mz_rho+1)
ab_rho_log2 = log2(ab_rho+1)
on_rho_log2 = log2(on_rho+1)

# transpose so that tissues are row.names and genes are headers
mz_sws1_log2_t <- t(mz_sws1_log2)
row.names(mz_sws1_log2_t) <- c("Brain","Eye","Heart","Kidney","Muscle","Testis")
nb_sws1_log2_t <- t(nb_sws1_log2)
row.names(nb_sws1_log2_t) <- c("Brain","Eye","Heart","Kidney","Muscle","Testis")
mz_rho_log2_t <- t(mz_rho_log2)
row.names(mz_rho_log2_t) <- c("Brain","Eye","Heart","Kidney","Muscle","Testis")
ab_rho_log2_t <- t(ab_rho_log2)
row.names(ab_rho_log2_t) <- c("Brain","Eye","Heart","Kidney","Muscle","Testis")
on_rho_log2_t <- t(on_rho_log2)
row.names(on_rho_log2_t) <- c("Brain","Eye","Heart","Kidney","Muscle","Testis")

# gene-by-gene pearson correlation based on normalised and log2 transformed expression across 6 tissues
# 1 is total positive linear correlation, 0 is no linear correlation, and −1 is total negative linear correlation
mz_sws1_log2_t_cor <- cor(mz_sws1_log2_t,method="pearson")
nb_sws1_log2_t_cor <- cor(nb_sws1_log2_t,method="pearson")
mz_rho_log2_t_cor <- cor(mz_rho_log2_t,method="pearson")
ab_rho_log2_t_cor <- cor(ab_rho_log2_t,method="pearson")
on_rho_log2_t_cor <- cor(on_rho_log2_t,method="pearson")

# library(Hmisc) # looking at p-vals below
# mz_sws1_log2_t_corP <- rcorr(mz_sws1_log2_t,type="pearson")
# nb_sws1_log2_t_corP <- rcorr(nb_sws1_log2_t,type="pearson")

## The below would not work since two data points is not enough
# # Since the genes and networks you are testing for correlation are mainly brain and eye specific, test for correlation using only those two tissues
# # Ideally (and if we had multiple individuals) we would just use eye expression
# mz_sws1_log2_t_BE <- mz_sws1_log2_t[1:2,]
# nb_sws1_log2_t_BE <- nb_sws1_log2_t[1:2,]
# mz_sws1_log2_t_BE_cor <- cor(mz_sws1_log2_t_BE,method="pearson")
# nb_sws1_log2_t_BE_cor <- cor(nb_sws1_log2_t_BE,method="pearson")

# Plotting correlation heatmap
# correlation matrix has redundant information. We’ll use the functions below to set half of it to NA.
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
# reorder the correlation matrix according to the correlation coefficient. This is useful to identify the hidden pattern in the matrix. hclust for hierarchical clustering order
# Helper function to reorder the correlation matrix
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


#### NOTE: as of 05/08/2019, M. zebra of below will need amending for colouring and bolding certain genes.
# this is because I added genes we did not predict (NR2C2 > sws1 and GATA2A > rho) to look at whether their correlation was less (it was not!)
# Kept the old plots instead - no new plots have been generated

# Mz sws1
mz_sws1_log2_t_cor <- reorder_cormat(mz_sws1_log2_t_cor) # reorder the correlation matrix to identify patterns
mz_sws1_log2_t_corUT <- get_upper_tri(mz_sws1_log2_t_cor) # take the upper triangle only since rest is redundant
melted_mz_sws1_log2_t_corUT <- melt(mz_sws1_log2_t_corUT) # melt for plotting

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R4eA.Mz_sws1_6tissue_pearson-correlation.tiff', units="in", width=12, height=8, res=300)
ggplot(data = melted_mz_sws1_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 7, hjust = 1, face = c(rep("plain", 30),"bold",rep("plain", 28)),color = c(rep("black", 30),"red",rep("black", 28))), # highlight sws1 - need to count how many pos. along it is 
        axis.text.y = element_text(size=7, face = c(rep("plain", 30),"bold",rep("plain", 28)),color = c(rep("black", 30),"red",rep("black", 28))), # highlight sws1 - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "M. zebra sws1 TF-TG network nodes") + labs (y = "M. zebra sws1 TF-TG network nodes") +
  coord_fixed()
dev.off()

# Nb sws1
nb_sws1_log2_t_cor <- reorder_cormat(nb_sws1_log2_t_cor) # reorder the correlation matrix to identify patterns
nb_sws1_log2_t_corUT <- get_upper_tri(nb_sws1_log2_t_cor) # take the upper triangle only since rest is redundant
melted_nb_sws1_log2_t_corUT <- melt(nb_sws1_log2_t_corUT) # melt for plotting

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R4eB.Nb_sws1_6tissue_pearson-correlation.tiff', units="in", width=12, height=8, res=300)
ggplot(data = melted_nb_sws1_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +2
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10.5, hjust = 1, face = c(rep("plain", 12),"bold",rep("plain", 11),"bold","plain"),color = c(rep("black", 12),"#009E73",rep("black", 11),"red","black")), # highlight sws1 - need to count how many pos. along it is 
        axis.text.y = element_text(size=10.5, face = c(rep("plain", 12),"bold",rep("plain", 11),"bold","plain"),color = c(rep("black", 12),"#009E73",rep("black", 11),"red","black")), # highlight sws1 - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "N. brichardi sws1 TF-TG network nodes") + labs (y = "N. brichardi sws1 TF-TG network nodes") +
  coord_fixed()
dev.off()

# Mz rho
mz_rho_log2_t_cor <- reorder_cormat(mz_rho_log2_t_cor) # reorder the correlation matrix to identify patterns
mz_rho_log2_t_corUT <- get_upper_tri(mz_rho_log2_t_cor) # take the upper triangle only since rest is redundant
melted_mz_rho_log2_t_corUT <- melt(mz_rho_log2_t_corUT) # melt for plotting

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R4eC.Mz_rho_6tissue_pearson-correlation.tiff', units="in", width=12, height=8, res=300)
ggplot(data = melted_mz_rho_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1, face = c(rep("plain", 2),"bold",rep("plain", 9),"bold",rep("plain",26)),color = c(rep("black", 2),"red",rep("black", 9),"#009E73",rep("black",26))), # highlight rho - need to count how many pos. along it is 
        axis.text.y = element_text(size=9, face = c(rep("plain", 2),"bold",rep("plain", 9),"bold",rep("plain",26)),color = c(rep("black", 2),"red",rep("black", 9),"#009E73",rep("black",26))), # highlight rho - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "M. zebra rho TF-TG network nodes") + labs (y = "M. zebra rho TF-TG network nodes") +
  coord_fixed()
dev.off()

# Ab rho
ab_rho_log2_t_cor <- reorder_cormat(ab_rho_log2_t_cor) # reorder the correlation matrix to identify patterns
ab_rho_log2_t_corUT <- get_upper_tri(ab_rho_log2_t_cor) # take the upper triangle only since rest is redundant
melted_ab_rho_log2_t_corUT <- melt(ab_rho_log2_t_corUT) # melt for plotting

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R4eD.Ab_rho_6tissue_pearson-correlation.tiff', units="in", width=12, height=8, res=300)
ggplot(data = melted_ab_rho_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1, face = c("plain","bold",rep("plain", 5),"bold",rep("plain",18),"bold",rep("plain",4)),color = c("black","#009E73",rep("black", 5),"red",rep("black",18),"#009E73",rep("black",4))), # highlight rho - need to count how many pos. along it is 
        axis.text.y = element_text(size=9, face = c("plain","bold",rep("plain", 5),"bold",rep("plain",18),"bold",rep("plain",4)),color = c("black","#009E73",rep("black", 5),"red",rep("black",18),"#009E73",rep("black",4))), # highlight rho - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "A. burtoni rho TF-TG network nodes") + labs (y = "A. burtoni rho TF-TG network nodes") +
  coord_fixed()
dev.off()

# On rho
on_rho_log2_t_cor <- reorder_cormat(on_rho_log2_t_cor) # reorder the correlation matrix to identify patterns
on_rho_log2_t_corUT <- get_upper_tri(on_rho_log2_t_cor) # take the upper triangle only since rest is redundant
melted_on_rho_log2_t_corUT <- melt(on_rho_log2_t_corUT) # melt for plotting

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R4eE.On_rho_6tissue_pearson-correlation.tiff', units="in", width=12, height=8, res=300)
ggplot(data = melted_on_rho_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lon", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1, face = c(rep("plain", 26),"bold","plain","bold","plain","bold"),color = c(rep("black", 26),"red","black","#009E73","black","#009E73")), # highlight rho - need to count how many pos. along it is 
        axis.text.y = element_text(size=10, face = c(rep("plain", 26),"bold","plain","bold","plain","bold"),color = c(rep("black", 26),"red","black","#009E73","black","#009E73")), # highlight rho - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "O. niloticus rho TF-TG network nodes") + labs (y = "O. niloticus rho TF-TG network nodes") +
  coord_fixed()
dev.off()


### Create one sheet with all plots combined and one legend

p1 <- ggplot(data = melted_mz_sws1_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 7, hjust = 1, face = c(rep("plain", 30),"bold",rep("plain", 28)),color = c(rep("black", 30),"red",rep("black", 28))), # highlight sws1 - need to count how many pos. along it is 
        axis.text.y = element_text(size=7, face = c(rep("plain", 30),"bold",rep("plain", 28)),color = c(rep("black", 30),"red",rep("black", 28))), # highlight sws1 - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "M. zebra sws1 TF-TG network nodes") + labs (y = "M. zebra sws1 TF-TG network nodes") +
  coord_fixed()

p2 <- ggplot(data = melted_nb_sws1_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10.5, hjust = 1, face = c(rep("plain", 12),"bold",rep("plain", 11),"bold","plain"),color = c(rep("black", 12),"#009E73",rep("black", 11),"red","black")), # highlight sws1 - need to count how many pos. along it is 
        axis.text.y = element_text(size=10.5, face = c(rep("plain", 12),"bold",rep("plain", 11),"bold","plain"),color = c(rep("black", 12),"#009E73",rep("black", 11),"red","black")), # highlight sws1 - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "N. brichardi sws1 TF-TG network nodes") + labs (y = "N. brichardi sws1 TF-TG network nodes") +
  coord_fixed()

p3 <- ggplot(data = melted_mz_rho_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1, face = c(rep("plain", 2),"bold",rep("plain", 9),"bold",rep("plain",26)),color = c(rep("black", 2),"red",rep("black", 9),"#009E73",rep("black",26))), # highlight rho - need to count how many pos. along it is 
        axis.text.y = element_text(size=9, face = c(rep("plain", 2),"bold",rep("plain", 9),"bold",rep("plain",26)),color = c(rep("black", 2),"red",rep("black", 9),"#009E73",rep("black",26))), # highlight rho - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "M. zebra rho TF-TG network nodes") + labs (y = "M. zebra rho TF-TG network nodes") +
  coord_fixed()

p4 <- ggplot(data = melted_ab_rho_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1, face = c("plain","bold",rep("plain", 5),"bold",rep("plain",18),"bold",rep("plain",4)),color = c("black","#009E73",rep("black", 5),"red",rep("black",18),"#009E73",rep("black",4))), # highlight rho - need to count how many pos. along it is 
        axis.text.y = element_text(size=9, face = c("plain","bold",rep("plain", 5),"bold",rep("plain",18),"bold",rep("plain",4)),color = c("black","#009E73",rep("black", 5),"red",rep("black",18),"#009E73",rep("black",4))), # highlight rho - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "A. burtoni rho TF-TG network nodes") + labs (y = "A. burtoni rho TF-TG network nodes") +
  coord_fixed()

p5 <- ggplot(data = melted_on_rho_log2_t_corUT, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lon", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1, face = c(rep("plain", 26),"bold","plain","bold","plain","bold"),color = c(rep("black", 26),"red","black","#009E73","black","#009E73")), # highlight rho - need to count how many pos. along it is 
        axis.text.y = element_text(size=10, face = c(rep("plain", 26),"bold","plain","bold","plain","bold"),color = c(rep("black", 26),"red","black","#009E73","black","#009E73")), # highlight rho - need to count how many pos. along it is
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 16),
        legend.key.size = unit(1,"cm")) +
  labs (x = "O. niloticus rho TF-TG network nodes") + labs (y = "O. niloticus rho TF-TG network nodes") +
  coord_fixed()
  
# extract legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(p1)

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R4e.sws1_rho_6tissue_pearson-correlation.tiff', units="in", width=24, height=16, res=300)
grid.arrange(p1 + theme(legend.position="none"),p2 + theme(legend.position="none"),p3 + theme(legend.position="none"),p4 + theme(legend.position="none"),p5 + theme(legend.position="none"),mylegend,nrow=2)
dev.off()



