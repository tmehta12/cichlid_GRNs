setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/")

library(reshape2)
library(ggplot2)
library(grid)

############################################################################################
# Promoter evolutionary rates taken by using O. niloticus 1st exon (+/- 100bp) as reference to BLAT against other genomes
# Then using OGIDs to ensure coherent orthogrouping with BLAT hits to take as start of 1st exon, using the original Brawand annotations for gene end and those that were not coherent with OGIDs
# Then pulled out the cichlid promoter sequences based on new cds annotation, aligned with MAFFT and then ran PAML on these aligned sequences
# Fourfold sites extracted in normal manner, using annotation, aligned and then ran PAML
############################################################################################

#### ##### ##### ##### ##### #### ##### ##### ##### #####

# Fig. 2aA - Evolutionary rate in promoter and fourfold degenerate sites at each ancestral node

MzPn_evolrate = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc1_4fold_Prom-BLAT.evolrate.txt2")
colnames(MzPn_evolrate) <- c('ortho', 'Fourfold', 'Promoter')
MzPn_evolrate <- cbind(MzPn_evolrate, Species = 'Anc1 - M. zebra/P. nyererei')
MzPn_evolrate_l <- melt(MzPn_evolrate,id.vars='Species', measure.vars=c('Fourfold','Promoter'))
wilcox.test(value ~ variable, data=MzPn_evolrate_l, paired=FALSE) # W = 7225000, p-value < 2.2e-16 - significant difference between promoter and fourfold

MzPnAb_evolrate = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc2_4fold_Prom-BLAT.evolrate.txt2")
colnames(MzPnAb_evolrate) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAb_evolrate <- cbind(MzPnAb_evolrate, Species = 'Anc2 - M. zebra/P. nyererei/A. burtoni')
MzPnAb_evolrate_l <- melt(MzPnAb_evolrate,id.vars='Species', measure.vars=c('Fourfold','Promoter'))
wilcox.test(value ~ variable, data=MzPnAb_evolrate_l, paired=FALSE) # W = 8862400, p-value < 2.2e-16 - significant difference between promoter and fourfold

MzPnAbNb_evolrate = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc3_4fold_Prom-BLAT.evolrate.txt2")
colnames(MzPnAbNb_evolrate) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAbNb_evolrate <- cbind(MzPnAbNb_evolrate, Species = 'Anc3 - M. zebra/P. nyererei/A. burtoni/N. brichardi')
MzPnAbNb_evolrate_l <- melt(MzPnAbNb_evolrate,id.vars='Species', measure.vars=c('Fourfold','Promoter'))
wilcox.test(value ~ variable, data=MzPnAbNb_evolrate_l, paired=FALSE) # W = 9276600, p-value < 2.2e-16 - significant difference between promoter and fourfold

# collate all
anc_evolrate <- rbind(MzPn_evolrate_l,MzPnAb_evolrate_l,MzPnAbNb_evolrate_l)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/FigR2aA_anc_evolrateplot.tiff', units="in", width=10, height=10, res=300)
ggplot(anc_evolrate, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species) + labs(title="Evolutionary rate at fourfold degenerate sites and promoter regions in the ancestral nodes of the five cichlids") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
dev.off()

# Fig. 2aB
# Plot the correlation between fourfold and promoter evolutionary rate to look at the orthgroups at each extreme
anc_evolrate2 <- rbind(MzPn_evolrate,MzPnAb_evolrate,MzPnAbNb_evolrate)

# take the 30 (if possible within the log2(promoter)>-5 (0.0313) and log2(fourfold)<-10(0.001) boundaries) topmost low fourfold and highest promoter evolutionary rate orthogroups
MzPn_evolrate$foldiff <- MzPn_evolrate$Promoter/MzPn_evolrate$Fourfold
MzPn_evolrate_peaks <- MzPn_evolrate[with(MzPn_evolrate, order(-foldiff)),][1:12,]
# MzPn_evolrate_peaks <- MzPn_evolrate[with(MzPn_evolrate, order(Fourfold, -Promoter)),][1:12,] # only 12 in boundary
MzPnAb_evolrate$foldiff <- MzPnAb_evolrate$Promoter/MzPnAb_evolrate$Fourfold
MzPnAb_evolrate_peaks <- MzPnAb_evolrate[with(MzPnAb_evolrate, order(-foldiff)),][1:30,]
# MzPnAb_evolrate_peaks <- MzPnAb_evolrate[with(MzPnAb_evolrate, order(Fourfold, -Promoter)),][1:10,]
MzPnAbNb_evolrate$foldiff <- MzPnAbNb_evolrate$Promoter/MzPnAbNb_evolrate$Fourfold
MzPnAbNb_evolrate_peaks <- MzPnAbNb_evolrate[with(MzPnAbNb_evolrate, order(-foldiff)),][c(1:30,167),] # also add foxd2 (TF) which is further down
# MzPnAbNb_evolrate_peaks <- MzPnAbNb_evolrate[with(MzPnAbNb_evolrate, order(Fourfold, -Promoter)),][1:10,]

# Get total numbers that are outliers
sum(MzPn_evolrate$Promoter > 0.0313 & MzPn_evolrate$Fourfold < 0.001, na.rm=TRUE) #Anc1 - 12
sum(MzPnAb_evolrate$Promoter > 0.0313 & MzPnAb_evolrate$Fourfold < 0.001, na.rm=TRUE) #Anc2 - 45
sum(MzPnAbNb_evolrate$Promoter > 0.0313 & MzPnAbNb_evolrate$Fourfold < 0.001, na.rm=TRUE) #Anc3 - 351

anc_evolrate_peaks <- rbind(MzPn_evolrate_peaks,MzPnAb_evolrate_peaks,MzPnAbNb_evolrate_peaks)
# match the orthogroups to gene name
OGIDs <- read.delim(file = "/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5", header = F, row.names = 1, stringsAsFactors=FALSE) # Load in the OGIDs file
anc_evolrate_peaks2 <- (OGIDs[match(anc_evolrate_peaks[,1],row.names(OGIDs)),(c(9,11,14))]) 
# manually amend some gene names according to other columns or BLAST
anc_evolrate_peaks2$V10[4] <- "ddb_G0267958"
anc_evolrate_peaks2$V10[7] <- "nrxn1"
anc_evolrate_peaks2$V10[8] <- "sp3b"
anc_evolrate_peaks2$V10[9] <- "ebf1"
anc_evolrate_peaks2$V10[10] <-"adamtsl5"
anc_evolrate_peaks2$V10[21] <- "csl2"
anc_evolrate_peaks2$V10[25] <- "hnrnpabb"
anc_evolrate_peaks2$V10[30] <- "pi4kb"
anc_evolrate_peaks2$V10[32] <- "trim35"
anc_evolrate_peaks2$V10[38] <- "fabp2"
anc_evolrate_peaks2$V10[39] <- "mpv17l"
anc_evolrate_peaks2$V10[40] <- "stbd1"
anc_evolrate_peaks2$V10[41] <- "snx6"
anc_evolrate_peaks2$V10[51] <- "melk"
anc_evolrate_peaks2$V10[56] <- "ptprn2"
anc_evolrate_peaks2$V10[57] <- "mad2l2"
anc_evolrate_peaks2$V10[58] <- "wdr77"
anc_evolrate_peaks2$V10[60] <- "gorasp1"
anc_evolrate_peaks2$V10[62] <- "clip1"
anc_evolrate_peaks2$V10[64] <- "lrp12"
anc_evolrate_peaks2$V10[69] <- "lpar1"

# cbind V10 column
anc_evolrate_peaks3 <- cbind(anc_evolrate_peaks,anc_evolrate_peaks2[,1])
colnames(anc_evolrate_peaks3) <- c("ortho","Fourfold","Promoter","Species","foldiff","gene")

# add the anc_evolrateplot_peaks labels to the plot and offset their labelling
anc_evolrateplot_correlation <- ggplot(anc_evolrate2, aes(x = log2(Fourfold), y = log2(Promoter))) + geom_point() +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species) + labs(title="Promoter vs fourfold evolutionary rate correlation at ancestral nodes") + 
  labs (x = "log2 (Fourfold evolutionary rate)") + labs (y = "log2 (Promoter evolutionary rate)") + theme(legend.position="none") + geom_hline(yintercept = -5,linetype="dashed",color="grey") + geom_vline(xintercept = -10,linetype="dashed",color="grey")

library(ggrepel)
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/FigR2aB_anc_evolrateplot_correlation.tiff', units="in", width=15, height=15, res=300)
anc_evolrateplot_correlation + geom_text_repel(data = anc_evolrate_peaks3, 
                                               mapping = aes(label=gene), 
                                               size = 3,
                                               fontface = 'bold', 
                                               color = 'dark grey',
                                               box.padding = unit(0.5, "lines"),
                                               point.padding = unit(0.5, "lines"))
dev.off()

# subset the genes that fall within log2(promoter)>-5 (0.0313) and log2(fourfold)<-10(0.001) in each node
anc_evolrate2_Anc1peak <- subset(anc_evolrate2, Promoter > 0.0313 & Fourfold < 0.001 & Species=="Anc1 - M. zebra/P. nyererei") # 12 genes
anc_evolrate2_Anc2peak <- subset(anc_evolrate2, Promoter > 0.0313 & Fourfold < 0.001 & Species=="Anc2 - M. zebra/P. nyererei/A. burtoni") # 45 genes
anc_evolrate2_Anc3peak <- subset(anc_evolrate2, Promoter > 0.0313 & Fourfold < 0.001 & Species=="Anc3 - M. zebra/P. nyererei/A. burtoni/N. brichardi") # 351 genes
anc_evolrate2_Anc123peak <- rbind(anc_evolrate2_Anc1peak,anc_evolrate2_Anc2peak,anc_evolrate2_Anc3peak)
anc_evolrate2_Anc123peak2 <- (OGIDs[match(anc_evolrate2_Anc123peak[,1],row.names(OGIDs)),(c(9,11,14))]) 
anc_evolrate2_Anc123peak3 <- cbind(anc_evolrate2_Anc123peak,anc_evolrate2_Anc123peak2) # match and bind to gene IDs

#### ##### ##### ##### ##### #### ##### ##### ##### #####

# Fig. 2aC - Evolutionary rate in promoter and fourfold degenerate sites at each branch (+ switched and non-switched at each branch, further down)

Mz_evolrate = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Mz_4fold_Prom-BLAT.evolrate.txt2")
colnames(Mz_evolrate) <- c('ortho', 'Fourfold', 'Promoter')
Mz_evolrate <- cbind(Mz_evolrate, Species = 'M. zebra')
Mz_evolrate_l <- melt(Mz_evolrate,id.vars='Species', measure.vars=c('Fourfold','Promoter'))
wilcox.test(value ~ variable, data=Mz_evolrate_l, paired=FALSE) # W = 9436400, p-value < 2.2e-16 - significant difference between promoter and fourfold

Pn_evolrate = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Pn_4fold_Prom-BLAT.evolrate.txt2")
colnames(Pn_evolrate) <- c('ortho', 'Fourfold', 'Promoter')
Pn_evolrate <- cbind(Pn_evolrate, Species = 'P. nyererei')
Pn_evolrate_l <- melt(Pn_evolrate,id.vars='Species', measure.vars=c('Fourfold','Promoter'))
wilcox.test(value ~ variable, data=Pn_evolrate_l, paired=FALSE) # W = 9417100, p-value < 2.2e-16 - significant difference between promoter and fourfold

Ab_evolrate = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Ab_4fold_Prom-BLAT.evolrate.txt2")
colnames(Ab_evolrate) <- c('ortho', 'Fourfold', 'Promoter')
Ab_evolrate <- cbind(Ab_evolrate, Species = 'A. burtoni')
Ab_evolrate_l <- melt(Ab_evolrate,id.vars='Species', measure.vars=c('Fourfold','Promoter'))
wilcox.test(value ~ variable, data=Ab_evolrate_l, paired=FALSE) # W = 10050000, p-value = 7.838e-07 - significant difference between promoter and fourfold

Nb_evolrate = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Nb_4fold_Prom-BLAT.evolrate.txt2")
colnames(Nb_evolrate) <- c('ortho', 'Fourfold', 'Promoter')
Nb_evolrate <- cbind(Nb_evolrate, Species = 'N. brichardi')
Nb_evolrate_l <- melt(Nb_evolrate,id.vars='Species', measure.vars=c('Fourfold','Promoter'))
wilcox.test(value ~ variable, data=Nb_evolrate_l, paired=FALSE) # W = 9999300, p-value = 1.047e-07 - significant difference between promoter and fourfold

On_evolrate = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/On_4fold_Prom-BLAT.evolrate.txt2")
colnames(On_evolrate) <- c('ortho', 'Fourfold', 'Promoter')
On_evolrate <- cbind(On_evolrate, Species = 'O. niloticus')
On_evolrate_l <- melt(On_evolrate,id.vars='Species', measure.vars=c('Fourfold','Promoter'))
wilcox.test(value ~ variable, data=On_evolrate_l, paired=FALSE) # W = 9691300, p-value = 1.152e-14 - significant difference between promoter and fourfold

# collate all
branch_evolrate <- rbind(Mz_evolrate_l,Pn_evolrate_l,Ab_evolrate_l,Nb_evolrate_l,On_evolrate_l)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/FigR2aC_branch_evolrateplot.tiff', units="in", width=10, height=10, res=300)
ggplot(branch_evolrate, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species) + labs(title="Evolutionary rate at fourfold degenerate sites and promoter regions in the five cichlids") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
dev.off()

# Fig. 2aD
# Plot the correlation between fourfold and promoter evolutionary rate in each species to look at the orthgroups at each extreme
branch_evolrate2 <- rbind(Mz_evolrate,Pn_evolrate,Ab_evolrate,Nb_evolrate,On_evolrate)

# take the 30 (if possible within the log2(promoter)>-5 (0.0313) and log2(fourfold)<-10(0.001) boundaries) topmost low fourfold and highest promoter evolutionary rate orthogroups
Mz_evolrate$foldiff <- Mz_evolrate$Promoter/Mz_evolrate$Fourfold
Mz_evolrate_peaks <- Mz_evolrate[with(Mz_evolrate, order(-foldiff)),][1:30,]
# Mz_evolrate_peaks <- Mz_evolrate[with(Mz_evolrate, order(Fourfold, -Promoter)),][1:10,]
Pn_evolrate$foldiff <- Pn_evolrate$Promoter/Pn_evolrate$Fourfold
Pn_evolrate_peaks <- Pn_evolrate[with(Pn_evolrate, order(-foldiff)),][1:29,]
# Pn_evolrate_peaks <- Pn_evolrate[with(Pn_evolrate, order(Fourfold, -Promoter)),][1:10,]
Ab_evolrate$foldiff <- Ab_evolrate$Promoter/Ab_evolrate$Fourfold
Ab_evolrate_peaks <- Ab_evolrate[with(Ab_evolrate, order(-foldiff)),][1:30,]
# Ab_evolrate_peaks <- Ab_evolrate[with(Ab_evolrate, order(Fourfold, -Promoter)),][1:10,]
Nb_evolrate$foldiff <- Nb_evolrate$Promoter/Nb_evolrate$Fourfold
Nb_evolrate_peaks <- Nb_evolrate[with(Nb_evolrate, order(-foldiff)),][1:30,]
# Nb_evolrate_peaks <- Nb_evolrate[with(Nb_evolrate, order(Fourfold, -Promoter)),][1:10,]
On_evolrate$foldiff <- On_evolrate$Promoter/On_evolrate$Fourfold
On_evolrate_peaks <- On_evolrate[with(On_evolrate, order(-foldiff)),][1:30,]
# On_evolrate_peaks <- On_evolrate[with(On_evolrate, order(Fourfold, -Promoter)),][1:10,]

# Get total numbers that are outliers
sum(Mz_evolrate$Promoter > 0.0313 & Mz_evolrate$Fourfold < 0.001, na.rm=TRUE) # Mz - 56
sum(Pn_evolrate$Promoter > 0.0313 & Pn_evolrate$Fourfold < 0.001, na.rm=TRUE) # Pn - 29
sum(Ab_evolrate$Promoter > 0.0313 & Ab_evolrate$Fourfold < 0.001, na.rm=TRUE) # Ab - 32
sum(Nb_evolrate$Promoter > 0.0313 & Nb_evolrate$Fourfold < 0.001, na.rm=TRUE) # Nb - 49
sum(On_evolrate$Promoter > 0.0313 & On_evolrate$Fourfold < 0.001, na.rm=TRUE) # On - 352

branch_evolrate_peaks <- rbind(Mz_evolrate_peaks,Pn_evolrate_peaks,Ab_evolrate_peaks,Nb_evolrate_peaks,On_evolrate_peaks)
# match the orthogroups to gene name
branch_evolrate_peaks2 <- (OGIDs[match(branch_evolrate_peaks[,1],row.names(OGIDs)),(c(9,11,14))])
# manually amend some gene names according to other columns or BLAST
branch_evolrate_peaks2$gene <- ifelse(branch_evolrate_peaks2$V10=="NULL",yes = branch_evolrate_peaks2$V12, no = branch_evolrate_peaks2$V10) # create an extra column that collates
branch_evolrate_peaks2$gene[19] <- "unc.prot.LOC100699166"
branch_evolrate_peaks2$gene[20] <- "mypop"
branch_evolrate_peaks2$gene[21] <- "gorasp1"
branch_evolrate_peaks2$gene[35] <- "ppt2" 
branch_evolrate_peaks2$gene[37] <- "sapcd1"
branch_evolrate_peaks2$gene[42] <- "cartpt"
branch_evolrate_peaks2$gene[50] <- "prss2" 
branch_evolrate_peaks2$gene[52] <-"ttf2" 
branch_evolrate_peaks2$gene[55] <- "cx32.2"
branch_evolrate_peaks2$gene[71] <- "nrxn1"
branch_evolrate_peaks2$gene[80] <- "stx2"
branch_evolrate_peaks2$gene[98] <- "scn2a" 
branch_evolrate_peaks2$gene[100] <- "nsg2" 
branch_evolrate_peaks2$gene[102] <- "gpr45"
branch_evolrate_peaks2$gene[106] <- "hla-dp1"
branch_evolrate_peaks2$gene[119] <- "kcng"
branch_evolrate_peaks2$gene[126] <- "akap2" 
branch_evolrate_peaks2$gene[127] <- "ccdc85a" 
branch_evolrate_peaks2$gene[130] <- "rpl22"
branch_evolrate_peaks2$gene[132] <- "pwwp2a"
branch_evolrate_peaks2$gene[139] <- "fndc3a"
branch_evolrate_peaks2$gene[141] <- "unc.prot.LOC102079744" 
branch_evolrate_peaks2$gene[146] <- "sncaip"

# cbind gene column
branch_evolrate_peaks3 <- cbind(branch_evolrate_peaks,branch_evolrate_peaks2[,4])
colnames(branch_evolrate_peaks3) <- c("ortho","Fourfold","Promoter","Species","foldiff","gene")

# add the branch_evolrateplot_peaks labels to the plot and offset their labelling
branch_evolrateplot_correlation <- ggplot(branch_evolrate2, aes(x = log2(Fourfold), y = log2(Promoter))) + geom_point() +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species) + labs(title="Promoter vs fourfold evolutionary rate correlation at branches") + 
  labs (x = "log2 (Fourfold evolutionary rate)") + labs (y = "log2 (Promoter evolutionary rate)") + theme(legend.position="none") + geom_hline(yintercept = -5,linetype="dashed",color="grey") + geom_vline(xintercept = -10,linetype="dashed",color="grey")

library(ggrepel)
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/FigR2aD_branch_evolrateplot_correlation.tiff', units="in", width=15, height=15, res=300)
branch_evolrateplot_correlation + geom_text_repel(data = branch_evolrate_peaks3, 
                                                  mapping = aes(label=gene), 
                                                  size = 3,
                                                  fontface = 'bold', 
                                                  color = 'dark grey',
                                                  box.padding = unit(0.5, "lines"),
                                                  point.padding = unit(0.5, "lines"))
dev.off()

# subset the genes that fall within log2(promoter)>-5 (0.0313) and log2(fourfold)<-10(0.001) in each node
branch_evolrate2_Mzpeak <- subset(branch_evolrate2, Promoter > 0.0313 & Fourfold < 0.001 & Species=="M. zebra") # 56 genes
branch_evolrate2_Pnpeak <- subset(branch_evolrate2, Promoter > 0.0313 & Fourfold < 0.001 & Species=="P. nyererei") # 29 genes
branch_evolrate2_Abpeak <- subset(branch_evolrate2, Promoter > 0.0313 & Fourfold < 0.001 & Species=="A. burtoni") # 32 genes
branch_evolrate2_Nbpeak <- subset(branch_evolrate2, Promoter > 0.0313 & Fourfold < 0.001 & Species=="N. brichardi") # 49 genes
branch_evolrate2_Onpeak <- subset(branch_evolrate2, Promoter > 0.0313 & Fourfold < 0.001 & Species=="O. niloticus") # 352 genes

#### ##### ##### ##### ##### #### ##### ##### ##### #####

# Fig. 2aE - Evolutionary rate in promoters and fourfold degenerate sites of co-expressed 1:1 orthologous cichlid genes at each ancestral node

# Pull out all the files you require from the directory
base_dir_evol <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/"
module_Ancevolrate_files <- list.files(path=base_dir_evol, pattern = "\\CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2$")
for (i in 1:length(module_Ancevolrate_files)){
  tmp = read.delim(file = paste0(base_dir_evol,module_Ancevolrate_files[i]), header = F)
  assign(module_Ancevolrate_files[i], tmp)
}


# Create a loop to run on all module files (run plots separately)
module_Ancdflist <- list(`0-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`1-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`2-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`3-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`4-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`5-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`6-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`7-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`8-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`9-CollatedAnc_4fold_Prom-BLAT.evolrate.txt2_modules.2`)

for(i in 1:length(module_Ancdflist))
{
  colnames(module_Ancdflist[[i]]) <- c('ortho', 'Fourfold', 'Promoter', 'GeneID', 'Module', 'Species')
  module_Ancdflist[[i]] = melt(module_Ancdflist[[i]],id.vars='Species', measure.vars=c('Fourfold','Promoter'))
  module_Ancdflist[[i]]$Species_f = factor(module_Ancdflist[[i]]$Species, levels=c('Anc1 - M. zebra/P. nyererei','Anc2 - M. zebra/P. nyererei/A. burtoni','Anc3 - M. zebra/P. nyererei/A. burtoni/N. brichardi')) # fix the order of Species according to phylo for facet_grid
}

# create plots for each
Module0_Ancevolrateplot <- ggplot(module_Ancdflist[[1]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 0 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module1_Ancevolrateplot <- ggplot(module_Ancdflist[[2]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 1 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module2_Ancevolrateplot <- ggplot(module_Ancdflist[[3]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 2 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module3_Ancevolrateplot <- ggplot(module_Ancdflist[[4]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 3 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module4_Ancevolrateplot <- ggplot(module_Ancdflist[[5]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 4 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module5_Ancevolrateplot <- ggplot(module_Ancdflist[[6]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 5 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module6_Ancevolrateplot <- ggplot(module_Ancdflist[[7]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 6 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module7_Ancevolrateplot <- ggplot(module_Ancdflist[[8]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 7 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module8_Ancevolrateplot <- ggplot(module_Ancdflist[[9]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 8 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module9_Ancevolrateplot <- ggplot(module_Ancdflist[[10]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 9 ancestral nodes evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/FigR2aE_module_anc_evolrateplot.tiff', units="in", width=20, height=20, res=300)
multiplot(Module0_Ancevolrateplot,Module1_Ancevolrateplot,Module2_Ancevolrateplot,Module3_Ancevolrateplot,Module4_Ancevolrateplot,Module5_Ancevolrateplot,Module6_Ancevolrateplot,Module7_Ancevolrateplot,Module8_Ancevolrateplot,Module9_Ancevolrateplot, cols = 2) # quick multiplot
dev.off()


#### ##### ##### ##### #####

# Fig. 2aF - Evolutionary rate in promoter and fourfold degenerate sites in each module for each branch of the five species
# pre-prepared files before loading here in NetworkReconstruction_v5.sh script

# Pull out all the files you require from the directory
base_dir_evol <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/"
module_evolrate_files <- list.files(path=base_dir_evol, pattern = "\\Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2$")
for (i in 1:length(module_evolrate_files)){
  tmp = read.delim(file = paste0(base_dir_evol,module_evolrate_files[i]), header = F)
  assign(module_evolrate_files[i], tmp)
}

# Module 0 - this is for running individually
# colnames(`0-Collated_evolRate.out_modules.2`) <- c('ortho', 'Fourfold', 'Promoter', 'GeneID', 'Module', 'Species')
# `0-Collated_evolRate.out_modules.2_l` <- melt(`0-Collated_evolRate.out_modules.2`,id.vars='Species', measure.vars=c('Fourfold','Promoter'))
# `0-Collated_evolRate.out_modules.2_l`$Species_f = factor(`0-Collated_evolRate.out_modules.2_l`$Species, levels=c('M. zebra','P. nyererei','A. burtoni','N. brichardi','O. niloticus')) # fix the order of Species according to phylo for facet_grid
# Module0_evolrateplot <- ggplot(`0-Collated_evolRate.out_modules.2_l`, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() +
#   theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 0 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
#   labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none") 

# Create a loop to run on all module files (run plots separately)
module_dflist <- list(`0-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`1-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`2-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`3-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`4-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`5-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`6-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`7-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`8-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`,`9-Collated_4fold_Prom-BLAT.evolrate.txt2_modules.2`)

for(i in 1:length(module_dflist))
{
  colnames(module_dflist[[i]]) <- c('ortho', 'Fourfold', 'Promoter', 'GeneID', 'Module', 'Species')
  module_dflist[[i]] = melt(module_dflist[[i]],id.vars='Species', measure.vars=c('Fourfold','Promoter'))
  module_dflist[[i]]$Species_f = factor(module_dflist[[i]]$Species, levels=c('M. zebra','P. nyererei','A. burtoni','N. brichardi','O. niloticus')) # fix the order of Species according to phylo for facet_grid
}

# create plots for each
Module0_evolrateplot <- ggplot(module_dflist[[1]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 0 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module1_evolrateplot <- ggplot(module_dflist[[2]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 1 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module2_evolrateplot <- ggplot(module_dflist[[3]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 2 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module3_evolrateplot <- ggplot(module_dflist[[4]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 3 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module4_evolrateplot <- ggplot(module_dflist[[5]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 4 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module5_evolrateplot <- ggplot(module_dflist[[6]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 5 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module6_evolrateplot <- ggplot(module_dflist[[7]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 6 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module7_evolrateplot <- ggplot(module_dflist[[8]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 7 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module8_evolrateplot <- ggplot(module_dflist[[9]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 8 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
Module9_evolrateplot <- ggplot(module_dflist[[10]], aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species_f) + labs(title="Module 9 evolutionary rate") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/FigR2aF_module_branch_evolrateplot.tiff', units="in", width=20, height=20, res=300)
multiplot(Module0_evolrateplot,Module1_evolrateplot,Module2_evolrateplot,Module3_evolrateplot,Module4_evolrateplot,Module5_evolrateplot,Module6_evolrateplot,Module7_evolrateplot,Module8_evolrateplot,Module9_evolrateplot, cols = 2) # quick multiplot
dev.off()

#### ##### ##### ##### ##### #### ##### ##### ##### #####

# Fig. 2aG - Evolutionary rate in promoter regions of switching and non-switched 1:1 orthologous cichlid genes at each ancestral node

### You should test the switches derived from Anc, found here - THESE ARE THE ONES TO USE:
# /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/prediction_mode_results/results/
# e.g. Anc3vs4_clusterassign_noOG-allassigned_StateChange.txt

MzPn_evolrate_sw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc1_4fold_Prom-BLAT.evolrate.txt2.Ancsw")
colnames(MzPn_evolrate_sw) <- c('ortho', 'Fourfold', 'Promoter')
MzPn_evolrate_sw <- cbind(MzPn_evolrate_sw, Species = 'Anc1 - M. zebra/P. nyererei', type = 'switch')
MzPn_evolrate_sw_l <- melt(MzPn_evolrate_sw)
wilcox.test(value ~ variable, data=MzPn_evolrate_sw_l, paired=FALSE) # W = 6679.5, p-value = 1.884e-07 - significant difference between promoter and fourfold

MzPn_evolrate_nsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc1_4fold_Prom-BLAT.evolrate.txt2.Ancnsw")
colnames(MzPn_evolrate_nsw) <- c('ortho', 'Fourfold', 'Promoter')
MzPn_evolrate_nsw <- cbind(MzPn_evolrate_nsw, Species = 'Anc1 - M. zebra/P. nyererei', type = 'noswitch')
MzPn_evolrate_nsw_l <- melt(MzPn_evolrate_nsw)
wilcox.test(value ~ variable, data=MzPn_evolrate_nsw_l, paired=FALSE) # W = 6792500, p-value < 2.2e-16 - significant difference between promoter and fourfold

MzPnAb_evolrate_sw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc2_4fold_Prom-BLAT.evolrate.txt2.Ancsw")
colnames(MzPnAb_evolrate_sw) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAb_evolrate_sw <- cbind(MzPnAb_evolrate_sw, Species = 'Anc2 - M. zebra/P. nyererei/A. burtoni', type = 'switch')
MzPnAb_evolrate_sw_l <- melt(MzPnAb_evolrate_sw)
wilcox.test(value ~ variable, data=MzPnAb_evolrate_sw_l, paired=FALSE) # W = 9090, p-value = 0.1039 - NO significant difference between promoter and fourfold

MzPnAb_evolrate_nsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc2_4fold_Prom-BLAT.evolrate.txt2.Ancnsw")
colnames(MzPnAb_evolrate_nsw) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAb_evolrate_nsw <- cbind(MzPnAb_evolrate_nsw, Species = 'Anc2 - M. zebra/P. nyererei/A. burtoni', type = 'noswitch')
MzPnAb_evolrate_nsw_l <- melt(MzPnAb_evolrate_nsw)
wilcox.test(value ~ variable, data=MzPnAb_evolrate_nsw_l, paired=FALSE) # W = 8303000, p-value < 2.2e-16 - significant difference between promoter and fourfold

MzPnAbNb_evolrate_sw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc3_4fold_Prom-BLAT.evolrate.txt2.Ancsw")
colnames(MzPnAbNb_evolrate_sw) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAbNb_evolrate_sw <- cbind(MzPnAbNb_evolrate_sw, Species = 'Anc3 - M. zebra/P. nyererei/A. burtoni/N. brichardi', type = 'switch')
MzPnAbNb_evolrate_sw_l <- melt(MzPnAbNb_evolrate_sw)
wilcox.test(value ~ variable, data=MzPnAbNb_evolrate_sw_l, paired=FALSE) # W = 44208, p-value = 0.0007895 - significant difference between promoter and fourfold

MzPnAbNb_evolrate_nsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc3_4fold_Prom-BLAT.evolrate.txt2.Ancnsw")
colnames(MzPnAbNb_evolrate_nsw) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAbNb_evolrate_nsw <- cbind(MzPnAbNb_evolrate_nsw, Species = 'Anc3 - M. zebra/P. nyererei/A. burtoni/N. brichardi', type = 'noswitch')
MzPnAbNb_evolrate_nsw_l <- melt(MzPnAbNb_evolrate_nsw)
wilcox.test(value ~ variable, data=MzPnAbNb_evolrate_nsw_l, paired=FALSE) # W = 8038700, p-value < 2.2e-16 - significant difference between promoter and fourfold

# collate all to plot together
anc_sw_nsw_evolrate <- rbind(MzPn_evolrate_sw_l,MzPn_evolrate_nsw_l,MzPnAb_evolrate_sw_l,MzPnAb_evolrate_nsw_l,MzPnAbNb_evolrate_sw_l,MzPnAbNb_evolrate_nsw_l)

# for plotting, remove fourfold sites
anc_sw_nsw_evolrate2 <- anc_sw_nsw_evolrate[ grep("Fourfold", anc_sw_nsw_evolrate$variable, invert = TRUE) , ] 

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/FigR2aG_anc_sw_nsw_promevolrateplot.tiff', units="in", width=10, height=10, res=300)
# ggplot(anc_sw_nsw_evolrate, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() +
#   theme_bw() + theme(panel.grid = element_blank()) + facet_grid(. ~Species + type) + labs(title="Evolutionary rate of switching and non-switching genes at fourfold degenerate sites and promoter regions in the ancestral nodes of the five cichlids") + guides(fill=guide_legend(title="Genomic region")) + 
#   labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
ggplot(anc_sw_nsw_evolrate2, aes(x = factor(type), y = log2(value), fill=factor(type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(. ~Species) + labs(title="Evolutionary rate of switching and non-switching genes at promoter regions in the ancestral nodes of the five cichlids") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
dev.off()

### Test - looking at evolutionary rate associated with switching/non-switching genes
## These are consistent switches using the species e.g. Mz and Pn same switch vs rest etc.
# Pn & Mz same switch vs Ab, Nb and On
# OGIDS.txt5-clusterassign-1to1_MzPnvsAbNbOn #89
# Pn, Mz & Ab same switch vs Nb and On
# OGIDS.txt5-clusterassign-1to1_MzPnAbvsNbOn #172
# Pn, Mz, Ab & Nb same switch vs On
# OGIDS.txt5-clusterassign-1to1_MzPnAbNbvsOn #983

MzPn_evolrate_Spsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc1_4fold_Prom-BLAT.evolrate.txt2.SpAncsw")
colnames(MzPn_evolrate_Spsw) <- c('ortho', 'Fourfold', 'Promoter')
MzPn_evolrate_Spsw <- cbind(MzPn_evolrate_Spsw, Species = 'Anc1 - M. zebra/P. nyererei', type = 'switch')
MzPn_evolrate_Spsw_l <- melt(MzPn_evolrate_Spsw)
wilcox.test(value ~ variable, data=MzPn_evolrate_Spsw_l, paired=FALSE) # W = 1126, p-value = 5.873e-05 - significant difference between promoter and fourfold

MzPn_evolrate_Spnsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc1_4fold_Prom-BLAT.evolrate.txt2.SpAncnsw")
colnames(MzPn_evolrate_Spnsw) <- c('ortho', 'Fourfold', 'Promoter')
MzPn_evolrate_Spnsw <- cbind(MzPn_evolrate_Spnsw, Species = 'Anc1 - M. zebra/P. nyererei', type = 'noswitch')
MzPn_evolrate_Spnsw_l <- melt(MzPn_evolrate_Spnsw)
wilcox.test(value ~ variable, data=MzPn_evolrate_Spnsw_l, paired=FALSE) # W = 7045700, p-value < 2.2e-16 - significant difference between promoter and fourfold

MzPnAb_evolrate_Spsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc2_4fold_Prom-BLAT.evolrate.txt2.SpAncsw")
colnames(MzPnAb_evolrate_Spsw) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAb_evolrate_Spsw <- cbind(MzPnAb_evolrate_Spsw, Species = 'Anc2 - M. zebra/P. nyererei/A. burtoni', type = 'switch')
MzPnAb_evolrate_Spsw_l <- melt(MzPnAb_evolrate_Spsw)
wilcox.test(value ~ variable, data=MzPnAb_evolrate_Spsw_l, paired=FALSE) # W = 4675, p-value = 0.0004842 - significant difference between promoter and fourfold

MzPnAb_evolrate_Spnsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc2_4fold_Prom-BLAT.evolrate.txt2.SpAncnsw")
colnames(MzPnAb_evolrate_Spnsw) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAb_evolrate_Spnsw <- cbind(MzPnAb_evolrate_Spnsw, Species = 'Anc2 - M. zebra/P. nyererei/A. burtoni', type = 'noswitch')
MzPnAb_evolrate_Spnsw_l <- melt(MzPnAb_evolrate_Spnsw)
wilcox.test(value ~ variable, data=MzPnAb_evolrate_Spnsw_l, paired=FALSE) # W = 8459600, p-value < 2.2e-16 - significant difference between promoter and fourfold

MzPnAbNb_evolrate_Spsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc3_4fold_Prom-BLAT.evolrate.txt2.SpAncsw")
colnames(MzPnAbNb_evolrate_Spsw) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAbNb_evolrate_Spsw <- cbind(MzPnAbNb_evolrate_Spsw, Species = 'Anc3 - M. zebra/P. nyererei/A. burtoni/N. brichardi', type = 'switch')
MzPnAbNb_evolrate_Spsw_l <- melt(MzPnAbNb_evolrate_Spsw)
wilcox.test(value ~ variable, data=MzPnAbNb_evolrate_Spsw_l, paired=FALSE) # W = 192930, p-value = 0.001611 - significant difference between promoter and fourfold

MzPnAbNb_evolrate_Spnsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Anc3_4fold_Prom-BLAT.evolrate.txt2.SpAncnsw")
colnames(MzPnAbNb_evolrate_Spnsw) <- c('ortho', 'Fourfold', 'Promoter')
MzPnAbNb_evolrate_Spnsw <- cbind(MzPnAbNb_evolrate_Spnsw, Species = 'Anc3 - M. zebra/P. nyererei/A. burtoni/N. brichardi', type = 'noswitch')
MzPnAbNb_evolrate_Spnsw_l <- melt(MzPnAbNb_evolrate_Spnsw)
wilcox.test(value ~ variable, data=MzPnAbNb_evolrate_Spnsw_l, paired=FALSE) # W = 6794400, p-value < 2.2e-16 - significant difference between promoter and fourfold

# collate all to plot together
anc_Spsw_Spnsw_evolrate <- rbind(MzPn_evolrate_Spsw_l,MzPn_evolrate_Spnsw_l,MzPnAb_evolrate_Spsw_l,MzPnAb_evolrate_Spnsw_l,MzPnAbNb_evolrate_Spsw_l,MzPnAbNb_evolrate_Spnsw_l)

# for plotting, remove fourfold sites
anc_Spsw_Spnsw_evolrate2 <- anc_Spsw_Spnsw_evolrate[ grep("Fourfold", anc_Spsw_Spnsw_evolrate$variable, invert = TRUE) , ] 

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/anc_Spsw_Spnsw_evolrateplot.tiff', units="in", width=10, height=10, res=300)
# ggplot(anc_Spsw_Spnsw_evolrate, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() +
#   theme_bw() + theme(panel.grid = element_blank()) + facet_grid(. ~Species + type) + labs(title="Evolutionary rate of switching and non-switching genes at fourfold degenerate sites and promoter regions in the ancestral nodes of the five cichlids") + guides(fill=guide_legend(title="Genomic region")) + 
#   labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
ggplot(anc_Spsw_Spnsw_evolrate2, aes(x = factor(type), y = log2(value), fill=factor(type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(. ~Species) + labs(title="Evolutionary rate of switching and non-switching genes at fourfold degenerate sites and promoter regions in the ancestral nodes (branch species same switch) of the five cichlids") + guides(fill=guide_legend(title="Genomic region")) + 
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
dev.off()

#### ##### ##### ##### ##### #### ##### ##### ##### #####

# Fig. 2aH - Evolutionary rate in promoter regions of switching and non-switched 1:1 orthologous cichlid genes at each branch

Mz_evolrate_sw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Mz_4fold_Prom-BLAT.evolrate.txt2.sw")
colnames(Mz_evolrate_sw) <- c('ortho', 'Fourfold', 'Promoter')
Mz_evolrate_sw <- cbind(Mz_evolrate_sw, Species = 'M. zebra', Type ='switch')
Mz_evolrate_sw_l <- melt(Mz_evolrate_sw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

Mz_evolrate_nsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Mz_4fold_Prom-BLAT.evolrate.txt2.nsw")
colnames(Mz_evolrate_nsw) <- c('ortho', 'Fourfold', 'Promoter')
Mz_evolrate_nsw <- cbind(Mz_evolrate_nsw, Species = 'M. zebra', Type ='noswitch')
Mz_evolrate_nsw_l <- melt(Mz_evolrate_nsw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

# Individual plot
# Mz_evolrate_sw_nsw <- rbind(Mz_evolrate_sw_l,Mz_evolrate_nsw_l)
# ggplot(Mz_evolrate_sw_nsw, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() +
#   theme_bw() + theme(panel.grid = element_blank()) + facet_grid(. ~Species + Type) + labs(title="Evolutionary rate of switching and non-switching genes at fourfold degenerate sites and promoter regions in M. zebra") + guides(fill=guide_legend(title="Genomic region")) + 
#   labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")

Pn_evolrate_sw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Pn_4fold_Prom-BLAT.evolrate.txt2.sw")
colnames(Pn_evolrate_sw) <- c('ortho', 'Fourfold', 'Promoter')
Pn_evolrate_sw <- cbind(Pn_evolrate_sw, Species = 'P. nyererei', Type ='switch')
Pn_evolrate_sw_l <- melt(Pn_evolrate_sw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

Pn_evolrate_nsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Pn_4fold_Prom-BLAT.evolrate.txt2.nsw")
colnames(Pn_evolrate_nsw) <- c('ortho', 'Fourfold', 'Promoter')
Pn_evolrate_nsw <- cbind(Pn_evolrate_nsw, Species = 'P. nyererei', Type ='noswitch')
Pn_evolrate_nsw_l <- melt(Pn_evolrate_nsw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

Ab_evolrate_sw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Ab_4fold_Prom-BLAT.evolrate.txt2.sw")
colnames(Ab_evolrate_sw) <- c('ortho', 'Fourfold', 'Promoter')
Ab_evolrate_sw <- cbind(Ab_evolrate_sw, Species = 'A. burtoni', Type ='switch')
Ab_evolrate_sw_l <- melt(Ab_evolrate_sw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

Ab_evolrate_nsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Ab_4fold_Prom-BLAT.evolrate.txt2.nsw")
colnames(Ab_evolrate_nsw) <- c('ortho', 'Fourfold', 'Promoter')
Ab_evolrate_nsw <- cbind(Ab_evolrate_nsw, Species = 'A. burtoni', Type ='noswitch')
Ab_evolrate_nsw_l <- melt(Ab_evolrate_nsw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

Nb_evolrate_sw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Nb_4fold_Prom-BLAT.evolrate.txt2.sw")
colnames(Nb_evolrate_sw) <- c('ortho', 'Fourfold', 'Promoter')
Nb_evolrate_sw <- cbind(Nb_evolrate_sw, Species = 'N. brichardi', Type ='switch')
Nb_evolrate_sw_l <- melt(Nb_evolrate_sw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

Nb_evolrate_nsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/Nb_4fold_Prom-BLAT.evolrate.txt2.nsw")
colnames(Nb_evolrate_nsw) <- c('ortho', 'Fourfold', 'Promoter')
Nb_evolrate_nsw <- cbind(Nb_evolrate_nsw, Species = 'N. brichardi', Type ='noswitch')
Nb_evolrate_nsw_l <- melt(Nb_evolrate_nsw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

On_evolrate_sw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/On_4fold_Prom-BLAT.evolrate.txt2.sw")
colnames(On_evolrate_sw) <- c('ortho', 'Fourfold', 'Promoter')
On_evolrate_sw <- cbind(On_evolrate_sw, Species = 'O. niloticus', Type ='switch')
On_evolrate_sw_l <- melt(On_evolrate_sw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

On_evolrate_nsw = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_220917/On_4fold_Prom-BLAT.evolrate.txt2.nsw")
colnames(On_evolrate_nsw) <- c('ortho', 'Fourfold', 'Promoter')
On_evolrate_nsw <- cbind(On_evolrate_nsw, Species = 'O. niloticus', Type ='noswitch')
On_evolrate_nsw_l <- melt(On_evolrate_nsw,id.vars=c('Species','Type'), measure.vars=c('Fourfold','Promoter'))

# Collate all to plot together
All_evolrate_sw_nsw_l <- rbind(Mz_evolrate_sw_l,Mz_evolrate_nsw_l,Pn_evolrate_sw_l,Pn_evolrate_nsw_l,Ab_evolrate_sw_l,Ab_evolrate_nsw_l,Nb_evolrate_sw_l,Nb_evolrate_nsw_l,On_evolrate_sw_l,On_evolrate_nsw_l)

# for plotting, remove fourfold sites
All_evolrate_sw_nsw_l2 <- All_evolrate_sw_nsw_l[ grep("Fourfold", All_evolrate_sw_nsw_l$variable, invert = TRUE) , ] 

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/Branch_evolrate_sw_nsw.tiff', units="in", width=25, height=15, res=300)
# ggplot(All_evolrate_sw_nsw_l, aes(x = factor(variable), y = log2(value), fill=factor(variable))) + geom_boxplot() +
#   theme_bw() + theme(panel.grid = element_blank()) + facet_grid(. ~Species + Type) + labs(title="Evolutionary rate of switching and non-switching genes at fourfold degenerate sites and promoter regions in five species") + guides(fill=guide_legend(title="Genomic region")) +
#   labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
ggplot(All_evolrate_sw_nsw_l2, aes(x = factor(Type), y = log2(value), fill=factor(Type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(. ~Species) + labs(title="Evolutionary rate of switching and non-switching genes at promoter regions in five species") + guides(fill=guide_legend(title="Genomic region")) +
  labs (x = "Genomic region") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")
dev.off()

#### ##### ##### ##### ##### #### ##### ##### ##### #####

# Fig.R2b - Evolutionary rate in promoter regions at each branch of conserved or divergently expressed genes (fourfold rates removed at the point of plotting)

# Pull out all the files from each directory
base_dir_br <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/br/"
br_comm_files <- list.files(path=base_dir_br, pattern = "\\.txt4$")
for (i in 1:length(br_comm_files)){
  tmp = read.delim(file = paste0(base_dir_br,br_comm_files[i]), header = T)
  assign(br_comm_files[i], tmp)
}

base_dir_ey <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/ey/"
ey_comm_files <- list.files(path=base_dir_ey, pattern = "\\.txt4$")
for (i in 1:length(ey_comm_files)){
  tmp = read.delim(file = paste0(base_dir_ey,ey_comm_files[i]), header = T)
  assign(ey_comm_files[i], tmp)
}

base_dir_ht <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/ht/"
ht_comm_files <- list.files(path=base_dir_ht, pattern = "\\.txt4$")
for (i in 1:length(ht_comm_files)){
  tmp = read.delim(file = paste0(base_dir_ht,ht_comm_files[i]), header = T)
  assign(ht_comm_files[i], tmp)
}

base_dir_kd <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/kd/"
kd_comm_files <- list.files(path=base_dir_kd, pattern = "\\.txt4$")
for (i in 1:length(kd_comm_files)){
  tmp = read.delim(file = paste0(base_dir_kd,kd_comm_files[i]), header = T)
  assign(kd_comm_files[i], tmp)
}

base_dir_ms <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/ms/"
ms_comm_files <- list.files(path=base_dir_ms, pattern = "\\.txt4$")
for (i in 1:length(ms_comm_files)){
  tmp = read.delim(file = paste0(base_dir_ms,ms_comm_files[i]), header = T)
  assign(ms_comm_files[i], tmp)
}

base_dir_ts <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/ts/"
ts_comm_files <- list.files(path=base_dir_ts, pattern = "\\.txt4$")
for (i in 1:length(ts_comm_files)){
  tmp = read.delim(file = paste0(base_dir_ts,ts_comm_files[i]), header = T)
  assign(ts_comm_files[i], tmp)
}

base_dir_br <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/br/"
br_non_comm_files <- list.files(path=base_dir_br, pattern = "\\.txt3$")
for (i in 1:length(br_non_comm_files)){
  tmp = read.delim(file = paste0(base_dir_br,br_non_comm_files[i]), header = T)
  assign(br_non_comm_files[i], tmp)
}

base_dir_ey <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/ey/"
ey_non_comm_files <- list.files(path=base_dir_ey, pattern = "\\.txt3$")
for (i in 1:length(ey_non_comm_files)){
  tmp = read.delim(file = paste0(base_dir_ey,ey_non_comm_files[i]), header = T)
  assign(ey_non_comm_files[i], tmp)
}

base_dir_ht <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/ht/"
ht_non_comm_files <- list.files(path=base_dir_ht, pattern = "\\.txt3$")
for (i in 1:length(ht_non_comm_files)){
  tmp = read.delim(file = paste0(base_dir_ht,ht_non_comm_files[i]), header = T)
  assign(ht_non_comm_files[i], tmp)
}

base_dir_kd <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/kd/"
kd_non_comm_files <- list.files(path=base_dir_kd, pattern = "\\.txt3$")
for (i in 1:length(kd_non_comm_files)){
  tmp = read.delim(file = paste0(base_dir_kd,kd_non_comm_files[i]), header = T)
  assign(kd_non_comm_files[i], tmp)
}

base_dir_ms <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/ms/"
ms_non_comm_files <- list.files(path=base_dir_ms, pattern = "\\.txt3$")
for (i in 1:length(ms_non_comm_files)){
  tmp = read.delim(file = paste0(base_dir_ms,ms_non_comm_files[i]), header = T)
  assign(ms_non_comm_files[i], tmp)
}

base_dir_ts <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Arboretum_SingleTissue_New/Results/ts/"
ts_non_comm_files <- list.files(path=base_dir_ts, pattern = "\\.txt3$")
for (i in 1:length(ts_non_comm_files)){
  tmp = read.delim(file = paste0(base_dir_ts,ts_non_comm_files[i]), header = T)
  assign(ts_non_comm_files[i], tmp)
}

# Subset the evolrate rate tables of each branch using the OGIDs in the comm and noncomm files, per tissue and per species

# Brain
br_comm_evolrate_Mz <- (Mz_evolrate[match(`comm-br-expr-mod_a_Mz.txt4`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
br_comm_evolrate_Mz <- na.omit(br_comm_evolrate_Mz) # removes rows with NA in them
br_noncomm_evolrate_Mz <- (Mz_evolrate[match(`noncomm-br-expr-mod_a_Mz.txt3`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
br_noncomm_evolrate_Mz <- na.omit(br_noncomm_evolrate_Mz) # removes rows with NA in them

br_comm_evolrate_Mz_l <- melt(br_comm_evolrate_Mz)
br_comm_evolrate_Mz_l <- cbind(br_comm_evolrate_Mz_l, type = 'conserved expr')
ggplot(br_comm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_noncomm_evolrate_Mz_l <- melt(br_noncomm_evolrate_Mz)
br_noncomm_evolrate_Mz_l <- cbind(br_noncomm_evolrate_Mz_l, type = 'divergent expr')
ggplot(br_noncomm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_comm_and_noncomm_evolrate_Mz_l <- rbind(br_comm_evolrate_Mz_l,br_noncomm_evolrate_Mz_l) # combine conserved and divergent to plot
br_comm_and_noncomm_evolrate_Mz_l2 <- subset(br_comm_and_noncomm_evolrate_Mz_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(br_comm_and_noncomm_evolrate_Mz_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_comm_evolrate_Pn <- (Pn_evolrate[match(`comm-br-expr-mod_a_Pn.txt4`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
br_comm_evolrate_Pn <- na.omit(br_comm_evolrate_Pn) # removes rows with NA in them
br_noncomm_evolrate_Pn <- (Pn_evolrate[match(`noncomm-br-expr-mod_a_Pn.txt3`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
br_noncomm_evolrate_Pn <- na.omit(br_noncomm_evolrate_Pn) # removes rows with NA in them

br_comm_evolrate_Pn_l <- melt(br_comm_evolrate_Pn)
br_comm_evolrate_Pn_l <- cbind(br_comm_evolrate_Pn_l, type = 'conserved expr')
ggplot(br_comm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_noncomm_evolrate_Pn_l <- melt(br_noncomm_evolrate_Pn)
br_noncomm_evolrate_Pn_l <- cbind(br_noncomm_evolrate_Pn_l, type = 'divergent expr')
ggplot(br_noncomm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_comm_and_noncomm_evolrate_Pn_l <- rbind(br_comm_evolrate_Pn_l,br_noncomm_evolrate_Pn_l) # combine conserved and divergent to plot
br_comm_and_noncomm_evolrate_Pn_l2 <- subset(br_comm_and_noncomm_evolrate_Pn_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(br_comm_and_noncomm_evolrate_Pn_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


br_comm_evolrate_Ab <- (Ab_evolrate[match(`comm-br-expr-mod_a_Ab.txt4`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
br_comm_evolrate_Ab <- na.omit(br_comm_evolrate_Ab) # removes rows with NA in them
br_noncomm_evolrate_Ab <- (Ab_evolrate[match(`noncomm-br-expr-mod_a_Ab.txt3`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
br_noncomm_evolrate_Ab <- na.omit(br_noncomm_evolrate_Ab) # removes rows with NA in them

br_comm_evolrate_Ab_l <- melt(br_comm_evolrate_Ab)
br_comm_evolrate_Ab_l <- cbind(br_comm_evolrate_Ab_l, type = 'conserved expr')
ggplot(br_comm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_noncomm_evolrate_Ab_l <- melt(br_noncomm_evolrate_Ab)
br_noncomm_evolrate_Ab_l <- cbind(br_noncomm_evolrate_Ab_l, type = 'divergent expr')
ggplot(br_noncomm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_comm_and_noncomm_evolrate_Ab_l <- rbind(br_comm_evolrate_Ab_l,br_noncomm_evolrate_Ab_l) # combine conserved and divergent to plot
br_comm_and_noncomm_evolrate_Ab_l2 <- subset(br_comm_and_noncomm_evolrate_Ab_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(br_comm_and_noncomm_evolrate_Ab_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


br_comm_evolrate_Nb <- (Nb_evolrate[match(`comm-br-expr-mod_a_Nb.txt4`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
br_comm_evolrate_Nb <- na.omit(br_comm_evolrate_Nb) # removes rows with NA in them
br_noncomm_evolrate_Nb <- (Nb_evolrate[match(`noncomm-br-expr-mod_a_Nb.txt3`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
br_noncomm_evolrate_Nb <- na.omit(br_noncomm_evolrate_Nb) # removes rows with NA in them

br_comm_evolrate_Nb_l <- melt(br_comm_evolrate_Nb)
br_comm_evolrate_Nb_l <- cbind(br_comm_evolrate_Nb_l, type = 'conserved expr')
ggplot(br_comm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_noncomm_evolrate_Nb_l <- melt(br_noncomm_evolrate_Nb)
br_noncomm_evolrate_Nb_l <- cbind(br_noncomm_evolrate_Nb_l, type = 'divergent expr')
ggplot(br_noncomm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_comm_and_noncomm_evolrate_Nb_l <- rbind(br_comm_evolrate_Nb_l,br_noncomm_evolrate_Nb_l) # combine conserved and divergent to plot
br_comm_and_noncomm_evolrate_Nb_l2 <- subset(br_comm_and_noncomm_evolrate_Nb_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(br_comm_and_noncomm_evolrate_Nb_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


br_comm_evolrate_On <- (On_evolrate[match(`comm-br-expr-mod_a_On.txt4`[,1],On_evolrate[,1]),(c(1,2,3,4))])
br_comm_evolrate_On <- na.omit(br_comm_evolrate_On) # removes rows with NA in them
br_noncomm_evolrate_On <- (On_evolrate[match(`noncomm-br-expr-mod_a_On.txt3`[,1],On_evolrate[,1]),(c(1,2,3,4))])
br_noncomm_evolrate_On <- na.omit(br_noncomm_evolrate_On) # removes rows with NA in them

br_comm_evolrate_On_l <- melt(br_comm_evolrate_On)
br_comm_evolrate_On_l <- cbind(br_comm_evolrate_On_l, type = 'conserved expr')
ggplot(br_comm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_noncomm_evolrate_On_l <- melt(br_noncomm_evolrate_On)
br_noncomm_evolrate_On_l <- cbind(br_noncomm_evolrate_On_l, type = 'divergent expr')
ggplot(br_noncomm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

br_comm_and_noncomm_evolrate_On_l <- rbind(br_comm_evolrate_On_l,br_noncomm_evolrate_On_l) # combine conserved and divergent to plot
br_comm_and_noncomm_evolrate_On_l2 <- subset(br_comm_and_noncomm_evolrate_On_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(br_comm_and_noncomm_evolrate_On_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

# add column for species in each and combine all species dataframes for the tissue to do a combined plot
br_comm_and_noncomm_evolrate_Mz_l3 <- cbind(br_comm_and_noncomm_evolrate_Mz_l2, species = 'M. zebra')
br_comm_and_noncomm_evolrate_Pn_l3 <- cbind(br_comm_and_noncomm_evolrate_Pn_l2, species = 'P. nyererei')
br_comm_and_noncomm_evolrate_Ab_l3 <- cbind(br_comm_and_noncomm_evolrate_Ab_l2, species = 'A. burtoni')
br_comm_and_noncomm_evolrate_Nb_l3 <- cbind(br_comm_and_noncomm_evolrate_Nb_l2, species = 'N. brichardi')
br_comm_and_noncomm_evolrate_On_l3 <- cbind(br_comm_and_noncomm_evolrate_On_l2, species = 'O. niloticus')
br_comm_and_noncomm_evolrate_combined <- rbind(br_comm_and_noncomm_evolrate_Mz_l3,br_comm_and_noncomm_evolrate_Pn_l3,br_comm_and_noncomm_evolrate_Ab_l3,br_comm_and_noncomm_evolrate_Nb_l3,br_comm_and_noncomm_evolrate_On_l3)

br_comm_and_noncomm_evolrateplot <- ggplot(br_comm_and_noncomm_evolrate_combined, aes(x = type, y = log2(value), fill=factor(type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~species) + labs(title="Brain") + guides(fill=guide_legend(title="Gene expression range")) + 
  labs (x = "Gene expression range") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")

# "Boxpot of evolutionary rate in promoter regions of genes with conserved or divergent gene expression in the brain of five cichlid species"

# Eye
ey_comm_evolrate_Mz <- (Mz_evolrate[match(`comm-ey-expr-mod_a_Mz.txt4`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
ey_comm_evolrate_Mz <- na.omit(ey_comm_evolrate_Mz) # removes rows with NA in them
ey_noncomm_evolrate_Mz <- (Mz_evolrate[match(`noncomm-ey-expr-mod_a_Mz.txt3`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
ey_noncomm_evolrate_Mz <- na.omit(ey_noncomm_evolrate_Mz) # removes rows with NA in them

ey_comm_evolrate_Mz_l <- melt(ey_comm_evolrate_Mz)
ey_comm_evolrate_Mz_l <- cbind(ey_comm_evolrate_Mz_l, type = 'conserved expr')
ggplot(ey_comm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_noncomm_evolrate_Mz_l <- melt(ey_noncomm_evolrate_Mz)
ey_noncomm_evolrate_Mz_l <- cbind(ey_noncomm_evolrate_Mz_l, type = 'divergent expr')
ggplot(ey_noncomm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_comm_and_noncomm_evolrate_Mz_l <- rbind(ey_comm_evolrate_Mz_l,ey_noncomm_evolrate_Mz_l) # combine conserved and divergent to plot
ey_comm_and_noncomm_evolrate_Mz_l2 <- subset(ey_comm_and_noncomm_evolrate_Mz_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ey_comm_and_noncomm_evolrate_Mz_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_comm_evolrate_Pn <- (Pn_evolrate[match(`comm-ey-expr-mod_a_Pn.txt4`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
ey_comm_evolrate_Pn <- na.omit(ey_comm_evolrate_Pn) # removes rows with NA in them
ey_noncomm_evolrate_Pn <- (Pn_evolrate[match(`noncomm-ey-expr-mod_a_Pn.txt3`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
ey_noncomm_evolrate_Pn <- na.omit(ey_noncomm_evolrate_Pn) # removes rows with NA in them

ey_comm_evolrate_Pn_l <- melt(ey_comm_evolrate_Pn)
ey_comm_evolrate_Pn_l <- cbind(ey_comm_evolrate_Pn_l, type = 'conserved expr')
ggplot(ey_comm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_noncomm_evolrate_Pn_l <- melt(ey_noncomm_evolrate_Pn)
ey_noncomm_evolrate_Pn_l <- cbind(ey_noncomm_evolrate_Pn_l, type = 'divergent expr')
ggplot(ey_noncomm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_comm_and_noncomm_evolrate_Pn_l <- rbind(ey_comm_evolrate_Pn_l,ey_noncomm_evolrate_Pn_l) # combine conserved and divergent to plot
ey_comm_and_noncomm_evolrate_Pn_l2 <- subset(ey_comm_and_noncomm_evolrate_Pn_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ey_comm_and_noncomm_evolrate_Pn_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ey_comm_evolrate_Ab <- (Ab_evolrate[match(`comm-ey-expr-mod_a_Ab.txt4`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
ey_comm_evolrate_Ab <- na.omit(ey_comm_evolrate_Ab) # removes rows with NA in them
ey_noncomm_evolrate_Ab <- (Ab_evolrate[match(`noncomm-ey-expr-mod_a_Ab.txt3`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
ey_noncomm_evolrate_Ab <- na.omit(ey_noncomm_evolrate_Ab) # removes rows with NA in them

ey_comm_evolrate_Ab_l <- melt(ey_comm_evolrate_Ab)
ey_comm_evolrate_Ab_l <- cbind(ey_comm_evolrate_Ab_l, type = 'conserved expr')
ggplot(ey_comm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_noncomm_evolrate_Ab_l <- melt(ey_noncomm_evolrate_Ab)
ey_noncomm_evolrate_Ab_l <- cbind(ey_noncomm_evolrate_Ab_l, type = 'divergent expr')
ggplot(ey_noncomm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_comm_and_noncomm_evolrate_Ab_l <- rbind(ey_comm_evolrate_Ab_l,ey_noncomm_evolrate_Ab_l) # combine conserved and divergent to plot
ey_comm_and_noncomm_evolrate_Ab_l2 <- subset(ey_comm_and_noncomm_evolrate_Ab_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ey_comm_and_noncomm_evolrate_Ab_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ey_comm_evolrate_Nb <- (Nb_evolrate[match(`comm-ey-expr-mod_a_Nb.txt4`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
ey_comm_evolrate_Nb <- na.omit(ey_comm_evolrate_Nb) # removes rows with NA in them
ey_noncomm_evolrate_Nb <- (Nb_evolrate[match(`noncomm-ey-expr-mod_a_Nb.txt3`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
ey_noncomm_evolrate_Nb <- na.omit(ey_noncomm_evolrate_Nb) # removes rows with NA in them

ey_comm_evolrate_Nb_l <- melt(ey_comm_evolrate_Nb)
ey_comm_evolrate_Nb_l <- cbind(ey_comm_evolrate_Nb_l, type = 'conserved expr')
ggplot(ey_comm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_noncomm_evolrate_Nb_l <- melt(ey_noncomm_evolrate_Nb)
ey_noncomm_evolrate_Nb_l <- cbind(ey_noncomm_evolrate_Nb_l, type = 'divergent expr')
ggplot(ey_noncomm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_comm_and_noncomm_evolrate_Nb_l <- rbind(ey_comm_evolrate_Nb_l,ey_noncomm_evolrate_Nb_l) # combine conserved and divergent to plot
ey_comm_and_noncomm_evolrate_Nb_l2 <- subset(ey_comm_and_noncomm_evolrate_Nb_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ey_comm_and_noncomm_evolrate_Nb_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ey_comm_evolrate_On <- (On_evolrate[match(`comm-ey-expr-mod_a_On.txt4`[,1],On_evolrate[,1]),(c(1,2,3,4))])
ey_comm_evolrate_On <- na.omit(ey_comm_evolrate_On) # removes rows with NA in them
ey_noncomm_evolrate_On <- (On_evolrate[match(`noncomm-ey-expr-mod_a_On.txt3`[,1],On_evolrate[,1]),(c(1,2,3,4))])
ey_noncomm_evolrate_On <- na.omit(ey_noncomm_evolrate_On) # removes rows with NA in them

ey_comm_evolrate_On_l <- melt(ey_comm_evolrate_On)
ey_comm_evolrate_On_l <- cbind(ey_comm_evolrate_On_l, type = 'conserved expr')
ggplot(ey_comm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_noncomm_evolrate_On_l <- melt(ey_noncomm_evolrate_On)
ey_noncomm_evolrate_On_l <- cbind(ey_noncomm_evolrate_On_l, type = 'divergent expr')
ggplot(ey_noncomm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ey_comm_and_noncomm_evolrate_On_l <- rbind(ey_comm_evolrate_On_l,ey_noncomm_evolrate_On_l) # combine conserved and divergent to plot
ey_comm_and_noncomm_evolrate_On_l2 <- subset(ey_comm_and_noncomm_evolrate_On_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ey_comm_and_noncomm_evolrate_On_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

# add column for species in each and combine all species dataframes for the tissue to do a combined plot
ey_comm_and_noncomm_evolrate_Mz_l3 <- cbind(ey_comm_and_noncomm_evolrate_Mz_l2, species = 'M. zebra')
ey_comm_and_noncomm_evolrate_Pn_l3 <- cbind(ey_comm_and_noncomm_evolrate_Pn_l2, species = 'P. nyererei')
ey_comm_and_noncomm_evolrate_Ab_l3 <- cbind(ey_comm_and_noncomm_evolrate_Ab_l2, species = 'A. burtoni')
ey_comm_and_noncomm_evolrate_Nb_l3 <- cbind(ey_comm_and_noncomm_evolrate_Nb_l2, species = 'N. brichardi')
ey_comm_and_noncomm_evolrate_On_l3 <- cbind(ey_comm_and_noncomm_evolrate_On_l2, species = 'O. niloticus')
ey_comm_and_noncomm_evolrate_combined <- rbind(ey_comm_and_noncomm_evolrate_Mz_l3,ey_comm_and_noncomm_evolrate_Pn_l3,ey_comm_and_noncomm_evolrate_Ab_l3,ey_comm_and_noncomm_evolrate_Nb_l3,ey_comm_and_noncomm_evolrate_On_l3)

ey_comm_and_noncomm_evolrateplot <- ggplot(ey_comm_and_noncomm_evolrate_combined, aes(x = type, y = log2(value), fill=factor(type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~species) + labs(title="Eye") + guides(fill=guide_legend(title="Gene expression range")) + 
  labs (x = "Gene expression range") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")

# "Boxpot of evolutionary rate in promoter regions of genes with conserved or divergent gene expression in the eye of five cichlid species"

# Heart
ht_comm_evolrate_Mz <- (Mz_evolrate[match(`comm-ht-expr-mod_a_Mz.txt4`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
ht_comm_evolrate_Mz <- na.omit(ht_comm_evolrate_Mz) # removes rows with NA in them
ht_noncomm_evolrate_Mz <- (Mz_evolrate[match(`noncomm-ht-expr-mod_a_Mz.txt3`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
ht_noncomm_evolrate_Mz <- na.omit(ht_noncomm_evolrate_Mz) # removes rows with NA in them

ht_comm_evolrate_Mz_l <- melt(ht_comm_evolrate_Mz)
ht_comm_evolrate_Mz_l <- cbind(ht_comm_evolrate_Mz_l, type = 'conserved expr')
ggplot(ht_comm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_noncomm_evolrate_Mz_l <- melt(ht_noncomm_evolrate_Mz)
ht_noncomm_evolrate_Mz_l <- cbind(ht_noncomm_evolrate_Mz_l, type = 'divergent expr')
ggplot(ht_noncomm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_comm_and_noncomm_evolrate_Mz_l <- rbind(ht_comm_evolrate_Mz_l,ht_noncomm_evolrate_Mz_l) # combine conserved and divergent to plot
ht_comm_and_noncomm_evolrate_Mz_l2 <- subset(ht_comm_and_noncomm_evolrate_Mz_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ht_comm_and_noncomm_evolrate_Mz_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_comm_evolrate_Pn <- (Pn_evolrate[match(`comm-ht-expr-mod_a_Pn.txt4`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
ht_comm_evolrate_Pn <- na.omit(ht_comm_evolrate_Pn) # removes rows with NA in them
ht_noncomm_evolrate_Pn <- (Pn_evolrate[match(`noncomm-ht-expr-mod_a_Pn.txt3`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
ht_noncomm_evolrate_Pn <- na.omit(ht_noncomm_evolrate_Pn) # removes rows with NA in them

ht_comm_evolrate_Pn_l <- melt(ht_comm_evolrate_Pn)
ht_comm_evolrate_Pn_l <- cbind(ht_comm_evolrate_Pn_l, type = 'conserved expr')
ggplot(ht_comm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_noncomm_evolrate_Pn_l <- melt(ht_noncomm_evolrate_Pn)
ht_noncomm_evolrate_Pn_l <- cbind(ht_noncomm_evolrate_Pn_l, type = 'divergent expr')
ggplot(ht_noncomm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_comm_and_noncomm_evolrate_Pn_l <- rbind(ht_comm_evolrate_Pn_l,ht_noncomm_evolrate_Pn_l) # combine conserved and divergent to plot
ht_comm_and_noncomm_evolrate_Pn_l2 <- subset(ht_comm_and_noncomm_evolrate_Pn_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ht_comm_and_noncomm_evolrate_Pn_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ht_comm_evolrate_Ab <- (Ab_evolrate[match(`comm-ht-expr-mod_a_Ab.txt4`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
ht_comm_evolrate_Ab <- na.omit(ht_comm_evolrate_Ab) # removes rows with NA in them
ht_noncomm_evolrate_Ab <- (Ab_evolrate[match(`noncomm-ht-expr-mod_a_Ab.txt3`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
ht_noncomm_evolrate_Ab <- na.omit(ht_noncomm_evolrate_Ab) # removes rows with NA in them

ht_comm_evolrate_Ab_l <- melt(ht_comm_evolrate_Ab)
ht_comm_evolrate_Ab_l <- cbind(ht_comm_evolrate_Ab_l, type = 'conserved expr')
ggplot(ht_comm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_noncomm_evolrate_Ab_l <- melt(ht_noncomm_evolrate_Ab)
ht_noncomm_evolrate_Ab_l <- cbind(ht_noncomm_evolrate_Ab_l, type = 'divergent expr')
ggplot(ht_noncomm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_comm_and_noncomm_evolrate_Ab_l <- rbind(ht_comm_evolrate_Ab_l,ht_noncomm_evolrate_Ab_l) # combine conserved and divergent to plot
ht_comm_and_noncomm_evolrate_Ab_l2 <- subset(ht_comm_and_noncomm_evolrate_Ab_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ht_comm_and_noncomm_evolrate_Ab_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ht_comm_evolrate_Nb <- (Nb_evolrate[match(`comm-ht-expr-mod_a_Nb.txt4`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
ht_comm_evolrate_Nb <- na.omit(ht_comm_evolrate_Nb) # removes rows with NA in them
ht_noncomm_evolrate_Nb <- (Nb_evolrate[match(`noncomm-ht-expr-mod_a_Nb.txt3`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
ht_noncomm_evolrate_Nb <- na.omit(ht_noncomm_evolrate_Nb) # removes rows with NA in them

ht_comm_evolrate_Nb_l <- melt(ht_comm_evolrate_Nb)
ht_comm_evolrate_Nb_l <- cbind(ht_comm_evolrate_Nb_l, type = 'conserved expr')
ggplot(ht_comm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_noncomm_evolrate_Nb_l <- melt(ht_noncomm_evolrate_Nb)
ht_noncomm_evolrate_Nb_l <- cbind(ht_noncomm_evolrate_Nb_l, type = 'divergent expr')
ggplot(ht_noncomm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_comm_and_noncomm_evolrate_Nb_l <- rbind(ht_comm_evolrate_Nb_l,ht_noncomm_evolrate_Nb_l) # combine conserved and divergent to plot
ht_comm_and_noncomm_evolrate_Nb_l2 <- subset(ht_comm_and_noncomm_evolrate_Nb_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ht_comm_and_noncomm_evolrate_Nb_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ht_comm_evolrate_On <- (On_evolrate[match(`comm-ht-expr-mod_a_On.txt4`[,1],On_evolrate[,1]),(c(1,2,3,4))])
ht_comm_evolrate_On <- na.omit(ht_comm_evolrate_On) # removes rows with NA in them
ht_noncomm_evolrate_On <- (On_evolrate[match(`noncomm-ht-expr-mod_a_On.txt3`[,1],On_evolrate[,1]),(c(1,2,3,4))])
ht_noncomm_evolrate_On <- na.omit(ht_noncomm_evolrate_On) # removes rows with NA in them

ht_comm_evolrate_On_l <- melt(ht_comm_evolrate_On)
ht_comm_evolrate_On_l <- cbind(ht_comm_evolrate_On_l, type = 'conserved expr')
ggplot(ht_comm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_noncomm_evolrate_On_l <- melt(ht_noncomm_evolrate_On)
ht_noncomm_evolrate_On_l <- cbind(ht_noncomm_evolrate_On_l, type = 'divergent expr')
ggplot(ht_noncomm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ht_comm_and_noncomm_evolrate_On_l <- rbind(ht_comm_evolrate_On_l,ht_noncomm_evolrate_On_l) # combine conserved and divergent to plot
ht_comm_and_noncomm_evolrate_On_l2 <- subset(ht_comm_and_noncomm_evolrate_On_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ht_comm_and_noncomm_evolrate_On_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

# add column for species in each and combine all species dataframes for the tissue to do a combined plot
ht_comm_and_noncomm_evolrate_Mz_l3 <- cbind(ht_comm_and_noncomm_evolrate_Mz_l2, species = 'M. zebra')
ht_comm_and_noncomm_evolrate_Pn_l3 <- cbind(ht_comm_and_noncomm_evolrate_Pn_l2, species = 'P. nyererei')
ht_comm_and_noncomm_evolrate_Ab_l3 <- cbind(ht_comm_and_noncomm_evolrate_Ab_l2, species = 'A. burtoni')
ht_comm_and_noncomm_evolrate_Nb_l3 <- cbind(ht_comm_and_noncomm_evolrate_Nb_l2, species = 'N. brichardi')
ht_comm_and_noncomm_evolrate_On_l3 <- cbind(ht_comm_and_noncomm_evolrate_On_l2, species = 'O. niloticus')
ht_comm_and_noncomm_evolrate_combined <- rbind(ht_comm_and_noncomm_evolrate_Mz_l3,ht_comm_and_noncomm_evolrate_Pn_l3,ht_comm_and_noncomm_evolrate_Ab_l3,ht_comm_and_noncomm_evolrate_Nb_l3,ht_comm_and_noncomm_evolrate_On_l3)

ht_comm_and_noncomm_evolrateplot <- ggplot(ht_comm_and_noncomm_evolrate_combined, aes(x = type, y = log2(value), fill=factor(type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~species) + labs(title="Heart") + guides(fill=guide_legend(title="Gene expression range")) + 
  labs (x = "Gene expression range") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")

# "Boxpot of evolutionary rate in promoter regions of genes with conserved or divergent gene expression in the heart of five cichlid species"

# Kidney
kd_comm_evolrate_Mz <- (Mz_evolrate[match(`comm-kd-expr-mod_a_Mz.txt4`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
kd_comm_evolrate_Mz <- na.omit(kd_comm_evolrate_Mz) # removes rows with NA in them
kd_noncomm_evolrate_Mz <- (Mz_evolrate[match(`noncomm-kd-expr-mod_a_Mz.txt3`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
kd_noncomm_evolrate_Mz <- na.omit(kd_noncomm_evolrate_Mz) # removes rows with NA in them

kd_comm_evolrate_Mz_l <- melt(kd_comm_evolrate_Mz)
kd_comm_evolrate_Mz_l <- cbind(kd_comm_evolrate_Mz_l, type = 'conserved expr')
ggplot(kd_comm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_noncomm_evolrate_Mz_l <- melt(kd_noncomm_evolrate_Mz)
kd_noncomm_evolrate_Mz_l <- cbind(kd_noncomm_evolrate_Mz_l, type = 'divergent expr')
ggplot(kd_noncomm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_comm_and_noncomm_evolrate_Mz_l <- rbind(kd_comm_evolrate_Mz_l,kd_noncomm_evolrate_Mz_l) # combine conserved and divergent to plot
kd_comm_and_noncomm_evolrate_Mz_l2 <- subset(kd_comm_and_noncomm_evolrate_Mz_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(kd_comm_and_noncomm_evolrate_Mz_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_comm_evolrate_Pn <- (Pn_evolrate[match(`comm-kd-expr-mod_a_Pn.txt4`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
kd_comm_evolrate_Pn <- na.omit(kd_comm_evolrate_Pn) # removes rows with NA in them
kd_noncomm_evolrate_Pn <- (Pn_evolrate[match(`noncomm-kd-expr-mod_a_Pn.txt3`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
kd_noncomm_evolrate_Pn <- na.omit(kd_noncomm_evolrate_Pn) # removes rows with NA in them

kd_comm_evolrate_Pn_l <- melt(kd_comm_evolrate_Pn)
kd_comm_evolrate_Pn_l <- cbind(kd_comm_evolrate_Pn_l, type = 'conserved expr')
ggplot(kd_comm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_noncomm_evolrate_Pn_l <- melt(kd_noncomm_evolrate_Pn)
kd_noncomm_evolrate_Pn_l <- cbind(kd_noncomm_evolrate_Pn_l, type = 'divergent expr')
ggplot(kd_noncomm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_comm_and_noncomm_evolrate_Pn_l <- rbind(kd_comm_evolrate_Pn_l,kd_noncomm_evolrate_Pn_l) # combine conserved and divergent to plot
kd_comm_and_noncomm_evolrate_Pn_l2 <- subset(kd_comm_and_noncomm_evolrate_Pn_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(kd_comm_and_noncomm_evolrate_Pn_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


kd_comm_evolrate_Ab <- (Ab_evolrate[match(`comm-kd-expr-mod_a_Ab.txt4`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
kd_comm_evolrate_Ab <- na.omit(kd_comm_evolrate_Ab) # removes rows with NA in them
kd_noncomm_evolrate_Ab <- (Ab_evolrate[match(`noncomm-kd-expr-mod_a_Ab.txt3`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
kd_noncomm_evolrate_Ab <- na.omit(kd_noncomm_evolrate_Ab) # removes rows with NA in them

kd_comm_evolrate_Ab_l <- melt(kd_comm_evolrate_Ab)
kd_comm_evolrate_Ab_l <- cbind(kd_comm_evolrate_Ab_l, type = 'conserved expr')
ggplot(kd_comm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_noncomm_evolrate_Ab_l <- melt(kd_noncomm_evolrate_Ab)
kd_noncomm_evolrate_Ab_l <- cbind(kd_noncomm_evolrate_Ab_l, type = 'divergent expr')
ggplot(kd_noncomm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_comm_and_noncomm_evolrate_Ab_l <- rbind(kd_comm_evolrate_Ab_l,kd_noncomm_evolrate_Ab_l) # combine conserved and divergent to plot
kd_comm_and_noncomm_evolrate_Ab_l2 <- subset(kd_comm_and_noncomm_evolrate_Ab_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(kd_comm_and_noncomm_evolrate_Ab_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


kd_comm_evolrate_Nb <- (Nb_evolrate[match(`comm-kd-expr-mod_a_Nb.txt4`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
kd_comm_evolrate_Nb <- na.omit(kd_comm_evolrate_Nb) # removes rows with NA in them
kd_noncomm_evolrate_Nb <- (Nb_evolrate[match(`noncomm-kd-expr-mod_a_Nb.txt3`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
kd_noncomm_evolrate_Nb <- na.omit(kd_noncomm_evolrate_Nb) # removes rows with NA in them

kd_comm_evolrate_Nb_l <- melt(kd_comm_evolrate_Nb)
kd_comm_evolrate_Nb_l <- cbind(kd_comm_evolrate_Nb_l, type = 'conserved expr')
ggplot(kd_comm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_noncomm_evolrate_Nb_l <- melt(kd_noncomm_evolrate_Nb)
kd_noncomm_evolrate_Nb_l <- cbind(kd_noncomm_evolrate_Nb_l, type = 'divergent expr')
ggplot(kd_noncomm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_comm_and_noncomm_evolrate_Nb_l <- rbind(kd_comm_evolrate_Nb_l,kd_noncomm_evolrate_Nb_l) # combine conserved and divergent to plot
kd_comm_and_noncomm_evolrate_Nb_l2 <- subset(kd_comm_and_noncomm_evolrate_Nb_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(kd_comm_and_noncomm_evolrate_Nb_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


kd_comm_evolrate_On <- (On_evolrate[match(`comm-kd-expr-mod_a_On.txt4`[,1],On_evolrate[,1]),(c(1,2,3,4))])
kd_comm_evolrate_On <- na.omit(kd_comm_evolrate_On) # removes rows with NA in them
kd_noncomm_evolrate_On <- (On_evolrate[match(`noncomm-kd-expr-mod_a_On.txt3`[,1],On_evolrate[,1]),(c(1,2,3,4))])
kd_noncomm_evolrate_On <- na.omit(kd_noncomm_evolrate_On) # removes rows with NA in them

kd_comm_evolrate_On_l <- melt(kd_comm_evolrate_On)
kd_comm_evolrate_On_l <- cbind(kd_comm_evolrate_On_l, type = 'conserved expr')
ggplot(kd_comm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_noncomm_evolrate_On_l <- melt(kd_noncomm_evolrate_On)
kd_noncomm_evolrate_On_l <- cbind(kd_noncomm_evolrate_On_l, type = 'divergent expr')
ggplot(kd_noncomm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

kd_comm_and_noncomm_evolrate_On_l <- rbind(kd_comm_evolrate_On_l,kd_noncomm_evolrate_On_l) # combine conserved and divergent to plot
kd_comm_and_noncomm_evolrate_On_l2 <- subset(kd_comm_and_noncomm_evolrate_On_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(kd_comm_and_noncomm_evolrate_On_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

# add column for species in each and combine all species dataframes for the tissue to do a combined plot
kd_comm_and_noncomm_evolrate_Mz_l3 <- cbind(kd_comm_and_noncomm_evolrate_Mz_l2, species = 'M. zebra')
kd_comm_and_noncomm_evolrate_Pn_l3 <- cbind(kd_comm_and_noncomm_evolrate_Pn_l2, species = 'P. nyererei')
kd_comm_and_noncomm_evolrate_Ab_l3 <- cbind(kd_comm_and_noncomm_evolrate_Ab_l2, species = 'A. burtoni')
kd_comm_and_noncomm_evolrate_Nb_l3 <- cbind(kd_comm_and_noncomm_evolrate_Nb_l2, species = 'N. brichardi')
kd_comm_and_noncomm_evolrate_On_l3 <- cbind(kd_comm_and_noncomm_evolrate_On_l2, species = 'O. niloticus')
kd_comm_and_noncomm_evolrate_combined <- rbind(kd_comm_and_noncomm_evolrate_Mz_l3,kd_comm_and_noncomm_evolrate_Pn_l3,kd_comm_and_noncomm_evolrate_Ab_l3,kd_comm_and_noncomm_evolrate_Nb_l3,kd_comm_and_noncomm_evolrate_On_l3)

kd_comm_and_noncomm_evolrateplot <- ggplot(kd_comm_and_noncomm_evolrate_combined, aes(x = type, y = log2(value), fill=factor(type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~species) + labs(title="Kidney") + guides(fill=guide_legend(title="Gene expression range")) + 
  labs (x = "Gene expression range") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")

# "Boxpot of evolutionary rate in promoter regions of genes with conserved or divergent gene expression in the kidney of five cichlid species"

# Muscle
ms_comm_evolrate_Mz <- (Mz_evolrate[match(`comm-ms-expr-mod_a_Mz.txt4`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
ms_comm_evolrate_Mz <- na.omit(ms_comm_evolrate_Mz) # removes rows with NA in them
ms_noncomm_evolrate_Mz <- (Mz_evolrate[match(`noncomm-ms-expr-mod_a_Mz.txt3`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
ms_noncomm_evolrate_Mz <- na.omit(ms_noncomm_evolrate_Mz) # removes rows with NA in them

ms_comm_evolrate_Mz_l <- melt(ms_comm_evolrate_Mz)
ms_comm_evolrate_Mz_l <- cbind(ms_comm_evolrate_Mz_l, type = 'conserved expr')
ggplot(ms_comm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_noncomm_evolrate_Mz_l <- melt(ms_noncomm_evolrate_Mz)
ms_noncomm_evolrate_Mz_l <- cbind(ms_noncomm_evolrate_Mz_l, type = 'divergent expr')
ggplot(ms_noncomm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_comm_and_noncomm_evolrate_Mz_l <- rbind(ms_comm_evolrate_Mz_l,ms_noncomm_evolrate_Mz_l) # combine conserved and divergent to plot
ms_comm_and_noncomm_evolrate_Mz_l2 <- subset(ms_comm_and_noncomm_evolrate_Mz_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ms_comm_and_noncomm_evolrate_Mz_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_comm_evolrate_Pn <- (Pn_evolrate[match(`comm-ms-expr-mod_a_Pn.txt4`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
ms_comm_evolrate_Pn <- na.omit(ms_comm_evolrate_Pn) # removes rows with NA in them
ms_noncomm_evolrate_Pn <- (Pn_evolrate[match(`noncomm-ms-expr-mod_a_Pn.txt3`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
ms_noncomm_evolrate_Pn <- na.omit(ms_noncomm_evolrate_Pn) # removes rows with NA in them

ms_comm_evolrate_Pn_l <- melt(ms_comm_evolrate_Pn)
ms_comm_evolrate_Pn_l <- cbind(ms_comm_evolrate_Pn_l, type = 'conserved expr')
ggplot(ms_comm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_noncomm_evolrate_Pn_l <- melt(ms_noncomm_evolrate_Pn)
ms_noncomm_evolrate_Pn_l <- cbind(ms_noncomm_evolrate_Pn_l, type = 'divergent expr')
ggplot(ms_noncomm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_comm_and_noncomm_evolrate_Pn_l <- rbind(ms_comm_evolrate_Pn_l,ms_noncomm_evolrate_Pn_l) # combine conserved and divergent to plot
ms_comm_and_noncomm_evolrate_Pn_l2 <- subset(ms_comm_and_noncomm_evolrate_Pn_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ms_comm_and_noncomm_evolrate_Pn_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ms_comm_evolrate_Ab <- (Ab_evolrate[match(`comm-ms-expr-mod_a_Ab.txt4`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
ms_comm_evolrate_Ab <- na.omit(ms_comm_evolrate_Ab) # removes rows with NA in them
ms_noncomm_evolrate_Ab <- (Ab_evolrate[match(`noncomm-ms-expr-mod_a_Ab.txt3`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
ms_noncomm_evolrate_Ab <- na.omit(ms_noncomm_evolrate_Ab) # removes rows with NA in them

ms_comm_evolrate_Ab_l <- melt(ms_comm_evolrate_Ab)
ms_comm_evolrate_Ab_l <- cbind(ms_comm_evolrate_Ab_l, type = 'conserved expr')
ggplot(ms_comm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_noncomm_evolrate_Ab_l <- melt(ms_noncomm_evolrate_Ab)
ms_noncomm_evolrate_Ab_l <- cbind(ms_noncomm_evolrate_Ab_l, type = 'divergent expr')
ggplot(ms_noncomm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_comm_and_noncomm_evolrate_Ab_l <- rbind(ms_comm_evolrate_Ab_l,ms_noncomm_evolrate_Ab_l) # combine conserved and divergent to plot
ms_comm_and_noncomm_evolrate_Ab_l2 <- subset(ms_comm_and_noncomm_evolrate_Ab_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ms_comm_and_noncomm_evolrate_Ab_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ms_comm_evolrate_Nb <- (Nb_evolrate[match(`comm-ms-expr-mod_a_Nb.txt4`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
ms_comm_evolrate_Nb <- na.omit(ms_comm_evolrate_Nb) # removes rows with NA in them
ms_noncomm_evolrate_Nb <- (Nb_evolrate[match(`noncomm-ms-expr-mod_a_Nb.txt3`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
ms_noncomm_evolrate_Nb <- na.omit(ms_noncomm_evolrate_Nb) # removes rows with NA in them

ms_comm_evolrate_Nb_l <- melt(ms_comm_evolrate_Nb)
ms_comm_evolrate_Nb_l <- cbind(ms_comm_evolrate_Nb_l, type = 'conserved expr')
ggplot(ms_comm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_noncomm_evolrate_Nb_l <- melt(ms_noncomm_evolrate_Nb)
ms_noncomm_evolrate_Nb_l <- cbind(ms_noncomm_evolrate_Nb_l, type = 'divergent expr')
ggplot(ms_noncomm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_comm_and_noncomm_evolrate_Nb_l <- rbind(ms_comm_evolrate_Nb_l,ms_noncomm_evolrate_Nb_l) # combine conserved and divergent to plot
ms_comm_and_noncomm_evolrate_Nb_l2 <- subset(ms_comm_and_noncomm_evolrate_Nb_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ms_comm_and_noncomm_evolrate_Nb_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ms_comm_evolrate_On <- (On_evolrate[match(`comm-ms-expr-mod_a_On.txt4`[,1],On_evolrate[,1]),(c(1,2,3,4))])
ms_comm_evolrate_On <- na.omit(ms_comm_evolrate_On) # removes rows with NA in them
ms_noncomm_evolrate_On <- (On_evolrate[match(`noncomm-ms-expr-mod_a_On.txt3`[,1],On_evolrate[,1]),(c(1,2,3,4))])
ms_noncomm_evolrate_On <- na.omit(ms_noncomm_evolrate_On) # removes rows with NA in them

ms_comm_evolrate_On_l <- melt(ms_comm_evolrate_On)
ms_comm_evolrate_On_l <- cbind(ms_comm_evolrate_On_l, type = 'conserved expr')
ggplot(ms_comm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_noncomm_evolrate_On_l <- melt(ms_noncomm_evolrate_On)
ms_noncomm_evolrate_On_l <- cbind(ms_noncomm_evolrate_On_l, type = 'divergent expr')
ggplot(ms_noncomm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ms_comm_and_noncomm_evolrate_On_l <- rbind(ms_comm_evolrate_On_l,ms_noncomm_evolrate_On_l) # combine conserved and divergent to plot
ms_comm_and_noncomm_evolrate_On_l2 <- subset(ms_comm_and_noncomm_evolrate_On_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ms_comm_and_noncomm_evolrate_On_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

# add column for species in each and combine all species dataframes for the tissue to do a combined plot
ms_comm_and_noncomm_evolrate_Mz_l3 <- cbind(ms_comm_and_noncomm_evolrate_Mz_l2, species = 'M. zebra')
ms_comm_and_noncomm_evolrate_Pn_l3 <- cbind(ms_comm_and_noncomm_evolrate_Pn_l2, species = 'P. nyererei')
ms_comm_and_noncomm_evolrate_Ab_l3 <- cbind(ms_comm_and_noncomm_evolrate_Ab_l2, species = 'A. burtoni')
ms_comm_and_noncomm_evolrate_Nb_l3 <- cbind(ms_comm_and_noncomm_evolrate_Nb_l2, species = 'N. brichardi')
ms_comm_and_noncomm_evolrate_On_l3 <- cbind(ms_comm_and_noncomm_evolrate_On_l2, species = 'O. niloticus')
ms_comm_and_noncomm_evolrate_combined <- rbind(ms_comm_and_noncomm_evolrate_Mz_l3,ms_comm_and_noncomm_evolrate_Pn_l3,ms_comm_and_noncomm_evolrate_Ab_l3,ms_comm_and_noncomm_evolrate_Nb_l3,ms_comm_and_noncomm_evolrate_On_l3)

ms_comm_and_noncomm_evolrateplot <- ggplot(ms_comm_and_noncomm_evolrate_combined, aes(x = type, y = log2(value), fill=factor(type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~species) + labs(title="Muscle") + guides(fill=guide_legend(title="Gene expression range")) + 
  labs (x = "Gene expression range") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")

# "Boxpot of evolutionary rate in promoter regions of genes with conserved or divergent gene expression in the muscle of five cichlid species"

# Testis
ts_comm_evolrate_Mz <- (Mz_evolrate[match(`comm-ts-expr-mod_a_Mz.txt4`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
ts_comm_evolrate_Mz <- na.omit(ts_comm_evolrate_Mz) # removes rows with NA in them
ts_noncomm_evolrate_Mz <- (Mz_evolrate[match(`noncomm-ts-expr-mod_a_Mz.txt3`[,1],Mz_evolrate[,1]),(c(1,2,3,4))])
ts_noncomm_evolrate_Mz <- na.omit(ts_noncomm_evolrate_Mz) # removes rows with NA in them

ts_comm_evolrate_Mz_l <- melt(ts_comm_evolrate_Mz)
ts_comm_evolrate_Mz_l <- cbind(ts_comm_evolrate_Mz_l, type = 'conserved expr')
ggplot(ts_comm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_noncomm_evolrate_Mz_l <- melt(ts_noncomm_evolrate_Mz)
ts_noncomm_evolrate_Mz_l <- cbind(ts_noncomm_evolrate_Mz_l, type = 'divergent expr')
ggplot(ts_noncomm_evolrate_Mz_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_comm_and_noncomm_evolrate_Mz_l <- rbind(ts_comm_evolrate_Mz_l,ts_noncomm_evolrate_Mz_l) # combine conserved and divergent to plot
ts_comm_and_noncomm_evolrate_Mz_l2 <- subset(ts_comm_and_noncomm_evolrate_Mz_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ts_comm_and_noncomm_evolrate_Mz_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_comm_evolrate_Pn <- (Pn_evolrate[match(`comm-ts-expr-mod_a_Pn.txt4`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
ts_comm_evolrate_Pn <- na.omit(ts_comm_evolrate_Pn) # removes rows with NA in them
ts_noncomm_evolrate_Pn <- (Pn_evolrate[match(`noncomm-ts-expr-mod_a_Pn.txt3`[,1],Pn_evolrate[,1]),(c(1,2,3,4))])
ts_noncomm_evolrate_Pn <- na.omit(ts_noncomm_evolrate_Pn) # removes rows with NA in them

ts_comm_evolrate_Pn_l <- melt(ts_comm_evolrate_Pn)
ts_comm_evolrate_Pn_l <- cbind(ts_comm_evolrate_Pn_l, type = 'conserved expr')
ggplot(ts_comm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_noncomm_evolrate_Pn_l <- melt(ts_noncomm_evolrate_Pn)
ts_noncomm_evolrate_Pn_l <- cbind(ts_noncomm_evolrate_Pn_l, type = 'divergent expr')
ggplot(ts_noncomm_evolrate_Pn_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_comm_and_noncomm_evolrate_Pn_l <- rbind(ts_comm_evolrate_Pn_l,ts_noncomm_evolrate_Pn_l) # combine conserved and divergent to plot
ts_comm_and_noncomm_evolrate_Pn_l2 <- subset(ts_comm_and_noncomm_evolrate_Pn_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ts_comm_and_noncomm_evolrate_Pn_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ts_comm_evolrate_Ab <- (Ab_evolrate[match(`comm-ts-expr-mod_a_Ab.txt4`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
ts_comm_evolrate_Ab <- na.omit(ts_comm_evolrate_Ab) # removes rows with NA in them
ts_noncomm_evolrate_Ab <- (Ab_evolrate[match(`noncomm-ts-expr-mod_a_Ab.txt3`[,1],Ab_evolrate[,1]),(c(1,2,3,4))])
ts_noncomm_evolrate_Ab <- na.omit(ts_noncomm_evolrate_Ab) # removes rows with NA in them

ts_comm_evolrate_Ab_l <- melt(ts_comm_evolrate_Ab)
ts_comm_evolrate_Ab_l <- cbind(ts_comm_evolrate_Ab_l, type = 'conserved expr')
ggplot(ts_comm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_noncomm_evolrate_Ab_l <- melt(ts_noncomm_evolrate_Ab)
ts_noncomm_evolrate_Ab_l <- cbind(ts_noncomm_evolrate_Ab_l, type = 'divergent expr')
ggplot(ts_noncomm_evolrate_Ab_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_comm_and_noncomm_evolrate_Ab_l <- rbind(ts_comm_evolrate_Ab_l,ts_noncomm_evolrate_Ab_l) # combine conserved and divergent to plot
ts_comm_and_noncomm_evolrate_Ab_l2 <- subset(ts_comm_and_noncomm_evolrate_Ab_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ts_comm_and_noncomm_evolrate_Ab_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ts_comm_evolrate_Nb <- (Nb_evolrate[match(`comm-ts-expr-mod_a_Nb.txt4`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
ts_comm_evolrate_Nb <- na.omit(ts_comm_evolrate_Nb) # removes rows with NA in them
ts_noncomm_evolrate_Nb <- (Nb_evolrate[match(`noncomm-ts-expr-mod_a_Nb.txt3`[,1],Nb_evolrate[,1]),(c(1,2,3,4))])
ts_noncomm_evolrate_Nb <- na.omit(ts_noncomm_evolrate_Nb) # removes rows with NA in them

ts_comm_evolrate_Nb_l <- melt(ts_comm_evolrate_Nb)
ts_comm_evolrate_Nb_l <- cbind(ts_comm_evolrate_Nb_l, type = 'conserved expr')
ggplot(ts_comm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_noncomm_evolrate_Nb_l <- melt(ts_noncomm_evolrate_Nb)
ts_noncomm_evolrate_Nb_l <- cbind(ts_noncomm_evolrate_Nb_l, type = 'divergent expr')
ggplot(ts_noncomm_evolrate_Nb_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_comm_and_noncomm_evolrate_Nb_l <- rbind(ts_comm_evolrate_Nb_l,ts_noncomm_evolrate_Nb_l) # combine conserved and divergent to plot
ts_comm_and_noncomm_evolrate_Nb_l2 <- subset(ts_comm_and_noncomm_evolrate_Nb_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ts_comm_and_noncomm_evolrate_Nb_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())


ts_comm_evolrate_On <- (On_evolrate[match(`comm-ts-expr-mod_a_On.txt4`[,1],On_evolrate[,1]),(c(1,2,3,4))])
ts_comm_evolrate_On <- na.omit(ts_comm_evolrate_On) # removes rows with NA in them
ts_noncomm_evolrate_On <- (On_evolrate[match(`noncomm-ts-expr-mod_a_On.txt3`[,1],On_evolrate[,1]),(c(1,2,3,4))])
ts_noncomm_evolrate_On <- na.omit(ts_noncomm_evolrate_On) # removes rows with NA in them

ts_comm_evolrate_On_l <- melt(ts_comm_evolrate_On)
ts_comm_evolrate_On_l <- cbind(ts_comm_evolrate_On_l, type = 'conserved expr')
ggplot(ts_comm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_noncomm_evolrate_On_l <- melt(ts_noncomm_evolrate_On)
ts_noncomm_evolrate_On_l <- cbind(ts_noncomm_evolrate_On_l, type = 'divergent expr')
ggplot(ts_noncomm_evolrate_On_l, aes(x = variable, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

ts_comm_and_noncomm_evolrate_On_l <- rbind(ts_comm_evolrate_On_l,ts_noncomm_evolrate_On_l) # combine conserved and divergent to plot
ts_comm_and_noncomm_evolrate_On_l2 <- subset(ts_comm_and_noncomm_evolrate_On_l, variable!="Fourfold") # Remove rows with fourfold
ggplot(ts_comm_and_noncomm_evolrate_On_l2, aes(x = type, y = log2(value))) + geom_boxplot() +
  theme_bw() + theme(panel.grid = element_blank())

# add column for species in each and combine all species dataframes for the tissue to do a combined plot
ts_comm_and_noncomm_evolrate_Mz_l3 <- cbind(ts_comm_and_noncomm_evolrate_Mz_l2, species = 'M. zebra')
ts_comm_and_noncomm_evolrate_Pn_l3 <- cbind(ts_comm_and_noncomm_evolrate_Pn_l2, species = 'P. nyererei')
ts_comm_and_noncomm_evolrate_Ab_l3 <- cbind(ts_comm_and_noncomm_evolrate_Ab_l2, species = 'A. burtoni')
ts_comm_and_noncomm_evolrate_Nb_l3 <- cbind(ts_comm_and_noncomm_evolrate_Nb_l2, species = 'N. brichardi')
ts_comm_and_noncomm_evolrate_On_l3 <- cbind(ts_comm_and_noncomm_evolrate_On_l2, species = 'O. niloticus')
ts_comm_and_noncomm_evolrate_combined <- rbind(ts_comm_and_noncomm_evolrate_Mz_l3,ts_comm_and_noncomm_evolrate_Pn_l3,ts_comm_and_noncomm_evolrate_Ab_l3,ts_comm_and_noncomm_evolrate_Nb_l3,ts_comm_and_noncomm_evolrate_On_l3)

ts_comm_and_noncomm_evolrateplot <- ggplot(ts_comm_and_noncomm_evolrate_combined, aes(x = type, y = log2(value), fill=factor(type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~species) + labs(title="Testis") + guides(fill=guide_legend(title="Gene expression range")) + 
  labs (x = "Gene expression range") + labs (y = "log2 (evolutionary rate)") + theme(legend.position="none")

# "Boxpot of evolutionary rate in promoter regions of genes with conserved or divergent gene expression in the testis of five cichlid species"

# multiplot of conserved/divergent expression in each tissue > cons_div_expr_evolplot.png
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/cons_div_expr_evolplot.tiff', units="in", width=20, height=20, res=300)
multiplot(br_comm_and_noncomm_evolrateplot, ey_comm_and_noncomm_evolrateplot, ht_comm_and_noncomm_evolrateplot, kd_comm_and_noncomm_evolrateplot, ms_comm_and_noncomm_evolrateplot, ts_comm_and_noncomm_evolrateplot, cols = 2)
dev.off()

### Rank genes based on fold diff (prom vs fourfold evol rate) and are divergently expressed in each tissue and species
# Load in the OGIDs to cbind the genes
OGIDs <- read.delim(file = "/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5", header = F, row.names = 1)

# create the per-species and per tissue files
br_noncomm_evolrate_Mz$tissue <- "Brain"
br_noncomm_evolrate_Pn$tissue <- "Brain"
br_noncomm_evolrate_Ab$tissue <- "Brain"
br_noncomm_evolrate_Nb$tissue <- "Brain"
br_noncomm_evolrate_On$tissue <- "Brain"
ey_noncomm_evolrate_Mz$tissue <- "Eye"
ey_noncomm_evolrate_Pn$tissue <- "Eye"
ey_noncomm_evolrate_Ab$tissue <- "Eye"
ey_noncomm_evolrate_Nb$tissue <- "Eye"
ey_noncomm_evolrate_On$tissue <- "Eye"
ht_noncomm_evolrate_Mz$tissue <- "Heart"
ht_noncomm_evolrate_Pn$tissue <- "Heart"
ht_noncomm_evolrate_Ab$tissue <- "Heart"
ht_noncomm_evolrate_Nb$tissue <- "Heart"
ht_noncomm_evolrate_On$tissue <- "Heart"
kd_noncomm_evolrate_Mz$tissue <- "Kidney"
kd_noncomm_evolrate_Pn$tissue <- "Kidney"
kd_noncomm_evolrate_Ab$tissue <- "Kidney"
kd_noncomm_evolrate_Nb$tissue <- "Kidney"
kd_noncomm_evolrate_On$tissue <- "Kidney"
ms_noncomm_evolrate_Mz$tissue <- "Muscle"
ms_noncomm_evolrate_Pn$tissue <- "Muscle"
ms_noncomm_evolrate_Ab$tissue <- "Muscle"
ms_noncomm_evolrate_Nb$tissue <- "Muscle"
ms_noncomm_evolrate_On$tissue <- "Muscle"
ts_noncomm_evolrate_Mz$tissue <- "Testis"
ts_noncomm_evolrate_Pn$tissue <- "Testis"
ts_noncomm_evolrate_Ab$tissue <- "Testis"
ts_noncomm_evolrate_Nb$tissue <- "Testis"
ts_noncomm_evolrate_On$tissue <- "Testis"


tissue_noncomm_evolrate_sp <- list(br_noncomm_evolrate_Mz,br_noncomm_evolrate_Pn,br_noncomm_evolrate_Ab,br_noncomm_evolrate_Nb,br_noncomm_evolrate_On,ey_noncomm_evolrate_Mz,ey_noncomm_evolrate_Pn,ey_noncomm_evolrate_Ab,ey_noncomm_evolrate_Nb,ey_noncomm_evolrate_On,ht_noncomm_evolrate_Mz,ht_noncomm_evolrate_Pn,ht_noncomm_evolrate_Ab,ht_noncomm_evolrate_Nb,ht_noncomm_evolrate_On,kd_noncomm_evolrate_Mz,kd_noncomm_evolrate_Pn,kd_noncomm_evolrate_Ab,kd_noncomm_evolrate_Nb,kd_noncomm_evolrate_On,ms_noncomm_evolrate_Mz,ms_noncomm_evolrate_Pn,ms_noncomm_evolrate_Ab,ms_noncomm_evolrate_Nb,ms_noncomm_evolrate_On,ts_noncomm_evolrate_Mz,ts_noncomm_evolrate_Pn,ts_noncomm_evolrate_Ab,ts_noncomm_evolrate_Nb,ts_noncomm_evolrate_On)

for(i in 1:length(tissue_noncomm_evolrate_sp))
{
  tissue_noncomm_evolrate_sp[[i]]$foldiff <- tissue_noncomm_evolrate_sp[[i]]$Promoter/tissue_noncomm_evolrate_sp[[i]]$Fourfold # add in fold diff between promoter and fourfold
  tissue_noncomm_evolrate_sp[[i]] <- tissue_noncomm_evolrate_sp[[i]][with(tissue_noncomm_evolrate_sp[[i]], order(-foldiff)),] # order based on fold diff
}

# create combined files for each tissue, match to gene symbols
br_noncomm_evolrate_combined <- do.call("rbind",tissue_noncomm_evolrate_sp[c(1,2,3,4,5)]) # rbind each tissue element from list
br_noncomm_evolrate_combined2 <- (OGIDs[match(br_noncomm_evolrate_combined$ortho,rownames(OGIDs)),(c(9,11,14))])
colnames(br_noncomm_evolrate_combined2) <- c('Ga_gene','Dr_gene','Hs_gene')
br_noncomm_evolrate_combined3 <- cbind(br_noncomm_evolrate_combined,br_noncomm_evolrate_combined2)
write.table(br_noncomm_evolrate_combined3,file="/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Br_noncomm_rankedevolrate.txt",col.names = T,row.names = F,quote = F,sep = '\t')

ey_noncomm_evolrate_combined <- do.call("rbind",tissue_noncomm_evolrate_sp[c(6,7,8,9,10)]) # rbind each tissue element from list
ey_noncomm_evolrate_combined2 <- (OGIDs[match(ey_noncomm_evolrate_combined$ortho,rownames(OGIDs)),(c(9,11,14))])
colnames(ey_noncomm_evolrate_combined2) <- c('Ga_gene','Dr_gene','Hs_gene')
ey_noncomm_evolrate_combined3 <- cbind(ey_noncomm_evolrate_combined,ey_noncomm_evolrate_combined2)
write.table(ey_noncomm_evolrate_combined3,file="/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/ey_noncomm_rankedevolrate.txt",col.names = T,row.names = F,quote = F,sep = '\t')

ht_noncomm_evolrate_combined <- do.call("rbind",tissue_noncomm_evolrate_sp[c(11,12,13,14,15)]) # rbind each tissue element from list
ht_noncomm_evolrate_combined2 <- (OGIDs[match(ht_noncomm_evolrate_combined$ortho,rownames(OGIDs)),(c(9,11,14))])
colnames(ht_noncomm_evolrate_combined2) <- c('Ga_gene','Dr_gene','Hs_gene')
ht_noncomm_evolrate_combined3 <- cbind(ht_noncomm_evolrate_combined,ht_noncomm_evolrate_combined2)
write.table(ht_noncomm_evolrate_combined3,file="/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/ht_noncomm_rankedevolrate.txt",col.names = T,row.names = F,quote = F,sep = '\t')

kd_noncomm_evolrate_combined <- do.call("rbind",tissue_noncomm_evolrate_sp[c(16,17,18,19,20)]) # rbind each tissue element from list
kd_noncomm_evolrate_combined2 <- (OGIDs[match(kd_noncomm_evolrate_combined$ortho,rownames(OGIDs)),(c(9,11,14))])
colnames(kd_noncomm_evolrate_combined2) <- c('Ga_gene','Dr_gene','Hs_gene')
kd_noncomm_evolrate_combined3 <- cbind(kd_noncomm_evolrate_combined,kd_noncomm_evolrate_combined2)
write.table(kd_noncomm_evolrate_combined3,file="/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/kd_noncomm_rankedevolrate.txt",col.names = T,row.names = F,quote = F,sep = '\t')

ms_noncomm_evolrate_combined <- do.call("rbind",tissue_noncomm_evolrate_sp[c(21,22,23,24,25)]) # rbind each tissue element from list
ms_noncomm_evolrate_combined2 <- (OGIDs[match(ms_noncomm_evolrate_combined$ortho,rownames(OGIDs)),(c(9,11,14))])
colnames(ms_noncomm_evolrate_combined2) <- c('Ga_gene','Dr_gene','Hs_gene')
ms_noncomm_evolrate_combined3 <- cbind(ms_noncomm_evolrate_combined,ms_noncomm_evolrate_combined2)
write.table(ms_noncomm_evolrate_combined3,file="/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/ms_noncomm_rankedevolrate.txt",col.names = T,row.names = F,quote = F,sep = '\t')

ts_noncomm_evolrate_combined <- do.call("rbind",tissue_noncomm_evolrate_sp[c(26,27,28,29,30)]) # rbind each tissue element from list
ts_noncomm_evolrate_combined2 <- (OGIDs[match(ts_noncomm_evolrate_combined$ortho,rownames(OGIDs)),(c(9,11,14))])
colnames(ts_noncomm_evolrate_combined2) <- c('Ga_gene','Dr_gene','Hs_gene')
ts_noncomm_evolrate_combined3 <- cbind(ts_noncomm_evolrate_combined,ts_noncomm_evolrate_combined2)
write.table(ts_noncomm_evolrate_combined3,file="/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/ts_noncomm_rankedevolrate.txt",col.names = T,row.names = F,quote = F,sep = '\t')

#### ##### ##### ##### ##### #### ##### ##### ##### #####

# FigS-R2a - Boxplot of evolutionary rate in promoter regions of 1:1 orthologous cichlid genes with intermediate or narrow (tissue-specific) expression range (tau)

# Plotting evolutionary rate at promoter region to tau to compare tissue specific and intermediate expression (might be better than above conserved/diverged however, will not be tissue-specific)

Ab_tau_finalOGID <- read.table(file="~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Ab_tau_final.txt2")
Mz_tau_finalOGID <- read.table(file="~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Mz_tau_final.txt2")
Pn_tau_finalOGID <- read.table(file="~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Pn_tau_final.txt2")
Nb_tau_finalOGID <- read.table(file="~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/Nb_tau_final.txt2")
On_tau_finalOGID <- read.table(file="~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/On_tau_final.txt2")

Ab_tau_finalOGID$type[Ab_tau_finalOGID$V2 < 0.5] <- "Broad"
Ab_tau_finalOGID$type[Ab_tau_finalOGID$V2 > 0.5 & Ab_tau_finalOGID$V2 < 0.9] <- "Intermediate"
Ab_tau_finalOGID$type[Ab_tau_finalOGID$V2 > 0.9] <- "Narrow"
Ab_tau_finalevolrate <- (Ab_evolrate[match(Ab_tau_finalOGID[,3],Ab_evolrate[,1]),(c(1,2,3,4))])
Ab_tau_finalOGID2 <- cbind(Ab_tau_finalOGID,Ab_tau_finalevolrate[,c(2,3,4)])
Ab_tau_finalOGID3 <- na.omit(Ab_tau_finalOGID2)

Mz_tau_finalOGID$type[Mz_tau_finalOGID$V2 < 0.5] <- "Broad"
Mz_tau_finalOGID$type[Mz_tau_finalOGID$V2 > 0.5 & Mz_tau_finalOGID$V2 < 0.9] <- "Intermediate"
Mz_tau_finalOGID$type[Mz_tau_finalOGID$V2 > 0.9] <- "Narrow"
Mz_tau_finalevolrate <- (Mz_evolrate[match(Mz_tau_finalOGID[,3],Mz_evolrate[,1]),(c(1,2,3,4))])
Mz_tau_finalOGID2 <- cbind(Mz_tau_finalOGID,Mz_tau_finalevolrate[,c(2,3,4)])
Mz_tau_finalOGID3 <- na.omit(Mz_tau_finalOGID2)

Pn_tau_finalOGID$type[Pn_tau_finalOGID$V2 < 0.5] <- "Broad"
Pn_tau_finalOGID$type[Pn_tau_finalOGID$V2 > 0.5 & Pn_tau_finalOGID$V2 < 0.9] <- "Intermediate"
Pn_tau_finalOGID$type[Pn_tau_finalOGID$V2 > 0.9] <- "Narrow"
Pn_tau_finalevolrate <- (Pn_evolrate[match(Pn_tau_finalOGID[,3],Pn_evolrate[,1]),(c(1,2,3,4))])
Pn_tau_finalOGID2 <- cbind(Pn_tau_finalOGID,Pn_tau_finalevolrate[,c(2,3,4)])
Pn_tau_finalOGID3 <- na.omit(Pn_tau_finalOGID2)

Nb_tau_finalOGID$type[Nb_tau_finalOGID$V2 < 0.5] <- "Broad"
Nb_tau_finalOGID$type[Nb_tau_finalOGID$V2 > 0.5 & Nb_tau_finalOGID$V2 < 0.9] <- "Intermediate"
Nb_tau_finalOGID$type[Nb_tau_finalOGID$V2 > 0.9] <- "Narrow"
Nb_tau_finalevolrate <- (Nb_evolrate[match(Nb_tau_finalOGID[,3],Nb_evolrate[,1]),(c(1,2,3,4))])
Nb_tau_finalOGID2 <- cbind(Nb_tau_finalOGID,Nb_tau_finalevolrate[,c(2,3,4)])
Nb_tau_finalOGID3 <- na.omit(Nb_tau_finalOGID2)

On_tau_finalOGID$type[On_tau_finalOGID$V2 < 0.5] <- "Broad"
On_tau_finalOGID$type[On_tau_finalOGID$V2 > 0.5 & On_tau_finalOGID$V2 < 0.9] <- "Intermediate"
On_tau_finalOGID$type[On_tau_finalOGID$V2 > 0.9] <- "Narrow"
On_tau_finalevolrate <- (On_evolrate[match(On_tau_finalOGID[,3],On_evolrate[,1]),(c(1,2,3,4))])
On_tau_finalOGID2 <- cbind(On_tau_finalOGID,On_tau_finalevolrate[,c(2,3,4)])
On_tau_finalOGID3 <- na.omit(On_tau_finalOGID2)

Tau_finalOGID3_combined <- rbind(Mz_tau_finalOGID3,Pn_tau_finalOGID3,Ab_tau_finalOGID3,Nb_tau_finalOGID3,On_tau_finalOGID3)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/tau_vs_promoter_evolrate.tiff', units="in", width=15, height=15, res=300)
ggplot(Tau_finalOGID3_combined, aes(x = type, y = log2(Promoter), fill=factor(type))) + geom_boxplot() + stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species) + labs(title="Tau vs Promoter evolutionary rate") + guides(fill=guide_legend(title="Tau expression range")) + 
  labs (x = "Tau expression range") + labs (y = "log2 (Promoter evolutionary rate)") + theme(legend.position="none")
dev.off()

# Plot the correlation between tau and promoter evolutionary rate, highlighting the outliers

branch_evolrate2_Mzpeak <- subset(branch_evolrate2, Promoter > 0.0313 & Fourfold < 0.001 & Species=="M. zebra")

# only take the top 30 with a high promoter:fourfold fold diff and tau>0.9 and Promoter > 0.0313
Mz_tau_finalOGID3$foldiff <- Mz_tau_finalOGID3$Promoter/Mz_tau_finalOGID3$Fourfold
Mz_tau_finalOGID3a <- subset(Mz_tau_finalOGID3, type=="Narrow" & Promoter > 0.0313)
Mz_tau_finalOGID3_peaks <- Mz_tau_finalOGID3a[with(Mz_tau_finalOGID3a, order(-foldiff)),][1:30,]
 
Pn_tau_finalOGID3$foldiff <- Pn_tau_finalOGID3$Promoter/Pn_tau_finalOGID3$Fourfold
Pn_tau_finalOGID3a <- subset(Pn_tau_finalOGID3, type=="Narrow" & Promoter > 0.0313)
Pn_tau_finalOGID3_peaks <- Pn_tau_finalOGID3a[with(Pn_tau_finalOGID3a, order(-foldiff)),][1:28,]

Ab_tau_finalOGID3$foldiff <- Ab_tau_finalOGID3$Promoter/Ab_tau_finalOGID3$Fourfold
Ab_tau_finalOGID3a <- subset(Ab_tau_finalOGID3, type=="Narrow" & Promoter > 0.0313)
Ab_tau_finalOGID3_peaks <- Ab_tau_finalOGID3a[with(Ab_tau_finalOGID3a, order(-foldiff)),][1:30,]

Nb_tau_finalOGID3$foldiff <- Nb_tau_finalOGID3$Promoter/Nb_tau_finalOGID3$Fourfold
Nb_tau_finalOGID3a <- subset(Nb_tau_finalOGID3, type=="Narrow" & Promoter > 0.0313)
Nb_tau_finalOGID3_peaks <- Nb_tau_finalOGID3a[with(Nb_tau_finalOGID3a, order(-foldiff)),][1:30,]

On_tau_finalOGID3$foldiff <- On_tau_finalOGID3$Promoter/On_tau_finalOGID3$Fourfold
On_tau_finalOGID3a <- subset(On_tau_finalOGID3, type=="Narrow" & Promoter > 0.0313)
On_tau_finalOGID3_peaks <- On_tau_finalOGID3a[with(On_tau_finalOGID3a, order(-foldiff)),][1:30,]

Tau_finalOGID3_combined_peaks <- rbind(Mz_tau_finalOGID3_peaks,Pn_tau_finalOGID3_peaks,Ab_tau_finalOGID3_peaks,Nb_tau_finalOGID3_peaks,On_tau_finalOGID3_peaks)

# match the orthogroups to gene name
Tau_finalOGID3_combined_peaks2 <- (OGIDs[match(Tau_finalOGID3_combined_peaks[,3],row.names(OGIDs)),(c(9,11,14))])
# amend some gene names according to other columns - as the df has factors levels, we need to state this so the factor is outputted, instead of the level 
Tau_finalOGID3_combined_peaks2$gene <- ifelse(Tau_finalOGID3_combined_peaks2$V10=="NULL",yes = levels(Tau_finalOGID3_combined_peaks2$V12)[Tau_finalOGID3_combined_peaks2$V12], no = levels(Tau_finalOGID3_combined_peaks2$V10)[Tau_finalOGID3_combined_peaks2$V10]) # create an extra column that collates

# cbind gene column
Tau_finalOGID3_combined_peaks3 <- cbind(Tau_finalOGID3_combined_peaks,Tau_finalOGID3_combined_peaks2[,4])
colnames(Tau_finalOGID3_combined_peaks3) <- c("cichlid.gene","V2","ortho","type","Fourfold","Promoter","Species","foldiff","gene")

# add the anc_evolrateplot_peaks labels to the plot and offset their labelling
tau_evolrateplot_correlation <- ggplot(Tau_finalOGID3_combined, aes(x = V2, y = log2(Promoter))) + geom_point() +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species) + labs(title="Tau vs Promoter evolutionary rate correlation") + guides(fill=guide_legend(title="Tau expression range")) + 
  labs (x = "Tau expression range") + labs (y = "log2 (Promoter evolutionary rate)") + theme(legend.position="none") + geom_hline(yintercept = -5,linetype="dashed",color="grey") + geom_vline(xintercept = 0.9,linetype="dashed",color="grey") #+ geom_smooth(method=lm)

library(ggrepel)
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EvolutionaryRate/tau_vs_promoter_evolrate_correlation.tiff', units="in", width=15, height=15, res=300)
tau_evolrateplot_correlation + geom_text_repel(data = Tau_finalOGID3_combined_peaks3, 
                                               mapping = aes(label=gene), 
                                               size = 3,
                                               fontface = 'bold', 
                                               color = 'dark grey',
                                               box.padding = unit(0.5, "lines"),
                                               point.padding = unit(0.5, "lines"))
dev.off()


#### ##### ##### ##### ##### #### ##### ##### ##### #####

# FigS-XX - Promoter evolutionary rate phylo and GO enrichment

## Here, we will plot the GO enrichment - note, that when filtered based on q < 0.05, left with no results
# Instead, just plot -log10(p-val < 0.05)
# prepared all files for this from line 5783 in NetworkReconstruction_v4.sh script


# load in the data 
Files1 <- dir(file.path(base_dir_evol))
promsort_GOgrep <- glob2rx("*GOOUTPUT_details_filtered2.txt3") # create a grep pattern to select only specific files
promsort_GOgrep2 <- grep(promsort_GOgrep, Files1) # run the grep
Files2 <- Files1[promsort_GOgrep2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files2)){
  tmp = read.delim(file = paste0(base_dir_evol,Files2[i]),header = F)
  assign(Files2[i], tmp)
}

# these are for individual plots

Pn_topprom_GO <- ggplot(branch1_PN_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3, aes(x = reorder(V2, -V3), y = -log10(V3))) + geom_bar(stat = "identity", fill = "#435e89") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), nudge_y = 0.1) + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)")

Nb_topprom_GO <- ggplot(branch2_NB_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3, aes(x = reorder(V2, -V3), y = -log10(V3))) + geom_bar(stat = "identity", fill = "#435e89") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), nudge_y = 0.1) + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)")

Mz_topprom_GO <- ggplot(branch3_MZ_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3, aes(x = reorder(V2, -V3), y = -log10(V3))) + geom_bar(stat = "identity", fill = "#435e89") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), nudge_y = 0.1) + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)")

On_topprom_GO <- ggplot(branch4_ON_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3, aes(x = reorder(V2, -V3), y = -log10(V3))) + geom_bar(stat = "identity", fill = "#435e89") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), nudge_y = 0.1) + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)")

Ab_topprom_GO <- ggplot(branch5_AB_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3, aes(x = reorder(V2, -V3), y = -log10(V3))) + geom_bar(stat = "identity", fill = "#435e89") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), nudge_y = 0.1) + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)")

MzPn_topprom_GO <- ggplot(branch7_MzPn_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3, aes(x = reorder(V2, -V3), y = -log10(V3))) + geom_bar(stat = "identity", fill = "#435e89") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), nudge_y = 0.1) + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)")

# branch 8 has no enrichment

MzPnAbNb_topprom_GO <- ggplot(branch9_MzPnAbNb_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3, aes(x = reorder(V2, -V3), y = -log10(V3))) + geom_bar(stat = "identity", fill = "#435e89") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), nudge_y = 0.1) + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)")

# this is to create a combined plot that we will use
toppromGO_combined <- rbind(branch1_PN_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3,branch2_NB_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3,branch3_MZ_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3,branch4_ON_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3,branch5_AB_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3,branch7_MzPn_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3,branch9_MzPnAbNb_evolRate.out.promsort.promsort50.GOOUTPUT_details_filtered2.txt3)
toppromGO_combined2 <- gsub("_promsort50", "\\1", toppromGO_combined$V1)
toppromGO_combined$Species <- toppromGO_combined2
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/FigS-R2a_top-prom_evolGO.tiff', units="in", width=12, height=12, res=300)
ggplot(toppromGO_combined, aes(x = reorder(V2, -V3), y = -log10(V3), fill=Species)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_fill_brewer(palette="Set3") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=3, group='V1') + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)")
dev.off()


##### multiplot function #####
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# create plots called p1, p2 and p3 then run with:
# multiplot(p1,p2,p3,cols=3)
##### ##### ##### ##### #####
