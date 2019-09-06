### Plotting GAT enrichment of CNEs and aCNEs across 10kb promoter regions

library(ggplot2)
library(hexbin)
library(reshape2)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/")
base_dir <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/7.CNEs/GATenrichment/"

# Load data (pre-prepared on cluster - q-val < 0.2 (FDR) filtered only)
Files <- dir(file.path(base_dir))
GATgrep <- glob2rx("gatnormed_*.txt") # create a grep pattern to select only specific files
GATgrep2 <- grep(GATgrep, Files) # run the grep
Files2 <- Files[GATgrep2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files2)){
  tmp = read.delim(file = paste0(base_dir,Files2[i]),sep='\t',header = F)
  assign(Files2[i], tmp)
}

# For individual plots and carrying out Mann-Whitney U to test difference between aCNE/hCNE switch/no switch
colnames(gatnormed_Mz.txt) <- c("annotation","fold","log2(fold)","p-val","q-val","type","species")
colnames(gatnormed_Pn.txt) <- c("annotation","fold","log2(fold)","p-val","q-val","type","species")
colnames(gatnormed_Ab.txt) <- c("annotation","fold","log2(fold)","p-val","q-val","type","species")
colnames(gatnormed_Nb.txt) <- c("annotation","fold","log2(fold)","p-val","q-val","type","species")
colnames(gatnormed_On.txt) <- c("annotation","fold","log2(fold)","p-val","q-val","type","species")
Mz_sw_positions <- c('Mz_aCNE_sw','Mz_CNE_sw','Mz_aCNE_nsw','Mz_CNE_nsw') # alter positioning on plot
Pn_sw_positions <- c('Pn_aCNE_sw','Pn_CNE_sw','Pn_aCNE_nsw','Pn_CNE_nsw') # alter positioning on plot
Ab_sw_positions <- c('Ab_aCNE_sw','Ab_CNE_sw','Ab_aCNE_nsw','Ab_CNE_nsw') # alter positioning on plot
Nb_sw_positions <- c('Nb_aCNE_sw','Nb_CNE_sw','Nb_aCNE_nsw','Nb_CNE_nsw') # alter positioning on plot
On_sw_positions <- c('On_aCNE_sw','On_CNE_sw','On_aCNE_nsw','On_CNE_nsw') # alter positioning on plot

ggplot(gatnormed_Mz.txt, aes(x = type, y = log2(fold), fill=factor(type))) + geom_point() +
  scale_x_discrete(limits = Mz_sw_positions) + theme(legend.position="none")
# Wilcoxon rank sum test (Mann-Whitney U-test) with continuity correction
# Mz_sw_grep <- glob2rx("Mz*_sw") # create wildcard grep pattern
# Mz_nosw_grep <- glob2rx("Mz*_nsw") # create wildcard grep pattern
# Mz_sw <- subset(gatnormed_Mz.txt, grepl(Mz_sw_grep, gatnormed_Mz.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# Mz_nosw <- subset(gatnormed_Mz.txt, grepl(Mz_nosw_grep, gatnormed_Mz.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# test the difference between aCNE/CNE sw and aCNE/CNE nsw
# wilcox.test(fold ~ type, data=Mz_sw, paired=FALSE) # W = 0, p-value = 1
# wilcox.test(fold ~ type, data=Mz_nosw, paired=FALSE) # W = 1, p-value = 1
# test the difference between aCNE sw and nsw, CNE sw and nsw
# Mz_aCNE_grep <- glob2rx("Mz_aCNE*") # create wildcard grep pattern
# Mz_CNE_grep <- glob2rx("Mz_CNE*") # create wildcard grep pattern
# Mz_aCNE <- subset(gatnormed_Mz.txt, grepl(Mz_aCNE_grep, gatnormed_Mz.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# Mz_CNE <- subset(gatnormed_Mz.txt, grepl(Mz_CNE_grep, gatnormed_Mz.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# wilcox.test(fold ~ type, data=Mz_aCNE, paired=FALSE) # W = 0, p-value = 1
# wilcox.test(fold ~ type, data=Mz_CNE, paired=FALSE) # W = 0, p-value = 1
# If the p-value is greater than 0.05, then we can accept the hypothesis H0 of statistical equality of the means of two groups.
# null hypothesis: the same (pval>0.05), alternative: there is a difference (pval<0.05) 
# Result:
# There is no difference between the two - the R 'W' output is the same as a 'U' output (sum rank)

ggplot(gatnormed_Pn.txt, aes(x = type, y = log2(fold), fill=factor(type))) + geom_point() +
  scale_x_discrete(limits = Pn_sw_positions) + theme(legend.position="none")
# Pn_sw_grep <- glob2rx("Pn*_sw") # create wildcard grep pattern
# Pn_nosw_grep <- glob2rx("Pn*_nsw") # create wildcard grep pattern
# Pn_sw <- subset(gatnormed_Pn.txt, grepl(Pn_sw_grep, gatnormed_Pn.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# Pn_nosw <- subset(gatnormed_Pn.txt, grepl(Pn_nosw_grep, gatnormed_Pn.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# wilcox.test(fold ~ type, data=Pn_sw, paired=FALSE) # W = 1, p-value = 1
# wilcox.test(fold ~ type, data=Pn_nosw, paired=FALSE) # W = 0, p-value = 1


ggplot(gatnormed_Ab.txt, aes(x = type, y = log2(fold), fill=factor(type))) + geom_point() +
  scale_x_discrete(limits = Ab_sw_positions) + theme(legend.position="none")
# Ab_sw_grep <- glob2rx("Ab*_sw") # create wildcard grep pattern
# Ab_nosw_grep <- glob2rx("Ab*_nsw") # create wildcard grep pattern
# Ab_sw <- subset(gatnormed_Ab.txt, grepl(Ab_sw_grep, gatnormed_Ab.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# Ab_nosw <- subset(gatnormed_Ab.txt, grepl(Ab_nosw_grep, gatnormed_Ab.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# wilcox.test(fold ~ type, data=Ab_sw, paired=FALSE) # W = 873, p-value = 0.0001813
# wilcox.test(fold ~ type, data=Ab_nosw, paired=FALSE) # NA - no data

ggplot(gatnormed_Nb.txt, aes(x = type, y = log2(fold), fill=factor(type))) + geom_point() +
  scale_x_discrete(limits = Nb_sw_positions) + theme(legend.position="none")
# Nb_sw_grep <- glob2rx("Nb*_sw") # create wildcard grep pattern
# Nb_nosw_grep <- glob2rx("Nb*_nsw") # create wildcard grep pattern
# Nb_sw <- subset(gatnormed_Nb.txt, grepl(Nb_sw_grep, gatnormed_Nb.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# Nb_nosw <- subset(gatnormed_Nb.txt, grepl(Nb_nosw_grep, gatnormed_Nb.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# wilcox.test(fold ~ type, data=Nb_sw, paired=FALSE) # W = 480, p-value = 0.003106
# wilcox.test(fold ~ type, data=Nb_nosw, paired=FALSE) # W = 1897, p-value = 3.787e-05

ggplot(gatnormed_On.txt, aes(x = type, y = log2(fold), fill=factor(type))) + geom_point() +
  scale_x_discrete(limits = On_sw_positions) + theme(legend.position="none")
# On_sw_grep <- glob2rx("On*_sw") # create wildcard grep pattern
# On_nosw_grep <- glob2rx("On*_nsw") # create wildcard grep pattern
# On_sw <- subset(gatnormed_On.txt, grepl(On_sw_grep, gatnormed_On.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# On_nosw <- subset(gatnormed_On.txt, grepl(On_nosw_grep, gatnormed_On.txt[[6]]), drop = TRUE) # pull out rows using grep pattern
# wilcox.test(fold ~ type, data=On_sw, paired=FALSE) # W = 1350, p-value = 4.971e-07
# wilcox.test(fold ~ type, data=On_nosw, paired=FALSE) # W = 5824, p-value < 2.2e-16

# sw range: xx
# nsw range: xx

# Generate a combined plot of the GAT enrichment - NOTE: 'Ab_aCNE_nosw','Ab_CNE_nosw' had no enrichment
GAT_combined <- rbind(gatnormed_Mz.txt,gatnormed_Pn.txt,gatnormed_Ab.txt,gatnormed_Nb.txt,gatnormed_On.txt)
colnames(GAT_combined) <- c("annotation","fold","log2(fold)","p-val","q-val","type","species")

# alter positioning on plot - removed 'Ab_aCNE_nosw','Ab_CNE_nosw' since they had no enrichment
# need to use this method instead of scale_x_discrete(limits = positions) defining positions as this does not work with the facet grid
#GAT_combined$type <- factor(GAT_combined$type, levels=c('Mz_aCNE_sw','Mz_CNE_sw','Mz_aCNE_nosw','Mz_CNE_nosw','Pn_aCNE_sw','Pn_CNE_sw','Pn_aCNE_nosw','Pn_CNE_nosw','Ab_aCNE_sw','Ab_CNE_sw','Nb_aCNE_sw','Nb_CNE_sw','Nb_aCNE_nosw','Nb_CNE_nosw','On_aCNE_sw','On_CNE_sw','On_aCNE_nosw','On_CNE_nosw'))
GAT_combined$type <- factor(GAT_combined$type, levels=c('Mz_aCNE_sw','Mz_CNE_sw','Mz_aCNE_nsw','Mz_CNE_nsw','Pn_aCNE_sw','Pn_CNE_sw','Pn_aCNE_nsw','Pn_CNE_nsw','Ab_aCNE_sw','Ab_CNE_sw','Ab_aCNE_nsw','Ab_CNE_nsw','Nb_aCNE_sw','Nb_CNE_sw','Nb_aCNE_nsw','Nb_CNE_nsw','On_aCNE_sw','On_CNE_sw','On_aCNE_nsw','On_CNE_nsw'))

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/CNEs/FigS-R2h_CNE_enrichment_state_change.tiff', units="in", width=18, height=12, res=300)
ggplot(GAT_combined, aes(x = type, y = log2(fold), fill=factor(type))) + geom_point() +
  facet_grid(.~species,scales="free_x", space="free_x") + 
  labs(title="GAT enrichment/association of CNEs/aCNEs in switching/non-switching promoter sequences") + 
  guides(fill=guide_legend(title="Type")) + labs (x = "Species, regulatory element and type of promoter overlap") + 
  labs (y = "log2 (fold enrichment)") + theme(legend.position="none")
dev.off()
