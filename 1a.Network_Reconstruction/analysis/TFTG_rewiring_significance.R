### Testing the significance of candidate gene re-wiring scores (DyNet degree-corrected) vs non-candidates

library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(dplyr)
library(broom)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/May2019_ReviewerComments/1.3.candidate_rewiring/")

TFTG1to1rw = read.table("TFTG1to1_DyNet_rewiring_CandsMarked_MoreThanMean.txt", header=TRUE)
TFTGallrw = read.table("TFTGalledges_DyNet_rewiring_CandsMarked_MoreThanMean.txt", header=TRUE)

# Directional test - A one-sided t-test or kstest would be fine.
#  if p-value > 0.05 accept the null hypothesis and no difference between both data sets i.e., both data sets are equal.
#  if p-value < 0.05 reject the null hypothesis and significant difference between both data sets i.e., both data sets are not equal.

# NOTE: for the ks.test - two-sample case alternative="greater" includes distributions for which x is stochastically smaller than y (the CDF of x lies above and hence to the left of that for y)
# This is in contrast to t.test or wilcox.test and thus, for your test, use alternative="less"

# Testing significant difference in degree-corrected DyNet rewiring score of (using a Kolmogorov–Smirnov (KS) test)
# X Cand greater than (>) Y Non-Cand - but to the right of the distribution so alt="less"
# A. all re-wired candidates vs all re-wired genes

TFTG1to1rwCandX <- subset(TFTG1to1rw, Candidate==1)
TFTG1to1rwNonCandY <- subset(TFTG1to1rw, Candidate==0)
ks.test(TFTG1to1rwCandX$Dn.Score_.degree_corrected, TFTG1to1rwNonCandY$Dn.Score_.degree_corrected, alternative = "less") # p-value = 0.9852 - no significant difference between 31/90 candidates and the other 1to1 re-wired genes

TFTGallrwCandX <- subset(TFTGallrw, Candidate==1)
TFTGallrwNonCandY <- subset(TFTGallrw, Candidate==0)
ks.test(TFTGallrwCandX$Dn.Score_.degree_corrected, TFTGallrwNonCandY$Dn.Score_.degree_corrected, alternative = "less") # p-value = 0.0002021 - significant difference between 89/90 candidates and the other re-wired genes

# X Cand and More than Mean greater than (>) Y Cand Less than Mean - but to the right of the distribution so alt="less"
# B. all re-wired candidates with rewiring scores above the mean (0.14 and 0.23) vs all re-wired genes

TFTG1to1rwMoreThanMeanX <- subset(TFTG1to1rw, MoreThanMean==1)
TFTG1to1rwNonMoreThanMeanY <- subset(TFTG1to1rw, MoreThanMean==0)
ks.test(TFTG1to1rwMoreThanMeanX$Dn.Score_.degree_corrected, TFTG1to1rwNonMoreThanMeanY$Dn.Score_.degree_corrected, alternative = "less") # p-value = 0.000639 - significant difference between 9/31 candidates and the other 1to1 re-wired genes

TFTGallrwMoreThanMeanX <- subset(TFTGallrw, MoreThanMean==1)
TFTGallrwNonMoreThanMeanY <- subset(TFTGallrw, MoreThanMean==0)
ks.test(TFTGallrwMoreThanMeanX$Dn.Score_.degree_corrected, TFTGallrwNonMoreThanMeanY$Dn.Score_.degree_corrected, alternative = "less") # p-value = 5.976e-14 - no significant difference between 60/89 candidates and the other re-wired genes


# Mann-Whitney Test Intepretation (DO NOT USE THIS AS IT IS NOT A DIRECTIONAL TEST)
#  if p-value > 0.05 accept the null hypothesis and no difference between both data sets i.e., both data sets are equal.
#  if p-value < 0.05 reject the null hypothesis and significant difference between both data sets i.e., both data sets are not equal.

# Testing significant difference in degree-corrected DyNet rewiring score of 

# A. all re-wired candidates vs all re-wired genes
wilcox.test(Dn.Score_.degree_corrected. ~ Candidate, data=TFTG1to1rw, paired=FALSE) # Mann-Whitney p-value = 0.002633 - significant difference between 31/90 candidates and the other 1to1 re-wired genes
wilcox.test(Dn.Score_.degree_corrected. ~ Candidate, data=TFTGallrw, paired=FALSE) # Mann-Whitney p-value = 0.01564 - no significant difference between 89/90 candidates and the other re-wired genes

# B. all re-wired candidates with rewiring scores above the mean (0.14 and 0.23) vs all re-wired genes
wilcox.test(Dn.Score_.degree_corrected. ~ MoreThanMean, data=TFTG1to1rw, paired=FALSE) # Mann-Whitney p-value = 0.001365 - significant difference between 9/31 candidates and the other 1to1 re-wired genes
wilcox.test(Dn.Score_.degree_corrected. ~ MoreThanMean, data=TFTGallrw, paired=FALSE) # Mann-Whitney p-value = 8.031e-15 - significant difference between 60/89 candidates and the other re-wired genes

mean(TFTG1to1rw$Dn.Score_.degree_corrected.) # mean = 0.1657865 = 0.17
mean(TFTGallrw$Dn.Score_.degree_corrected.) # mean = 0.2274561 = 0.23

# Plotting new Fig. 2c and Fig. 2f to include candidates and their significance on box plots
TFTG1to1rw_scores <- melt(TFTG1to1rw,id.vars=c('Candidate','Gene_Symbol','MoreThanMean'), measure.vars=c('DyNet_Rewiring_.Dn.score.','Dn.Score_.degree_corrected.'))
TFTGallrw_scores <- melt(TFTGallrw,id.vars=c('Candidate','Gene_Symbol','MoreThanMean'), measure.vars=c('DyNet_Rewiring_.Dn.score.','Dn.Score_.degree_corrected.'))
TFTG1to1rw_scores_DC <- melt(TFTG1to1rw,id.vars=c('Candidate','Gene_Symbol','MoreThanMean'), measure.vars=('Dn.Score_.degree_corrected.'))
TFTGallrw_scores_DC <- melt(TFTGallrw,id.vars=c('Candidate','Gene_Symbol','MoreThanMean'), measure.vars=('Dn.Score_.degree_corrected.'))


# Change 0 and 1 from Candidates to label facets
TFTG1to1rw_scores$Candidate[TFTG1to1rw_scores$Candidate == "0"] <- "Non-candidate"
TFTG1to1rw_scores$Candidate[TFTG1to1rw_scores$Candidate == "1"] <- "Candidate"
TFTG1to1rw_scores_DC$Candidate[TFTG1to1rw_scores_DC$Candidate == "0"] <- "Non-candidate"
TFTG1to1rw_scores_DC$Candidate[TFTG1to1rw_scores_DC$Candidate == "1"] <- "Candidate"
TFTGallrw_scores$Candidate[TFTGallrw_scores$Candidate == "0"] <- "Non-candidate"
TFTGallrw_scores$Candidate[TFTGallrw_scores$Candidate == "1"] <- "Candidate"
TFTGallrw_scores_DC$Candidate[TFTGallrw_scores_DC$Candidate == "0"] <- "Non-candidate"
TFTGallrw_scores_DC$Candidate[TFTGallrw_scores_DC$Candidate == "1"] <- "Candidate"

# Change the variables to be distinct in each acse of 1to1 and All dfs
TFTG1to1rw_scores_DC$variable <- "1-to-1 TF-TG Dn rewiring score (degree-corrected)"
TFTGallrw_scores_DC$variable <- "All TF-TG Dn rewiring score (degree-corrected)"

# rbind 1to1 and all degree-corrected scores only
TFTGrw_DC <- rbind(TFTG1to1rw_scores_DC,TFTGallrw_scores_DC)

# demarcating particular gene examples (as per mention in manuscript + supp. text) to include as labels on plot
TFTGrw_DC$CandExample <- ifelse(TFTGrw_DC$Gene_Symbol == "rh2" & TFTGrw_DC$MoreThanMean == 1 | TFTGrw_DC$Gene_Symbol == "opn1sw1" & TFTGrw_DC$MoreThanMean == 1 | TFTGrw_DC$Gene_Symbol == "opn1sw2" & TFTGrw_DC$MoreThanMean == 1 | TFTGrw_DC$Gene_Symbol == "rho" & TFTGrw_DC$MoreThanMean == 1 | TFTGrw_DC$Gene_Symbol == "gdf10b" & TFTGrw_DC$MoreThanMean == 1 | TFTGrw_DC$Gene_Symbol == "actr1b" & TFTGrw_DC$MoreThanMean == 1 | TFTGrw_DC$Gene_Symbol == "pax6a" & TFTGrw_DC$MoreThanMean == 1 | TFTGrw_DC$Gene_Symbol == "draxin" & TFTGrw_DC$MoreThanMean == 1 | TFTGrw_DC$Gene_Symbol == "cntn3a.1" & TFTGrw_DC$MoreThanMean == 1 | TFTGrw_DC$Gene_Symbol == "irx1b" & TFTGrw_DC$MoreThanMean == 1, 1, 0) # only demarcate the opsins (rh2, opn1sw1 > sws1, opn1sw2 > sws2, rho all more than mean), gdf10b, actr1b, draxin, cntn3a.1 > cntn4, pax6a, irx1b
sum(TFTGrw_DC$CandExample == 1) # sanity check that 11 (actually 10examples but gdf10b appears in both sets) have been demarcated - OK!
TFTGrw_DC$Gene_Symbol <- gsub("opn1sw1", "sws1", TFTGrw_DC$Gene_Symbol) # rename some gene symbols - opn1sw1 > sws1
TFTGrw_DC$Gene_Symbol <- gsub("opn1sw2", "sws2", TFTGrw_DC$Gene_Symbol) # rename some gene symbols - opn1sw2 > sws2
TFTGrw_DC$Gene_Symbol <- gsub("cntn3a.1", "cntn4", TFTGrw_DC$Gene_Symbol) # rename some gene symbols - cntn3a.1 > cntn4

# demarcating candidates that also have a value more than the mean (this is to demarcate those dots on the plot with different colour)
TFTGrw_DC$CandMoreThanMean <- ifelse(TFTGrw_DC$Candidate == "Candidate" & TFTGrw_DC$MoreThanMean == 1, 1, 0) # assign a 1 for candidate and higher than mean, and 0 for any gene lower than mean
sum(TFTGrw_DC$CandMoreThanMean == 1) # sanity check that 69 have been demarcated - OK!

# Violin plot of Dn-rewiring score (degree-corrected) from 1to1 and all TF-TG edges with example candidate genes demarcated
tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/Fig2c.Violinplot_Dn_degree-corrected_rewiring_1to1_allgenes_and_Candidates.tiff', units="in", width=13, height=8, res=300)
ggplot(TFTGrw_DC, aes(x = factor(variable), y = value, label = Gene_Symbol, fill=factor(variable))) + geom_violin() +
  theme_bw() + theme(panel.grid = element_blank()) + guides(fill=guide_legend(title="Source of score")) + 
  labs (x = "Transcription Factor (TF) - Target Gene (TG) network edges") + labs (y = expression(DyNet~(D[n])~rewiring~("degree"-corrected)~score)) + theme(legend.position="none") + ylim(-7, 15) +
  scale_y_continuous(breaks = seq(0, 0.3, by=0.025), limits=c(0,0.3)) +
  scale_x_discrete(labels=c('6,844 1-to-1 orthologs', '14,590 1-to-1 and non-1-to-1 orthologs')) +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black")) +
  geom_text_repel(
    data          = subset(TFTGrw_DC, CandExample == 1),   
    nudge_y       = 36 - subset(TFTGrw_DC, CandExample == 1)$CandExample,
    segment.size  = 0.5,
    segment.color = "grey30",
    direction     = "x",
    box.padding = 0.5,
    size = 6,
    point.padding = 1.2,
    colour = "#CC3300"
  ) + 
  geom_point(size = ifelse(TFTGrw_DC$CandMoreThanMean == 1, 4, 0.8), color = ifelse(TFTGrw_DC$CandMoreThanMean == 1, "orange", "grey20")) +
  stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 5, fill= "white") +
  scale_fill_manual(values=c("#009E73", "#0072B2")) # scale_fill_brewer(palette="Dark2")
dev.off()
# + labs(title= expression(Violin~plot~DyNet~(D[n])~rewiring~("degree"-corrected)~score~from~1-to-1~and~all~TF-TG~edges~with~demarcated~candidate~genes)) 

# Figure legend:
# Violin plots of DyNet (Dn) rewiring (degree-corrected) score from 6,844 1-to-1 orthologs in 215,810 TF-TG network edges (green, left violin) and 14,590 1-to-1 and non-1-to-1 orthologs in 843,168 TF-TG network edges (blue, right violin). Mean rewiring score shown within each violin plot (white diamond). Degree-corrected rewiring score shown for non-candidate genes (grey dots through centre) and candidate morphogenetic trait genes (orange dots) with rewiring scores higher than the mean, and selected candidate examples are demarcated within.


## 3.3 Significance of “more predicted unique TF regulators of sws1 in M. zebra (38 TFs) as compared to N. brichardi (6 TFs) (Fig. 3a)”
# ranked the edges based on q-value and ordering in terms of significance. 
# Carry out Mann-Whitney test of significant difference in unique regulator motif prediction q-vals.

MzNb_sws1Uniq = read.table("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/May2019_ReviewerComments/3.3.sws1_significance/MzNb-Edge_Attributes_Collated4c.coexpr_promONLY_FDRrank.sws1_ONLY_uniqRANKED.txt", header=FALSE)

# Mann-Whitney Test Intepretation:
#  if p-value > 0.05 accept the null hypothesis and no difference between both data sets i.e., both data sets are equal.
#  if p-value < 0.05 reject the null hypothesis and significant difference between both data sets i.e., both data sets are not equal.

# Testing significant difference in unique edges of Mz and Nb sws1
wilcox.test(V7 ~ V9, data=MzNb_sws1Uniq, paired=FALSE) # Mann-Whitney p-value = 0.1731 - no significant difference between q-value unique regulatory edges of Mz and Nb sws1 networks

Mz_sws1Uniq <- MzNb_sws1Uniq[MzNb_sws1Uniq$V9 == "Mz",]
Mz_sws1Uniq2 <- Mz_sws1Uniq[!is.na(Mz_sws1Uniq$V7), ]
Nb_sws1Uniq <- MzNb_sws1Uniq[MzNb_sws1Uniq$V9 == "Nb",]
Nb_sws1Uniq2 <- Nb_sws1Uniq[!is.na(Nb_sws1Uniq$V7), ]
                     
mean(Mz_sws1Uniq2$V7) # mean = 0.02348811 = 0.02
mean(Nb_sws1Uniq2$V7) # mean = 0.03040333 = 0.03

MzNb_sws1Uniq_qval <- melt(MzNb_sws1Uniq,id.vars=c('V1','V9'), measure.vars='V7')
MzNb_sws1Uniq_qval3 <- distinct(MzNb_sws1Uniq_qval,V1, .keep_all= TRUE) # Remove duplicate rows of the dataframe based on ID
# MzNb_sws1Uniq_qval3 <- MzNb_sws1Uniq_qval2[!is.na(MzNb_sws1Uniq_qval2$value), ]
Mz_sws1Uniq_qval3 <- subset(MzNb_sws1Uniq_qval3, V9 == "Mz")
Nb_sws1Uniq_qval3 <- subset(MzNb_sws1Uniq_qval3, V9 == "Nb")
mean(Mz_sws1Uniq_qval3$value, na.rm=TRUE) # mean = 0.01404464 = 0.01
mean(Nb_sws1Uniq_qval3$value, na.rm=TRUE) # mean = 0.026766 = 0.03
wilcox.test(value ~ V9, data=MzNb_sws1Uniq_qval3, paired=FALSE) # Mann-Whitney p-value = 0.0741 - no significant difference between q-value unique regulatory edges of Mz and Nb sws1 networks

# demarcate your examples of Mz unique that you included in main text - nr2c2 (OG3997_0:OG11366_0); rxrb (OG1991_0:OG11366_0); atrx (OG4786_1:OG11366_0)
MzNb_sws1Uniq_qval3$GeneSymbol[MzNb_sws1Uniq_qval3$V1 == "OG3997_0:OG11366_0"] <- "NR2C2"
MzNb_sws1Uniq_qval3$GeneSymbol[MzNb_sws1Uniq_qval3$V1 == "OG1991_0:OG11366_0"] <- "RXRB"
MzNb_sws1Uniq_qval3$GeneSymbol[MzNb_sws1Uniq_qval3$V1 == "OG4786_1:OG11366_0"] <- "ATRX"
MzNb_sws1Uniq_qval3$MzCandExample <- ifelse(MzNb_sws1Uniq_qval3$V9 == "Mz" & MzNb_sws1Uniq_qval3$GeneSymbol != "NA", 1, 0)
MzNb_sws1Uniq_qval3$NbCandExample <- ifelse(MzNb_sws1Uniq_qval3$V9 == "Nb" & MzNb_sws1Uniq_qval3$GeneSymbol != "NA", 1, 0)
MzNb_sws1Uniq_qval3$Mean[MzNb_sws1Uniq_qval3$V9 == "Mz"] <- "0.01404464" # since you need to resize and recolour based on two separate means, add a mean column
MzNb_sws1Uniq_qval3$Mean[MzNb_sws1Uniq_qval3$V9 == "Nb"] <- "0.026766" # since you need to resize and recolour based on two separate means, add a mean column

# Unique TF-sws1 edges

# Violin plot of the q-val of TF-sws1 edges
tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/3.3.Violinplot_MzNb_uniqueTF-sws1_edges.tiff', units="in", width=38, height=22, res=200)
ggplot(MzNb_sws1Uniq_qval3, aes(x = factor(V9), y = value, label = GeneSymbol, fill=factor(V9))) + geom_violin() +
  theme_bw() + theme(panel.grid = element_blank()) + guides(fill=guide_legend(title="Source of score")) + 
  labs (x = expression(bold(Unique~TF-bolditalic("sws1")~edges))) + labs (y = "Significance of edge (FDR < 0.05)") + theme(legend.position="none") + ylim(-7, 15) +
  scale_y_continuous(breaks = seq(0, 0.05, by=0.005), limits=c(0,0.05)) +
  scale_x_discrete(labels=expression(italic(N.~brichardi), italic(M.~zebra))) +
  theme(axis.text.x = element_text(size=68),
        axis.text.y = element_text(size=54),  
        axis.title.x = element_text(size=86, margin = margin(t = 20, r = 0, b = 0, l = 0), face="bold"),
        axis.title.y = element_text(size=86, margin = margin(t = 0, r = 20, b = 0, l = 0), face="bold")) +
  geom_text_repel(
    data          = subset(MzNb_sws1Uniq_qval3, MzCandExample == 1),   
    nudge_x       = 15 - subset(MzNb_sws1Uniq_qval3, MzCandExample == 1)$MzCandExample,
    segment.size  = 0.5,
    segment.color = "grey30",
    direction     = "x",
    box.padding = 0.5,
    size = 20,
    point.padding = 1.2,
    colour = "#CC3300"
  ) +
  geom_text_repel(
    data          = subset(MzNb_sws1Uniq_qval3, NbCandExample == 1),   
    nudge_x       = -5 - subset(MzNb_sws1Uniq_qval3, NbCandExample == 1)$NbCandExample,
    segment.size  = 0.5,
    segment.color = "grey30",
    direction     = "x",
    box.padding = 0.5,
    size = 20,
    point.padding = 1.2,
    colour = "#CC3300"
  ) +
  scale_fill_manual(values=c("#009E73", "#0072B2")) +
  geom_point(size = ifelse(MzNb_sws1Uniq_qval3$value < MzNb_sws1Uniq_qval3$Mean, 24, 6), color = ifelse(MzNb_sws1Uniq_qval3$value < MzNb_sws1Uniq_qval3$Mean, "orange", "grey20")) +
  stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 14, fill= "white")
dev.off()

sum(MzNb_sws1Uniq_qval3$V9 == "Mz" & MzNb_sws1Uniq_qval3$value < MzNb_sws1Uniq_qval3$Mean, na.rm=TRUE) # 18 TFs < mean
sum(MzNb_sws1Uniq_qval3$V9 == "Nb" & MzNb_sws1Uniq_qval3$value < MzNb_sws1Uniq_qval3$Mean, na.rm=TRUE) # 2 TFs < mean
