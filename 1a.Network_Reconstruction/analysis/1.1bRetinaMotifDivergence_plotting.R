### 1.1b. Testing the divergence of retina/lens-related TFs motifs in eye-expressed genes

library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(dplyr)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/May2019_ReviewerComments/1.1b.RetinaMotifDivergence/")

# create a scatter plot of the enrichment of motifs in module 1 gene promoters VS the eye expression level in each species, colouring based on whether module 1 or not

Cons = read.table("Mz_Pn_Ab_Nb_On_Conserved.Final.txt", header=TRUE, row.names=NULL)
NotCons = read.table("Mz_Pn_Ab_Nb_On_NotConserved.Final.txt", header=TRUE, sep = "\t", row.names=NULL)

# Conserved TF motifs in module 1 gene promoters
p1 <- ggplot(Cons, aes(x=Mz_eye_expr_log.x.1., y=Mz_qval, color=factor(Mz_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=Cons$TF, size=2.5, vjust = -1) +
  labs(x = "M. zebra log(x+1) eye expression", y = "M. zebra motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p1a <- p1 + xlim(0,6.3) + ylim(0.047,0) 

p2 <- ggplot(Cons, aes(x=Pn_eye_expr_log.x.1., y=Pn_qval, color=factor(Pn_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=Cons$TF, size=2.5, vjust = -1) +
  labs(x = "P. nyererei log(x+1) eye expression", y = "P. nyererei motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p2a <- p2 + xlim(0,6.3) + ylim(0.047,0)

p3 <- ggplot(Cons, aes(x=Ab_eye_expr_log.x.1., y=Ab_qval, color=factor(Ab_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=Cons$TF, size=2.5, vjust = -1) +
  labs(x = "A. burtoni log(x+1) eye expression", y = "A. burtoni motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p3a <- p3 + xlim(0,6.3) + ylim(0.047,0)

p4 <- ggplot(Cons, aes(x=Nb_eye_expr_log.x.1., y=Nb_qval, color=factor(Nb_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=Cons$TF, size=2.5, vjust = -1) +
  labs(x = "N. brichardi log(x+1) eye expression", y = "N. brichardi motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p4a <- p4 + xlim(0,6.3) + ylim(0.047,0)

p5 <- ggplot(Cons, aes(x=On_eye_expr_log.x.1., y=On_qval, color=factor(On_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=Cons$TF, size=2.5, vjust = -1) +
  labs(x = "O. niloticus log(x+1) eye expression", y = "O. niloticus motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p5a <- p5 + xlim(0,6.3) + ylim(0.047,0)

# # extract legend function
# g_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
# 
# myConslegend <- g_legend(p1) - can't use this since not all legends are the same

# tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R1f.retinalTF_Cons-enrichment_expr_5sp.tiff', units="in", width=24, height=16, res=300)
# grid.arrange(p1 + theme(legend.position="none"),p2 + theme(legend.position="none"),p3 + theme(legend.position="none"),p4 + theme(legend.position="none"),p5 + theme(legend.position="none"),myConslegend,nrow=1)
# dev.off()

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R1fA.retinalTF_Cons-enrichment_expr_5sp.tiff', units="in", width=40, height=6, res=300)
grid.arrange(p1a,p2a,p3a,p4a,p5a,nrow=1)
dev.off()

# Conserved in 5 species TF motif enrichment in module 1 gene promoters
# CRX is not intersting here
# Others that are include SOX13, MEOX1, NKX6-2




# Not Conserved TF motifs in module 1 gene promoters

p6 <- ggplot(NotCons, aes(x=Mz_eye_expr_log.x.1., y=Mz_qval, color=factor(Mz_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=NotCons$TF, size=2.5, vjust = -1) +
  labs(x = "M. zebra log(x+1) eye expression", y = "M. zebra motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p6a <- p6 + xlim(0,6.3) + ylim(0.05,0) 

p7 <- ggplot(NotCons, aes(x=Pn_eye_expr_log.x.1., y=Pn_qval, color=factor(Pn_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=NotCons$TF, size=2.5, vjust = -1) +
  labs(x = "P. nyererei log(x+1) eye expression", y = "P. nyererei motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p7a <- p7 + xlim(0,6.3) + ylim(0.05,0)

p8 <- ggplot(NotCons, aes(x=Ab_eye_expr_log.x.1., y=Ab_qval, color=factor(Ab_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=NotCons$TF, size=2.5, vjust = -1) +
  labs(x = "A. burtoni log(x+1) eye expression", y = "A. burtoni motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p8a <- p8 + xlim(0,6.3) + ylim(0.05,0)

p9 <- ggplot(NotCons, aes(x=Nb_eye_expr_log.x.1., y=Nb_qval, color=factor(Nb_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=NotCons$TF, size=2.5, vjust = -1) +
  labs(x = "N. brichardi log(x+1) eye expression", y = "N. brichardi motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p9a <- p9 + xlim(0,6.3) + ylim(0.05,0)

p10 <- ggplot(NotCons, aes(x=On_eye_expr_log.x.1., y=On_qval, color=factor(On_TF_module.assign))) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=NotCons$TF, size=2.5, vjust = -1) +
  labs(x = "O. niloticus log(x+1) eye expression", y = "O. niloticus motif enrichment (q-value < 0.05) in module 1 gene promoters", color = "Module assignment")
p10a <- p10 + xlim(0,6.3) + ylim(0.05,0)

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R1fB.retinalTF_NotCons-enrichment_expr_5sp.tiff', units="in", width=40, height=6, res=300)
grid.arrange(p6a,p7a,p8a,p9a,p10a,nrow=1)
dev.off()

# then, calculate the variance in enrichment and eye expression for each cons and not cons TF and plot
ConsMVar <- read.table("Mz_Pn_Ab_Nb_On_Conserved.FinalMotifVar.txt", header=TRUE, row.names=1)
ConsEVar <- read.table("Mz_Pn_Ab_Nb_On_Conserved.FinalExprVar.txt", header=TRUE, row.names=1)

NotConsMVar <- read.table("Mz_Pn_Ab_Nb_On_NotConserved.FinalMotifVar.txt", header=TRUE, row.names=1)
NotConsEVar <- read.table("Mz_Pn_Ab_Nb_On_NotConserved.FinalExprVar.txt", header=TRUE, row.names=1)

ConsMVar1 <- data.frame(apply(ConsMVar, 1, var, na.rm=TRUE))
ConsEVar1 <- data.frame(apply(ConsEVar, 1, var, na.rm=TRUE))
ConsEMVar <- cbind(ConsMVar1,ConsEVar1)
colnames(ConsEMVar) <- c("M-var", "E-var")

NotConsMVar1 <- data.frame(apply(NotConsMVar, 1, var, na.rm=TRUE))
NotConsEVar1 <- data.frame(apply(NotConsEVar, 1, var, na.rm=TRUE))
NotConsEMVar <- cbind(NotConsMVar1,NotConsEVar1)
colnames(NotConsEMVar) <- c("M-var", "E-var")

p11 <- ggplot(ConsEMVar, aes(x=`E-var`, y=`M-var`)) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=row.names(ConsEMVar), size=2.5, vjust = -1) +
  labs(x = "Variance of log(x+1) eye expression\nin all five species", y = "Variance of motif enrichment (q-value < 0.05)\nin module 1 gene promoters of all five species", color = "Module assignment") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) + 
  scale_y_reverse()

p12 <- ggplot(NotConsEMVar, aes(x=`E-var`, y=`M-var`)) +
  geom_point() + 
  theme_minimal() +
  geom_text(label=row.names(NotConsEMVar), size=2.5, vjust = -1) +
  labs(x = "Variance of log(x+1) eye expression\nin all five species", y = "Variance of motif enrichment (q-value < 0.05)\nin module 1 gene promoters of all five species", color = "Module assignment") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) + 
  scale_y_reverse()

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R1fC.retinalTF_Variance-enrichment_expr_5sp.tiff', units="in", width=16, height=6, res=300)
grid.arrange(p11,p12,nrow=1)
dev.off()

#### NOTE: As of 13/08/19, ran the analysis below but Sushmita suggested the above instead.


Mz = read.table("mz-JASPAR_fimo.switched.expr.map2", header=TRUE, sep = "\t")
Pn = read.table("pn-JASPAR_fimo.switched.expr.map2", header=TRUE, sep = "\t")
Ab = read.table("ab-JASPAR_fimo.switched.expr.map2", header=TRUE, sep = "\t")
Nb = read.table("nb-JASPAR_fimo.switched.expr.map2", header=TRUE, sep = "\t")
On = read.table("on-JASPAR_fimo.switched.expr.map2", header=TRUE, sep = "\t")
names(Mz)[21] <- "Co-expression"
names(Mz)[22] <- "Species"
names(Pn)[21] <- "Co-expression"
names(Pn)[22] <- "Species"
names(Ab)[21] <- "Co-expression"
names(Ab)[22] <- "Species"
names(Nb)[21] <- "Co-expression"
names(Nb)[22] <- "Species"
names(On)[21] <- "Co-expression"
names(On)[22] <- "Species"

# Mann-Whitney Test Intepretation:
#  if p-value > 0.05 accept the null hypothesis and no difference between both data sets i.e., both data sets are equal.
#  if p-value < 0.05 reject the null hypothesis and significant difference between both data sets i.e., both data sets are not equal.

# Testing significant difference in TF raw expression in TF-TG both eye-expressed (no state change) or TG eye-expressed and TF not
wilcox.test(as.numeric(geneA_Eye_RawExpr) ~ `Co-expression`, data=Mz, paired=FALSE) # Mann-Whitney p-value < 2.2e-16 - significant difference of M. zebra TF raw expression between state-changed and non-state-changed TFs in interactions with eye-expressed TGs 
wilcox.test(as.numeric(geneA_Eye_RawExpr) ~ `Co-expression`, data=Pn, paired=FALSE) # Mann-Whitney p-value = 1.438e-05 - significant difference of P. nyererei TF raw expression between state-changed and non-state-changed TFs in interactions with eye-expressed TGs 
wilcox.test(as.numeric(geneA_Eye_RawExpr) ~ `Co-expression`, data=Ab, paired=FALSE) # Mann-Whitney p-value = 0.4531 - no significant difference of A. burtoni TF raw expression between state-changed and non-state-changed TFs in interactions with eye-expressed TGs 
wilcox.test(as.numeric(geneA_Eye_RawExpr) ~ `Co-expression`, data=Nb, paired=FALSE) # Mann-Whitney p-value = 8.257e-08 - significant difference of N. brichardi TF raw expression between state-changed and non-state-changed TFs in interactions with eye-expressed TGs 
wilcox.test(as.numeric(geneA_Eye_RawExpr) ~ `Co-expression`, data=On, paired=FALSE) # Mann-Whitney p-value < 2.2e-16 - significant difference of O. niloticus TF raw expression between state-changed and non-state-changed TFs in interactions with eye-expressed TGs 


# Violin plot of "TF raw expression in predicted TF-TG interactions of eye-expressed TGs"
# x-axis "TF-TG eye co-expressed (no state-change)" and "TG eye-expressed and TF state-changed", main label "TF-TG interactions of eye-expressed TGs"; y-axis "raw expression value"

# create one unified table
All <- rbind(Mz,Pn,Ab,Nb,On)

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/1.1b.Expression_divergence_of_retinalTFs.tiff', units="in", width=13, height=8, res=300)
ggplot(All, aes(x = factor(`Co-expression`), y = log2(as.numeric(geneA_Eye_RawExpr)), fill=factor(`Co-expression`))) + geom_violin() + 
  stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~Species) + 
  # labs(title="Decay of TFBSs and TF raw expression in predicted TF-TG interactions of eye-expressed TGs") + 
  labs (x = "TF-TG interactions of eye-expressed TGs") + labs (y = "log2 (TF raw expression value)") + theme(legend.position="none") +
  scale_fill_manual(values=c("#009E73", "#0072B2")) +
  scale_x_discrete(labels=c('TF-TG eye co-expressed\n(no state-change)', 'TG eye-expressed\nand TF state-changed')) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=16),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black"))
dev.off()

# In each species (except On), there is a higher mean log2(TF raw expression value) of TF-TG eye co-expressed (no state-change) than TG eye-expressed and TF state-changed away from eye-expressed (module 1)

### Ran further analyses to test state-changed TFs along the phylogeny e.g. On vs NbAbPnMz; Nb vs AbPnMz; Ab vs PnMz and their change in raw expression
# created a column that is 'geneA_OGID:geneB_OGID' then awk match (adding string for match or not), then do if-else match on modules and add string for match or not
# 1. For each state-changed TFs e.g. in On, characterize how many orthologous TF-TGs exist in (A) module 1 (of MzPnAbNb) vs (B) TF(state-changed On)-TGs(module1)
# 2. If (A) is higher than (B) then gain of sites in module 1
# 3. If (A) is lower than (B) then a loss of sites in module 1 (decay)
# 4. Difference in TF raw expression (species TF minus ancestral species TF e.g. Mz vs On) $7-$28: negative value indicates a loss of expression with/without state change, positive indicates a gain of expression with/without state-change

# Violin plot of "Decay of TFBSs and TF raw expression in predicted TF-TG interactions of eye-expressed TGs"
# x-axis  "TF-TG interactions of module 1 (eye-expressed) TGs"; y-axis "Difference in log2(TF raw expression value) of species against O. niloticus"

MzPnAbNb_On_SC_Expr = read.table("MzPnAbNb-fromONandOrth-JASPAR_fimo.expr.map2", header=FALSE, sep = "\t")
# A = Orthologous TF-TG, TF-TG module 1 (eye co-expressed) and TF state-changed from O. niloticus > For On it is "TF-TG module 1 (eye co-expressed)"
# B = Orthologous TF-TG, TG in module 1 (eye-expressed) and TF not state-changed to module 1 from O. niloticus > For On it is "TG in module 1 (eye-expressed) and TF not in module 1"
MzPnAbNb_On_SC_Expr$V45 = factor(MzPnAbNb_On_SC_Expr$V43, levels=c('M. zebra','P. nyererei','A. burtoni','N. brichardi','O. niloticus'))

Mz_On_SC_Expr <- subset(MzPnAbNb_On_SC_Expr, V43 == "M. zebra" | V43 == "O. niloticus")
Mz_On_SC_Expr <- Mz_On_SC_Expr[!(Mz_On_SC_Expr$V43 == "M. zebra" & Mz_On_SC_Expr$V42 == "B" ),] 
wilcox.test(as.numeric(V7) ~ V43, data=Mz_On_SC_Expr, paired=FALSE) # Mann-Whitney p-value = 6.546e-05 - significant difference of raw expression value of orthologous TFs (state-changed) from O. niloticus non module 1 to M. zebra module 1
Pn_On_SC_Expr <- subset(MzPnAbNb_On_SC_Expr, V43 == "P. nyererei" | V43 == "O. niloticus")
Pn_On_SC_Expr <- Pn_On_SC_Expr[!(Pn_On_SC_Expr$V43 == "P. nyererei" & Pn_On_SC_Expr$V42 == "B" ),] 
wilcox.test(as.numeric(V7) ~ V43, data=Pn_On_SC_Expr, paired=FALSE) # Mann-Whitney p-value < 2.2e-16 - significant difference of raw expression value of orthologous TFs (state-changed) from O. niloticus non module 1 to P. nyererei module 1
Ab_On_SC_Expr <- subset(MzPnAbNb_On_SC_Expr, V43 == "A. burtoni" | V43 == "O. niloticus")
wilcox.test(as.numeric(V7) ~ V43, data=Ab_On_SC_Expr, paired=FALSE) # Mann-Whitney p-value < 2.2e-16 - significant difference of raw expression value of orthologous TFs (state-changed) from O. niloticus non module 1 to A. burtoni module 1
Nb_On_SC_Expr <- subset(MzPnAbNb_On_SC_Expr, V43 == "N. brichardi" | V43 == "O. niloticus")
wilcox.test(as.numeric(V7) ~ V43, data=Nb_On_SC_Expr, paired=FALSE) # Mann-Whitney p-value = 0.0144 - significant difference of raw expression value of orthologous TFs (state-changed) from O. niloticus non module 1 to N. brichardi module 1

tiff('/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Nature_submission/Reviewer_comments/New_Figs/FigS-R1f.State-changing_of_retinalTFs_and_expression_fromON.tiff', units="in", width=13, height=8, res=300)
ggplot(MzPnAbNb_On_SC_Expr, aes(x = factor(V42), y = log2(as.numeric(V7)), fill=factor(V42))) + geom_violin() + 
  stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~V45) + 
  # labs(title="Decay of TFBSs and TF raw expression in predicted TF-TG interactions of eye-expressed TGs") + 
  # guides(fill=guide_legend(title="State-changing")) + 
  labs (x = "TF-TG interactions of module 1 (eye-expressed) TGs") + labs (y = "log2(TF eye expression value)") + theme(legend.position="none") +
  scale_fill_manual(values=c("#009E73", "#0072B2")) +
  # scale_x_discrete(labels=c('TF-TG eye co-expressed\n(no state-change)', 'TG eye-expressed\nand TF state-changed')) +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20, colour = "black", face="italic"))
dev.off()

## Conclusions:
# The On file represents all the edges and TFs from which the orthologous TFs in other species have switched their module assignment
# In Mz and Ab, a large proportion (912/1524 = 51% and 1161/1524 = 65%) of TFs and their TFBSs have state-changed (from O. niloticus) and are now found as regulators of module 1 (eye-expressed genes)
# Examples include CRX(>sws1, >otx5, >wnt9b), TBX4(>) and TFE3(>)
# CRX is retinal, TBX4 is fin bud, TFE3 is immune system activation
# In Pn and Nb, the proportion is lower (28% and 25%)
# Examples include CRX(>valopb)
# Scenario B (Orthologous TF-TG, TG in module 1 (eye-expressed) and TF not state-changed to module 1 from O. niloticus > For On it is "TG in module 1 (eye-expressed) and TF not in module 1")
  # exists 159 (Mz) and 405 (Pn) and 0 times in Ab and Nb - this means there are very few instances of orthologous TF-TGs where the TG is in module 1 and switched TFs (from On) are moving to another module
# This means that for orthologous TF-TG interactions, TFs favour switching to module 1 > this is supported by a higher mean TF expression value in Mz and Pn (A vs B)
# Using the Mann-Whitney test, we show that there is a significant change in eye expression value of orthologous TFs (state-changed) from O. niloticus non module 1 to other four species module 1
# This means that there is a signficant change in TF eye expression value associated with state-changes to module 1 from other modules in O. niloticus
# Furthermore, the mean expression value for TFs switched to module 1 (scenario A) in Mz, Pn and Ab is higher than the expression of non-module 1 TFs in orthologous TF-TG interactions of O. niloticus
# This indicating that eye-expressed genes are potentially under the regulatory control of state-changed TFs from O. niloticus    


# The below plot shows the difference in log2(TF raw expression value) but is not a meaningful plot so do not use
# ggplot(MzPnAbNb_On_SC_Expr, aes(x = factor(V42), y = as.numeric(V44), fill=factor(V42))) + geom_violin() + 
#   stat_summary(fun.y = "mean", geom = "point", shape= 23, size= 3, fill= "white") +
#   theme_bw() + theme(panel.grid = element_blank()) + facet_grid(.~V45) + 
#   # labs(title="Decay of TFBSs and TF raw expression in predicted TF-TG interactions of eye-expressed TGs") + 
#   # guides(fill=guide_legend(title="State-changing")) + 
#   labs (x = "TF-TG interactions of module 1 (eye-expressed) TGs") + labs (y = "Difference in log2(TF raw expression value) of species against O. niloticus") + theme(legend.position="none") +
#   scale_fill_manual(values=c("#009E73", "#0072B2")) +
#   # scale_x_discrete(labels=c('TF-TG eye co-expressed\n(no state-change)', 'TG eye-expressed\nand TF state-changed')) +
#   theme(axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16),  
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=20),
#         strip.text.x = element_text(size = 20, colour = "black"))


