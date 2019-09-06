### Tissue specificity of gene expression across modules
library(ggplot2)
library(hexbin)
library(reshape2)


# tau is a measure for tissue specificity of gene expression as described by Yanai and colleagues³.
# Usage: Summarise tissue replicates using mean() or similar. For it to function correctly, NA's must be 0. For a gene (rows) by tissue (columns) matrix, apply then as:
# ³Yanai et al. 2004: Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification. Bioinformatics 21:650-659.

# First, load in the module gene expression (as they are log transformed and normalised)
setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/")
base_dir2 <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/"

# load expression data - this is of the modules with orthgroups mapped
# Here, the data has been log transformed as: log(x+1), where x is the raw FPKM, and "log" is the natural log, and then each gene has been normalized to have mean zero.
Files4 <- dir(file.path(base_dir2))
modexprgrepOGID <- glob2rx("*-exprtab_a.txt") # create a grep pattern to select only specific files
modexprgrepOGID2 <- grep(modexprgrepOGID, Files4) # run the grep
Files5 <- Files4[modexprgrepOGID2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files5)){
  tmp = read.delim(file = paste0(base_dir2,Files5[i]),header = T, row.names = 1)
  assign(Files5[i], tmp)
}

# change colnames
colnames(`Ab-exprtab_a.txt`) <- c('Br','Ey','Ht','Kd','Ms','Ts')
colnames(`Mz-exprtab_a.txt`) <- c('Br','Ey','Ht','Kd','Ms','Ts')
colnames(`Pn-exprtab_a.txt`) <- c('Br','Ey','Ht','Kd','Ms','Ts')
colnames(`Nb-exprtab_a.txt`) <- c('Br','Ey','Ht','Kd','Ms','Ts')
colnames(`On-exprtab_a.txt`) <- c('Br','Ey','Ht','Kd','Ms','Ts')

# create a copy of dataframes to work on
Ab_exprtab_a_tau1 <- `Ab-exprtab_a.txt`
Mz_exprtab_a_tau1 <- `Mz-exprtab_a.txt`
Pn_exprtab_a_tau1 <- `Pn-exprtab_a.txt`
Nb_exprtab_a_tau1 <- `Nb-exprtab_a.txt`
On_exprtab_a_tau1 <- `On-exprtab_a.txt`

# tau does not handle values <0, so for all genes with expression <1, set to O
Ab_exprtab_a_tau1[Ab_exprtab_a_tau1 < 1] <- 0
Mz_exprtab_a_tau1[Mz_exprtab_a_tau1 < 1] <- 0
Pn_exprtab_a_tau1[Pn_exprtab_a_tau1 < 1] <- 0
Nb_exprtab_a_tau1[Nb_exprtab_a_tau1 < 1] <- 0
On_exprtab_a_tau1[On_exprtab_a_tau1 < 1] <- 0


# Change all NA to O
Ab_exprtab_a_tau1[is.na(Ab_exprtab_a_tau1)] <- 0
Mz_exprtab_a_tau1[is.na(Mz_exprtab_a_tau1)] <- 0
Pn_exprtab_a_tau1[is.na(Pn_exprtab_a_tau1)] <- 0
Nb_exprtab_a_tau1[is.na(Nb_exprtab_a_tau1)] <- 0
On_exprtab_a_tau1[is.na(On_exprtab_a_tau1)] <- 0

# store tau function
tau<-function(x){
  if(any(is.na(x))) stop('NA\'s need to be 0.')
  if(any(x<0)) stop('Negative input values not permitted. Maybe data is log transformed?')
  t<-sum(1-x/max(x))/(length(x)-1)
}

# calculate tau (breadth of expression) across all tissues for each gene - usage 'tau1 <- apply(dataset, 1, tau)'
Ab_exprtab_a_tau2 <- apply(Ab_exprtab_a_tau1, 1, tau)
Mz_exprtab_a_tau2 <- apply(Mz_exprtab_a_tau1, 1, tau)
Pn_exprtab_a_tau2 <- apply(Pn_exprtab_a_tau1, 1, tau)
Nb_exprtab_a_tau2 <- apply(Nb_exprtab_a_tau1, 1, tau)
On_exprtab_a_tau2 <- apply(On_exprtab_a_tau1, 1, tau)

# The values of τ vary from 0 to 1; ubiquitous or broad expr (τ ≤ 0.5); intermediate expr (0.5 < τ < 0.9); and tissue-specific or narrow expr (τ ≥ 0.9)

# Plotting - how does expression pattern correlate with co-expression clustering (switching vs non-switched) and protein coding evolution?

# A. split the τ results of each species into switching and non-switching genes (only 1:1 orthologs)

# read in the switching gene files
Files4 <- dir(file.path(base_dir2))
switchgenes <- glob2rx("OGIDS.txt5-clusterassign-1to1_*vsall") # create a grep pattern to select only specific files
switchgenes2 <- grep(switchgenes, Files4) # run the grep
Files6 <- Files4[switchgenes2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files6)){
  tmp = read.delim(file = paste0(base_dir2,Files6[i]),header = F, row.names = 1)
  assign(Files6[i], tmp)
}

# read in the no-switching gene files
Files4 <- dir(file.path(base_dir2))
noswitchgenes <- glob2rx("OGIDS.txt5-clusterassign-1to1_*vsallns") # create a grep pattern to select only specific files
noswitchgenes2 <- grep(noswitchgenes, Files4) # run the grep
Files7 <- Files4[noswitchgenes2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files7)){
  tmp = read.delim(file = paste0(base_dir2,Files7[i]),header = F, row.names = 1)
  assign(Files7[i], tmp)
}

# pull out the genes only in each switch or noswitch list
OGIDS1to1Absw <- `OGIDS.txt5-clusterassign-1to1_Abvsall`[,3]
OGIDS1to1Mzsw <- `OGIDS.txt5-clusterassign-1to1_Mzvsall`[,1]
OGIDS1to1Pnsw <- `OGIDS.txt5-clusterassign-1to1_Pnvsall`[,2]
OGIDS1to1Nbsw <- `OGIDS.txt5-clusterassign-1to1_Nbvsall`[,4]
OGIDS1to1Onsw <- `OGIDS.txt5-clusterassign-1to1_Onvsall`[,5]
OGIDS1to1Abns <- `OGIDS.txt5-clusterassign-1to1_Abvsallns`[,3]
OGIDS1to1Mzns <- `OGIDS.txt5-clusterassign-1to1_Mzvsallns`[,1]
OGIDS1to1Pnns <- `OGIDS.txt5-clusterassign-1to1_Pnvsallns`[,2]
OGIDS1to1Nbns <- `OGIDS.txt5-clusterassign-1to1_Nbvsallns`[,4]
OGIDS1to1Onns <- `OGIDS.txt5-clusterassign-1to1_Onvsallns`[,5]

# subset the tau results according to switch or noswitch
AbTau.nam <- names(Ab_exprtab_a_tau2)
AbTau.nam.sw <- AbTau.nam %in% OGIDS1to1Absw
AbTau.nam.sw2 <- Ab_exprtab_a_tau2[AbTau.nam.sw] # switched
AbTau.nam.ns <- AbTau.nam %in% OGIDS1to1Abns
AbTau.nam.ns2 <- Ab_exprtab_a_tau2[AbTau.nam.ns] # no switch

MzTau.nam <- names(Mz_exprtab_a_tau2)
MzTau.nam.sw <- MzTau.nam %in% OGIDS1to1Mzsw
MzTau.nam.sw2 <- Mz_exprtab_a_tau2[MzTau.nam.sw] # switched
MzTau.nam.ns <- MzTau.nam %in% OGIDS1to1Mzns
MzTau.nam.ns2 <- Mz_exprtab_a_tau2[MzTau.nam.ns] # no switch

PnTau.nam <- names(Pn_exprtab_a_tau2)
PnTau.nam.sw <- PnTau.nam %in% OGIDS1to1Pnsw
PnTau.nam.sw2 <- Pn_exprtab_a_tau2[PnTau.nam.sw] # switched
PnTau.nam.ns <- PnTau.nam %in% OGIDS1to1Pnns
PnTau.nam.ns2 <- Pn_exprtab_a_tau2[PnTau.nam.ns] # no switch

NbTau.nam <- names(Nb_exprtab_a_tau2)
NbTau.nam.sw <- NbTau.nam %in% OGIDS1to1Nbsw
NbTau.nam.sw2 <- Nb_exprtab_a_tau2[NbTau.nam.sw] # switched
NbTau.nam.ns <- NbTau.nam %in% OGIDS1to1Nbns
NbTau.nam.ns2 <- Nb_exprtab_a_tau2[NbTau.nam.ns] # no switch

OnTau.nam <- names(On_exprtab_a_tau2)
OnTau.nam.sw <- OnTau.nam %in% OGIDS1to1Onsw
OnTau.nam.sw2 <- On_exprtab_a_tau2[OnTau.nam.sw] # switched
OnTau.nam.ns <- OnTau.nam %in% OGIDS1to1Onns
OnTau.nam.ns2 <- On_exprtab_a_tau2[OnTau.nam.ns] # no switch

# B. plot the switching vs non-switching against tau first

AbTau.sw <- as.data.frame(AbTau.nam.sw2)
AbTau.sw <- cbind(AbTau.sw,type='switch')
colnames(AbTau.sw) <- c('tau','type')
AbTau.ns <- as.data.frame(AbTau.nam.ns2)
AbTau.ns <- cbind(AbTau.ns,type='noswitch')
colnames(AbTau.ns) <- c('tau','type')
AbTau.sw_ns <- rbind(AbTau.sw,AbTau.ns) # create a single dataframe containing tau values and the switching/no switching property
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/A.tau_sw_ns/Ab_tau_sw_ns.tiff', units="in", width=7, height=7, res=300)
ggplot(AbTau.sw_ns, aes(x = type, y = tau,fill=type)) + geom_violin() +
  theme_bw() + theme(panel.grid = element_blank()) + labs (x = "A. burtoni module switching/no switching") + labs (y = "Tissue specificity of gene expression (tau, τ)") + theme(legend.position="none") +
  geom_hline(yintercept = 0.9,linetype="dashed") + annotate("text", x=2.59, y=0.91, hjust=1, label="tissue-specific (τ ≥ 0.9)", color = "black") +
  geom_hline(yintercept = 0.5,linetype="dashed") + annotate("text", x=2.59, y=0.51, hjust=1, label="intermediate expression (0.5 < τ < 0.9)", color = "black") +
  annotate("text", x=2.59, y=0.49, hjust=1, label="broad expression (τ ≤ 0.5)", color = "black")
dev.off()

MzTau.sw <- as.data.frame(MzTau.nam.sw2)
MzTau.sw <- cbind(MzTau.sw,type='switch')
colnames(MzTau.sw) <- c('tau','type')
MzTau.ns <- as.data.frame(MzTau.nam.ns2)
MzTau.ns <- cbind(MzTau.ns,type='noswitch')
colnames(MzTau.ns) <- c('tau','type')
MzTau.sw_ns <- rbind(MzTau.sw,MzTau.ns) # create a single dataframe containing tau values and the switching/no switching property
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/A.tau_sw_ns/Mz_tau_sw_ns.tiff', units="in", width=7, height=7, res=300)
ggplot(MzTau.sw_ns, aes(x = type, y = tau,fill=type)) + geom_violin() +
  theme_bw() + theme(panel.grid = element_blank()) + labs (x = "M. zebra module switching/no switching") + labs (y = "Tissue specificity of gene expression (tau, τ)") + theme(legend.position="none") +
  geom_hline(yintercept = 0.9,linetype="dashed") + annotate("text", x=2.59, y=0.91, hjust=1, label="tissue-specific (τ ≥ 0.9)", color = "black") +
  geom_hline(yintercept = 0.5,linetype="dashed") + annotate("text", x=2.59, y=0.51, hjust=1, label="intermediate expression (0.5 < τ < 0.9)", color = "black") +
  annotate("text", x=2.59, y=0.49, hjust=1, label="broad expression (τ ≤ 0.5)", color = "black")
dev.off()

PnTau.sw <- as.data.frame(PnTau.nam.sw2)
PnTau.sw <- cbind(PnTau.sw,type='switch')
colnames(PnTau.sw) <- c('tau','type')
PnTau.ns <- as.data.frame(PnTau.nam.ns2)
PnTau.ns <- cbind(PnTau.ns,type='noswitch')
colnames(PnTau.ns) <- c('tau','type')
PnTau.sw_ns <- rbind(PnTau.sw,PnTau.ns) # create a single dataframe containing tau values and the switching/no switching property
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/A.tau_sw_ns/Pn_tau_sw_ns.tiff', units="in", width=7, height=7, res=300)
ggplot(PnTau.sw_ns, aes(x = type, y = tau,fill=type)) + geom_violin() +
  theme_bw() + theme(panel.grid = element_blank()) + labs (x = "P. nyererei module switching/no switching") + labs (y = "Tissue specificity of gene expression (tau, τ)") + theme(legend.position="none") +
  geom_hline(yintercept = 0.9,linetype="dashed") + annotate("text", x=2.59, y=0.91, hjust=1, label="tissue-specific (τ ≥ 0.9)", color = "black") +
  geom_hline(yintercept = 0.5,linetype="dashed") + annotate("text", x=2.59, y=0.51, hjust=1, label="intermediate expression (0.5 < τ < 0.9)", color = "black") +
  annotate("text", x=2.59, y=0.49, hjust=1, label="broad expression (τ ≤ 0.5)", color = "black")
dev.off()

NbTau.sw <- as.data.frame(NbTau.nam.sw2)
NbTau.sw <- cbind(NbTau.sw,type='switch')
colnames(NbTau.sw) <- c('tau','type')
NbTau.ns <- as.data.frame(NbTau.nam.ns2)
NbTau.ns <- cbind(NbTau.ns,type='noswitch')
colnames(NbTau.ns) <- c('tau','type')
NbTau.sw_ns <- rbind(NbTau.sw,NbTau.ns) # create a single dataframe containing tau values and the switching/no switching property
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/A.tau_sw_ns/Nb_tau_sw_ns.tiff', units="in", width=7, height=7, res=300)
ggplot(NbTau.sw_ns, aes(x = type, y = tau,fill=type)) + geom_violin() +
  theme_bw() + theme(panel.grid = element_blank()) + labs (x = "N. brichardi module switching/no switching") + labs (y = "Tissue specificity of gene expression (tau, τ)") + theme(legend.position="none") +
  geom_hline(yintercept = 0.9,linetype="dashed") + annotate("text", x=2.59, y=0.91, hjust=1, label="tissue-specific (τ ≥ 0.9)", color = "black") +
  geom_hline(yintercept = 0.5,linetype="dashed") + annotate("text", x=2.59, y=0.51, hjust=1, label="intermediate expression (0.5 < τ < 0.9)", color = "black") +
  annotate("text", x=2.59, y=0.49, hjust=1, label="broad expression (τ ≤ 0.5)", color = "black")
dev.off()

OnTau.sw <- as.data.frame(OnTau.nam.sw2)
OnTau.sw <- cbind(OnTau.sw,type='switch')
colnames(OnTau.sw) <- c('tau','type')
OnTau.ns <- as.data.frame(OnTau.nam.ns2)
OnTau.ns <- cbind(OnTau.ns,type='noswitch')
colnames(OnTau.ns) <- c('tau','type')
OnTau.sw_ns <- rbind(OnTau.sw,OnTau.ns) # create a single dataframe containing tau values and the switching/no switching property
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/A.tau_sw_ns/On_tau_sw_ns.tiff', units="in", width=7, height=7, res=300)
ggplot(OnTau.sw_ns, aes(x = type, y = tau,fill=type)) + geom_violin() +
  theme_bw() + theme(panel.grid = element_blank()) + labs (x = "O. niloticus module switching/no switching") + labs (y = "Tissue specificity of gene expression (tau, τ)") + theme(legend.position="none") +
  geom_hline(yintercept = 0.9,linetype="dashed") + annotate("text", x=2.59, y=0.91, hjust=1, label="tissue-specific (τ ≥ 0.9)", color = "black") +
  geom_hline(yintercept = 0.5,linetype="dashed") + annotate("text", x=2.59, y=0.51, hjust=1, label="intermediate expression (0.5 < τ < 0.9)", color = "black") +
  annotate("text", x=2.59, y=0.49, hjust=1, label="broad expression (τ ≤ 0.5)", color = "black")
dev.off()

# Do a combined plot
MzTau.sw_ns2 <- cbind(MzTau.sw_ns,Species='M. zebra')
PnTau.sw_ns2 <- cbind(PnTau.sw_ns,Species='P. nyererei')
AbTau.sw_ns2 <- cbind(AbTau.sw_ns,Species='A. burtoni')
NbTau.sw_ns2 <- cbind(NbTau.sw_ns,Species='N. brichardi')
OnTau.sw_ns2 <- cbind(OnTau.sw_ns,Species='O. niloticus')
combined_Tau.sw_ns <- rbind(MzTau.sw_ns2,PnTau.sw_ns2,AbTau.sw_ns2,NbTau.sw_ns2,OnTau.sw_ns2)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/A.tau_sw_ns/FigR1e_combined_tau_sw_ns.tiff', units="in", width=25, height=15, res=300)
ggplot(combined_Tau.sw_ns, aes(x = type, y = tau,fill=type)) + geom_violin() +
  theme_bw() + theme(panel.grid = element_blank()) + facet_grid(. ~Species) + labs (x = "Module state-change/no state-change") + labs (y = "Tissue specificity of gene expression (tau, τ)") + theme(legend.position="none") +
  geom_hline(yintercept = 0.9,linetype="dashed") + annotate("text", x=2.59, y=0.91, hjust=1, label="tissue-specific (τ ≥ 0.9)", color = "black") +
  geom_hline(yintercept = 0.5,linetype="dashed") + annotate("text", x=2.59, y=0.51, hjust=1, label="intermediate expression (0.5 < τ < 0.9)", color = "black") +
  annotate("text", x=2.59, y=0.49, hjust=1, label="broad expression (τ ≤ 0.5)", color = "black")
dev.off()

# C. Plot the tau of switching vs non-switching against omega (dN/dS) [TO DO]


# D. Break the above down into class of gene (TFs and duplicates)


# E. Plot tau against evolutionary rate at promoter regions (to see if this pattern is better than using conserved/divergent expressed from single-tissue clustering)
# Will do this in another script 'evolRate_figs.R' so here, just write out the tau data frames
Ab_tau_final <- as.data.frame(Ab_exprtab_a_tau2)
Mz_tau_final <- as.data.frame(Mz_exprtab_a_tau2)
Pn_tau_final <- as.data.frame(Pn_exprtab_a_tau2)
Nb_tau_final <- as.data.frame(Nb_exprtab_a_tau2)
On_tau_final <- as.data.frame(On_exprtab_a_tau2)

write.table(Ab_tau_final,file="Ab_tau_final.txt",row.names = T, col.names = F,quote=F,sep = "\t")
write.table(Mz_tau_final,file="Mz_tau_final.txt",row.names = T, col.names = F,quote=F,sep = "\t")
write.table(Pn_tau_final,file="Pn_tau_final.txt",row.names = T, col.names = F,quote=F,sep = "\t")
write.table(Nb_tau_final,file="Nb_tau_final.txt",row.names = T, col.names = F,quote=F,sep = "\t")
write.table(On_tau_final,file="On_tau_final.txt",row.names = T, col.names = F,quote=F,sep = "\t")

######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~

### Calculating expression divergence in tissues between species - using Euclidian distance

# Test on pairwise comparison
# The distance between samples i and j can be written as
# dist(i,j)=(Yi−Yj)⊤(Yi−Yj)
# with Yi and Yj columns i and j. This result can be very convenient in practice as computations can be made much faster using matrix multiplication.

# In this formula, the expression data xi and yi are subtracted directly from each other.
# We should therefore make sure that the expression data are properly normalized when using the Euclidean distance, for example by converting the measured gene expression levels to log-ratios (done when we use the module expression data)
# Unlike the correlation-based distance functions, the Euclidean distance takes the magnitude of the expression data into account. It, therefore, preserves more information about the data and may be preferable.

# EUCLIDIAN DISTANCE NEEDS TO BE BASED ON 1:1 ORTHOLOGS ONLY (USE ROW NAMES AS ORTHOGROUP) - DONE
# Using the module expression data as this has been log transformed (log(x+1)) and normalised so that each species shows mean 0 - DONE
# This methodology (log transformation and normalisation) only stretches and shifts expression values, it does not alter the shape of the distribution.
# Plot the Distribution of expression level for each species (y axis, Density 0.0 to 0.5; x-axis, Species expression -3 to +3) > are they normally distributed? - Yes, DONE
# Then calculate the euclidian distance of each orthogroup between the species for the different tissues
# These can be plotted, species pairwise against evolutionary rate (of the whole tree for each orthogroup at promoters and dN/dS?) - You can sum the distances and maybe log2 to plot
# Multispecies expression divergence can then be a sum of all the values for each orthgroup, plotted against the (overall rate of that orthgroup in the whole tree?)


# First, load in the module gene expression (as they are log transformed and normalised)
base_dir2 <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/"

# load expression data - this is of the modules, 1to1 only (6844 orthgroups)
Files4 <- dir(file.path(base_dir2))
modexprgrepOGID1to1 <- glob2rx("*-exprtab_a-OGID3-1to1.txt") # create a grep pattern to select only specific files
modexprgrepOGID1to1_2 <- grep(modexprgrepOGID1to1, Files4) # run the grep
Files7 <- Files4[modexprgrepOGID1to1_2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files7)){
  tmp = read.delim(file = paste0(base_dir2,Files7[i]),header = F, row.names = 1)
  assign(Files7[i], tmp)
}

# change colnames
colnames(`Ab-exprtab_a-OGID3-1to1.txt`) <- c('Ab_gene', 'Ab_br', 'Ab_ey','Ab_ht','Ab_kd','Ab_ms','Ab_ts')
colnames(`Mz-exprtab_a-OGID3-1to1.txt`) <- c('Mz_gene', 'Mz_br', 'Mz_ey','Mz_ht','Mz_kd','Mz_ms','Mz_ts')
colnames(`Pn-exprtab_a-OGID3-1to1.txt`) <- c('Pn_gene', 'Pn_br', 'Pn_ey','Pn_ht','Pn_kd','Pn_ms','Pn_ts')
colnames(`Nb-exprtab_a-OGID3-1to1.txt`) <- c('Nb_gene', 'Nb_br', 'Nb_ey','Nb_ht','Nb_kd','Nb_ms','Nb_ts')
colnames(`On-exprtab_a-OGID3-1to1.txt`) <- c('On_gene', 'On_br', 'On_ey','On_ht','On_kd','On_ms','On_ts')

# Remove genename for downstream analysis
Ab_exprtab_a_OGID3_1to1.txt2 <- `Ab-exprtab_a-OGID3-1to1.txt`[,-1]
Mz_exprtab_a_OGID3_1to1.txt2 <- `Mz-exprtab_a-OGID3-1to1.txt`[,-1]
Pn_exprtab_a_OGID3_1to1.txt2 <- `Pn-exprtab_a-OGID3-1to1.txt`[,-1]
Nb_exprtab_a_OGID3_1to1.txt2 <- `Nb-exprtab_a-OGID3-1to1.txt`[,-1]
On_exprtab_a_OGID3_1to1.txt2 <- `On-exprtab_a-OGID3-1to1.txt`[,-1]

# cbind all data and then separate by tissue
OGIDS1to1oudata<-cbind.data.frame(`Mz-exprtab_a-OGID3-1to1.txt`,`Pn-exprtab_a-OGID3-1to1.txt`,`Ab-exprtab_a-OGID3-1to1.txt`,`Nb-exprtab_a-OGID3-1to1.txt`,`On-exprtab_a-OGID3-1to1.txt`) # since the order is the same we can just cbind
Br_cols <- seq(2, 35, by = 7) # increment by 7 from col2 to select tissue columns
Br_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Br_cols])
Ey_cols <- seq(3, 35, by = 7) # increment by 7 from col3 to select tissue columns
Ey_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Ey_cols])
Ht_cols <- seq(4, 35, by = 7) # increment by 7 from col4 to select tissue columns
Ht_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Ht_cols])
Kd_cols <- seq(5, 35, by = 7) # increment by 7 from col5 to select tissue columns
Kd_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Kd_cols])
Ms_cols <- seq(6, 35, by = 7) # increment by 7 from col6 to select tissue columns
Ms_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Ms_cols])
Ts_cols <- seq(7, 35, by = 7) # increment by 7 from col7 to select tissue columns
Ts_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Ts_cols])

colnames(Br_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Ey_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Ht_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Kd_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Ms_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Ts_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')


## Create density plots for the distribution of expression level in each species

# A. burtoni
Ab_exprtab_a_OGID3_1to1_l <- melt(as.matrix(Ab_exprtab_a_OGID3_1to1.txt2)) # since we want to capture the rownames, melt as matrix
# Basic histogram with black outline, white fill
ggplot(Ab_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(Ab_exprtab_a_OGID3_1to1_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Normaldist_expr/Aburtoni_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(Ab_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(A.~burtoni)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous A. burtoni genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.

# M. zebra
Mz_exprtab_a_OGID3_1to1_l <- melt(as.matrix(Mz_exprtab_a_OGID3_1to1.txt2)) # since we want to capture the rownames, melt as matrix
# Basic histogram with black outline, white fill
ggplot(Mz_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(Mz_exprtab_a_OGID3_1to1_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Normaldist_expr/Mzebra_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(Mz_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(M.~zebra)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous M. zebra genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.


# P. nyererei
Pn_exprtab_a_OGID3_1to1_l <- melt(as.matrix(Pn_exprtab_a_OGID3_1to1.txt2)) # since we want to capture the rownames, melt as matrix
# Basic histogram with black outline, white fill
ggplot(Pn_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(Pn_exprtab_a_OGID3_1to1_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Normaldist_expr/Pnyererei_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(Pn_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(P.~nyererei)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous P. nyererei genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.

# N. brichardi
Nb_exprtab_a_OGID3_1to1_l <- melt(as.matrix(Nb_exprtab_a_OGID3_1to1.txt2)) # since we want to capture the rownames, melt as matrix
# Basic histogram with black outline, white fill
ggplot(Nb_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(Nb_exprtab_a_OGID3_1to1_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Normaldist_expr/Nbrichardi_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(Nb_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(N.~brichardi)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous N. brichardi genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.

# O. niloticus
On_exprtab_a_OGID3_1to1_l <- melt(as.matrix(On_exprtab_a_OGID3_1to1.txt2)) # since we want to capture the rownames, melt as matrix
# Basic histogram with black outline, white fill
ggplot(On_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(On_exprtab_a_OGID3_1to1_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/Normaldist_expr/Oniloticus_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(On_exprtab_a_OGID3_1to1_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(O.~niloticus)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous O. niloticus genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.

## Calculate the euclidian distance of each orthogroup between the species means for the different tissues

# cbind all the species expression files with tissue colheaders to compute pairwise distances
All_exprtab_a_OGID3_1to1 <- cbind(Mz_exprtab_a_OGID3_1to1.txt2,Pn_exprtab_a_OGID3_1to1.txt2,Ab_exprtab_a_OGID3_1to1.txt2,Nb_exprtab_a_OGID3_1to1.txt2,On_exprtab_a_OGID3_1to1.txt2)

# Now to compute all the distances at once, we have the function dist.
# If we were interested in the distance between samples, we would transpose the matrix d <- dist(t(All_exprtab_a_OGID3_1to1))

# 1. We want to compute pairwise distance across species for each tissue

# load in the whole tree evolutionary rate at promoter regions - this is for plotting
Tree_prom_evolrate <- read.delim(file = "/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/New_090817/WholeTree_4fold_Prom-MAFextract.evolrate.txt2", header = F, row.names = 1)
colnames(Tree_prom_evolrate) <- c('4fold', 'prom')

# Load in the OGIDs to cbind the genes
OGIDs <- read.delim(file = "/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5", header = F, row.names = 1)

## creating plots of expression divergence vs evolutionary rate

# Brain
Br_exprtab_a_OGID3_1to1cols <- seq(1, 30, by = 6) # increment by 6 from col1 to select tissue columns
Br_exprtab_a_OGID3_1to1 <- as.matrix(All_exprtab_a_OGID3_1to1[,Br_exprtab_a_OGID3_1to1cols])

## To run dist on each row save each row in a list
Br_exprtab_a_OGID3_1to1.list <- setNames(split(Br_exprtab_a_OGID3_1to1, seq(nrow(Br_exprtab_a_OGID3_1to1))), rownames(Br_exprtab_a_OGID3_1to1))
Br_exprtab_a_OGID3_1to1.list.ED <- t(sapply(Br_exprtab_a_OGID3_1to1.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Br_exprtab_a_OGID3_1to1[1,])
colnames(Br_exprtab_a_OGID3_1to1.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Br_exprtab_a_OGID3_1to1.list.ED[is.na(Br_exprtab_a_OGID3_1to1.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Br_exprtab_a_OGID3_1to1.list.ED_rowsums <- as.data.frame(rowSums(Br_exprtab_a_OGID3_1to1.list.ED))
colnames(Br_exprtab_a_OGID3_1to1.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Br_exprtab_a_OGID3_1to1.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2 <- cbind(Br_exprtab_a_OGID3_1to1.list.ED,Br_exprtab_a_OGID3_1to1.list.ED_rowsums,Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate)
Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3 <- na.omit(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5 <- cbind(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3,Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale - do a y-intercept line for evol rate > log2(0.2) 
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/brain_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Brain expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = -2.321928,linetype="dashed") + geom_vline(xintercept = 1.807355,linetype="dashed")
dev.off()
# In the plot there are a few genes with high divergence (log2(4) = 2) && high evolutionary rate in promoter region (log2(0.2) = -2.3)
subset(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & ED.rowsum < 8 & prom > 0.2)[,c(11,13,14,15,16)]



# Eye
Ey_exprtab_a_OGID3_1to1cols <- seq(2, 30, by = 6) # increment by 6 from col1 to select tissue columns
Ey_exprtab_a_OGID3_1to1 <- as.matrix(All_exprtab_a_OGID3_1to1[,Ey_exprtab_a_OGID3_1to1cols])

## To run dist on each row save each row in a list
Ey_exprtab_a_OGID3_1to1.list <- setNames(split(Ey_exprtab_a_OGID3_1to1, seq(nrow(Ey_exprtab_a_OGID3_1to1))), rownames(Ey_exprtab_a_OGID3_1to1))
Ey_exprtab_a_OGID3_1to1.list.ED <- t(sapply(Ey_exprtab_a_OGID3_1to1.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Ey_exprtab_a_OGID3_1to1[1,])
colnames(Ey_exprtab_a_OGID3_1to1.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Ey_exprtab_a_OGID3_1to1.list.ED[is.na(Ey_exprtab_a_OGID3_1to1.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Ey_exprtab_a_OGID3_1to1.list.ED_rowsums <- as.data.frame(rowSums(Ey_exprtab_a_OGID3_1to1.list.ED))
colnames(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2 <- cbind(Ey_exprtab_a_OGID3_1to1.list.ED,Ey_exprtab_a_OGID3_1to1.list.ED_rowsums,Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate)
Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3 <- na.omit(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5 <- cbind(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3,Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/eye_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Eye expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = -2.321928,linetype="dashed") + geom_vline(xintercept = 1.807355,linetype="dashed")
dev.off()


# Heart
Ht_exprtab_a_OGID3_1to1cols <- seq(3, 30, by = 6) # increment by 6 from col1 to select tissue columns
Ht_exprtab_a_OGID3_1to1 <- as.matrix(All_exprtab_a_OGID3_1to1[,Ht_exprtab_a_OGID3_1to1cols])

## To run dist on each row save each row in a list
Ht_exprtab_a_OGID3_1to1.list <- setNames(split(Ht_exprtab_a_OGID3_1to1, seq(nrow(Ht_exprtab_a_OGID3_1to1))), rownames(Ht_exprtab_a_OGID3_1to1))
Ht_exprtab_a_OGID3_1to1.list.ED <- t(sapply(Ht_exprtab_a_OGID3_1to1.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Ht_exprtab_a_OGID3_1to1[1,])
colnames(Ht_exprtab_a_OGID3_1to1.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Ht_exprtab_a_OGID3_1to1.list.ED[is.na(Ht_exprtab_a_OGID3_1to1.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Ht_exprtab_a_OGID3_1to1.list.ED_rowsums <- as.data.frame(rowSums(Ht_exprtab_a_OGID3_1to1.list.ED))
colnames(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2 <- cbind(Ht_exprtab_a_OGID3_1to1.list.ED,Ht_exprtab_a_OGID3_1to1.list.ED_rowsums,Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate)
Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3 <- na.omit(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5 <- cbind(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3,Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/heart_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Heart expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = -2.321928,linetype="dashed") + geom_vline(xintercept = 1.807355,linetype="dashed")
dev.off()


# Kidney
Kd_exprtab_a_OGID3_1to1cols <- seq(4, 30, by = 6) # increment by 6 from col1 to select tissue columns
Kd_exprtab_a_OGID3_1to1 <- as.matrix(All_exprtab_a_OGID3_1to1[,Kd_exprtab_a_OGID3_1to1cols])

## To run dist on each row save each row in a list
Kd_exprtab_a_OGID3_1to1.list <- setNames(split(Kd_exprtab_a_OGID3_1to1, seq(nrow(Kd_exprtab_a_OGID3_1to1))), rownames(Kd_exprtab_a_OGID3_1to1))
Kd_exprtab_a_OGID3_1to1.list.ED <- t(sapply(Kd_exprtab_a_OGID3_1to1.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Kd_exprtab_a_OGID3_1to1[1,])
colnames(Kd_exprtab_a_OGID3_1to1.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Kd_exprtab_a_OGID3_1to1.list.ED[is.na(Kd_exprtab_a_OGID3_1to1.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Kd_exprtab_a_OGID3_1to1.list.ED_rowsums <- as.data.frame(rowSums(Kd_exprtab_a_OGID3_1to1.list.ED))
colnames(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2 <- cbind(Kd_exprtab_a_OGID3_1to1.list.ED,Kd_exprtab_a_OGID3_1to1.list.ED_rowsums,Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate)
Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3 <- na.omit(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5 <- cbind(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3,Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/kidney_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Kidney expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = -2.321928,linetype="dashed") + geom_vline(xintercept = 1.807355,linetype="dashed")
dev.off()



# Muscle
Ms_exprtab_a_OGID3_1to1cols <- seq(5, 30, by = 6) # increment by 6 from col1 to select tissue columns
Ms_exprtab_a_OGID3_1to1 <- as.matrix(All_exprtab_a_OGID3_1to1[,Ms_exprtab_a_OGID3_1to1cols])

## To run dist on each row save each row in a list
Ms_exprtab_a_OGID3_1to1.list <- setNames(split(Ms_exprtab_a_OGID3_1to1, seq(nrow(Ms_exprtab_a_OGID3_1to1))), rownames(Ms_exprtab_a_OGID3_1to1))
Ms_exprtab_a_OGID3_1to1.list.ED <- t(sapply(Ms_exprtab_a_OGID3_1to1.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Ms_exprtab_a_OGID3_1to1[1,])
colnames(Ms_exprtab_a_OGID3_1to1.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Ms_exprtab_a_OGID3_1to1.list.ED[is.na(Ms_exprtab_a_OGID3_1to1.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Ms_exprtab_a_OGID3_1to1.list.ED_rowsums <- as.data.frame(rowSums(Ms_exprtab_a_OGID3_1to1.list.ED))
colnames(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2 <- cbind(Ms_exprtab_a_OGID3_1to1.list.ED,Ms_exprtab_a_OGID3_1to1.list.ED_rowsums,Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate)
Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3 <- na.omit(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5 <- cbind(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3,Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/muscle_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Muscle expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = -2.321928,linetype="dashed") + geom_vline(xintercept = 1.807355,linetype="dashed")
dev.off()



# Testis
Ts_exprtab_a_OGID3_1to1cols <- seq(6, 30, by = 6) # increment by 6 from col1 to select tissue columns
Ts_exprtab_a_OGID3_1to1 <- as.matrix(All_exprtab_a_OGID3_1to1[,Ts_exprtab_a_OGID3_1to1cols])

## To run dist on each row save each row in a list
Ts_exprtab_a_OGID3_1to1.list <- setNames(split(Ts_exprtab_a_OGID3_1to1, seq(nrow(Ts_exprtab_a_OGID3_1to1))), rownames(Ts_exprtab_a_OGID3_1to1))
Ts_exprtab_a_OGID3_1to1.list.ED <- t(sapply(Ts_exprtab_a_OGID3_1to1.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Ts_exprtab_a_OGID3_1to1[1,])
colnames(Ts_exprtab_a_OGID3_1to1.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Ts_exprtab_a_OGID3_1to1.list.ED[is.na(Ts_exprtab_a_OGID3_1to1.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Ts_exprtab_a_OGID3_1to1.list.ED_rowsums <- as.data.frame(rowSums(Ts_exprtab_a_OGID3_1to1.list.ED))
colnames(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2 <- cbind(Ts_exprtab_a_OGID3_1to1.list.ED,Ts_exprtab_a_OGID3_1to1.list.ED_rowsums,Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate)
Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3 <- na.omit(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5 <- cbind(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate3,Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/testis_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Testis expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = -2.321928,linetype="dashed") + geom_vline(xintercept = 1.807355,linetype="dashed")
dev.off()

## visualising a distance matrix - this is just a test: not the correct way of doing it. Your matrix needs to be rows and columns of equal ID with pairwise comaprisons
library(factoextra)
Ts_exprtab_a_OGID3_1to1.list.EDTEST <- t(sapply(Ts_exprtab_a_OGID3_1to1.list, get_dist, method="euclidian"))
fviz_dist(Ts_exprtab_a_OGID3_1to1.list.EDTEST,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# the number of orthgroups in four quadrants (total of 5861 OGs after filtering)
  # A - top left: increased evol rate, but low divergence (ED.rowsum < 3.5 & prom > 22.7)
  # B - top right: increased evol rate, high divergence (ED.rowsum > 3.5 & prom > 22.7)
  # C - bottom right: low to intermediate evol rate, high divergence (ED.rowsum > 3.5 & prom < 22.7)
  # D - bottom left: low to intermediate evol rate, low divergence (ED.rowsum < 3.5 & prom < 22.7)
# Looked for some zebrafish housekeeping genes however none are in the 1to1 set e.g. gapdh, bactin2, acta1b

BrA <- nrow(subset(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)) # A - 1329
BrB <- nrow(subset(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)) # B - 1260
BrC <- nrow(subset(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)) # C - 1635
BrD <- nrow(subset(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)) # D - 1637

EyA <- nrow(subset(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)) # A - 958
EyB <- nrow(subset(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)) # B - 1631
EyC <- nrow(subset(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)) # C - 2011
EyD <- nrow(subset(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)) # D - 1261

HtA <- nrow(subset(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)) # A - 1136
HtB <- nrow(subset(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)) # B - 1453
HtC <- nrow(subset(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)) # C - 1807
HtD <- nrow(subset(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)) # D - 1465

KdA <- nrow(subset(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)) # A - 956
KdB <- nrow(subset(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)) # B - 1633
KdC <- nrow(subset(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)) # C - 2092
KdD <- nrow(subset(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)) # D - 1180

MsA <- nrow(subset(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)) # A - 957
MsB <- nrow(subset(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)) # B - 1632
MsC <- nrow(subset(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)) # C - 1986
MsD <- nrow(subset(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)) # D - 1286

TsA <- nrow(subset(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)) # A - 622
TsB <- nrow(subset(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)) # B - 1967
TsC <- nrow(subset(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)) # C - 2553
TsD <- nrow(subset(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)) # D - 719

# Plot the above numbers
EDquadrantsBr <- data.frame(BrA,BrB,BrC,BrD)
colnames(EDquadrantsBr) <- c('A','B','C','D')
EDquadrantsEy <- data.frame(EyA,EyB,EyC,EyD)
colnames(EDquadrantsEy) <- c('A','B','C','D')
EDquadrantsHt <- data.frame(HtA,HtB,HtC,HtD)
colnames(EDquadrantsHt) <- c('A','B','C','D')
EDquadrantsKd <- data.frame(KdA,KdB,KdC,KdD)
colnames(EDquadrantsKd) <- c('A','B','C','D')
EDquadrantsMs <- data.frame(MsA,MsB,MsC,MsD)
colnames(EDquadrantsMs) <- c('A','B','C','D')
EDquadrantsTs <- data.frame(TsA,TsB,TsC,TsD)
colnames(EDquadrantsTs) <- c('A','B','C','D')
EDquadrant <- rbind(EDquadrantsBr,EDquadrantsEy,EDquadrantsHt,EDquadrantsKd,EDquadrantsMs,EDquadrantsTs)
rownames(EDquadrant) <- c('Br','Ey','Ht','Kd','Ms','Ts')
EDquadrant2 <- melt(EDquadrant)
EDquadrant2$variable <- c('Br A', 'Ey A', 'Ht A','Kd A','Ms A','Ts A','Br B', 'Ey B', 'Ht B','Kd B','Ms B','Ts B','Br C', 'Ey C', 'Ht C','Kd C','Ms C','Ts C','Br D', 'Ey D', 'Ht D','Kd D','Ms D','Ts D')
EDquadrant2$type <- rep(c('A','B','C','D'), each=6)

positions <- c('Br A', 'Ey A', 'Ht A','Kd A','Ms A','Ts A','Br B', 'Ey B', 'Ht B','Kd B','Ms B','Ts B','Br C', 'Ey C', 'Ht C','Kd C','Ms C','Ts C','Br D', 'Ey D', 'Ht D','Kd D','Ms D','Ts D') # this alters the ordering of the bars on the plot
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/Euclidian_quadrants.tiff', units="in", width=12, height=12, res=300)
ggplot(EDquadrant2, aes(x = variable,y = value)) + scale_fill_brewer(palette="Pastel1") +
         geom_bar(stat="identity",aes(fill = type)) + scale_x_discrete(limits = positions) +
  labs (x = "Quadrant of expression divergence and evolutionary rate in tissue") + labs (y = "Number of 1:1 orthogroups") + guides(fill=guide_legend(title="Quadrant")) + labs(title="Quadrant of expression divergence as Euclidian distance in tissues of 1:1 orthologous genes and their evolutionary rate in promoter regions") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


### At this point, you should run a GO enrichment of the orthgroups with large expression difference and increased evol rate.
## In each plot, added lines for x = 2 and y = 5, pull out genes in that top right quadrant
Br_highdistevol <- subset(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
Ey_highdistevol <- subset(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)] # in this is rxrb, rarg
Ht_highdistevol <- subset(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
Kd_highdistevol <- subset(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
Ms_highdistevol <- subset(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
Ts_highdistevol <- subset(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
write.table(Br_highdistevol,file="Br_highdistevol.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ey_highdistevol,file="Ey_highdistevol.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ht_highdistevol,file="Ht_highdistevol.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Kd_highdistevol,file="Kd_highdistevol.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ms_highdistevol,file="Ms_highdistevol.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ts_highdistevol,file="Ts_highdistevol.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
# See line 7357 onwards in NetworkReconstruction_v4.sh for this - in summary, no statistically significant (q <0.05) enrichment however, we can use this output as mapping to GO terms instead

## Then, run a go enrichment for all the other quadrants
# Quadrant A
Br_highevol_lowdiv <- subset(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
Ey_highevol_lowdiv <- subset(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
Ht_highevol_lowdiv <- subset(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
Kd_highevol_lowdiv <- subset(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
Ms_highevol_lowdiv <- subset(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
Ts_highevol_lowdiv <- subset(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom > 22.7)[,c(11,13,14,15,16,12)]
write.table(Br_highevol_lowdiv,file="Br_highevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ey_highevol_lowdiv,file="Ey_highevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ht_highevol_lowdiv,file="Ht_highevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Kd_highevol_lowdiv,file="Kd_highevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ms_highevol_lowdiv,file="Ms_highevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ts_highevol_lowdiv,file="Ts_highevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)

# Quadrant C
Br_lowevol_highdiv <- subset(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Ey_lowevol_highdiv <- subset(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Ht_lowevol_highdiv <- subset(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Kd_lowevol_highdiv <- subset(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Ms_lowevol_highdiv <- subset(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Ts_lowevol_highdiv <- subset(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum > 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
write.table(Br_lowevol_highdiv,file="Br_lowevol_highdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ey_lowevol_highdiv,file="Ey_lowevol_highdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ht_lowevol_highdiv,file="Ht_lowevol_highdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Kd_lowevol_highdiv,file="Kd_lowevol_highdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ms_lowevol_highdiv,file="Ms_lowevol_highdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ts_lowevol_highdiv,file="Ts_lowevol_highdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)

# Quadrant D
Br_lowevol_lowdiv <- subset(Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Ey_lowevol_lowdiv <- subset(Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Ht_lowevol_lowdiv <- subset(Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Kd_lowevol_lowdiv <- subset(Kd_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Ms_lowevol_lowdiv <- subset(Ms_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
Ts_lowevol_lowdiv <- subset(Ts_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5, ED.rowsum < 3.5 & prom < 22.7)[,c(11,13,14,15,16,12)]
write.table(Br_lowevol_lowdiv,file="Br_lowevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ey_lowevol_lowdiv,file="Ey_lowevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ht_lowevol_lowdiv,file="Ht_lowevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Kd_lowevol_lowdiv,file="Kd_lowevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ms_lowevol_lowdiv,file="Ms_lowevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)
write.table(Ts_lowevol_lowdiv,file="Ts_lowevol_lowdiv.txt",row.names = TRUE,col.names = TRUE,sep = '\t',quote=FALSE)

## Check the average fold increase in promoter evolutionary rate (as comparison to 4fold, neutral) and average expression divergence for genes in each quadrant
# Select the quadrant boundaries by similar figures in fold increase and divergence for each e.g. similar average evol rate in quadrant A vs B, C vs D
mean(Br_highevol_lowdiv[,2])/mean(Br_highevol_lowdiv[,6]) # A - 270
mean(Br_highevol_lowdiv[,1]) # A - 2
mean(Br_highdistevol[,2])/mean(Br_highdistevol[,6]) # B - 282
mean(Br_highdistevol[,1]) # B - 5
mean(Br_lowevol_highdiv[,2])/mean(Br_lowevol_highdiv[,6]) # C - 5
mean(Br_lowevol_highdiv[,1]) # C - 5
mean(Br_lowevol_lowdiv[,2])/mean(Br_lowevol_lowdiv[,6]) # D - 8
mean(Br_lowevol_lowdiv[,1]) # D - 2

mean(Ey_highevol_lowdiv[,2])/mean(Ey_highevol_lowdiv[,6]) # A - 223
mean(Ey_highevol_lowdiv[,1]) # A - 2
mean(Ey_highdistevol[,2])/mean(Ey_highdistevol[,6]) # B - 319
mean(Ey_highdistevol[,1]) # B - 6
mean(Ey_lowevol_highdiv[,2])/mean(Ey_lowevol_highdiv[,6]) # C - 6
mean(Ey_lowevol_highdiv[,1]) # C - 6
mean(Ey_lowevol_lowdiv[,2])/mean(Ey_lowevol_lowdiv[,6]) # D - 6
mean(Ey_lowevol_lowdiv[,1]) # D - 2

mean(Ht_highevol_lowdiv[,2])/mean(Ht_highevol_lowdiv[,6]) # A - 249
mean(Ht_highevol_lowdiv[,1]) # A - 2
mean(Ht_highdistevol[,2])/mean(Ht_highdistevol[,6]) # B - 300
mean(Ht_highdistevol[,1]) # B - 6
mean(Ht_lowevol_highdiv[,2])/mean(Ht_lowevol_highdiv[,6]) # C - 7
mean(Ht_lowevol_highdiv[,1]) # C - 6
mean(Ht_lowevol_lowdiv[,2])/mean(Ht_lowevol_lowdiv[,6]) # D - 5
mean(Ht_lowevol_lowdiv[,1]) # D - 2

mean(Kd_highevol_lowdiv[,2])/mean(Kd_highevol_lowdiv[,6]) # A - 303
mean(Kd_highevol_lowdiv[,1]) # A - 2
mean(Kd_highdistevol[,2])/mean(Kd_highdistevol[,6]) # B - 262
mean(Kd_highdistevol[,1]) # B - 7
mean(Kd_lowevol_highdiv[,2])/mean(Kd_lowevol_highdiv[,6]) # C - 5
mean(Kd_lowevol_highdiv[,1]) # C - 7
mean(Kd_lowevol_lowdiv[,2])/mean(Kd_lowevol_lowdiv[,6]) # D - 9
mean(Kd_lowevol_lowdiv[,1]) # D - 2

mean(Ms_highevol_lowdiv[,2])/mean(Ms_highevol_lowdiv[,6]) # A - 445
mean(Ms_highevol_lowdiv[,1]) # A - 2
mean(Ms_highdistevol[,2])/mean(Ms_highdistevol[,6]) # B - 226
mean(Ms_highdistevol[,1]) # B - 6
mean(Ms_lowevol_highdiv[,2])/mean(Ms_lowevol_highdiv[,6]) # C - 6
mean(Ms_lowevol_highdiv[,1]) # C - 6
mean(Ms_lowevol_lowdiv[,2])/mean(Ms_lowevol_lowdiv[,6]) # D - 7
mean(Ms_lowevol_lowdiv[,1]) # D - 2

mean(Ts_highevol_lowdiv[,2])/mean(Ts_highevol_lowdiv[,6]) # A - 191
mean(Ts_highevol_lowdiv[,1]) # A - 3
mean(Ts_highdistevol[,2])/mean(Ts_highdistevol[,6]) # B - 319
mean(Ts_highdistevol[,1]) # B - 6
mean(Ts_lowevol_highdiv[,2])/mean(Ts_lowevol_highdiv[,6]) # C - 7
mean(Ts_lowevol_highdiv[,1]) # C - 6
mean(Ts_lowevol_lowdiv[,2])/mean(Ts_lowevol_lowdiv[,6]) # D - 5
mean(Ts_lowevol_lowdiv[,1]) # D - 3

# Calculate average fold change in evol rate (from 4fold to promoter) for each quadrant
(((mean(Br_highevol_lowdiv[,2])/mean(Br_highevol_lowdiv[,6]))+
  (mean(Ey_highevol_lowdiv[,2])/mean(Ey_highevol_lowdiv[,6]))+
  (mean(Ht_highevol_lowdiv[,2])/mean(Ht_highevol_lowdiv[,6]))+
  (mean(Kd_highevol_lowdiv[,2])/mean(Kd_highevol_lowdiv[,6]))+
  (mean(Ms_highevol_lowdiv[,2])/mean(Ms_highevol_lowdiv[,6]))+
  (mean(Ts_highevol_lowdiv[,2])/mean(Ts_highevol_lowdiv[,6]))
)/6) # A - 280.2963

(((mean(Br_highdistevol[,2])/mean(Br_highdistevol[,6]))+
    (mean(Ey_highdistevol[,2])/mean(Ey_highdistevol[,6]))+
    (mean(Ht_highdistevol[,2])/mean(Ht_highdistevol[,6]))+
    (mean(Kd_highdistevol[,2])/mean(Kd_highdistevol[,6]))+
    (mean(Ms_highdistevol[,2])/mean(Ms_highdistevol[,6]))+
    (mean(Ts_highdistevol[,2])/mean(Ts_highdistevol[,6]))
  )/6) # B - 284.6816

(((mean(Br_lowevol_highdiv[,2])/mean(Br_lowevol_highdiv[,6]))+
  (mean(Ey_lowevol_highdiv[,2])/mean(Ey_lowevol_highdiv[,6]))+
  (mean(Ht_lowevol_highdiv[,2])/mean(Ht_lowevol_highdiv[,6]))+
  (mean(Kd_lowevol_highdiv[,2])/mean(Kd_lowevol_highdiv[,6]))+
  (mean(Ms_lowevol_highdiv[,2])/mean(Ms_lowevol_highdiv[,6]))+
  (mean(Ts_lowevol_highdiv[,2])/mean(Ts_lowevol_highdiv[,6]))
)/6) # C - 6.017961

(((mean(Br_lowevol_lowdiv[,2])/mean(Br_lowevol_lowdiv[,6]))+
    (mean(Ey_lowevol_lowdiv[,2])/mean(Ey_lowevol_lowdiv[,6]))+
    (mean(Ht_lowevol_lowdiv[,2])/mean(Ht_lowevol_lowdiv[,6]))+
    (mean(Kd_lowevol_lowdiv[,2])/mean(Kd_lowevol_lowdiv[,6]))+
    (mean(Ms_lowevol_lowdiv[,2])/mean(Ms_lowevol_lowdiv[,6]))+
    (mean(Ts_lowevol_lowdiv[,2])/mean(Ts_lowevol_lowdiv[,6]))
)/6) # D - 6.927445

# Calculate average euclidian distance to check simialrity
((mean(Br_highevol_lowdiv[,1]) + mean(Ey_highevol_lowdiv[,1]) + mean(Ht_highevol_lowdiv[,1]) + mean(Kd_highevol_lowdiv[,1]) + mean(Ms_highevol_lowdiv[,1]) + mean(Br_highevol_lowdiv[,1]) + mean(Ey_highevol_lowdiv[,1]) + mean(Ht_highevol_lowdiv[,1]) + mean(Br_highevol_lowdiv[,1]) + mean(Ts_highevol_lowdiv[,1]))/6) # A - 3.987657
((mean(Br_highdistevol[,1]) + mean(Ey_highdistevol[,1]) + mean(Ht_highdistevol[,1]) + mean(Kd_highdistevol[,1]) + mean(Ms_highdistevol[,1]) + mean(Br_highdistevol[,1]) + mean(Ey_highdistevol[,1]) + mean(Ht_highdistevol[,1]) + mean(Br_highdistevol[,1]) + mean(Ts_highdistevol[,1]))/6) # B - 9.848171
((mean(Br_lowevol_highdiv[,1]) + mean(Ey_lowevol_highdiv[,1]) + mean(Ht_lowevol_highdiv[,1]) + mean(Kd_lowevol_highdiv[,1]) + mean(Ms_lowevol_highdiv[,1]) + mean(Br_lowevol_highdiv[,1]) + mean(Ey_lowevol_highdiv[,1]) + mean(Ht_lowevol_highdiv[,1]) + mean(Br_lowevol_highdiv[,1]) + mean(Ts_lowevol_highdiv[,1]))/6) # C - 9.922247
((mean(Br_lowevol_lowdiv[,1]) + mean(Ey_lowevol_lowdiv[,1]) + mean(Ht_lowevol_lowdiv[,1]) + mean(Kd_lowevol_lowdiv[,1]) + mean(Ms_lowevol_lowdiv[,1]) + mean(Br_lowevol_lowdiv[,1]) + mean(Ey_lowevol_lowdiv[,1]) + mean(Ht_lowevol_lowdiv[,1]) + mean(Br_lowevol_lowdiv[,1]) + mean(Ts_lowevol_lowdiv[,1]))/6) # D - 4.006146

### Look at the R-squared for correlation between expression divergence and sequence conservation (as evolutionary rate)
# R-squared being the statistical measure of how close the data are to the fitted regression line - expect this to be low!
brain.lm = lm(ED.rowsum ~ prom, data=Br_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5)
summary(brain.lm)$r.squared # 9.945546e-06
eye.lm = lm(ED.rowsum ~ prom, data=Ey_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5)
summary(eye.lm)$r.squared # 2.848292e-06
heart.lm = lm(ED.rowsum ~ prom, data=Ht_exprtab_a_OGID3_1to1.list.ED_rowsums.evolrate5)
summary(heart.lm)$r.squared # 6.040135e-05

### FigS-R2e – GO enrichment of 1:1 orthogroups in each quadrant of evolutionary rate in promoter region against multispecies Euclidian distance of expression divergence in six tissues

## Here, we will plot the GO enrichment - note, that when filtered based on q < 0.05, left with no results
# Instead, just plot -log10(p-val < 0.05)
# prepared all files for this from line 7575 in NetworkReconstruction_v4.sh script

# load in the data
Files8 <- dir(file.path(base_dir2))
quad_GOgrep <- glob2rx("*GOOUTPUT_details_filtered2.txt3") # create a grep pattern to select only specific files
quad_GOgrep2 <- grep(quad_GOgrep, Files8) # run the grep
Files9 <- Files8[quad_GOgrep2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files9)){
  tmp = read.delim(file = paste0(base_dir2,Files9[i]),header = F)
  assign(Files9[i], tmp)
}

# this is to create combined plots for each tissue
BrquadGO_combined <- rbind(Br_highevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3,Br_highdistevol.genes.GOOUTPUT_details_filtered2.txt3,Br_lowevol_highdiv.genes.GOOUTPUT_details_filtered2.txt3,Br_lowevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3)
BrquadGO_combined$V1 <- gsub("highdistevol", "B - highevol/highdiv", BrquadGO_combined$V1) ## rename highdistevol to highevol_highdiv
BrquadGO_combined$V1 <- gsub("highevol_lowdiv","A - highevol/lowdiv", BrquadGO_combined$V1)
BrquadGO_combined$V1 <- gsub("lowevol_highdiv","C - lowevol/highdiv", BrquadGO_combined$V1)
BrquadGO_combined$V1 <- gsub("lowevol_lowdiv","D - lowevol/lowdiv", BrquadGO_combined$V1)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/FigS-R2eA_Br_quadrantGO.tiff', units="in", width=12, height=12, res=300)
ggplot(BrquadGO_combined, aes(x = reorder(V2, -V3), y = -log10(V3), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_fill_brewer(palette="Set3") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=3, group='V1') + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)") + guides(fill=guide_legend(title="Quadrant"))
dev.off()

EyquadGO_combined <- rbind(Ey_highevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3,Ey_highdistevol.genes.GOOUTPUT_details_filtered2.txt3,Ey_lowevol_highdiv.genes.GOOUTPUT_details_filtered2.txt3,Ey_lowevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3)
EyquadGO_combined$V1 <- gsub("highdistevol", "B - highevol/highdiv", EyquadGO_combined$V1) ## rename highdistevol to highevol_highdiv
EyquadGO_combined$V1 <- gsub("highevol_lowdiv","A - highevol/lowdiv", EyquadGO_combined$V1)
EyquadGO_combined$V1 <- gsub("lowevol_highdiv","C - lowevol/highdiv", EyquadGO_combined$V1)
EyquadGO_combined$V1 <- gsub("lowevol_lowdiv","D - lowevol/lowdiv", EyquadGO_combined$V1)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/FigS-R2eB_Ey_quadrantGO.tiff', units="in", width=12, height=12, res=300)
ggplot(EyquadGO_combined, aes(x = reorder(V2, -V3), y = -log10(V3), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_fill_brewer(palette="Set3") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=3, group='V1') + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)") + guides(fill=guide_legend(title="Quadrant"))
dev.off()

HtquadGO_combined <- rbind(Ht_highevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3,Ht_highdistevol.genes.GOOUTPUT_details_filtered2.txt3,Ht_lowevol_highdiv.genes.GOOUTPUT_details_filtered2.txt3,Ht_lowevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3)
HtquadGO_combined$V1 <- gsub("highdistevol", "B - highevol/highdiv", HtquadGO_combined$V1) ## rename highdistevol to highevol_highdiv
HtquadGO_combined$V1 <- gsub("highevol_lowdiv","A - highevol/lowdiv", HtquadGO_combined$V1)
HtquadGO_combined$V1 <- gsub("lowevol_highdiv","C - lowevol/highdiv", HtquadGO_combined$V1)
HtquadGO_combined$V1 <- gsub("lowevol_lowdiv","D - lowevol/lowdiv", HtquadGO_combined$V1)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/FigS-R2eC_Ht_quadrantGO.tiff', units="in", width=12, height=12, res=300)
ggplot(HtquadGO_combined, aes(x = reorder(V2, -V3), y = -log10(V3), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_fill_brewer(palette="Set3") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=3, group='V1') + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)") + guides(fill=guide_legend(title="Quadrant"))
dev.off()

KdquadGO_combined <- rbind(Kd_highevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3,Kd_highdistevol.genes.GOOUTPUT_details_filtered2.txt3,Kd_lowevol_highdiv.genes.GOOUTPUT_details_filtered2.txt3,Kd_lowevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3)
KdquadGO_combined$V1 <- gsub("highdistevol", "B - highevol/highdiv", KdquadGO_combined$V1) ## rename highdistevol to highevol_highdiv
KdquadGO_combined$V1 <- gsub("highevol_lowdiv","A - highevol/lowdiv", KdquadGO_combined$V1)
KdquadGO_combined$V1 <- gsub("lowevol_highdiv","C - lowevol/highdiv", KdquadGO_combined$V1)
KdquadGO_combined$V1 <- gsub("lowevol_lowdiv","D - lowevol/lowdiv", KdquadGO_combined$V1)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/FigS-R2eD_Kd_quadrantGO.tiff', units="in", width=12, height=12, res=300)
ggplot(KdquadGO_combined, aes(x = reorder(V2, -V3), y = -log10(V3), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_fill_brewer(palette="Set3") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=3, group='V1') + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)") + guides(fill=guide_legend(title="Quadrant"))
dev.off()

MsquadGO_combined <- rbind(Ms_highevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3,Ms_highdistevol.genes.GOOUTPUT_details_filtered2.txt3,Ms_lowevol_highdiv.genes.GOOUTPUT_details_filtered2.txt3,Ms_lowevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3)
MsquadGO_combined$V1 <- gsub("highdistevol", "B - highevol/highdiv", MsquadGO_combined$V1) ## rename highdistevol to highevol_highdiv
MsquadGO_combined$V1 <- gsub("highevol_lowdiv","A - highevol/lowdiv", MsquadGO_combined$V1)
MsquadGO_combined$V1 <- gsub("lowevol_highdiv","C - lowevol/highdiv", MsquadGO_combined$V1)
MsquadGO_combined$V1 <- gsub("lowevol_lowdiv","D - lowevol/lowdiv", MsquadGO_combined$V1)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/FigS-R2eE_Ms_quadrantGO.tiff', units="in", width=12, height=12, res=300)
ggplot(MsquadGO_combined, aes(x = reorder(V2, -V3), y = -log10(V3), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_fill_brewer(palette="Set3") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=3, group='V1') + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)") + guides(fill=guide_legend(title="Quadrant"))
dev.off()

TsquadGO_combined <- rbind(Ts_highevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3,Ts_highdistevol.genes.GOOUTPUT_details_filtered2.txt3,Ts_lowevol_highdiv.genes.GOOUTPUT_details_filtered2.txt3,Ts_lowevol_lowdiv.genes.GOOUTPUT_details_filtered2.txt3)
TsquadGO_combined$V1 <- gsub("highdistevol", "B - highevol/highdiv", TsquadGO_combined$V1) ## rename highdistevol to highevol_highdiv
TsquadGO_combined$V1 <- gsub("highevol_lowdiv","A - highevol/lowdiv", TsquadGO_combined$V1)
TsquadGO_combined$V1 <- gsub("lowevol_highdiv","C - lowevol/highdiv", TsquadGO_combined$V1)
TsquadGO_combined$V1 <- gsub("lowevol_lowdiv","D - lowevol/lowdiv", TsquadGO_combined$V1)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/FigS-R2eF_Ts_quadrantGO.tiff', units="in", width=12, height=12, res=300)
ggplot(TsquadGO_combined, aes(x = reorder(V2, -V3), y = -log10(V3), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_fill_brewer(palette="Set3") +
  theme_bw() + theme(panel.grid = element_blank()) + coord_flip() + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=3, group='V1') + labs (x = "GO Term") + labs (y = "-log10 (P-value < 0.05)") + guides(fill=guide_legend(title="Quadrant"))
dev.off()


## Plot the number of module switched genes that overlap quadrants of ED of expr divergence and evol rate
# Created a table of the numbers on line 7964 onwards in NetworkReconstruction_v4.sh script
quadrantoverlapswitchgenes <- read.delim(file="~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/ED_quadrants/quadrantgenes_overlapstatechange.txt2", header = F, sep = "\t")
colnames(quadrantoverlapswitchgenes) <- c("variable","value","type")

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/EuclidianDist/FigS-R2f_Module_state_change_overlap_quadrants.tiff', units="in", width=12, height=12, res=300)
ggplot(quadrantoverlapswitchgenes, aes(x = variable,y = value)) + scale_fill_brewer(palette="Pastel1") +
  geom_bar(stat="identity",aes(fill = type)) +
  labs (x = "Tissue and quadrant of genes overlapping module state changes") + labs (y = "Number of 1:1 orthogroups") + guides(fill=guide_legend(title="Ancestral nodes of state change")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~

# Comparisons of total branch lengths of expression trees among the six tissues

# Above, I computed, separately for each 1:1 orthgroup, a pairwise matrix of Euclidian gene expression distance between all genes in a pair of species.
# Created 6 expression distance matrices, 1 for each tissue (6844 OGs x 10 species pairs)

# Here, we want to construct gene expression trees based on the above matrices using the neighbour-joining approach (implemented in the ape package)
# Then, use the branch lengths of trees as a measure of gene expression divergence
# Plot branch lengths of each tissue and do comparisons between species pairs
# We can then use the branch lengths instead, to plot against evolutionary rate

XX <- dist(Br_exprtab_a_OGID3_1to1)

t(sapply(Ms_exprtab_a_OGID3_1to1.list, dist))

XX <- dist(t(Ts_exprtab_a_OGID3_1to1))
XX2 <- as.matrix(XX)
XX3 <- nj(XX2)
XX4 <- function(x) nj(x)
XX5 <- boot.phylo(XX3, XX2, XX4)

plot(XX3, cex=1.5)
nodelabels(XX5, adj = c(1.2, 1), frame = "n", cex = 1)

par(mfcol = c(3, 2))






######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~


### From Jenny Chen talk at CSHL
# Rate of Expression evolution is not linear
# Plot expression difference in each tissue = residuals after regressing pairwise comparisons of gene expression
# This is all against O. niloticus
# Plot is squared expression difference (y-axis, 0 to 0.1) vs evolutionary time (substitutions per 100bp, 0 to 0.8)

# t is a matrix with 2 columns, each column is the species you are comparing.  Each row is a gene and the log10(tpm+0.01) value for that gene.
# If you have multiple samples for each tissue in a species, then just take the mean of those samples:

t = cbind(ref.mean, other.mean) # make matrix t
t = t[apply(t, 1, function(x) sum(is.na(x)))  == 0,] # remove things that are NA in both species
m = prcomp(t) # do a PCA
res = m$x[,2]  # take the residuals from the 2nd axis
if (m$rotation[2,2]<0) r = r*-1 # sign of 'r' tells you if "other.mean" is higher or lower than expected

# If you want to plot where the PCA axis is, you just have to center the 't' matrix since the prcomp does that automatically for you:

plot(ref.mean - mean(ref.mean, na.rm=T), other.mean - mean(other.mean, na.rm=T))
abline(a=0, b=m$rotation[2,1]/m$rotation[1,1], col="red", lwd=2)


######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~


### Fitting Ornstein Uhlenbeck (OU) to the data
## Use the 'OUCH' package

#install.packages("ouch",repos="http://kingaa.github.io/")
#install.packages("devtools")
library(devtools)
library(ouch)

# Running the following for each tissue, across the species

# First, load in the module gene expression (as they are log transformed and normalised)
base_dir2 <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/"

# load expression data - this is of the modules, 1to1 only (6844 orthgroups)
Files4 <- dir(file.path(base_dir2))
modexprgrep <- glob2rx("*-exprtab_a-OGID3-1to1.txt") # create a grep pattern to select only specific files
modexprgrep2 <- grep(modexprgrep, Files4) # run the grep
Files5 <- Files4[modexprgrep2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files5)){
  tmp = read.delim(file = paste0(base_dir2,Files5[i]),header = F, row.names = 1)
  assign(Files5[i], tmp)
}
# change colnames
colnames(`Ab-exprtab_a-OGID3-1to1.txt`) <- c('Ab_gene', 'Ab_br', 'Ab_ey','Ab_ht','Ab_kd','Ab_ms','Ab_ts')
colnames(`Mz-exprtab_a-OGID3-1to1.txt`) <- c('Mz_gene', 'Mz_br', 'Mz_ey','Mz_ht','Mz_kd','Mz_ms','Mz_ts')
colnames(`Pn-exprtab_a-OGID3-1to1.txt`) <- c('Pn_gene', 'Pn_br', 'Pn_ey','Pn_ht','Pn_kd','Pn_ms','Pn_ts')
colnames(`Nb-exprtab_a-OGID3-1to1.txt`) <- c('Nb_gene', 'Nb_br', 'Nb_ey','Nb_ht','Nb_kd','Nb_ms','Nb_ts')
colnames(`On-exprtab_a-OGID3-1to1.txt`) <- c('On_gene', 'On_br', 'On_ey','On_ht','On_kd','On_ms','On_ts')

# cbind all data and then separate by tissue
OGIDS1to1oudata<-cbind.data.frame(`Mz-exprtab_a-OGID3-1to1.txt`,`Pn-exprtab_a-OGID3-1to1.txt`,`Ab-exprtab_a-OGID3-1to1.txt`,`Nb-exprtab_a-OGID3-1to1.txt`,`On-exprtab_a-OGID3-1to1.txt`) # since the order is the same we can just cbind
Br_cols <- seq(2, 35, by = 7) # increment by 7 from col2 to select tissue columns
Br_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Br_cols])
Ey_cols <- seq(3, 35, by = 7) # increment by 7 from col3 to select tissue columns
Ey_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Ey_cols])
Ht_cols <- seq(4, 35, by = 7) # increment by 7 from col4 to select tissue columns
Ht_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Ht_cols])
Kd_cols <- seq(5, 35, by = 7) # increment by 7 from col5 to select tissue columns
Kd_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Kd_cols])
Ms_cols <- seq(6, 35, by = 7) # increment by 7 from col6 to select tissue columns
Ms_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Ms_cols])
Ts_cols <- seq(7, 35, by = 7) # increment by 7 from col7 to select tissue columns
Ts_OGIDS1to1oudata <- as.matrix(OGIDS1to1oudata[,Ts_cols])

colnames(Br_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Ey_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Ht_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Kd_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Ms_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')
colnames(Ts_OGIDS1to1oudata)<-c('Mzebra','Pnyererei','Aburtoni','Nbrichardi','Oniloticus')

## Then, we need to prepare the phylo tree with branch lengths that will be used for fitting the model
library(ape)
# load cichlid tree in (nexus format)
cichlidtree<-read.nexus("/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Gene_expression_data/OU_fromJenny/cichlidtree.nex")
# Get tree into OUCH format using the ape2ouch function in ouch
cichlidtree2<-ape2ouch(cichlidtree)
# You can also plot this tree to see if it looks right.
plot(cichlidtree2)
# Now lets look at it - you can see some differences from the ape phylo class: the internal nodes are numbered first followed by the tips etc.
cichlidtree2
# To access the times column we need to use the syntax - note these are different to your nexus inputs
cichlidtree2@times
## Now we want to get our data into a format that ouch can use.  This is a common stumbling block for people when they try to use ouch.  I have modified the example that was sent around on the r-sig-phylo mailing list a couple of months ago written by Graham Slater and Aaron King.  I have also added several sections to the end to explain how to fit multiple optima.
# Convert the tree to a dataframe
cichlidtree2dataframe<-as(cichlidtree2,"data.frame")
# Add the anc names to labels column
levels(cichlidtree2dataframe$labels) <- c(levels(cichlidtree2dataframe$labels), "Anc4")
cichlidtree2dataframe$labels[1]<-"Anc4"
levels(cichlidtree2dataframe$labels) <- c(levels(cichlidtree2dataframe$labels), "Anc1")
cichlidtree2dataframe$labels[2]<-"Anc1"
levels(cichlidtree2dataframe$labels) <- c(levels(cichlidtree2dataframe$labels), "Anc2")
cichlidtree2dataframe$labels[3]<-"Anc2"
levels(cichlidtree2dataframe$labels) <- c(levels(cichlidtree2dataframe$labels), "Anc3")
cichlidtree2dataframe$labels[4]<-"Anc3"

## We will be running the OU model by assigning each of the nodes and tips to an optima (regimes) - all need to be added as a factor for processing
cichlidtree2dataframe$regimes.null = as.factor(rep("anc", 3)) # the null hypothesis of a single optima across all nodes and branches
cichlidtree2dataframe$regimes.haplo = as.factor(c(1,1,2,2,1,1,2,2,2)) # an optima of the haplochromines only
cichlidtree2dataframe$regimes.haplo_Nb = as.factor(c(1,2,2,2,1,2,2,2,2)) # an optima of the haplochromines and Nb
cichlidtree2dataframe$regimes.riverine = as.factor(c(2,1,2,1,2,1,2,1,1)) # an optima of the riverines
ouch.cichlid.regime <- cichlidtree2dataframe # rename


## First run a Brownian model with the brown function using the expression data.

# Remake the outree this time constructing it from the ground up from the ouch.cichlid.regime dataframe using the ouchtree command in OUCH.
# This creates a tree with the nodes in the same order as those in the ouch.cichlid.regime

curOuchTree = ouchtree(ouch.cichlid.regime$nodes, ouch.cichlid.regime$ancestors, ouch.cichlid.regime$times, labels = as.character(ouch.cichlid.regime$labels))
plot(curOuchTree) # check the plotting of the tree

# 1. Add four columns before the expression data for each tissue, filled with NA - these are for the Anc expression, then rearrange the column order
Br_OGIDS1to1oudata2 <- cbind(Anc4=NA, Anc1=NA, Anc2=NA, Anc3=NA, Br_OGIDS1to1oudata[,c('Oniloticus','Nbrichardi','Aburtoni','Pnyererei','Mzebra')])
Ey_OGIDS1to1oudata2 <- cbind(Anc4=NA, Anc1=NA, Anc2=NA, Anc3=NA, Ey_OGIDS1to1oudata[,c('Oniloticus','Nbrichardi','Aburtoni','Pnyererei','Mzebra')])
Ht_OGIDS1to1oudata2 <- cbind(Anc4=NA, Anc1=NA, Anc2=NA, Anc3=NA, Ht_OGIDS1to1oudata[,c('Oniloticus','Nbrichardi','Aburtoni','Pnyererei','Mzebra')])
Kd_OGIDS1to1oudata2 <- cbind(Anc4=NA, Anc1=NA, Anc2=NA, Anc3=NA, Kd_OGIDS1to1oudata[,c('Oniloticus','Nbrichardi','Aburtoni','Pnyererei','Mzebra')])
Ms_OGIDS1to1oudata2 <- cbind(Anc4=NA, Anc1=NA, Anc2=NA, Anc3=NA, Ms_OGIDS1to1oudata[,c('Oniloticus','Nbrichardi','Aburtoni','Pnyererei','Mzebra')])
Ts_OGIDS1to1oudata2 <- cbind(Anc4=NA, Anc1=NA, Anc2=NA, Anc3=NA, Ts_OGIDS1to1oudata[,c('Oniloticus','Nbrichardi','Aburtoni','Pnyererei','Mzebra')])
# 2. Transpose the dataframe of each tissue and cbind to the ouch.cichlid.regime df (per tissue)
Br.ouch.cichlid.regime <- cbind(ouch.cichlid.regime,t(Br_OGIDS1to1oudata2))
Ey.ouch.cichlid.regime <- cbind(ouch.cichlid.regime,t(Ey_OGIDS1to1oudata2))
Ht.ouch.cichlid.regime <- cbind(ouch.cichlid.regime,t(Ht_OGIDS1to1oudata2))
Kd.ouch.cichlid.regime <- cbind(ouch.cichlid.regime,t(Kd_OGIDS1to1oudata2))
Ms.ouch.cichlid.regime <- cbind(ouch.cichlid.regime,t(Ms_OGIDS1to1oudata2))
Ts.ouch.cichlid.regime <- cbind(ouch.cichlid.regime,t(Ts_OGIDS1to1oudata2))
# 3. Save this as a list for each tissue, but remove the regimes etc. so that we can run brown and hansen
Br.ouch.cichlid.regime.list <- as.list(Br.ouch.cichlid.regime[,9:6852]) # create a list of only the expression data
Ey.ouch.cichlid.regime.list <- as.list(Ey.ouch.cichlid.regime[,9:6852]) # create a list of only the expression data
Ht.ouch.cichlid.regime.list <- as.list(Ht.ouch.cichlid.regime[,9:6852]) # create a list of only the expression data
Kd.ouch.cichlid.regime.list <- as.list(Kd.ouch.cichlid.regime[,9:6852]) # create a list of only the expression data
Ms.ouch.cichlid.regime.list <- as.list(Ms.ouch.cichlid.regime[,9:6852]) # create a list of only the expression data
Ts.ouch.cichlid.regime.list <- as.list(Ts.ouch.cichlid.regime[,9:6852]) # create a list of only the expression data

# 4. Assign node number as 'names' attributes to each row of the list using a for loop
list.names <- names(Br.ouch.cichlid.regime.list)
for (i in 1:length(list.names)) {
  names(Br.ouch.cichlid.regime.list[[i]]) <- Br.ouch.cichlid.regime$nodes # or could do rep(1:9)
}
list.names <- names(Ey.ouch.cichlid.regime.list)
for (i in 1:length(list.names)) {
  names(Ey.ouch.cichlid.regime.list[[i]]) <- Ey.ouch.cichlid.regime$nodes # or could do rep(1:9)
}
list.names <- names(Ht.ouch.cichlid.regime.list)
for (i in 1:length(list.names)) {
  names(Ht.ouch.cichlid.regime.list[[i]]) <- Ht.ouch.cichlid.regime$nodes # or could do rep(1:9)
}
list.names <- names(Kd.ouch.cichlid.regime.list)
for (i in 1:length(list.names)) {
  names(Kd.ouch.cichlid.regime.list[[i]]) <- Kd.ouch.cichlid.regime$nodes # or could do rep(1:9)
}
list.names <- names(Ms.ouch.cichlid.regime.list)
for (i in 1:length(list.names)) {
  names(Ms.ouch.cichlid.regime.list[[i]]) <- Ms.ouch.cichlid.regime$nodes # or could do rep(1:9)
}
list.names <- names(Ts.ouch.cichlid.regime.list)
for (i in 1:length(list.names)) {
  names(Ts.ouch.cichlid.regime.list[[i]]) <- Ts.ouch.cichlid.regime$nodes # or could do rep(1:9)
}

# 5. Run brownian motion and hansen (OU) for each row and regime, storing in a list using a for loop - you can access each result using list[[#]], then access each element e.g. loglik by list[[#]]@loglik

# hansen - Fits the Ornstein-Uhlenbeck-based Hansen model to data
# sqrt.alpha, sigma - used to initialize the optimization algorithm.  The selection strength matrix alpha and the random drift variance-covariance matrix sigma-sq are parameterized by their matrix square roots.  Specifically, these initial guesses are each packed into lower-triangular matrices (column by column).  The product of this matrix with its transpose is the alpha or sigma-sq matrix.
# fit - If fit=TRUE, then the likelihood will be maximized. Iffit=FALSE, the likelihood will be evaluated at the specified values ofsqrt.alpha and sigma; the optima theta will be returned as well.

# if you get unsuccessful convergence problems then there are a few things that you can do
# run more iterations - control=(list(maxit00)
# increase the relative tolerance (reltol) parameter. By increasing the tolerance (making the value smaller) it sometimes helps with convergence - tested and works with 1e-1 max (no errors of unsuccessful convergence, wich happens with 1e-2)
# change the moethod from optim to method="subplex"

# brown(data=Br.ouch.cichlid.regime.list["OG10000_0"], tree=curOuchTree) # this runs on a single OG
Br_brown <- vector("list",6844) # store all the brownian motion results as vectors in a list of length 6844
Br_OUhaplo <- vector("list",6844) # store all the OU haplo results as vectors in a list of length 6844
Br_OUhaploNb <- vector("list",6844) # store all the OU haplo+Nb results as vectors in a list of length 6844
Br_OUriverine <- vector("list",6844) # store all the OU riverine results as vectors in a list of length 6844
for (i in 1:length(Br.ouch.cichlid.regime.list)) {
  Br_brown[[i]] <- brown(data=Br.ouch.cichlid.regime.list[i], tree=curOuchTree)
  Br_OUhaplo[[i]] <-hansen(data=Br.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Br_OUhaploNb[[i]] <-hansen(data=Br.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo_Nb"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Br_OUriverine[[i]] <-hansen(data=Br.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.riverine"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
}

Ey_brown <- vector("list",6844) # store all the brownian motion results as vectors in a list of length 6844
Ey_OUhaplo <- vector("list",6844) # store all the OU haplo results as vectors in a list of length 6844
Ey_OUhaploNb <- vector("list",6844) # store all the OU haplo+Nb results as vectors in a list of length 6844
Ey_OUriverine <- vector("list",6844) # store all the OU riverine results as vectors in a list of length 6844
for (i in 1:length(Ey.ouch.cichlid.regime.list)) {
  Ey_brown[[i]] <- brown(data=Ey.ouch.cichlid.regime.list[i], tree=curOuchTree)
  Ey_OUhaplo[[i]] <-hansen(data=Ey.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Ey_OUhaploNb[[i]] <-hansen(data=Ey.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo_Nb"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Ey_OUriverine[[i]] <-hansen(data=Ey.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.riverine"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
}

Ht_brown <- vector("list",6844) # store all the brownian motion results as vectors in a list of length 6844
Ht_OUhaplo <- vector("list",6844) # store all the OU haplo results as vectors in a list of length 6844
Ht_OUhaploNb <- vector("list",6844) # store all the OU haplo+Nb results as vectors in a list of length 6844
Ht_OUriverine <- vector("list",6844) # store all the OU riverine results as vectors in a list of length 6844
for (i in 1:length(Ht.ouch.cichlid.regime.list)) {
  Ht_brown[[i]] <- brown(data=Ht.ouch.cichlid.regime.list[i], tree=curOuchTree)
  Ht_OUhaplo[[i]] <-hansen(data=Ht.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Ht_OUhaploNb[[i]] <-hansen(data=Ht.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo_Nb"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Ht_OUriverine[[i]] <-hansen(data=Ht.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.riverine"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
}

Kd_brown <- vector("list",6844) # store all the brownian motion results as vectors in a list of length 6844
Kd_OUhaplo <- vector("list",6844) # store all the OU haplo results as vectors in a list of length 6844
Kd_OUhaploNb <- vector("list",6844) # store all the OU haplo+Nb results as vectors in a list of length 6844
Kd_OUriverine <- vector("list",6844) # store all the OU riverine results as vectors in a list of length 6844
for (i in 1:length(Kd.ouch.cichlid.regime.list)) {
  Kd_brown[[i]] <- brown(data=Kd.ouch.cichlid.regime.list[i], tree=curOuchTree)
  Kd_OUhaplo[[i]] <-hansen(data=Kd.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Kd_OUhaploNb[[i]] <-hansen(data=Kd.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo_Nb"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Kd_OUriverine[[i]] <-hansen(data=Kd.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.riverine"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
}

Ms_brown <- vector("list",6844) # store all the brownian motion results as vectors in a list of length 6844
Ms_OUhaplo <- vector("list",6844) # store all the OU haplo results as vectors in a list of length 6844
Ms_OUhaploNb <- vector("list",6844) # store all the OU haplo+Nb results as vectors in a list of length 6844
Ms_OUriverine <- vector("list",6844) # store all the OU riverine results as vectors in a list of length 6844
for (i in 1:length(Ms.ouch.cichlid.regime.list)) {
  Ms_brown[[i]] <- brown(data=Ms.ouch.cichlid.regime.list[i], tree=curOuchTree)
  Ms_OUhaplo[[i]] <-hansen(data=Ms.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Ms_OUhaploNb[[i]] <-hansen(data=Ms.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo_Nb"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Ms_OUriverine[[i]] <-hansen(data=Ms.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.riverine"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
}

Ts_brown <- vector("list",6844) # store all the brownian motion results as vectors in a list of length 6844
Ts_OUhaplo <- vector("list",6844) # store all the OU haplo results as vectors in a list of length 6844
Ts_OUhaploNb <- vector("list",6844) # store all the OU haplo+Nb results as vectors in a list of length 6844
Ts_OUriverine <- vector("list",6844) # store all the OU riverine results as vectors in a list of length 6844
for (i in 1:length(Ts.ouch.cichlid.regime.list)) {
  Ts_brown[[i]] <- brown(data=Ts.ouch.cichlid.regime.list[i], tree=curOuchTree)
  Ts_OUhaplo[[i]] <-hansen(data=Ts.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Ts_OUhaploNb[[i]] <-hansen(data=Ts.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.haplo_Nb"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
  Ts_OUriverine[[i]] <-hansen(data=Ts.ouch.cichlid.regime.list[i],tree=curOuchTree, regimes=ouch.cichlid.regime["regimes.riverine"], sqrt.alpha=1, sigma=1, reltol=1e-1, fit=TRUE)
}

# 6. Having fit the OU model, you will get a sigma^2 and alpha value, that we can extract out

Br_OUhaploalpha <- vector("list",6844) # store the sqrt alpha values in this vector
for (i in 1:length(Br_OUhaplo)) {
  Br_OUhaploalpha[[i]] <- Br_OUhaplo[[i]]@sqrt.alpha # this pulls out the square root of alpha
  for (n in 1:length(list.names)) {
    names(Br_OUhaploalpha) <- list.names # add the OG ids to each object
  }
  Br_OUhaploalpha[[i]] <- Br_OUhaploalpha[[i]]^2 # we want just alpha so square each of the values in the list
}

Br_OUhaplosigmasq <- vector("list",6844) # store the sigma values in this vector
for (i in 1:length(Br_OUhaplo)) {
  Br_OUhaplosigmasq[[i]] <- Br_OUhaplo[[i]]@sigma # this pulls out the square root of alpha
  for (n in 1:length(list.names)) {
    names(Br_OUhaplosigmasq) <- list.names # add the OG ids to each object
  }
  Br_OUhaplosigmasq[[i]] <- Br_OUhaplosigmasq[[i]]^2 # we want just alpha so square each of the values in the list
}

# 7. Variance (Expression divergence) is calculated as sigma^2 / 2 * alpha
# Br_OUhaplosigmasq[[1]]/(Br_OUhaploalpha[[1]]*2) # this is the test, single run
Br_OUhaplo_exprdiv <- vector("list",6844)
for (i in 1:length(Br_OUhaplosigmasq)) {
  Br_OUhaplo_exprdiv[[i]] <- Br_OUhaplosigmasq[[i]]/(Br_OUhaploalpha[[i]]*2) # calculates the expression divergence
  for (n in 1:length(list.names)) {
    names(Br_OUhaplo_exprdiv) <- list.names # add the OG ids to each object
  }
}

# 8. For plotting, we do 1-divergence
Br_OUhaplo_exprdiv1 <- vector("list",6844)
for (i in 1:length(Br_OUhaplo_exprdiv)) {
  Br_OUhaplo_exprdiv1[[i]] <- 1-Br_OUhaplo_exprdiv[[i]] # does 1-divergence
  for (n in 1:length(list.names)) {
    names(Br_OUhaplo_exprdiv1) <- list.names # add the OG ids to each object
  }
}

XX <- as.data.frame(Br_OUhaplo_exprdiv1)
XX2 <- t(XX)
XX3 <- (Tree_prom_evolrate[match(rownames(XX2),rownames(Tree_prom_evolrate)),(c(1,2))])
XX4 <- cbind(XX2,XX3)
XX5 <- na.omit(XX4)
colnames(XX5) <- c('exprdiv','4fold','prom')

ggplot(XX5, aes(x = (exprdiv),y = log2(prom))) +
  geom_hex(bins=40) +
  labs (x = "1 - Expression Divergence") + labs (y = "log2 (evolutionary rate in promoter region)")

subset(XX5, exprdiv < -0) # find values with expression divergence < 0


######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~ ######## ~~~~~~~~~~~~~


# This is a function that you could just run in R (not part of 'ouch' package)
# ornstein_uhlenbeck <- function(T,n,nu,lambda,sigma,x0){
#   dw  <- rnorm(n, 0, sqrt(T/n))
#   dt  <- T/n
#   x <- c(x0)
#   for (i in 2:(n+1)) {
#     x[i]  <-  x[i-1] + lambda*(nu-x[i-1])*dt + sigma*dw[i-1] #nu is long run mean, lambda is mean reversion speed
#   }
#   return(x);
# }



########## ALL OF BELOW WILL BE INCORRECT, HERE I USED THE RAW EXPRESSION VALUES (FPKM), DID A LOG2 AND NORMALISED ACROSS THE WHOLE DATA FRAME - IT WOULD TAKE EACH OF THEM AS REPLICATES THEN
########## However, this will be useful for when you do have replicates!

### log2 transformation of raw FPKM expression data
### additional normalisation of log2 transformed FPKM


setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Gene_expression_data/fromChris_notlogtransf/expression_for_Tarang/")
base_dir <- "~/Documents/TGAC/Projects/Cichlid_GRNs/Gene_expression_data/fromChris_notlogtransf/expression_for_Tarang/"


# load expression data
Files <- dir(file.path(base_dir))
#if you get an error when reading files in then make sure that the files are UTF-8 format
for (i in 1:length(Files)){
  tmp = read.delim(file = paste0(base_dir,Files[i]), header = F, row.names = 1)
  assign(Files[i], tmp)
}

# log2 transformation of raw FPKM expression data
Ablog2 <- log2(Ab_expr_noloss_nz_names.geneexp)
Galog2 <- log2(Ga_expr_noloss_nz_names.geneexp)
Mzlog2 <- log2(Mz_expr_noloss_nz_names.geneexp)
Nblog2 <- log2(Nb_expr_noloss_nz_names.geneexp)
Onlog2 <- log2(On_expr_noloss_nz_names.geneexp)
Pnlog2 <- log2(Pn_expr_noloss_nz_names.geneexp)

# Change all -Inf (log2(0)) to NA for normalization below
Ablog2[Ablog2 == -Inf] <- NA
Galog2[Galog2 == -Inf] <- NA
Mzlog2[Mzlog2 == -Inf] <- NA
Nblog2[Nblog2 == -Inf] <- NA
Onlog2[Onlog2 == -Inf] <- NA
Pnlog2[Pnlog2 == -Inf] <- NA

# A function for normalising gene expression datasets in a comparative manner over different sequencing libraries and experiments. It builds on FPKM normalisation and needs FPKM normalised expression values as input.
# Usage: Input values must have NA's instead of 0 (or log2(0) = -Inf) and be log2 transformed. On a gene (rows) by library (columns) matrix, use then as in
# Note that this function also applies the cutoff suggested in the study¹ by removing genes below an expression value of log2(0.125) = -3.
# ¹Hart et al. 2013: Finding the active genes in deep RNA-seq gene expression studies. BMC Genomics 14:778.

# zFPKM function
z_fpkm<-function(i){
  if(all(i>0)) stop('Input not log2 transformed.')
  if(all(!is.na(i))) stop('0\'s need to be NA\'s.')
  my<-density(i,na.rm=T)$x[which.max(density(i,na.rm=T)$y)]
  U<-mean(i[i>my],na.rm=T)
  sigma<-(U-my)*(.5*pi)^.5
  z<-(i-my)/sigma
  z[z< -3]<-NA
  return(z)
}



# usage 'zFPKM_data1 <- apply(FPKM_data, 2, z_fpkm)'
XX <- Ablog2[,1:2]
XX2 <- apply(XX, 2, z_fpkm)

Ablog2_zFPKM <- apply(Ablog2, 2, z_fpkm)
Galog2_zFPKM <- apply(Galog2, 2, z_fpkm)
Mzlog2_zFPKM <- apply(Mzlog2, 2, z_fpkm)
Nblog2_zFPKM <- apply(Nblog2, 2, z_fpkm)
Onlog2_zFPKM <- apply(Onlog2, 2, z_fpkm)
Pnlog2_zFPKM <- apply(Pnlog2, 2, z_fpkm)

# Create a backup of these log2 transformed, normalised expression values
Ablog2_zFPKM_b <- apply(Ablog2, 2, z_fpkm)
Galog2_zFPKM_b <- apply(Galog2, 2, z_fpkm)
Mzlog2_zFPKM_b <- apply(Mzlog2, 2, z_fpkm)
Nblog2_zFPKM_b <- apply(Nblog2, 2, z_fpkm)
Onlog2_zFPKM_b <- apply(Onlog2, 2, z_fpkm)
Pnlog2_zFPKM_b <- apply(Pnlog2, 2, z_fpkm)

# tau is a measure for tissue specificity of gene expression as described by Yanai and colleagues³.
# Usage: Summarise tissue replicates using mean() or similar. For it to function correctly, NA's must be 0. For a gene (rows) by tissue (columns) matrix, apply then as:
# ³Yanai et al. 2004: Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification. Bioinformatics 21:650-659.

# tau does not handle values <0, so for all genes with expression <1, set to O
Ablog2_zFPKM[Ablog2_zFPKM < 1] <- 0
Galog2_zFPKM[Galog2_zFPKM < 1] <- 0
Mzlog2_zFPKM[Mzlog2_zFPKM < 1] <- 0
Nblog2_zFPKM[Nblog2_zFPKM < 1] <- 0
Onlog2_zFPKM[Onlog2_zFPKM < 1] <- 0
Pnlog2_zFPKM[Pnlog2_zFPKM < 1] <- 0

# Change all NA to O
Ablog2_zFPKM[is.na(Ablog2_zFPKM)] <- 0
Galog2_zFPKM[is.na(Galog2_zFPKM)] <- 0
Mzlog2_zFPKM[is.na(Mzlog2_zFPKM)] <- 0
Nblog2_zFPKM[is.na(Nblog2_zFPKM)] <- 0
Onlog2_zFPKM[is.na(Onlog2_zFPKM)] <- 0
Pnlog2_zFPKM[is.na(Pnlog2_zFPKM)] <- 0


# tau function
tau<-function(x){
  if(any(is.na(x))) stop('NA\'s need to be 0.')
  if(any(x<0)) stop('Negative input values not permitted. Maybe data is log transformed?')
  t<-sum(1-x/max(x))/(length(x)-1)
}

# usage 'tau1 <- apply(dataset, 1, tau)' - calculates tau (breadth of expression) across all tissues for each gene
Ablog2_zFPKM_tau <- apply(Ablog2_zFPKM, 1, tau)
Galog2_zFPKM_tau <- apply(Galog2_zFPKM, 1, tau)
Mzlog2_zFPKM_tau <- apply(Mzlog2_zFPKM, 1, tau)
Nblog2_zFPKM_tau <- apply(Nblog2_zFPKM, 1, tau)
Onlog2_zFPKM_tau <- apply(Onlog2_zFPKM, 1, tau)
Pnlog2_zFPKM_tau <- apply(Pnlog2_zFPKM, 1, tau)
# The values of τ vary from 0 to 1; ubiquitous or broad expr (τ ≤ 0.5); intermediate expr (0.5 < τ < 0.9); and tissue-specific or narrow expr (τ ≥ 0.9)






### Calculating expression divergence in tissues between species - using Euclidian distance

# Test on pairwise comparison
# The distance between samples i and j can be written as
# dist(i,j)=(Yi−Yj)⊤(Yi−Yj)
# with Yi and Yj columns i and j. This result can be very convenient in practice as computations can be made much faster using matrix multiplication.

# We can now use the formulas above to compute distance.
# Let’s compute distance between orthologs from two species (Mz and Ab), gene A mz.gene.s81.88 brain and gene B ab.gene.s268.38 brain, and then to the testis expression of each gene.
# In thoery, the same tissues should be closer to one another than different tissues from the same species

geneAb <- Mzlog2_zFPKM_b["mz.gene.s81.88",1]
geneBb <- Ablog2_zFPKM_b["ab.gene.s268.38",1]
geneAt <- Mzlog2_zFPKM_b["mz.gene.s81.88",6]
geneBt <- Ablog2_zFPKM_b["ab.gene.s268.38",6]

sqrt(sum((geneAb-geneBb)^2))
sqrt(sum((geneAb-geneAt)^2))
sqrt(sum((geneBb-geneBt)^2))

# faster computation
sqrt( crossprod(geneAb-geneBb) ) # 0.2303695
sqrt( crossprod(geneAb-geneAt) ) # 1.304473
sqrt( crossprod(geneBb-geneBt) ) # 0.8613025
# brain from two separate species is closer than expr of brain and testis in same species


# THIS NEEDS TO BE BASED ON 1:1 ORTHOLOGS ONLY (USE ROW NAMES AS ORTHOGROUP) - DONE
# log2 transform the data - DONE
# Normalize the data so that each species shows mean 0 and variance 1 - DONE
  # This methodology only stretches and shifts expression values, it does not alter the shape of the distribution.
# PLOT The Distribution of expression level for each species (y axis, Density 0.0 to 0.5; x-axis, Species expression -3 to +3) > are they normally distributed? - Yes, DONE
# Then calculate the euclidian distance of each orthogroup between the species for the different tissues
# These can be plotted, species pairwise against evolutionary rate (of the whole tree for each orthogroup at promoters and dN/dS?) - You can sum the distances and maybe log2 to plot
# Multispecies expression divergence can then be a sum of all the values for each orthgroup, plotted against the (overall rate of that orthgroup in the whole tree?)

# load expression data
Files2 <- dir(file.path(base_dir))
OGID1to1grep <- glob2rx("*OGIDs1to1*") # create a grep pattern to select only specific files
OGID1to1grep2 <- grep(OGID1to1grep, Files2) # run the grep
Files3 <- Files2[OGID1to1grep2] # subset the files from the folder to only select the ones you want using the grep
# load data in
for (i in 1:length(Files3)){
  tmp = read.delim(file = paste0(base_dir,Files3[i]), header = F, row.names = 1)
  assign(Files3[i], tmp)
}


# log2 transformation of raw FPKM expression data
AbOGIDs1to1log2 <- log2(`Ab_expr_noloss_nz_names-OGIDs1to1.geneexp`)
MzOGIDs1to1log2 <- log2(`Mz_expr_noloss_nz_names-OGIDs1to1.geneexp`)
NbOGIDs1to1log2 <- log2(`Nb_expr_noloss_nz_names-OGIDs1to1.geneexp`)
OnOGIDs1to1log2 <- log2(`On_expr_noloss_nz_names-OGIDs1to1.geneexp`)
PnOGIDs1to1log2 <- log2(`Pn_expr_noloss_nz_names-OGIDs1to1.geneexp`)

# Change all -Inf (log2(0)) to NA for normalization below
AbOGIDs1to1log2[AbOGIDs1to1log2 == -Inf] <- NA
MzOGIDs1to1log2[MzOGIDs1to1log2 == -Inf] <- NA
NbOGIDs1to1log2[NbOGIDs1to1log2 == -Inf] <- NA
OnOGIDs1to1log2[OnOGIDs1to1log2 == -Inf] <- NA
PnOGIDs1to1log2[PnOGIDs1to1log2 == -Inf] <- NA

# normalise the expression data using the function at the start of this script

####### NOTE - NEED TO CHECK THIS NORMALISATION, NEED TO TAKE MEAN AS ZERO INSTEAD AS NO REPLICATES - ALL BELOW INCORRECT THEN (19/06/2017)
####### INSTEAD YOU COULD INSTEAD USE THE MODULE EXPRESSION VALUES THAT ARE log(x+1) and then normalized to have mean zero (subtracted out the mean)


AbOGIDs1to1log2_zFPKM <- apply(AbOGIDs1to1log2, 2, z_fpkm)
MzOGIDs1to1log2_zFPKM <- apply(MzOGIDs1to1log2, 2, z_fpkm)
NbOGIDs1to1log2_zFPKM <- apply(NbOGIDs1to1log2, 2, z_fpkm)
OnOGIDs1to1log2_zFPKM <- apply(OnOGIDs1to1log2, 2, z_fpkm)
PnOGIDs1to1log2_zFPKM <- apply(PnOGIDs1to1log2, 2, z_fpkm)

# create density plots for the distribution of expression level in each species
library(reshape2)
library(ggplot2)

# A. burtoni
AbOGIDs1to1log2_zFPKM_l <- melt(AbOGIDs1to1log2_zFPKM)
# Basic histogram with black outline, white fill
ggplot(AbOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(AbOGIDs1to1log2_zFPKM_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('Aburtoni_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(AbOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(A.~burtoni)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous A. burtoni genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.

# M. zebra
MzOGIDs1to1log2_zFPKM_l <- melt(MzOGIDs1to1log2_zFPKM)
# Basic histogram with black outline, white fill
ggplot(MzOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(MzOGIDs1to1log2_zFPKM_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('Mzebra_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(MzOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(M.~zebra)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous M. zebra genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.

# P. nyererei
PnOGIDs1to1log2_zFPKM_l <- melt(PnOGIDs1to1log2_zFPKM)
# Basic histogram with black outline, white fill
ggplot(PnOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(PnOGIDs1to1log2_zFPKM_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('Pnyererei_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(PnOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(P.~nyererei)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous P. nyererei genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.

# N. brichardi
NbOGIDs1to1log2_zFPKM_l <- melt(NbOGIDs1to1log2_zFPKM)
# Basic histogram with black outline, white fill
ggplot(NbOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(NbOGIDs1to1log2_zFPKM_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('Nbrichardi_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(NbOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(N.~brichardi)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous N. brichardi genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.

# O. niloticus
OnOGIDs1to1log2_zFPKM_l <- melt(OnOGIDs1to1log2_zFPKM)
# Basic histogram with black outline, white fill
ggplot(OnOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# Density curve
ggplot(OnOGIDs1to1log2_zFPKM_l, aes(x=value)) + geom_density()
# Histogram overlaid with kernel density curve
tiff('Oniloticus_1to1exprrange.tiff', units="in", width=6, height=6, res=300) # this part will create a high-resolution tiff file
ggplot(OnOGIDs1to1log2_zFPKM_l, aes(x=value)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  labs (x = expression(italic(O.~niloticus)~gene~expression)) + labs (y= ("Density")) +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
dev.off()
# Caption - Distribution of expression level across 6,844 1:1 orthologous O. niloticus genes. Gene expression was determined by log probe intensity, and was normalized to have mean 0 and variance 1.

# Calculate the euclidian distance of each orthogroup between the species means for the different tissues

# Change colnames of all 1to1log2_zFPKM so that they are species and tissue-specific
colnames(AbOGIDs1to1log2_zFPKM) <- c('Ab_Br','Ab_Ey','Ab_Ht','Ab_Kd','Ab_Ms','Ab_Ts')
colnames(MzOGIDs1to1log2_zFPKM) <- c('Mz_Br','Mz_Ey','Mz_Ht','Mz_Kd','Mz_Ms','Mz_Ts')
colnames(NbOGIDs1to1log2_zFPKM) <- c('Nb_Br','Nb_Ey','Nb_Ht','Nb_Kd','Nb_Ms','Nb_Ts')
colnames(OnOGIDs1to1log2_zFPKM) <- c('On_Br','On_Ey','On_Ht','On_Kd','On_Ms','On_Ts')
colnames(PnOGIDs1to1log2_zFPKM) <- c('Pn_Br','Pn_Ey','Pn_Ht','Pn_Kd','Pn_Ms','Pn_Ts')


# cbind all above to compute pairwise distances
All_OGIDs1to1log2_zFPKM <- cbind(MzOGIDs1to1log2_zFPKM,PnOGIDs1to1log2_zFPKM,AbOGIDs1to1log2_zFPKM,NbOGIDs1to1log2_zFPKM,OnOGIDs1to1log2_zFPKM)

# In this formula, the expression data xi and yi are subtracted directly from each other.
# We should therefore make sure that the expression data are properly normalized when using the Euclidean distance, for example by converting the measured gene expression levels to log-ratios (done above)
# Unlike the correlation-based distance functions, the Euclidean distance takes the magnitude of the expression data into account. It, therefore, preserves more information about the data and may be preferable.
# Now to compute all the distances at once, we have the function dist.
# If we were interested in the distance between samples, we would transpose the matrix d <- dist(t(All_OGIDs1to1log2_zFPKM))

# 1. We want to compute pairwise distance across species for each tissue

# load in the whole tree evolutionary rate at promoter regions - this is for plotting
Tree_prom_evolrate <- read.delim(file = "/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Evolutionary_rate_analysis/wholeTree_4fold_prom_evolRates.sorted.out", header = F, row.names = 1)
colnames(Tree_prom_evolrate) <- c('4fold', 'prom')

# Load in the OGIDs to cbind the genes
OGIDs <- read.delim(file = "/Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Module_genesandexpr/OGIDS.txt5", header = F, row.names = 1)

# creating plots of expression divergence vs evolutionary rate
library(ggplot2)
library(hexbin)

# Brain
Br_OGIDs1to1log2_zFPKM <- cbind(All_OGIDs1to1log2_zFPKM[,"Mz_Br"],All_OGIDs1to1log2_zFPKM[,"Pn_Br"],All_OGIDs1to1log2_zFPKM[,"Ab_Br"],All_OGIDs1to1log2_zFPKM[,"Nb_Br"],All_OGIDs1to1log2_zFPKM[,"On_Br"])
colnames(Br_OGIDs1to1log2_zFPKM) <- c('Mz_Br','Pn_Br','Ab_Br','Nb_Br','On_Br')

## To run dist on each row save each row in a list
Br_OGIDs1to1log2_zFPKM.list <- setNames(split(Br_OGIDs1to1log2_zFPKM, seq(nrow(Br_OGIDs1to1log2_zFPKM))), rownames(Br_OGIDs1to1log2_zFPKM))
Br_OGIDs1to1log2_zFPKM.list.ED <- t(sapply(Br_OGIDs1to1log2_zFPKM.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Br_OGIDs1to1log2_zFPKM[1,])
colnames(Br_OGIDs1to1log2_zFPKM.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Br_OGIDs1to1log2_zFPKM.list.ED[is.na(Br_OGIDs1to1log2_zFPKM.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Br_OGIDs1to1log2_zFPKM.list.ED_rowsums <- as.data.frame(rowSums(Br_OGIDs1to1log2_zFPKM.list.ED))
colnames(Br_OGIDs1to1log2_zFPKM.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Br_OGIDs1to1log2_zFPKM.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2 <- cbind(Br_OGIDs1to1log2_zFPKM.list.ED,Br_OGIDs1to1log2_zFPKM.list.ED_rowsums,Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate)
Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3 <- na.omit(Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5 <- cbind(Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3,Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('brain_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Brain expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes")
dev.off()
# In the plot there are several genes with high divergence (log2(6.8) = ~2.5) && high evolutionary rate in promoter region (log2(140) = ~7.1)
  # PKNOX2 (x=2.7,y=7.2)
  # WNT10B


# Eye
Ey_OGIDs1to1log2_zFPKM <- cbind(All_OGIDs1to1log2_zFPKM[,"Mz_Ey"],All_OGIDs1to1log2_zFPKM[,"Pn_Ey"],All_OGIDs1to1log2_zFPKM[,"Ab_Ey"],All_OGIDs1to1log2_zFPKM[,"Nb_Ey"],All_OGIDs1to1log2_zFPKM[,"On_Ey"])
colnames(Ey_OGIDs1to1log2_zFPKM) <- c('Mz_Ey','Pn_Ey','Ab_Ey','Nb_Ey','On_Ey')

## To run dist on each row save each row in a list
Ey_OGIDs1to1log2_zFPKM.list <- setNames(split(Ey_OGIDs1to1log2_zFPKM, seq(nrow(Ey_OGIDs1to1log2_zFPKM))), rownames(Ey_OGIDs1to1log2_zFPKM))
Ey_OGIDs1to1log2_zFPKM.list.ED <- t(sapply(Ey_OGIDs1to1log2_zFPKM.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Ey_OGIDs1to1log2_zFPKM[1,])
colnames(Ey_OGIDs1to1log2_zFPKM.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Ey_OGIDs1to1log2_zFPKM.list.ED[is.na(Ey_OGIDs1to1log2_zFPKM.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums <- as.data.frame(rowSums(Ey_OGIDs1to1log2_zFPKM.list.ED))
colnames(Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2 <- cbind(Ey_OGIDs1to1log2_zFPKM.list.ED,Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums,Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate)
Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3 <- na.omit(Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5 <- cbind(Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3,Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('Eye_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Eye expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes")
dev.off()

# Heart
Ht_OGIDs1to1log2_zFPKM <- cbind(All_OGIDs1to1log2_zFPKM[,"Mz_Ht"],All_OGIDs1to1log2_zFPKM[,"Pn_Ht"],All_OGIDs1to1log2_zFPKM[,"Ab_Ht"],All_OGIDs1to1log2_zFPKM[,"Nb_Ht"],All_OGIDs1to1log2_zFPKM[,"On_Ht"])
colnames(Ht_OGIDs1to1log2_zFPKM) <- c('Mz_Ht','Pn_Ht','Ab_Ht','Nb_Ht','On_Ht')

## To run dist on each row save each row in a list
Ht_OGIDs1to1log2_zFPKM.list <- setNames(split(Ht_OGIDs1to1log2_zFPKM, seq(nrow(Ht_OGIDs1to1log2_zFPKM))), rownames(Ht_OGIDs1to1log2_zFPKM))
Ht_OGIDs1to1log2_zFPKM.list.ED <- t(sapply(Ht_OGIDs1to1log2_zFPKM.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Ht_OGIDs1to1log2_zFPKM[1,])
colnames(Ht_OGIDs1to1log2_zFPKM.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Ht_OGIDs1to1log2_zFPKM.list.ED[is.na(Ht_OGIDs1to1log2_zFPKM.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums <- as.data.frame(rowSums(Ht_OGIDs1to1log2_zFPKM.list.ED))
colnames(Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2 <- cbind(Ht_OGIDs1to1log2_zFPKM.list.ED,Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums,Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate)
Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3 <- na.omit(Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5 <- cbind(Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3,Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('Heart_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Heart expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes")
dev.off()


# Kidney
Kd_OGIDs1to1log2_zFPKM <- cbind(All_OGIDs1to1log2_zFPKM[,"Mz_Kd"],All_OGIDs1to1log2_zFPKM[,"Pn_Kd"],All_OGIDs1to1log2_zFPKM[,"Ab_Kd"],All_OGIDs1to1log2_zFPKM[,"Nb_Kd"],All_OGIDs1to1log2_zFPKM[,"On_Kd"])
colnames(Kd_OGIDs1to1log2_zFPKM) <- c('Mz_Kd','Pn_Kd','Ab_Kd','Nb_Kd','On_Kd')

## To run dist on each row save each row in a list
Kd_OGIDs1to1log2_zFPKM.list <- setNames(split(Kd_OGIDs1to1log2_zFPKM, seq(nrow(Kd_OGIDs1to1log2_zFPKM))), rownames(Kd_OGIDs1to1log2_zFPKM))
Kd_OGIDs1to1log2_zFPKM.list.ED <- t(sapply(Kd_OGIDs1to1log2_zFPKM.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Kd_OGIDs1to1log2_zFPKM[1,])
colnames(Kd_OGIDs1to1log2_zFPKM.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Kd_OGIDs1to1log2_zFPKM.list.ED[is.na(Kd_OGIDs1to1log2_zFPKM.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums <- as.data.frame(rowSums(Kd_OGIDs1to1log2_zFPKM.list.ED))
colnames(Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2 <- cbind(Kd_OGIDs1to1log2_zFPKM.list.ED,Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums,Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate)
Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3 <- na.omit(Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5 <- cbind(Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3,Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('Kidney_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Kd_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Kidney expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes")
dev.off()


# Muscle
Ms_OGIDs1to1log2_zFPKM <- cbind(All_OGIDs1to1log2_zFPKM[,"Mz_Ms"],All_OGIDs1to1log2_zFPKM[,"Pn_Ms"],All_OGIDs1to1log2_zFPKM[,"Ab_Ms"],All_OGIDs1to1log2_zFPKM[,"Nb_Ms"],All_OGIDs1to1log2_zFPKM[,"On_Ms"])
colnames(Ms_OGIDs1to1log2_zFPKM) <- c('Mz_Ms','Pn_Ms','Ab_Ms','Nb_Ms','On_Ms')

## To run dist on each row save each row in a list
Ms_OGIDs1to1log2_zFPKM.list <- setNames(split(Ms_OGIDs1to1log2_zFPKM, seq(nrow(Ms_OGIDs1to1log2_zFPKM))), rownames(Ms_OGIDs1to1log2_zFPKM))
Ms_OGIDs1to1log2_zFPKM.list.ED <- t(sapply(Ms_OGIDs1to1log2_zFPKM.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Ms_OGIDs1to1log2_zFPKM[1,])
colnames(Ms_OGIDs1to1log2_zFPKM.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Ms_OGIDs1to1log2_zFPKM.list.ED[is.na(Ms_OGIDs1to1log2_zFPKM.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums <- as.data.frame(rowSums(Ms_OGIDs1to1log2_zFPKM.list.ED))
colnames(Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2 <- cbind(Ms_OGIDs1to1log2_zFPKM.list.ED,Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums,Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate)
Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3 <- na.omit(Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5 <- cbind(Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3,Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('Muscle_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Ms_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Muscle expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes")
dev.off()


# Testis
Ts_OGIDs1to1log2_zFPKM <- cbind(All_OGIDs1to1log2_zFPKM[,"Mz_Ts"],All_OGIDs1to1log2_zFPKM[,"Pn_Ts"],All_OGIDs1to1log2_zFPKM[,"Ab_Ts"],All_OGIDs1to1log2_zFPKM[,"Nb_Ts"],All_OGIDs1to1log2_zFPKM[,"On_Ts"])
colnames(Ts_OGIDs1to1log2_zFPKM) <- c('Mz_Ts','Pn_Ts','Ab_Ts','Nb_Ts','On_Ts')

## To run dist on each row save each row in a list
Ts_OGIDs1to1log2_zFPKM.list <- setNames(split(Ts_OGIDs1to1log2_zFPKM, seq(nrow(Ts_OGIDs1to1log2_zFPKM))), rownames(Ts_OGIDs1to1log2_zFPKM))
Ts_OGIDs1to1log2_zFPKM.list.ED <- t(sapply(Ts_OGIDs1to1log2_zFPKM.list, dist)) # this has computed the Euclidian distance as you require, you can check first row with dist(Ts_OGIDs1to1log2_zFPKM[1,])
colnames(Ts_OGIDs1to1log2_zFPKM.list.ED) <- c('Mz-Pn','Mz-Ab','Mz-Nb','Mz-On','Pn-Ab','Pn-Nb','Pn-On','Ab-Nb','Ab-On','Nb-On') # colheaders are pairwise comparisons of Mz-Pn, Mz-Ab, Mz-Nb, Mz-On, Pn-Ab, Pn-Nb, Pn-On, Ab-Nb, Ab-On, Nb-On
# change all NA to 0
Ts_OGIDs1to1log2_zFPKM.list.ED[is.na(Ts_OGIDs1to1log2_zFPKM.list.ED)] <- 0

# for plotting against the evolutionary rate, sum each row and log2 the sum
Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums <- as.data.frame(rowSums(Ts_OGIDs1to1log2_zFPKM.list.ED))
colnames(Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums) <- 'ED.rowsum'

# match the evolutionary rate rows to the expression divergence (EDsum)
Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate <- (Tree_prom_evolrate[match(rownames(Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums),rownames(Tree_prom_evolrate)),(c(1,2))])

# Create a combined master table
Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2 <- cbind(Ts_OGIDs1to1log2_zFPKM.list.ED,Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums,Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate)
Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3 <- na.omit(Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate2) # remove rows with NA

# cbind the gene symbols
Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4 <- (OGIDs[match(rownames(Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3),rownames(OGIDs)),(c(9,11,14))])
colnames(Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4) <- c('Ga_gene','Dr_gene','Hs_gene')
Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5 <- cbind(Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate3,Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate4)

# In this plot, x-axis increased expression divergence moving left to right, y-axis increased evolutionary rate up the scale
tiff('Testis_euclidian_prom_evol.tiff', units="in", width=12, height=12, res=300)
ggplot(Ts_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5, aes(x = log2(ED.rowsum),y = log2(prom))) +
  geom_hex(bins=45) +
  labs (x = "log2 (sum of Euclidian distance)") + labs (y = "log2 (evolutionary rate in promoter region)") + labs(title="Testis expression divergence as Euclidian distance and evolutionary rate in 1:1 orthologous genes")
dev.off()



### Look at the R-squared for correlation between expression divergence and sequence conservation (as evolutionary rate)
# R-squared is a statistical measure of how close the data are to the fitted regression line.
brain.lm = lm(ED.rowsum ~ prom, data=Br_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5)
summary(brain.lm)$r.squared # 0.000512742
eye.lm = lm(ED.rowsum ~ prom, data=Ey_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5)
summary(eye.lm)$r.squared # 6.218064e-05
heart.lm = lm(ED.rowsum ~ prom, data=Ht_OGIDs1to1log2_zFPKM.list.ED_rowsums.evolrate5)
summary(heart.lm)$r.squared # 0.0009915897
