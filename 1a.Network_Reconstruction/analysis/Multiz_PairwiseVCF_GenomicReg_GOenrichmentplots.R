### Plotting GO enrichment of genic regions intersecting pairwise variants from multiz MSA

library(ggplot2)
library(hexbin)
library(reshape2)
library(gridExtra)
library(grid)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/")
base_dir <-("~/Documents/TGAC/Projects/Cichlid_GRNs/functional_variant_analysis/")

# load in the data 
Files1 <- dir(file.path(base_dir))
VCF_GOgrep <- glob2rx("*GOOUTPUT_details_filtered.txt3") # create a grep pattern to select only specific files
VCF_GOgrep2 <- grep(VCF_GOgrep, Files1) # run the grep
Files2 <- Files1[VCF_GOgrep2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files2)){
  tmp = read.delim(file = paste0(base_dir,Files2[i]),header = F)
  assign(Files2[i], tmp)
}

# these are for individual plots

PromTFBS_GO <- ggplot(`PromTFBS_collated-GOOUTPUT_details_filtered.txt3`, aes(x = reorder(V2, -V9), y = -log2(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="Promoter TFBSs") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

promTFBSnosw_GO <- ggplot(`promTFBSnosw_collated-GOOUTPUT_details_filtered.txt3`, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="No switched promoter TFBSs") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

ThreeUTRmiRNA_GO <- ggplot(`3UTRmiRNA_collated-GOOUTPUT_details_filtered.txt3`, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="3' UTR miRNA binding sites") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

CNE_GO <- ggplot(`CNE_collated-GOOUTPUT_details_filtered.txt3`, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="CNEs") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

PromNosw_GO <- ggplot(`PromNosw_collated-GOOUTPUT_details_filtered.txt3`, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="Non switched promoters") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

ThreeUTR_GO <- ggplot(`3UTR_collated-GOOUTPUT_details_filtered.txt3`, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="3' UTR") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

aCNE_GO <- ggplot(`aCNE_collated-GOOUTPUT_details_filtered.txt3`, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="aCNEs") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

Promsw_GO <- ggplot(`Promsw_collated-GOOUTPUT_details_filtered.txt3`, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="Switched promoters") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

# no enrichment for promTFBSsw
# promTFBSsw_GO <- ggplot(`promTFBSsw_collated-GOOUTPUT_details_filtered.txt3`, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
#   theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + theme(plot.title = element_text(size=28,face="bold")) + theme(axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold")) + theme(legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold")) + labs(title="Switched promoter TFBSs") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Pairwise species comparison")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')


# labs(title=expression(italic(A.~burtoni))) 

# as plots are all diferent size, do individual ones
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2la–GO_enrichment_SNPs_overlapping_PromNosw_log10.tiff', units="in", width=18, height=8, res=300)
PromNosw_GO
dev.off()
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2lb–GO_enrichment_SNPs_overlapping_Promsw_log10.tiff', units="in", width=12, height=6, res=300)
Promsw_GO
dev.off()
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2lc–GO_enrichment_SNPs_overlapping_PromTFBS_log10.tiff', units="in", width=18, height=8, res=300)
PromTFBS_GO
dev.off()
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2ld–GO_enrichment_SNPs_overlapping_promTFBSnosw_log10.tiff', units="in", width=18, height=8, res=300)
promTFBSnosw_GO
dev.off()
# no enrichment for promTFBSsw
# tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2le–GO_enrichment_SNPs_overlapping_promTFBSsw_log10.tiff', units="in", width=36, height=22, res=300)
# promTFBSsw_GO
# dev.off()
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2le–GO_enrichment_SNPs_overlapping_CNE_log10.tiff', units="in", width=28, height=14, res=300)
CNE_GO
dev.off()
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2lf–GO_enrichment_SNPs_overlapping_aCNE_log10.tiff', units="in", width=18, height=8, res=300)
aCNE_GO
dev.off()
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2lg–GO_enrichment_SNPs_overlapping_ThreeUTR_log10.tiff', units="in", width=14, height=6, res=300)
ThreeUTR_GO
dev.off()
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2lh–GO_enrichment_SNPs_overlapping_ThreeUTRmiRNA_log10.tiff', units="in", width=18, height=8, res=300)
ThreeUTRmiRNA_GO
dev.off()

# this is to create a combined plot - may not be best here!

plotlist <- list(PromNosw_GO,Promsw_GO,PromTFBS_GO,promTFBSnosw_GO,promTFBSsw_GO,CNE_GO,aCNE_GO,ThreeUTR_GO,ThreeUTRmiRNA_GO)

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2l–GO_enrichment_SNPs_overlapping_regulatory_regions_log10.tiff', units="in", width=72, height=56, res=300)
grid.arrange(grobs=plotlist[c(1:9)],ncol=3,nrow=3)
dev.off()

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2l–GO_enrichment_SNPs_overlapping_regulatory_regions_log10.tiff', units="in", width=36, height=22, res=300)
multiplot(PromNosw_GO,Promsw_GO,PromTFBS_GO,promTFBSnosw_GO,promTFBSsw_GO,CNE_GO,aCNE_GO,ThreeUTR_GO,ThreeUTRmiRNA_GO,cols=3)
dev.off()

# # this is to create a combined plot that we will use
# pdf('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2i_GOenrichment_miRNA_targets.pdf', width=36, height=16)
# multiplot(Mz_VCF_GO,Pn_VCF_GO,Ab_VCF_GO,Nb_VCF_GO,On_VCF_GO,cols=5)
# dev.off()

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