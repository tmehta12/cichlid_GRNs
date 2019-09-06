### Plotting miRNA targets GO enrichment

library(ggplot2)
library(hexbin)
library(reshape2)
library(grid)

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/6.miRNAs/GOOUTPUT/")
base_dir <-("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/6.miRNAs/GOOUTPUT/")

# load in the data 
Files1 <- dir(file.path(base_dir))
miRNA_GOgrep <- glob2rx("*GOOUTPUT_details_filtered.txt3") # create a grep pattern to select only specific files
miRNA_GOgrep2 <- grep(miRNA_GOgrep, Files1) # run the grep
Files2 <- Files1[miRNA_GOgrep2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files2)){
  tmp = read.delim(file = paste0(base_dir,Files2[i]),header = F)
  assign(Files2[i], tmp)
}

# these are for individual plots
Pn_topmiRNA_GO <- ggplot(pn_miRNAtop10targets_GOOUTPUT_details_filtered.txt3, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + labs(title=expression(italic(P.~nyererei))) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="miRNA")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

Nb_topmiRNA_GO <- ggplot(nb_miRNAtop10targets_GOOUTPUT_details_filtered.txt3, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + labs(title=expression(italic(N.~brichardi))) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="miRNA")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

Mz_topmiRNA_GO <- ggplot(mz_miRNAtop10targets_GOOUTPUT_details_filtered.txt3, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + labs(title=expression(italic(M.~zebra))) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="miRNA")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

On_topmiRNA_GO <- ggplot(on_miRNAtop10targets_GOOUTPUT_details_filtered.txt3, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + labs(title=expression(italic(O.~niloticus))) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="miRNA")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

Ab_topmiRNA_GO <- ggplot(ab_miRNAtop10targets_GOOUTPUT_details_filtered.txt3, aes(x = reorder(V2, -V9), y = -log10(V9), fill=V1)) + geom_bar(width=.7, position=position_dodge(), stat = "identity") + scale_colour_hue(l=40) +
  theme_bw() + labs (x = "GO Term") + labs (y = "log10 fold enrichment (FDR < 0.05)") + labs(title=expression(italic(A.~burtoni))) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="miRNA")) + theme(panel.grid = element_blank()) + coord_flip() # + geom_text(aes( label = V8), position=position_dodge(width=1), hjust = -0.5, size=1, group='V1')

# this is to create a combined plot that we will use
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/FigS-R2h_GOenrichment_miRNA_targets.tiff', units="in", width=36, height=16, res=300)
multiplot(Mz_topmiRNA_GO,Pn_topmiRNA_GO,Ab_topmiRNA_GO,Nb_topmiRNA_GO,On_topmiRNA_GO,cols=5)
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