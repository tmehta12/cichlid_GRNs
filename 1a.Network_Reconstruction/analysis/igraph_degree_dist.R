# The degree of a node refers to the number of links associated with a node.
# Degree can be measured as the links going in ("in degree"), out ("out degree"), or both.
# The degree() function takes a graph input and gives the degree of specified nodes.
# With the argument "v=V(graph)" you tell the function to give the degree of all nodes in the graph,
# while the "mode" argument specifies in, out, or both.

library("igraph")
# library("poweRlaw")
library("ggplot2")

setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/igraph/")
base_dir <-("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/")

# Load in the data
Files1 <- dir(file.path(base_dir))
Edge_grep <- glob2rx("*-Edge_Attributes_Collated4d_CUT.txt") # create a grep pattern to select only specific files
Edge_grep2 <- grep(Edge_grep, Files1) # run the grep
Files2 <- Files1[Edge_grep2] # subset the files from the folder to only select the ones you want using the grep

for (i in 1:length(Files2)){
  tmp = read.delim(file = paste0(base_dir,Files2[i]),header = T)
  assign(Files2[i], tmp)
}

# Convert edge lists to igraph graph object
Ab <- graph.data.frame(`Ab-Edge_Attributes_Collated4d_CUT.txt`,directed=TRUE)
Mz <- graph.data.frame(`Mz-Edge_Attributes_Collated4d_CUT.txt`,directed=TRUE)
Pn <- graph.data.frame(`Pn-Edge_Attributes_Collated4d_CUT.txt`,directed=TRUE)
Nb <- graph.data.frame(`Nb-Edge_Attributes_Collated4d_CUT.txt`,directed=TRUE)
On <- graph.data.frame(`On-Edge_Attributes_Collated4d_CUT.txt`,directed=TRUE)

# List of degrees
Ab.degrees <- degree(Ab)
Mz.degrees <- degree(Mz)
Pn.degrees <- degree(Pn)
Nb.degrees <- degree(Nb)
On.degrees <- degree(On)

# Count the frequencies of each degree
Ab.degree.histogram <- as.data.frame(table(Ab.degrees))
Mz.degree.histogram <- as.data.frame(table(Mz.degrees))
Pn.degree.histogram <- as.data.frame(table(Pn.degrees))
Nb.degree.histogram <- as.data.frame(table(Nb.degrees))
On.degree.histogram <- as.data.frame(table(On.degrees))

# Need to convert the first column to numbers, otherwise log-log will not work
Ab.degree.histogram[,1] <- as.numeric(Ab.degree.histogram[,1])
Mz.degree.histogram[,1] <- as.numeric(Mz.degree.histogram[,1])
Pn.degree.histogram[,1] <- as.numeric(Pn.degree.histogram[,1])
Nb.degree.histogram[,1] <- as.numeric(Nb.degree.histogram[,1])
On.degree.histogram[,1] <- as.numeric(On.degree.histogram[,1])

# Now, plot it!
tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3Ac–Ab_edge_degree_dist.tiff', units="in", width=10, height=10, res=300)
ggplot(Ab.degree.histogram, aes(x = Ab.degrees, y = Freq)) +
  geom_point() +
  scale_x_continuous("log10 degree\n(- nodes with this amount of connections)",
                     breaks = c(1, 3, 10, 30, 100, 300, 600),
                     trans = "log10") +
  scale_y_continuous("log10 frequency",
                     breaks = c(1, 3, 10, 30, 100, 300, 1000),
                     trans = "log10") +
  ggtitle(expression(paste(italic("A. burtoni"), " regulatory network edge degree distribution"))) +
theme_bw()
dev.off()

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3Aa–Mz_edge_degree_dist.tiff', units="in", width=10, height=10, res=300)
ggplot(Mz.degree.histogram, aes(x = Mz.degrees, y = Freq)) +
  geom_point() +
  scale_x_continuous("log10 degree\n(- nodes with this amount of connections)",
                     breaks = c(1, 3, 10, 30, 100, 300),
                     trans = "log10") +
  scale_y_continuous("log10 frequency",
                     breaks = c(1, 3, 10, 30, 100, 300, 1000),
                     trans = "log10") +
  ggtitle(expression(paste(italic("M. zebra"), " regulatory network edge degree distribution"))) +
theme_bw()
dev.off()

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3Ab–Pn_edge_degree_dist.tiff', units="in", width=10, height=10, res=300)
ggplot(Pn.degree.histogram, aes(x = Pn.degrees, y = Freq)) +
  geom_point() +
  scale_x_continuous("log10 degree\n(- nodes with this amount of connections)",
                     breaks = c(1, 3, 10, 30, 100, 300),
                     trans = "log10") +
  scale_y_continuous("log10 frequency",
                     breaks = c(1, 3, 10, 30, 100, 300, 1000),
                     trans = "log10") +
  ggtitle(expression(paste(italic("P. nyererei"), " regulatory network edge degree distribution"))) +
theme_bw()
dev.off()

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3Ad–Nb_edge_degree_dist.tiff', units="in", width=10, height=10, res=300)
ggplot(Nb.degree.histogram, aes(x = Nb.degrees, y = Freq)) +
  geom_point() +
  scale_x_continuous("log10 degree\n(- nodes with this amount of connections)",
                     breaks = c(1, 3, 10, 30, 100, 300),
                     trans = "log10") +
  scale_y_continuous("log10 frequency",
                     breaks = c(1, 3, 10, 30, 100, 300, 1000),
                     trans = "log10") +
  ggtitle(expression(paste(italic("N. brichardi"), " regulatory network edge degree distribution"))) +
theme_bw()
dev.off()

tiff('~/Documents/TGAC/Projects/Cichlid_GRNs/Manuscript/Figures/New/Supplementary/networks/FigS-R3Ae–On_edge_degree_dist.tiff', units="in", width=10, height=10, res=300)
ggplot(On.degree.histogram, aes(x = On.degrees, y = Freq)) +
  geom_point() +
  scale_x_continuous("log10 degree\n(- nodes with this amount of connections)",
                     breaks = c(1, 3, 10, 30, 100, 300),
                     trans = "log10") +
  scale_y_continuous("log10 frequency",
                     breaks = c(1, 3, 10, 30, 100, 300, 1000),
                     trans = "log10") +
  ggtitle(expression(paste(italic("O. niloticus"), " regulatory network edge degree distribution"))) +
theme_bw()
dev.off()

