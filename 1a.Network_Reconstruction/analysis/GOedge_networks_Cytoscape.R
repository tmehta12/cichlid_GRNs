########################################################################################################################################################################################################
#
# This script
# 1. loads several GO edge interactions into cytoscape
# 2. calculates network properties
# 3. implements these properties into networks and changes the style for visualisations

########################################################################################################################################################################################################

### 0. Install relevant packages
source("https://bioconductor.org/biocLite.R")

biocLite("igraph")
biocLite("RJSONIO")
biocLite("httr")

### 1. Basic setup
library(igraph)
library(RJSONIO)
library(httr)

port.number = 1234
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
print(base.url)

### 2. Loading Networks
# There are many ways to import network data into igraph. In this example, let’s try to load text data as Data Frame, and then convert it into igraph object.

# Set working directory
setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/")

# Load list of edges as Data Frame
network.df <- read.table("GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module0.txt_Final.sif")

# Convert it into igraph object
network <- graph.data.frame(network.df,directed=F)

# Remove duplicate edges & loops
g.tca <- simplify(network, remove.multiple=T, remove.loops=T)

# Name it
g.tca$name = "Ab_GO_Module0"

### 3. Convert igraph object into JSON
# basic data exchange format between Cytoscape and external tools is Cytoscape.js JSON. To send igraph object to Cytoscape, you need to convert it into JSON:
# This function will be published as a part of utility package, but not ready yet.
source('cytoscape_util.R')

# Convert it into Cytosccape.js JSON
cygraph <- toCytoscape(g.tca)

### 4. Visualize it in Cytoscape
# Now you are ready to visualize your data in Cytoscape.
send2cy(cygraph, 'default%20white', 'circular')
# If you want to visualize your network in igraph, obviously it is possible:
# plot(g.tca,vertex.label=V(g.tca)$name)
