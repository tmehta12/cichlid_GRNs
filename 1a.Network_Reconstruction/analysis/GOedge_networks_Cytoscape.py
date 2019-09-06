#!/usr/bin/env python

########################################################################################################################################################################################################
#
# This script IS INCOMPLETE
# CURRENTLY IT ONLY
# 1. loads several GO edge interactions into cytoscape
# then you need to manually do the other things > run NetworkAnalyzer
		# 1. map node size to 'degree': low values to small sizes
		# 2. edge size according to 'betweenness': low values to small sizes
		# 3. map node color to 'BetweennessCentrality': low values to bright colours
		# 4. map edge colour to 'NumberOfUnderlyingEdges': low values to bright colours
		# 5. remove duplicated edges and self-loops
		# 6. style > Node > label color > black
		# 7. style > Node > Border Width > 2.0
		# 8. style > Node > Border Paint > Grey
		# 9. circular/organic layout

########################################################################################################################################################################################################

### 0. Import the required packages

import requests
import json
from IPython.display import display
from IPython.display import Image

# Basic Setup
PORT_NUMBER = 1234
BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'
HEADERS = {'Content-Type': 'application/json'}

# Utility function to print result (JSON Printer)
def jp(data):
    print(json.dumps(data, indent=4))

### 1. Start from scratch: Delete current session
res = requests.delete(BASE + 'session')
# jp(res.json())

### 2. Load a network from file / URL
# URL Parameters
url_params = {
    'source': 'url',
    'collection': 'edgeGO'
}

# Array of data source.  URL of the file (remote or local)
network_files = [
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
#'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module3.txt_Final.sif', - empty
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module5.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module6.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module9.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
#'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module3.txt_Final.sif', - empty
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module5.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module6.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module9.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module3.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
#'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module5.txt_Final.sif', - empty
#'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module6.txt_Final.sif', - empty
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module9.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
#'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module3.txt_Final.sif', - empty
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module5.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module6.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module9.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module3.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module5.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module6.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module9.txt_Final.sif'
]

# Load network from URLs
res = requests.post(BASE + 'networks', params=url_params, data=json.dumps(network_files), headers=HEADERS)
# jp(res.json())

# suid = res.json()[0]['networkSUID'][0]

# Make a utility to get first SUID
# def get_suid(response):
#     return res.json()[0]['networkSUID'][0]

########################################################################################################################################################################################################

## test
res = requests.delete(BASE + 'session')
network_files = [
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
]
res = requests.post(BASE + 'networks', params=url_params, data=json.dumps(network_files), headers=HEADERS)
suid = res.json()[0]['networkSUID'][0]
# 0. calculate network properties
import networkx as nx
deg = nx.degree(res)
btw = nx.betweenness_centrality(res)

# 1. map node size to 'degree': low values to small sizes
nx.set_node_attributes(g, 'degree', deg)

# 2. edge size according to 'betweenness': low values to small sizes
nx.set_edge_attributes(g, 'betweenness', btw)
# 3. map node color to 'BetweennessCentrality': low values to bright colours
# 4. map edge colour to 'NumberOfUnderlyingEdges': low values to bright colours
# 5. remove duplicated edges and self-loops
# 6. style > Node > label color > black
# 7. style > Node > Border Width > 2.0
# 8. style > Node > Border Paint > Grey
# 9. organic layout



########################################################################################################################################################################################################

import os
import sys
from time import sleep
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from py2cytoscape import cyrest
import requests
import json
import networkx as nx

### 1. Define required functions

### Setup Cytoscape
# * Launch Cytoscape on your local machine.
# * Leave Cytoscape running in the background during the remainder of the tutorial.
# * Check cytoscape connection.

cytoscape=cyrest.cyclient()
cytoscape.version()


### 2. Create networks from local files

# HTTP Client for Python
import requests
# Standard JSON library
import json
# Basic Setup
PORT_NUMBER = 1234
# Specify your machine's URL (The IP address of the machine running Cytoscape and cyREST) if you use
# Docker or remote server for this notebook.
#IP = '192.168.1.1'
# If you are running both Notebook server and Cytoscape on a same machine, just use localhost
IP = 'localhost'
BASE = 'http://' + IP + ':' + str(PORT_NUMBER) + '/v1/'
# Header for posting data to the server as JSON
HEADERS = {'Content-Type': 'application/json'}

# The POST method is used to create new Cytoscape objects. For example,
# POST http://localhost:1234/v1/networks
# means create new network(s) by specified method. If you want to create networks from files on your machine or remote servers, all you need to do is create a list of file locations and post it to Cytoscape.

# a. Create a small utility function to create networks from list of local file URLs
def create_from_list(network_list):
    server_res = requests.post(BASE + 'networks?source=url&collection=Edge_GO_enrichment', data=json.dumps(network_list), headers=HEADERS)
    return server_res.json()

# b. create a list of the network files to load
network_files=[
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module3.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module5.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module6.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module9.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module3.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module5.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module6.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Mz-Edge_Attributes_Collated4d-Module9.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module3.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module5.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module6.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Nb-Edge_Attributes_Collated4d-Module9.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module3.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module5.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module6.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_On-Edge_Attributes_Collated4d-Module9.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module0.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module1.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module2.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module3.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module4.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module5.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module6.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module7.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module8.txt_Final.sif',
'file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Pn-Edge_Attributes_Collated4d-Module9.txt_Final.sif',
]

# c. create the networks
print(json.dumps(create_from_list(network_files), indent=4))


pathway_location = "file:///Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/Edge_Attributes/1.EdgeGOenrichment/GeneSymbols2_Ab-Edge_Attributes_Collated4d-Module0.txt_Final.sif"
res1 = requests.post(BASE + "networks?source=url", data=json.dumps([pathway_location]), headers=HEADERS)
result = json.loads(res1.content)
pathway_suid = result[0]["networkSUID"][0]
print("Pathway SUID = " + str(pathway_suid))




# network properties
cy = CyRestClient(ip='127.0.0.1', port=1234)
import networkx as nx
g = nx.scale_free_graph(500)
deg = nx.degree(g)
btw = nx.betweenness_centrality(g)

nx.set_node_attributes(g, 'degree', deg)
nx.set_node_attributes(g, 'betweenness', btw)

g_cy = cy.network.create_from_networkx(g)

cy.layout.apply(name='kamada-kawai', network=g_cy)

directed = cy.style.create('Directed')
cy.style.apply(directed, network=g_cy)

result = cy.edgebundling.apply(g_cy)

# png
network_png = g_cy.get_png()
from IPython.display import Image
Image(network_png)

# svg
network_svg = g_cy.get_svg()
from IPython.display import SVG
SVG(network_svg)

# pdf
network_pdf = g_cy.get_pdf()
f = open('scale_free_500.pdf', 'wb')
f.write(network_pdf)
f.close()

# javascript
import py2cytoscape.cytoscapejs as renderer
view = g_cy.get_first_view()
# style_for_widget = cy.style.get(my_yeast_style.get_name(), data_format='cytoscapejs')
renderer.render(view, 'Directed', background='radial-gradient(#FFFFFF 15%, #DDDDDD 105%)')
