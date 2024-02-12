library(RCX)
library(parallel)
library(data.table)
library(igraph)
setwd("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/Tools/NBS2_NetPropagation")

### Preparation folder
dir.create("./data/",showWarnings = F)
dir.create("./data/somaticAlterations/",showWarnings = F)
dir.create("./ppi_networks/",showWarnings = F)
dir.create("./results/",showWarnings = F)

### Open the cx network from the Ideker paper about PCNet
rcx = readCX("./ppi_networks/Parsimonious Composite Network (PCNet).cx")
# and convert into igraph
ig <- toIgraph(rcx)
save(ig,file="data/PCNet_igraph.RData")

### Filter networks based on genes present in PathwayCommons db and save the list of common genes
all_PC = fread("./data/PathwayCommons11.All.hgnc.txt")

genes_PC = unique(c(all_PC$PARTICIPANT_A,all_PC$PARTICIPANT_B))
commonGenes_PCnet = data.frame(gene = intersect(genes_PC,V(ig)$nodeName))
fwrite(commonGenes_PCnet,"./data/commonGenesPCnet.tsv",sep="\t")