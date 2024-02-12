library(RCX)
library(parallel)
library(data.table)
library(igraph)
setwd("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/Tools/NBS2_NetPropagation")

### Then Run in order the following scripts
# Run somaticAlterations.R


# Run data_processing.py

### Reduce PPI
listTumors = c('BLCA', 'BRCA', 'CESC', 'COADREAD', 'ESCA', 'GBM', 'HNSC', 'LGG', 'PRAD', 'SARC', 'STAD', 'TGCT', 'UCEC')
tumor = 'BRCA'

for(tumor in listTumors){
  ### PCNet
  folderSomaticAlt = "somaticMutationsProfiles_PCNet"
  PathwayCommon = fread(paste0("./data/",folderSomaticAlt,"/",tumor,"/",tumor,"_edge2features_all.txt"),data.table = F)
  signPC = paste0(PathwayCommon$V1,":",PathwayCommon$V2)
  
  df = as_long_data_frame(ig)
  signPCNet = paste0(df$from_nodeName,":",df$to_nodeName)
  signPCNet_rev = paste0(df$to_nodeName,":",df$from_nodeName)
  PathwayCommon_PCNet = PathwayCommon[signPC%in%c(signPCNet,signPCNet_rev),]
  fwrite(PathwayCommon_PCNet,paste0("./data/",folderSomaticAlt,"/",tumor,"/",tumor,"_edge2features_PCNet.txt"),sep="\t",col.names = F)
}

for(tumor in listTumors){
  ### MasterNet
  folderSomaticAlt = "somaticMutationsProfiles_MasterNet"
  PathwayCommon = fread(paste0("./data/",folderSomaticAlt,"/",tumor,"/",tumor,"_edge2features_all.txt"),data.table = F)
  signPC = paste0(PathwayCommon$V1,":",PathwayCommon$V2)
  
  df = as_long_data_frame(net)
  df$from_nodeName = nodes[df$from]
  df$to_nodeName = nodes[df$to]
  signMasterNet = c(paste0(df$from_nodeName,":",df$to_nodeName),paste0(df$to_nodeName,":",df$from_nodeName))
  PathwayCommon_MasterNet = PathwayCommon[signPC%in%signMasterNet,]
  fwrite(PathwayCommon_MasterNet,paste0("./data/",folderSomaticAlt,"/",tumor,"/",tumor,"_edge2features_MasterNet.txt"),sep="\t",col.names = F)
}
# Run NBS_tumorSpecific Opt for optimized the hyperparameters
# Run NBS_tumorSpecific
