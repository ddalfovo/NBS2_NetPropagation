library(data.table)
setwd("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/Tools/NBS2_NetPropagation")


### Reduce PPI
listTumors = c('BLCA', 'BRCA', 'CESC', 'COADREAD', 'ESCA', 'GBM', 'HNSC', 'LGG', 'PRAD', 'SARC', 'STAD', 'TGCT', 'UCEC')
# tumor = 'BRCA'

for(tumor in listTumors){
  ### PCNet
  folderSomaticAlt = "somaticAlterations"
  PathwayCommon = fread(paste0("./data/",folderSomaticAlt,"/",tumor,"/",tumor,"_edge2features_all.txt"),data.table = F)
  signPC = paste0(PathwayCommon$V1,":",PathwayCommon$V2)
  
  df = as_long_data_frame(ig)
  signPCNet = paste0(df$from_nodeName,":",df$to_nodeName)
  signPCNet_rev = paste0(df$to_nodeName,":",df$from_nodeName)
  PathwayCommon_PCNet = PathwayCommon[signPC%in%c(signPCNet,signPCNet_rev),]
  fwrite(PathwayCommon_PCNet,paste0("./data/",folderSomaticAlt,"/",tumor,"/",tumor,"_edge2features_PCNet.txt"),sep="\t",col.names = F)
}
