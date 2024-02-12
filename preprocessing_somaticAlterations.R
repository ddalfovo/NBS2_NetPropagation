library(data.table)
library(plyr)
library(readxl)
require(dplyr) 
require(tidyr)
library(parallel)
library(pbmcapply)
library(igraph)
library(tidymodels)
setwd("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/Tools/NBS2_NetPropagation")
# load("./data/cBioPortalAlterations.RData")

#############################################################################################
OG = fread("/shares/CIBIO-Storage/BCG/scratch1/Resources/geneLists/Oncogenes.txt",header=F)
TSG = fread("/shares/CIBIO-Storage/BCG/scratch1/Resources/geneLists/TumorSupp.txt",header=F)
cancerGenes = fread('/shares/CIBIO-Storage/BCG/scratch1/Resources/geneLists/CancerGenesList.txt', header = F)
colnames(OG) = colnames(cancerGenes) = colnames(TSG) = 'Hugo Symbol'

EUR = fread("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/GWAS_PGS/data/binary_ped/EUR_CancerCell2020/TCGA_Affy_EUR_20201124.fam")
# EUR

final_CN = fread(file = './data/somaticAlterations/gisticCN_PCNet.tsv',data.table = F)
final_CN = final_CN[gsub("-..$","",final_CN$SAMPLE)%in%gsub("-...-..$","",EUR$V2),]
final_MUT = fread(file = './data/somaticAlterations/mutations_PCNet.tsv',data.table = F)
final_MUT = final_MUT[gsub("-..$","",final_MUT$SAMPLE)%in%gsub("-...-..$","",EUR$V2),]
folderSomaticAlt = "somaticAlterations"

final_CN = final_CN[,-grep("SAMPLE.|DATASET.",colnames(final_CN))]
final_MUT = final_MUT[,-grep("SAMPLE.|DATASET.",colnames(final_MUT))]

tmp = as.data.frame(final_CN[,3:ncol(final_CN)])
res = mclapply(1:ncol(tmp),function(col){
  col.name = colnames(tmp)[col]
  col.values = tmp[,col]
  col.values[is.na(col.values)] <- 0
  col.values = replace(col.values,col.values%in%c(-1,1),c(0))
  # out.values = replace(col.values,col.values%in%c(-2,2),c(1))
  
  if(col.name%in%OG$`Hugo Symbol`&col.name%in%TSG$`Hugo Symbol`){
    out.values = replace(col.values,col.values%in%c(-2,2),c(1))
  } else if(col.name%in%TSG$`Hugo Symbol`) {
    out.values = replace(col.values,col.values%in%c(2),c(0))
    out.values = replace(out.values,col.values%in%c(-2),c(1))
  } else if(col.name%in%OG$`Hugo Symbol`) {
    out.values = replace(col.values,col.values%in%c(-2),c(0))
    out.values = replace(out.values,col.values%in%c(2),c(1))
  } else if(col.name%in%cancerGenes$`Hugo Symbol`&!(col.name%in%OG$`Hugo Symbol`|col.name%in%TSG$`Hugo Symbol`)){
    out.values = replace(col.values,col.values%in%c(-2,2),c(1))
  } else {
    out.values = replace(col.values,col.values%in%c(-2,2),c(1))
  }
  tmpNum = as.numeric(out.values)
  tmpNum[is.na(tmpNum)] <- 0
  tmpNum
},mc.cores = 50)
out_res = as.data.frame(do.call(cbind,res))
colnames(out_res) = colnames(tmp)
final_CN = cbind(final_CN[,1:2],out_res)
# sum(rowSums(final_CN[,-c(1:2)])>0)

#### MUT
tmp = as.data.frame(final_MUT[,3:ncol(final_MUT)])
res = mclapply(1:ncol(tmp),function(col){
  col.name = colnames(tmp)[col]
  col.values = tmp[,col]
  col.values[is.na(col.values)] <- 0
  
  as.numeric(as.character(col.values))
},mc.cores = 50)
out_res = (do.call(cbind,res))
colnames(out_res) = colnames(tmp)
final_MUT = cbind(final_MUT[,1:2],out_res)

# sum(rowSums(final_MUT[,-c(1:2)])>0)
final_CN_MUT = cbind(final_CN[,1:2],
                     (final_CN[,3:ncol(final_CN)])|(final_MUT[,3:ncol(final_MUT)]))


final_CN_MUT[final_CN_MUT==F] <- 0
final_CN_MUT[final_CN_MUT==T] <- 1

final_CN_MUT[,2] = gsub("-..$","",final_CN_MUT[,2])
final_CN_MUT = final_CN_MUT[!duplicated(final_CN_MUT[,2]),]
alterationsAll = final_CN_MUT[,-1]
labelAll = final_CN_MUT[,2:1]

colnames(alterationsAll)[1]="SAMPLE"


fwrite(alterationsAll,file=paste0("./data/",folderSomaticAlt,"/TCGA_alterations_all.tsv"),sep="\t")
fwrite(labelAll,file=paste0("./data/",folderSomaticAlt,"/TCGA_labels_all.tsv"),sep="\t")

# Split in one fold
idx = sample(1:nrow(labelAll))
# trainN = round(nrow(labelAll) * 2/3)
# valN = round(nrow(labelAll) * 1/3)
alterationsAll = alterationsAll[idx,]
labelAll = labelAll[idx,]

# num_groups = 3
# alterationsAll_splitted = alterationsAll %>% 
#   group_by((row_number()-1) %/% (n()/num_groups)) %>%
#   nest %>% pull(data)
# labelAll_splitted = labelAll %>% 
#   group_by((row_number()-1) %/% (n()/num_groups)) %>%
#   nest %>% pull(data)

## New split test
res = lapply(split(labelAll,labelAll$DATASET),function(df){
  set.seed(442)
  df_split <- initial_split(df, prop = 2/3)
  df_training <- training(df_split)
  df_validation <- testing(df_split)
  list(df_training,df_validation)
})

labelAll_training = data.frame(do.call(rbind,lapply(res,'[[',1)))
labelAll_validation = data.frame(do.call(rbind,lapply(res,'[[',2)))

fwrite(alterationsAll[match(labelAll_training$SAMPLE,alterationsAll$SAMPLE),],file=paste0("./data/",folderSomaticAlt,"/TCGA_alterations_training.tsv"),sep="\t")
fwrite(labelAll_training,file=paste0("./",folderSomaticAlt,"/TCGA_labels_training.tsv"),sep="\t",col.names = F)
fwrite(alterationsAll[match(labelAll_validation$SAMPLE,alterationsAll$SAMPLE),],file=paste0("./",folderSomaticAlt,"/TCGA_alterations_validation.tsv"),sep="\t")
fwrite(labelAll_validation,file=paste0("./",folderSomaticAlt,"/TCGA_labels_validation.tsv"),sep="\t",col.names = F)

# Make a single fold
# Separately for each tumor
tcga_sub = fread("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/TCGA/Cancer_subtypes_list.tsv")

# GBM recoded
gbm_lgg = data.frame(read_xlsx("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/TCGA/subtypes/GBM_LGG.xlsx",skip = 1))
gbm = gbm_lgg[gbm_lgg$Study=="Glioblastoma multiforme",]

tcga_sub[tcga_sub$`TCGA PanCanAtlas Cancer Type Acronym`=='GBM','Subtype'] <- gbm[match(tcga_sub[tcga_sub$`TCGA PanCanAtlas Cancer Type Acronym`=='GBM',COMMON],gbm$Case),"Original.Subtype"]

# PRAD
prad_def = fread('/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/TCGA/PRAD/TCGA_PRAD_subtypes.tsv')
tcga_sub$Subtype[tcga_sub$`Oncotree Code`=="PRAD"] <- "Other"
tcga_sub$Subtype[tcga_sub$COMMON%in%prad_def[prad_def$Subtype=="1-ERG",PATIENT_ID]] <- "ERG"
tcga_sub$Subtype[tcga_sub$COMMON%in%prad_def[prad_def$Subtype=="5-SPOP",PATIENT_ID]] <- "SPOP"


# tumorType = 'BRCA'
for(tumorType in unique(labelAll$DATASET)) {
  cat(tumorType,"\n")
  alterationsT = alterationsAll[labelAll$DATASET==tumorType,]
  labelT = labelAll[labelAll$DATASET==tumorType,]
  labelT$DATASET = tcga_sub[match(labelT$SAMPLE,tcga_sub$COMMON),Subtype]
  alterationsT = alterationsT[!is.na(labelT$DATASET),]
  labelT = labelT[!is.na(labelT$DATASET),]
  
  if(nrow(alterationsT)>0){
    dir.create(paste0("./",folderSomaticAlt,"/",tumorType))
    fwrite(alterationsT,file=paste0("./",folderSomaticAlt,"/",tumorType,"/",tumorType,"_alterations_all.tsv"),sep="\t")
    fwrite(labelT,file=paste0("./",folderSomaticAlt,"/",tumorType,"/",tumorType,"_labels_all.tsv"),sep="\t",col.names = F)
    
    idx = sample(1:nrow(labelT))
    
    res = lapply(split(labelT,labelT$DATASET),function(df){
      if(nrow(df)>9){
        set.seed(442)
        df_split <- initial_split(df, prop = 2/3)
        df_training <- training(df_split)
        df_validation <- testing(df_split)
        list(df_training,df_validation)
      }
    })
    
    labelT_training = data.frame(do.call(rbind,lapply(res,'[[',1)))
    labelT_validation = data.frame(do.call(rbind,lapply(res,'[[',2)))
    
    
    fwrite(alterationsT[match(labelT_training$SAMPLE,alterationsT$SAMPLE),],file=paste0("./",folderSomaticAlt,"/",tumorType,"/",tumorType,"_alterations_training.tsv"),sep="\t")
    fwrite(labelT_training,file=paste0("./",folderSomaticAlt,"/",tumorType,"/",tumorType,"_labels_training.tsv"),sep="\t",col.names = F)
    fwrite(alterationsT[match(labelT_validation$SAMPLE,alterationsT$SAMPLE),],file=paste0("./",folderSomaticAlt,"/",tumorType,"/",tumorType,"_alterations_validation.tsv"),sep="\t")
    fwrite(labelT_validation,file=paste0("./",folderSomaticAlt,"/",tumorType,"/",tumorType,"_labels_validation.tsv"),sep="\t",col.names = F)
    
    # Prepare the 3-fold cross to optimize the hyperparameters
    dir.create(paste0("./",folderSomaticAlt,"/",tumorType,"/foldOpt"))
    alterationsT_training = alterationsT[match(labelT_training$SAMPLE,alterationsT$SAMPLE),]
    
    
    num_groups = 3
    res = lapply(split(labelT_training,labelT_training$DATASET),function(df){
      set.seed(442)
      df_split <- initial_split(df, prop = 1/3)
      df_fold1 <- training(df_split)
      df_fold23 <- testing(df_split)
      set.seed(442)
      df_split <- initial_split(df_fold23, prop = 1/2)
      df_fold2 <- training(df_split)
      df_fold3 <- testing(df_split)
      list(df_fold1,df_fold2,df_fold3)
    })
    
    labelT_training_fold = list(data.frame(do.call(rbind,lapply(res,'[[',1))),
                                data.frame(do.call(rbind,lapply(res,'[[',2))),
                                data.frame(do.call(rbind,lapply(res,'[[',3))))
    
    for(i in 1:num_groups) {
      fwrite(alterationsT_training[match(data.frame(do.call(rbind,labelT_training_fold[-i]))$SAMPLE,alterationsT_training$SAMPLE),],file=paste0("./",folderSomaticAlt,"/",tumorType,"/foldOpt/",tumorType,"_alterations_fold",i,"_training.tsv"),sep="\t")
      fwrite(data.frame(do.call(rbind,labelT_training_fold[-i])),file=paste0("./",folderSomaticAlt,"/",tumorType,"/foldOpt/",tumorType,"_labels_fold",i,"_training.tsv"),sep="\t",col.names = F)
      fwrite(alterationsT_training[match(data.frame(do.call(rbind,labelT_training_fold[i]))$SAMPLE,alterationsT_training$SAMPLE),],file=paste0("./",folderSomaticAlt,"/",tumorType,"/foldOpt/",tumorType,"_alterations_fold",i,"_validation.tsv"),sep="\t")
      fwrite(data.frame(do.call(rbind,labelT_training_fold[i])),file=paste0("./",folderSomaticAlt,"/",tumorType,"/foldOpt/",tumorType,"_labels_fold",i,"_validation.tsv"),sep="\t",col.names = F)
    }
  }
}