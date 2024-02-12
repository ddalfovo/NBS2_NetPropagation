library(cgdsr)
# library(cBioPortalData)
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

load(file="data/PCNet_igraph.RData")
PathwayCommons_PCNet = fread("./data/commonGenesPCnet.tsv")
gene_mapping = fread("./data/Homo_sapiens.gene_info")

all_syn = strsplit(gene_mapping$Synonyms,"\\|")
gene_mapping_coding = gene_mapping[gene_mapping$type_of_gene=="protein-coding",]

convNameGenes = function(nodeName){
  res = mclapply(nodeName,function(name){
    idx = which(gene_mapping$Symbol==name)
    if(length(idx)==0){
      idx_syn = grep(name,gene_mapping$Synonyms)
      if(length(idx_syn)>0) {
        idx_res = which(sapply(strsplit(gene_mapping$Synonyms[idx_syn],"\\|"),function(search){
          any(search==name)
        }))
        if(length(idx_res)==1){
          gene_mapping[idx_syn[idx_res],]$Symbol
        } else {
          NA
        }
      } else {
        NA
      }
    } else if(length(idx)>1) {
      NA
    } else {
      name
    }
  },mc.cores = 50)
  return(res)
}
res = convNameGenes(V(ig)$nodeName)
nameConverted_PCnet = unlist(res)

PCNet_filtered = ig
PCNet_filtered = delete_vertices(PCNet_filtered,V(PCNet_filtered)$name[which(is.na(nameConverted_PCnet))])

save(PCNet_filtered, file="data/networks_filtered.RData")

###########################################################
### Processing the hallmark genes and use just this one ###
###########################################################
# hallmarksFile = readLines("./data/hallmarks.txt")
# idx = c(0,which(hallmarksFile==""),length(hallmarksFile)+1)
# pathways = lapply(1:(length(idx)-1),function(i){
#   hallmarksFile[(idx[i]+2):(idx[i+1]-1)]
# })
# names(pathways) = hallmarksFile[c(0,which(hallmarksFile==""))+1]
# hallmarks_genes = fread("data/hallmarks_genes.txt",header = F,col.names = "gene")

### ALL genes
hallmarks_genes = V(PCNet_filtered)$nodeName
######################
### GENERATE SOMATIC PROFILES given a data.frame of genes

somaticProfiles = function(hallmarks_genes, cores, chunks){
  mycgds = CGDS("http://www.cbioportal.org/")
  # cbio <- cBioPortal()
  
  TCGA_studies=getCancerStudies(mycgds)[,c(1,2)]
  TCGA_studies=TCGA_studies[grepl('TCGA, PanCancer Atlas', TCGA_studies$name),]
  
  num_groups = chunks
  hallmarks_genes_splitted = hallmarks_genes %>% 
    group_by((row_number()-1) %/% (n()/num_groups)) %>%
    nest %>% pull(data)
  
  chunk = pbmclapply(hallmarks_genes_splitted,function(hallmarks_genes_splitted_df){
    mylist_CN = mclapply(1:nrow(TCGA_studies), function(i){
      cat(i,'\n')
      mycancerstudy = TCGA_studies[i,1]
      # getGeneticProfiles(mycgds,mycancerstudy)
      mycaselist=paste0(mycancerstudy,'_all')
      mygeneticprofile = paste0(mycancerstudy,'_gistic')
      
      df=getProfileData(mycgds,hallmarks_genes_splitted_df$gene,mygeneticprofile,mycaselist)
      df$DATASET=toupper(strsplit(mycancerstudy, '_')[[1]][1])
      df$SAMPLE=gsub('\\.', '-', row.names(df))
      
      df 
    },mc.cores = 1)
    as.data.frame(rbindlist(mylist_CN))
  },mc.cores = cores)
  final_CN=as.data.frame(do.call(cbind,chunk))
  
  nc = ncol(final_CN)
  final_CN = final_CN[,c(nc-1,nc,1:(nc-2))]
  
  
  chunk = pbmclapply(hallmarks_genes_splitted,function(hallmarks_genes_splitted_df){
    mylist_MUT = mclapply(1:nrow(TCGA_studies), function(i){
      cat(i,'\n')
      mycancerstudy = TCGA_studies[i,1]
      getGeneticProfiles(mycgds,mycancerstudy)
      mycaselist=paste0(mycancerstudy,'_all')
      mygeneticprofile = paste0(mycancerstudy,'_mutations')
      
      df = as.matrix(getProfileData(mycgds,hallmarks_genes_splitted_df$gene,mygeneticprofile,mycaselist))
      tmp = df
      tmp[is.na(tmp)] <- 0
      tmp[tmp=='NaN'] <- 0
      tmp[tmp!='0'] <- 1
      tmp = as.data.frame(tmp)
      colnames(tmp) = colnames(df)
      rownames(tmp) = rownames(df)
      df = tmp
      df$DATASET=toupper(strsplit(mycancerstudy, '_')[[1]][1])
      df$SAMPLE=gsub('\\.', '-', row.names(df))
      
      df 
    }, mc.cores = 1)
    as.data.frame(rbindlist(mylist_MUT))
  },mc.cores = cores)
  final_MUT=as.data.frame(do.call(cbind,chunk))
  
  nc = ncol(final_MUT)
  final_MUT = final_MUT[,c(nc-1,nc,1:(nc-2))]
  return(list(final_MUT = final_MUT, final_CN = final_CN))
}

res_PCNet = somaticProfiles(PathwayCommons_PCNet, 10, 10)
fwrite(res_PCNet[['final_CN']],file = paste0('data/somaticAlterations/gisticCN.tsv'),sep='\t',col.names = T)
fwrite(res_PCNet[['final_MUT']],file = paste0('data/somaticAlterations/mutations.tsv'),sep='\t',col.names = T)

save.image("./data/somaticMutationsProfiles/cBioPortalAlterations.RData")
