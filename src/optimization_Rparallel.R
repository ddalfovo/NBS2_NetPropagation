library(parallel)
setwd("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/Tools/NBS2_NetPropagation")

alpha = c(0.1,0.3,0.5,0.7,0.9)
delta = c(0.01,0.032,0.1,0.32,1)
beta = c(2e-2,2e-3,2e-4,2e-5,2e-6)

rst_prob_fix = 0.5 # alpha
lam_fix = 1e-1 # delta
beta_loss_fix = 2e-4 # beta

comb = data.frame(cbind(alpha,delta,beta))
comb$delta = lam_fix
comb$beta = beta_loss_fix
alpha = lapply(1:nrow(comb),function(i)unname(unlist(comb[i,])))

comb = data.frame(cbind(alpha,delta,beta))
comb$alpha = rst_prob_fix
comb$beta = beta_loss_fix
delta = lapply(1:nrow(comb),function(i)unname(unlist(comb[i,])))

comb = data.frame(cbind(alpha,delta,beta))
comb$alpha = rst_prob_fix
comb$delta = lam_fix
beta = lapply(1:nrow(comb),function(i)unname(unlist(comb[i,])))

listOpt = c(alpha,delta,beta)
listOpt = unique(listOpt)

src = paste0(getwd(),"/src/")
out = paste0(getwd(),"/results/")

tumor = 'BRCA'
mclapply(listOpt,function(pSet){
  mclapply(c("1","2","3"),function(fold){
    system(paste0("python ",folder,"NBS_tumorSpec_foldOpt_parallel_Rparallel.py ",paste0(pSet,collapse=" ")," ",fold," ",out," ",tumor))
  },mc.cores = 1)
},mc.cores = 5)
