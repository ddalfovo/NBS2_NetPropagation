library(data.table)
setwd("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/Tools/NBS2_NetPropagation")

tumor = "BRCA"
list.pickle = list.files(paste0("./results/",tumor,"/pickles/"), pattern=paste0(tumor,"_pickle_"), full.names = T)

res = lapply(list.pickle,function(f){
  system(paste0("python src/optimization_GetAccuracy.py ",f, " ",tumor),intern = T)
})

# save.image(paste0("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/NBS/results/",tumor,"_pickles_results.RData"))
# load(paste0("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/NBS/results/",tumor,"_pickles_results.RData"))

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

r = lapply(res,function(el){
  elems = strsplit(gsub("^\\[","",gsub("\\]$","",el)),", ")
  iter = length(elems[[1]])
  acc = as.numeric(elems[[1]][iter])
  c(iter,acc)
})
names = strsplit(gsub(".pickle$","",basename(list.pickle)),"_")
names(names) = sapply(names,"[[",3)
names = lapply(names,function(n){
  n=n[4:6]
  as.numeric(n)
})
names.df = data.frame(do.call(rbind,names))
names.df$fold = names(names)
colnames(names.df) = c("alpha","delta","beta","fold")


optAlpha = lapply(alpha,function(el){
  idx = which(names.df$alpha == el[1]&names.df$delta == el[2]&names.df$beta == el[3])
  c(mean(sapply(r, "[[", 2)[idx]),length(sapply(r, "[[", 2)[idx]))
})

optDelta = lapply(delta,function(el){
  idx = which(names.df$alpha == el[1]&names.df$delta == el[2]&names.df$beta == el[3])
  c(mean(sapply(r, "[[", 2)[idx]),length(sapply(r, "[[", 2)[idx]))
})

optBeta = lapply(beta,function(el){
  idx = which(names.df$alpha == el[1]&names.df$delta == el[2]&names.df$beta == el[3])
  c(mean(sapply(r, "[[", 2)[idx]),length(sapply(r, "[[", 2)[idx]))
})

library(ggplot2)
pdf(paste0("./results/",tumor,"/results_fold.pdf"))
dfAlpha = data.frame(alpha = sapply(alpha,'[[',1),accAlpha = sapply(optAlpha,"[[",1))
ggplot(dfAlpha,aes(x=alpha,y=accAlpha)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks=dfAlpha$alpha) +
  theme_minimal() +
  ylab("accuracy") + xlab("alpha")


dfDelta = data.frame(delta = sapply(delta,'[[',2),accDelta = sapply(optDelta,"[[",1))
ggplot(dfDelta,aes(x=log10(delta),y=accDelta)) +
  geom_line() + geom_point() +
  scale_x_continuous(labels = dfDelta$delta, breaks = log10(dfDelta$delta)) +
  theme_minimal() +
  ylab("accuracy") + xlab("delta")

dfBeta = data.frame(beta = sapply(beta,'[[',3),accBeta = sapply(optBeta,"[[",1))
ggplot(dfBeta,aes(x=log10(beta),y=accBeta)) +
  geom_line() + geom_point() +
  scale_x_continuous(labels = dfBeta$beta,breaks = log10(dfBeta$beta)) +
  theme_minimal() +
  ylab("accuracy") + xlab("beta")
dev.off()
