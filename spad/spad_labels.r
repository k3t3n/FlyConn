## ----------------------------------------------------------------------------
##--Preamble
library(mclust)
rm(list=ls())

ACl = readRDS("../../Aclass_d11_t33_min100.rds")
ntwk = ACl$network
tresh = 0.95

spad = readRDS('/media/WDHDD/ckt/nat/data/nat_neuropil_distribution.rds')

spad_names = colnames(spad)
no.ntwk.clusters = max(ntwk)

spad_avg = colMeans(spad)
spad_prop = spad_avg/sum(spad_avg)
spad_ntwk_avgs = matrix(0,nrow=no.ntwk.clusters,ncol=ncol(spad))
spad_ntwk_prop = matrix(0,nrow=no.ntwk.clusters,ncol=ncol(spad))


for(i in 1:no.ntwk.clusters )
{
  ntwk_indicies = which(ntwk==i)
  spad_ntwk_avgs[i,] = colMeans(spad[ntwk_indicies,])
  spad_ntwk_prop[i,] = spad_ntwk_avgs[i,]/sum(spad_ntwk_avgs[i,])
}

spad_ntwk_norm = matrix(0,no.ntwk.clusters,length(spad_names))
for(j in 1:ncol(spad_ntwk_norm))
{
  spad_ntwk_norm[,j] = spad_ntwk_avgs[,j]/spad_avg[j]
}
  spad_ntwk_norm[which(is.infinite(spad_ntwk_norm))] = 0


max_avgs = apply(spad_ntwk_avgs, 2, max)
max_prop = apply(spad_ntwk_prop, 2, max)
max_norm_avgs = apply(spad_ntwk_norm, 2, max)


spad_labels_avgs = list()
spad_labels_prop = list()
spad_labels_norm = list()

for(i in 1:no.ntwk.clusters )
{
  temp  = NULL
  temp1 = NULL
  temp2 = NULL
  temp3 = NULL
  for(j in 1:length(spad_names))
  {

    if(spad_ntwk_avgs[i,j] >= (tresh*max_avgs[j]) && spad_ntwk_prop[i,j] < (tresh*max_prop[j])){
    # if(spad_ntwk_avgs[i,j] >= (15*spad_avg[j])){
            temp = c(temp,paste("(H)",spad_names[j],sep=""))
            temp = gsub(" ", "", temp, fixed = TRUE)
            # temp1 = c(temp1,spad_names[j])
    }
    if(spad_ntwk_avgs[i,j] < (tresh*max_avgs[j]) && spad_ntwk_prop[i,j] >= (tresh*max_prop[j])){
    # if(spad_ntwk_avgs[i,j] >= (15*spad_avg[j])){
            temp = c(temp,paste("(P)",spad_names[j],sep=""))
            temp = gsub(" ", "", temp, fixed = TRUE)
    }
    if(spad_ntwk_avgs[i,j] >= (tresh*max_avgs[j]) && spad_ntwk_prop[i,j] >= (tresh*max_prop[j])){
            temp = c(temp,paste("(H|P)",spad_names[j],sep=""))
            temp = gsub(" ", "", temp, fixed = TRUE)
    }
    # if(spad_ntwk_prop[i,j] >= (tresh*max_prop[j])){
    # # if(spad_ntwk_prop[i,j] >= (5*spad_prop[j])){
    #         temp2 = c(temp2,spad_names[j])
    # }
    if(spad_ntwk_norm[i,j] >= (tresh*max_norm_avgs[j])){
            temp3 = c(temp3,spad_names[j])
    }
  }#j
  if(is.null(temp)){
    spad_labels_avgs[i] = paste(" ")
  }else{
    spad_labels_avgs[[i]] = paste(temp,collapse=", ")
  }

  # if(is.null(spad_labels_prop[i])){
  #   spad_labels_prop[i] = paste(" ")
  # }else{
  #   spad_labels_prop[[i]] = paste(temp2,collapse=", ")
  # }

  if(is.null(temp3)){
    spad_labels_norm[i] = paste(" ")
  }else{
    spad_labels_norm[[i]] = paste(temp3,collapse=", ")
  }


}#i

ntwk_name = rep(NA,no.ntwk.clusters)
for(i in 1:no.ntwk.clusters)  ntwk_name[i] = paste('Ntwk',toString(i),sep='')

label_frame = data.frame(ntwk_name)
label_frame$spad_labels_avgs = unlist(spad_labels_avgs)
# label_frame$spad_labels_prop = unlist(spad_labels_prop)
# label_frame$spad_labels_norm = unlist(spad_labels_norm)

colnames(spad_ntwk_avgs) = spad_names; rownames(spad_ntwk_avgs) = label_frame$ntwk_name
colnames(spad_ntwk_prop) = spad_names; rownames(spad_ntwk_prop) = label_frame$ntwk_name
write.csv(spad_ntwk_avgs, "spad_ntwk_avgs.csv")
write.csv(spad_ntwk_prop, "spad_ntwk_prop.csv")
write.csv(label_frame,'spad_labels.csv')
