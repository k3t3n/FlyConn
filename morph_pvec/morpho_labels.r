## ----------------------------------------------------------------------------
## Preamble
library(factoextra)
library(dplyr)
rm(list=ls())
# graphics.off()


## ----------------------------------------------------------------------------
calc.morph.labels = TRUE
do.pvec.centerofmass = TRUE
create_new_csv = FALSE


morpho23 <- readRDS("morpho23feat.rds")
pvec104 = readRDS("pvec104feat.rds")

#load elbow finder
# source('/home/noob/WDHDD/clustering/exp/find_elbow.r')
#clean non-numeric and redundant vectors
morpho = morpho23[-c(1,2,5,6,12,19)]
pvec = pvec104[-c(1:4)]
MP = cbind(morpho,pvec) ##--main morpho_n_pvec

ACl = readRDS("../Aclass_d11_t33_min100.rds")
morpho_classification = ACl$morpho
ntwk = ACl$network
# ntwk  = readRDS("/home/noob/WDHDD/clustering/exp/ARGO/P05/MBHAC/ensemble_clustering/ntwk_d11.rds")
class_vec= readRDS('class_vec_morph.rds')


zero_vert = which(class_vec==0)
ntwk = ntwk[-zero_vert]
class_vec = class_vec[-zero_vert]
morpho_classification = morpho_classification[-zero_vert]

zero_vert = which(morpho_classification==0)
ntwk = ntwk[-zero_vert]
class_vec = class_vec[-zero_vert]
morpho_classification = morpho_classification[-zero_vert]

## ----------------------------------------------------------------------------
## Compare ratio of cluster-average to variable(feature)-average
if(calc.morph.labels)
{
  var_avg_all_neurons = colMeans(morpho)
  ratio = matrix(0, nrow=max(ntwk),ncol=ncol(morpho))
  var_avg = var_avg_all_neurons

  for(roi in 1:nrow(ratio))
  {
    roi_indicies = class_vec[which(ntwk==roi)]
    roi_names = morpho23$Neuron_name[roi_indicies]

    var_avg_per_cluster = colMeans(morpho[roi_indicies,])
    var_avg = rbind(var_avg, var_avg_per_cluster)

    ratio[roi,] = var_avg_per_cluster/var_avg_all_neurons
  }

  mo_name = names(morpho)

  mo_labels = vector(mode='list',length=nrow(ratio))

  col_maxs = apply(ratio, 2, max)
  col_mins = apply(ratio, 2, min)

  for(i in 1:nrow(ratio))
  {
    temp = NULL
    for(j in 1:length(mo_name))
    {
      if( ratio[i,j]>=1.8 )
      {
        temp = c(temp,paste("(H)",mo_name[j],sep=""))
        temp = gsub(" ", "", temp, fixed = TRUE)
      }
      if( ratio[i,j]<=0.55 )
      {
        temp = c(temp,paste("(L)",mo_name[j],sep=""))
        temp = gsub(" ", "", temp, fixed = TRUE)
      }
    }
    if(is.null(temp)){
      mo_labels[[i]] = paste(" ")
    } else {
      mo_labels[[i]] = paste(temp,collapse=", ")
    }
  }

  ntwk_name = rep(NA,nrow(ratio))
  for(i in 1:nrow(ratio)) ntwk_name[i] = paste('Class #',toString(i),sep='')
  label_frame = data.frame(ntwk_name)
  label_frame$morpho_label = unlist(mo_labels)

  colnames(var_avg) = mo_name
  rownames(var_avg) = c("All Neurons",ntwk_name)

  if(create_new_csv)
  {
    write.csv(label_frame,'morpho_labels.csv')
    write.csv(var_avg,'var_avg.csv')
  }

}



## ----------------------------------------------------------------------------
## Cumulative pvec, center of mass
calc.com <- function(x)
{
  temp = cumsum(as.numeric(x)/sum(x))
  ans = min(which(temp>=0.5))
}
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

if(do.pvec.centerofmass)
{
  com = apply(pvec,1,calc.com)  # center of mass
  median_com = rep(0, max(ntwk))
  mode_com = rep(0, max(ntwk))
  pvec_label = rep(NULL, max(ntwk))
  for(roi in 1:length(median_com))
  {
    roi_indicies = class_vec[which(ntwk==roi)]
    roi_names = morpho23$Neuron_name[roi_indicies]
    # if(roi==35) browser()
    median_com[roi] = median(com[roi_indicies])
    mode_com[roi] = getmode(com[roi_indicies])

    if (mode_com[roi]<=40) {
      pvec_label[roi] = 'proximal'
    } else if (mode_com[roi]>=60) {
      pvec_label[roi] = 'most-distant'
    } else if ((mode_com[roi]>=48) && (mode_com[roi]<=52) ) {
      pvec_label[roi] = 'well-balanced'
    }
  }

  ntwk_name = rep(NA,length(median_com))
  for(roi in 1:length(median_com))
    {
      ntwk_name[roi] = paste('Class #',toString(roi),sep='')
      if(is.na(pvec_label[roi])) pvec_label[roi] = paste(" ")
    }
  label_frame = data.frame(ntwk_name)
  label_frame$pvec_label = unlist(pvec_label)

  names(mode_com) = ntwk_name

  if(create_new_csv)
  {
    write.csv(label_frame,'pvec_labels.csv')
    write.csv(data.frame(mode_com),'pvec_mode_com.csv')
  }
}
