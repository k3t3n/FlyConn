## ----------------------------------------------------------------------------
## Preamble
library(factoextra)
library(corrplot)
library(dplyr)
rm(list=ls())
# graphics.off()

draw.plot=TRUE


## ----------------------------------------------------------------------------
morpho23 <- readRDS("morpho23feat.rds")
pvec104 = readRDS("pvec104feat.rds")

#load elbow finder
source('/home/noob/WDHDD/clustering/exp/find_elbow.r')
#clean non-numeric and redundant vectors
morpho = morpho23[-c(1,2,5,6,12)]
pvec = pvec104[-c(1:4)]

pca=prcomp(cbind(morpho,pvec),retx=TRUE,center=TRUE,scale=TRUE)
find_elbow(pca$sdev,verbose=FALSE)$IP
pca=prcomp(cbind(morpho,pvec),rank=d,retx=TRUE,center=TRUE,scale=TRUE)
var.explained = pca$sdev^2 /sum(pca$sdev^2)
var <- get_pca_var(pca)
ind <- get_pca_ind(pca)

ACl = readRDS("/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/data/Aclass_d11_t33_min100.rds")
morpho_classification = ACl$morpho

class_vec= readRDS('/media/WDHDD/clustering/morphology/class_vec_morph.rds')

roi = 4

roi_indicies = class_vec[which(morpho_classification==roi)]
roi_names = morpho23$Neuron_name[roi_indicies]


if(draw.plot)
{
  for(i in c(1,2,3,4,5))
  {
    for(j in 3)
    {
      p = fviz_pca_biplot(pca,
                  col.var = "contrib",
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                  axes=c(i,j),
                  select.var=list(name = c('coeff60','coeff61','coeff62','coeff63','Surface','Volume','Length')),
                  # select.var=list(name = c('coeff69','coeff70','coeff71','Depth','Surface','Volume','Length')),
                  select.ind=list(name = as.character(roi_indicies)),
                  # col.ind = grp,
                  label="var")
      x11()
      plot(p)
    }
  }
}
