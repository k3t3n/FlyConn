## ----------------------------------------------------------------------------
##--Preamble
library(mclust)
rm(list=ls())

spad = readRDS('/media/WDHDD/ckt/nat/data/nat_neuropil_distribution.rds')
pca=prcomp(spad,retx=TRUE,center=TRUE,scale=TRUE)

data = pca$x

mclust.options('subset'=nrow(data))

mc = Mclust(data, modelNames=c('VVI','VVV'), G=54)

saveRDS(mc, file='nat_spad_clustering_tau95_k54.rds')
