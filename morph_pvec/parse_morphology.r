## ----------------------------------------------------------------------------
## Preamble
library(dplyr)
rm(list=ls())
graphics.off()


## ----------------------------------------------------------------------------
create.new.morpho = FALSE
create.new.pvec = FALSE
do.PCA = TRUE
plot.svd = TRUE
pca.choice = 1
# 1. only morph
# 2. only pvec
# 3. combine both


## ----------------------------------------------------------------------------
if(create.new.morpho)
{

  morp <- read.csv('taiwan_female_nmo/morpho_1588104037847.csv')
  morp = unname(unlist(morp))

  num.morph.features = 23

  ##--assign column names
  df <- data.frame()
  names = as.character(morp[1:num.morph.features])
  for (k in names) df[[k]] <- as.character()



  for(j in 1:num.morph.features)
  {
    cat("\nProcessing column",j,"..")
    for(i in 1:(length(morp)/num.morph.features-1))
    {
          df[i,j] = as.character(morp[num.morph.features*i+j])
    }
  }

  for(j in 1:num.morph.features)
  {
    if(j!=2 && j!=5)
      df[[j]] = as.numeric(df[[j]])
  }
  saveRDS(df, file="morpho23feat.rds")
  morpho23 = df
} else {
  morpho23 <- readRDS("morpho23feat.rds")
}





## ----------------------------------------------------------------------------
##--Align Persistent Vectors
if(create.new.pvec)
{
  p = read.csv("taiwan_female_nmo/pvec_1588104055264.csv")

  num.neurons = length(p[[1]])
  #intialize single empty class vector
  class_vec = rep(0,num.neurons)
  #loop through all neurons
  m_index = NULL
  for(n in 1:num.neurons)
  {
    m_index = which(p$neuron_id==morpho23$neuron_id[n])
    class_vec[n] = m_index
  }

  ##--treating it as a matrix
  # pvec = matrix(data = NA, nrow = num.neurons, ncol =100)
  # colnames(pvec) = as.character(names(p)[5:104])
  #
  # for(j in 5:104)
  #     pvec[,(j-4)] =  p[[j]][class_vec]

  pvec104 = p
  for(j in 1:length(p))
      pvec104[[j]] =  p[[j]][class_vec]

  saveRDS(pvec104, file="pvec104feat.rds")
} else{
  pvec104 = readRDS("pvec104feat.rds")
}


if(do.PCA)
{
  #load elbow finder
  source('/home/noob/WDHDD/clustering/exp/find_elbow.r')
  #clean non-numeric and redundant vectors
  morpho = morpho23[-c(1,2,5,6,12)]
  pvec = pvec104[-c(1:4)]

  if(pca.choice==1)
  {
    pca=prcomp(morpho,retx=TRUE,center=TRUE,scale=TRUE)
    d = find_elbow(pca$sdev)$PL[2]
    saveRDS(pca$x[,1:d],file="morpho.rds")
  } else if(pca.choice==2)
  {
    pca=prcomp(pvec,retx=TRUE,center=TRUE,scale=TRUE)
    d = find_elbow(pca$sdev)$PL[2]
    saveRDS(pca$x[,1:d],file="pvec.rds")
  } else if(pca.choice==3)
  {
    pca=prcomp(cbind(morpho,pvec),retx=TRUE,center=TRUE,scale=TRUE)
    d = find_elbow(pca$sdev)$IP
    saveRDS(pca$x[,1:d],file="pvec_n_morpho.rds")
  } else cat("ERROR! Invalid choice")

  var.explained = pca$sdev^2 /sum(pca$sdev^2)


}


if(plot.svd)
{
  ##
  x11()
  plot(cumsum(var.explained),
        ylab="Cumulative Percentage Variance",
        xlab="Index")
  points(d,cumsum(var.explained)[d],col="red",pch=19);
  abline(v=d, col = "darkgray", lty = 3, lwd=2)
  abline(h=cumsum(var.explained)[d], col = "darkgray", lty = 3, lwd=2)
  if(pca.choice==3 || pca.choice==2 ){
    xtick = c(0,d,seq(20,length(pca$sdev),by=20))
    ytick = c(0,0.25,0.50,round(cumsum(var.explained)[d],2),1)
    axis(side = 1, at = xtick)
    axis(side = 2, at = ytick)
    grid()
  }

  ##
  x11()
  plot(pca$sdev,
        ylab="Magnitude of singular value",
        xlab="Index")
  points(d,pca$sdev[d],col="red",pch=19);
  abline(v=d, col = "darkgray", lty = 3, lwd=2)
  if(pca.choice==3 || pca.choice==2){
    xtick = c(0,d,seq(20,length(pca$sdev),by=20))
    axis(side = 1, at = xtick); grid()
  }
}
