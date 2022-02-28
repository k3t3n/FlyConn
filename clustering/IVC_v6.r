## ----------------------------------------------------------------------------
##--Preamble
library(tictoc,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(foreach,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(iterators,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(doParallel,lib.loc="/home/kmehta7/R/3.5.2/library/")
# library(tictoc)
# library(foreach)
# library(doParallel)
rm(list=ls())


getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}




## ----------------------------------------------------------------------------
##
calc.dist <- function(pi, pi_star)
{
  y               = clustering.labels(pi)
  num.of.vertices = length(y)
  num.trials      = length(pi)

  k_labels = unique(pi_star)
  if(length(which(k_labels==0))) k_labels = k_labels[-which(k_labels==0)]
  kappa = length(unique(pi_star))

  ##--Calculate partitions
  P <- vector(mode = "list", length = max(k_labels))
  for(k in k_labels)
  {
    P[[k]] = which(pi_star==k)
  }

  ##--Calculate center of each parition, for each clustering
  ##--Rows are
  ##--Columns are
  M <- matrix(rep(NA),num.trials,length(P))
  for(i in 1:num.trials)
  {
    for(k in k_labels)
    {
      M[i,k] = getmode(pi[[i]][P[[k]]])
    }
  }

  ##--create distance vector for each data point
  dist = matrix(NA,length(k_labels),num.of.vertices) #vector(mode = "list", length = length(k_labels))
  for(k in k_labels)
  {
    dist_vector = vector(mode = "list", length = num.of.vertices)
    for(j in 1:num.of.vertices)
    {
      dist_vector[[j]] = 1*(y[[j]] != M[,k])
    }
    dist[k,] = unlist(lapply(dist_vector,sum))
  }

  pi_new_ref_cluster = apply(dist,2,which.min)

  return(pi_new_ref_cluster)

} #end function




## ----------------------------------------------------------------------------
## Reshape clustering label list
clustering.labels <- function(pi, verbose=FALSE)
{
        num.of.vertices = length(pi[[1]])
        clust_lab = list()
        if(verbose) cat('\nInitializing list..')
        for (node_num in c(1:num.of.vertices))
        {
            clust_lab[[node_num]] = rep(0,length(pi))
        }
        if(verbose) cat('\nList Done!')
        for ( d in c(1:length(pi)) )
        {
            if(verbose) cat('\nWorking on run #',d)
            m = pi[[d]]
            for (node_num in c(1:length(m)))
            {
               clust_lab[[node_num]][d] = m[node_num]
            }
       }
       return(clust_lab)
} #end function





## ----------------------------------------------------------------------------
##
remove_miss_vertices  = FALSE
trial.range           = c(101,200) #which trials/graphs to choose
num.graphs            = 800
file.descrip          = "P05_MBHAC_4to100_d11"
input.dir             = "/home/kmehta7/exp/IVC/data/"


## auto generating file names
df.filename = paste(input.dir,"df_",file.descrip,"_g",num.graphs,".rds",sep="")
output.filename = paste("/scratch/kmehta7/IVC_",file.descrip,"_g",trial.range[1],"_g",trial.range[2],".rData",sep="")
num.cores = 40


##--The clusterings 'pi'
df = readRDS(df.filename)
num.trials = trial.range[2]-trial.range[1]+1 #length(df)-1
pi <- list()
for(i in 1:num.trials)
{
    pi[[i]] = df[,(i+trial.range[1])]
    if (!remove_miss_vertices)
        pi[[i]][which(is.na(pi[[i]]))]=0
}

# pi      = list()
# pi[[1]] = c(1,1,1,2,2,2,2,2,1)
# pi[[2]] = c(2,2,1,1,1,1,1,1,1)
# pi[[3]] = c(1,1,2,2,2,2,2,3,3)
# pi[[4]] = c(1,2,2,2,2,2,2,3,3)
# pi[[5]] = c(1,1,1,3,2,2,2,3,3)
# num.trials = length(pi)

## ----------------------------------------------------------------------------
##--Clustering assignment of each data point (across all clusterings) 'y_j'

pi_new = pi
# set.seed(3000)
shuffle.order = 1:num.trials #c(4, 2, 1, 5, 3)
cat("\nShuffle order is:\n",shuffle.order,"\n")
tic('Total Time')

registerDoParallel(cores=num.cores)
n=0

while (n < num.trials)
{
  output_vec <- foreach(j=1:num.cores) %dopar%
  {
    ref_cluster = shuffle.order[n+j]
    converged=FALSE
    it.cnt = 0
    while(!converged)
    {
      it.cnt = it.cnt+1; cat("\nIteration number",ref_cluster,".",it.cnt)
      ##--Set pi_star,reference clustering
      pi_star = pi_new[[ref_cluster]]
      pi_new[[ref_cluster]] = calc.dist(pi,pi_star)
      ##--check covergence
      if (all(pi_star == pi_new[[ref_cluster]])) converged=TRUE
    }#converge
    return(pi_new[[ref_cluster]])
  }#foreach

  for(i in 1:num.cores)
  {
    n = n+1
    pi_new[[n]] = output_vec[[i]]
  }
  # save(pi_new, file=output.filename)
}


save(pi, pi_new, file=output.filename)

cat("\n\n**********************************************************************")
cat("\n** Program terminated successfully ** "); toc()
cat("**********************************************************************\n")

##
