## ----------------------------------------------------------------------------
##--Preamble
library(tictoc)
library(foreach)
library(iterators)
library(doParallel)
# library(tictoc)
# library(foreach)
# library(doParallel)
rm(list=ls())


getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}



## ----------------------------------------------------------------------------
## Create a data frame of all clusterings
create.df <- function(all_classification, df.filename, num.graphs)
{
        # Aconn = readRDS("/media/WDHDD/clustering/exp/data/A_conn.rds")
        # neuron_name = attributes(Aconn)$Dimnames[[1]]

        neuron_name = readRDS("../data/neuron_names.rds")
        df=data.frame("neuron_name"=neuron_name)

        for(gr.index in 1:num.graphs)
        {
            cat("\nWorking on graph",gr.index)
            temp=rep(NA,length(neuron_name))
            gr.name = paste("../data/binary_adj_matrix/A50/A50_",toString(gr.index),".rds",sep="")
            Ap = readRDS(gr.name)
            for (v.index in 1:length(neuron_name))
            {
                 p.index = which( attributes(Ap)$Dimnames[[1]] == neuron_name[v.index] )
                 if(length(p.index)!=0)  temp[v.index] = all_classification[[gr.index]][p.index]
            }
            df = cbind(df,temp)
            names(df)[(gr.index+1)] = paste("G",toString(gr.index),sep="")
        }

        saveRDS(df,file=df.filename)
        return(df)
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




# ## ----------------------------------------------------------------------------
# ## Reshape clustering label list
# clustering.labels <- function(pi, verbose=FALSE)
# {
#         num.of.vertices = length(pi[[1]])
#         clust_lab = list()
#         if(verbose) cat('\nInitializing list..')
#         for (node_num in c(1:num.of.vertices))
#         {
#             clust_lab[[node_num]] = rep(0,length(pi))
#         }
#         if(verbose) cat('\nList Done!')
#         for ( d in c(1:length(pi)) )
#         {
#             if(verbose) cat('\nWorking on run #',d)
#             m = pi[[d]]
#             for (node_num in c(1:length(m)))
#             {
#                clust_lab[[node_num]][d] = m[node_num]
#             }
#        }
#        return(clust_lab)
# } #end function




## ----------------------------------------------------------------------------
## Co-association matrix
coass <- function(clust_lab)
{
        n = length(clust_lab)
        coass_matrix = array(rep(0,1),dim=c(n,n))


        for (i in 1:n)
        {
            for (j in 1:n)
            {
                coass_matrix[i,j] = sum( (clust_lab[[i]]==clust_lab[[j]])*1 )
            } #j

        } #i

        return(coass_matrix)
} #end function




## ----------------------------------------------------------------------------
##
put.togther <- function(CM, label = rep(0,n), sum_threshold,
                        verbose = FALSE, seed.num=0)
{
  n = nrow(CM)
  # label = rep(0,n)
  available_index=1

  if(seed.num)
  {
      set.seed(seed.num)
      order.n = sample(n)
  } else {order.n=1:n}

  for(i in order.n)
  {
      if(label[i]==0)
      {
        label[i] = available_index
        available_index = available_index + 1
      }

      same_as_i = which(CM[i,] >= sum_threshold)
      all_zero = same_as_i[which(label[same_as_i]==0)]
      label[all_zero] = label[i]
      matched_labels = unique(label[same_as_i])
      assignment_label = min(matched_labels)
      if(assignment_label==0) stop('\n\nERROR!! Check assignment_label\n')
# browser()
      for(j in matched_labels)
        label[which(label==j)] = assignment_label

  } #i

  if(verbose)
      print(sort(table(label),decreasing=TRUE))


  return(label)
} #end function







###############################################################################
###############################################################################
## user input
## ----------------------------------------------------------------------------
remove_miss_vertices  = FALSE
trial.range           = c(1,5) #which trials/graphs to choose
num.graphs            = 5
num.cores             = 5
create_dataframe      = TRUE
file.descrip          = "MBHAC_results"
input.dir             = "./output/"
tau = 0.95

## auto generating file names
df.filename = paste(input.dir,"df_",file.descrip,"_g",num.graphs,".rds",sep="")
IVC.filename = paste(input.dir,"IVC_",file.descrip,"_g",trial.range[1],"_g",trial.range[2],".rData",sep="")


load(paste(input.dir,file.descrip,".rData",sep=""))
if(create_dataframe)
{
    cat("\nCreating new dataframe")
    cat("\n------------------------\n")
    create.df(all_classification, df.filename, num.graphs)
}

browser()
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


save(pi, pi_new, file=IVC.filename)

cat("\n\n**********************************************************************")
cat("\n** Program terminated successfully ** "); toc()
cat("**********************************************************************\n")

##

browser()
## auto generating filenames
IVC.filename = paste("./output/IVC_",file.descrip,".rData",sep="")
CM.filename = paste("./output/CM_",file.descrip,".rds",sep="")
output.filename = paste("./output/cl_",file.descrip,".rds",sep="")

# load data
load(IVC.filename); rm(pi)
num.trials = length(pi_new)



final_classification_list <- vector(mode='list',length=num.trials)

num.vertices = length(pi_new[[1]]) ; cat('\n\nNumber of vertices:',num.vertices)



        cat('\n\nNumber of graphs/trials combined:',num.trials)

        clust_lab = clustering.labels(pi_new)

        if(file.exists(CM.filename))
        {
          CM = readRDS(CM.filename)
          cat('\n\nLoaded existing coassociation matrix!')
        } else{
          cat('\n\nCreating NEW coassociation matrix..')
          CM = coass(clust_lab)
          saveRDS(CM, file=CM.filename)
        }
        diag(CM) = num.trials

        cat('\n\nTau set at:',tau,'\n')

        for( sum_threshold in c( (num.trials): round(tau*num.trials) ) )
        {
            final_classification_list[[sum_threshold]] = put.togther(CM=CM,sum_threshold=sum_threshold)
        }

saveRDS(final_classification_list,file=output.filename)

cat("\n\n******************************************")
cat("\n** Done! Program terminated successfully **\n\n")
