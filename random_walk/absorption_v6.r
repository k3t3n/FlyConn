## ----------------------------------------------------------------------------
##--Preamble
library(igraph,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(coda,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(MCMCpack,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(mclust,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(Matrix,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(tictoc,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(foreach,lib.loc="/home/kmehta7/R/3.5.2/library/")
# library(Rmpi)
library(iterators,lib.loc="/home/kmehta7/R/3.5.2/library/")
# library(doMPI,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(doParallel,lib.loc="/home/kmehta7/R/3.5.2/library/")

rm(list=ls())

## ----------------------------------------------------------------------------
ell <- function(VP)
{
    EP = rep(VP, each=2)[-1]
    EP = EP[-length(EP)]
    trans.prob = E(gT)$weight[get.edge.ids(gT, EP)]
    cost       = E(gC)$weight[get.edge.ids(gC, EP)]
    sum(trans.prob*cost)
}




## ----------------------------------------------------------------------------
P   <- readRDS('/home/kmehta7/random_walk/absorption/data/PMAT_95.rds')
gps <- readRDS('/home/kmehta7/random_walk/absorption/data/cl_95.rds')
gps = gps[which(gps!=0)]
gps = unname(table(gps))

# P   <- rbind(c(0, 0.05,   0,    0,   0,     0,    0),
#              c(0, 0.10, 0.3,    0, 0.25,    0,    0),
#              c(0,    0,   0, 0.10,    0,    0,    0),
#              c(0,    0,   0,    0,    0,    0,    0),
#              c(0,    0,   0, 0.50,    0, 0.20,  0.4),
#              c(0,    0,   0, 0.15,    0,    0,    0),
#              c(0,    0,   0,    0,    0,    0,    0))
# gps <- c(1000, 250, 100, 500, 300, 400,150)

##--
set.seed(60)
num.cores = 28
path.num.max = 10^4
verbose = TRUE
do.sum = FALSE
binary.transition = TRUE
filename = '/scratch/kmehta7/absorption_v6_mean10e4_95.rData'


#intialize parameters
g        = graph_from_adjacency_matrix(P, mode = "directed", weighted = TRUE, diag = TRUE)
blk_size = gps/min(gps)
deg_out  = degree(g, v = V(g), mode = "out", loops = FALSE, normalized = FALSE)

#transition
T = P
T[which(T>0)]=1
diag(T) = 0
if(!binary.transition) for(ind in 1:nrow(T)) if(deg_out[ind]) T[ind,] = T[ind,]/deg_out[ind]
gT = graph_from_adjacency_matrix(T, mode = "directed", weighted = TRUE, diag = FALSE)

#cost
C = T
for(ind in which(t(C)>0))
{
  i = ceiling(ind/nrow(P))
  j = ind-(ncol(P)*(i-1))
  C[i,j] = 1/(P[i,j]*blk_size[i]*blk_size[j])
}
gC = graph_from_adjacency_matrix(C, mode = "directed", weighted = TRUE, diag = FALSE)

A = matrix(rep(NA),nrow(P),ncol(P)) #vector(mode="list",length=nrow(P))
D = matrix(rep(NA),nrow(P),ncol(P))
# A = NULL
# D = NULL

tic(" Total Time")
registerDoParallel(cores=num.cores)
n=0


for(source_node in 1:nrow(P))
{
  dfs_temp=dfs(g, root=source_node, neimode="out", unreachable = FALSE, order = TRUE, order.out = TRUE,
           father = FALSE, dist = FALSE, in.callback = NULL,
           out.callback = NULL, extra = NULL, rho = parent.frame())$order
  for(target_node in 1:ncol(P))
  {
    if(!(target_node%in%dfs_temp) || source_node==target_node)
    {
      A[source_node,target_node] = 0
      D[source_node,target_node] = 0
      tot=(source_node-1)*nrow(P) + target_node
      if(verbose){cat('\nDone with',tot,'..'); cat('NO PATHS!')}
    }
  }
}

A = as.vector(t(A))
D = as.vector(t(D))

tot.ind = which(is.na(A))


while (n <= num.cores*floor(length(tot.ind)/num.cores))#length(A))
{
    output_vec <- foreach(i=1:num.cores) %dopar%
    {
          #check tot
          if((n+i)>length(tot.ind))
          {
            pairwise_absorption = 0
            drift = 0
            out_list = list(pairwise_absorption=pairwise_absorption, drift=drift)
            return(out_list)
          }
          ##
          tot = tot.ind[n+i]
          source_node = ceiling(tot/nrow(P))
          target_node = tot-(ncol(P)*(source_node-1))
          path_list <- list()
          seednum = 400

          ##
            for(path.num in 1:path.num.max)
            {
              set.seed(seednum+path.num*100)
              #intialize
              path_vector = c(source_node)
              P_reduced=P
              g_reduced=g
              hit.target=FALSE



              while(hit.target == FALSE)
              {
                # print(path_vector)
                root=path_vector[length(path_vector)]
                dfs_temp=dfs(g_reduced, root=root, neimode="out", unreachable = FALSE, order = TRUE, order.out = TRUE,
                       father = FALSE, dist = FALSE, in.callback = NULL,
                       out.callback = NULL, extra = NULL, rho = parent.frame())$order

                #break condition
                if(!(target_node%in%dfs_temp))
                {
                  path_vector = NA
                  break
                }
                ##

                dfs_temp=dfs_temp[-which(is.na(dfs_temp))]
                P_reduced[-dfs_temp,]=0; P_reduced[,-dfs_temp]=0
                g_reduced = graph_from_adjacency_matrix(P_reduced, mode = "directed", weighted = TRUE, diag = FALSE)
                nz = which(P_reduced[root,]!=0)
                nz = setdiff(nz,path_vector)

                #break condition
                if(length(nz)==0)
                {
                  path_vector = NA
                  break
                }
                ##

                if (length(nz)==1){
                  nxt.vert = nz
                } else{
                  nxt.vert = sample(nz,1)
                }

                P_reduced[root,]=0; P_reduced[,root]=0
                path_vector = c(path_vector, nxt.vert)
                if (nxt.vert==target_node) hit.target=TRUE
              }
              path_list[[path.num]] = path_vector
            }
            no.of.hops=unlist(lapply(path_list,length))
            path_list = unique(path_list[which(no.of.hops>1)])

            # print(path_list)


          #--avg path length
          path_length_list = lapply(path_list, ell)
          if(do.sum){
            pairwise_absorption = sum(unlist(path_length_list))
          } else{
            pairwise_absorption = mean(unlist(path_length_list))
          }

          #--shortest path length
          shortest_path = min((unlist(path_length_list)))
          ##--driftiness
          drift = pairwise_absorption/shortest_path
          ##--clear garbage
          rm(path_list);rm(path_length_list);gc(verbose=FALSE)

          if(verbose){cat('\nDone with',tot,'..') }

          out_list = list(pairwise_absorption=pairwise_absorption, drift=drift)
          return(out_list)
    }#dopar

    for(j in 1:num.cores)
    {
        A[tot.ind[n+j]] = output_vec[[j]]$pairwise_absorption
        D[tot.ind[n+j]] = output_vec[[j]]$drift
    }
    save(A,D, file=filename)
    n = n+num.cores
}#while


A=matrix(A,nrow=nrow(P),ncol=ncol(P),byrow=TRUE)
D=matrix(D,nrow=nrow(P),ncol=ncol(P),byrow=TRUE)

A[is.nan(A)]=0
D[is.nan(D)]=0

  avg_out_absorption = rowSums(A)/(ncol(P)-1)
  avg_out_driftiness = rowSums(D)/(ncol(P)-1)
  cat('\n\nOut Absorption:\n',avg_out_absorption)
  cat('\n\nOut Driftiness:\n',avg_out_driftiness)

  avg_in_absorption = colSums(A)/(nrow(P)-1)
  avg_in_driftiness = colSums(D)/(nrow(P)-1)
  cat('\n\nIn Absorption:\n',avg_in_absorption)
  cat('\n\nIn Driftiness:\n',avg_in_driftiness)


cat('\n\n')

save(P, C, T, A, D,
     avg_out_absorption, avg_in_absorption,
     avg_out_driftiness, avg_in_driftiness,
     file=filename)

cat("\n\n**********************************************************************")
cat("\n** Program terminated successfully ** "); toc()
cat("**********************************************************************\n")
