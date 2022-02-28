
## ----------------------------------------------------------------------------
## Packages

library(igraph,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(coda,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(MCMCpack,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(mclust,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(Matrix,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(RSpectra,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(fpc,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(clusterCrit,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(withr,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(crayon,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(ggplot2,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(e1071,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(caret,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(tictoc,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(foreach,lib.loc="/home/kmehta7/R/3.5.2/library/")
# library(Rmpi)
library(iterators,lib.loc="/home/kmehta7/R/3.5.2/library/")
# library(doMPI,lib.loc="/home/kmehta7/R/3.5.2/library/")
library(doParallel,lib.loc="/home/kmehta7/R/3.5.2/library/")

cat("\014")     # clear screen
rm(list=ls())   # clean workspace
# if (!is.null(dev.list())) dev.off()   # close open x11 windows
setwd("/home/kmehta7/exp/P05/MBHAC/") # set working directory







## ----------------------------------------------------------------------------
## Estimate P

#g is a matrix of connectivity,
#(either original or the clustered ones, depending on what you feed in as "groups")
#estP should be the estimated P-matrix based on how the groups are assigned
estP <- function(g,groups) # This will take in groups=newClassification
  #notice that max(groups) may not be the number of groups!
  #Maybe this is what causes an error? The number should be the maximum as well.
{
    #group.labels <- as.numeric(dimnames(table(groups))[[1]])
    d = max(as.numeric(dimnames(table(groups))[[1]]))
    R <- matrix(0,nrow=d,ncol=d) #makes a matrix of dimension equal to the number of groups found
    #print(paste("The table of groups is"))
    #print(table(groups))
    for(i in 1:d){
      a<-which(groups==i)
      for(j in 1:d){
        b<-which(groups==j)
        if((length(a)>0)&&(length(b)>0)){
            if(i == j){
                R[i,j]<-sum(g[a,b])/(length(a)*(length(b)-1))
            }
            else{
                R[i,j]<-sum(g[a,b])/(length(a)*length(b))
            }
      }
    }#j
  }#i
  return(R)
}






### Adjacency Binary Matrix
binary_adj_matrix <- function(A_conn, A_prob, sd)
{
    set.seed(sd)

    #create the binary Adjacency Matrix
    u   = runif(length(A_prob), min = 0, max = 1) #sample from unifrom distribution
    Adj_entries = (A_prob>u)*1
    Adj = A_conn
    Adj@x = Adj_entries


    # remove all-zero zero rows and columns
    if ( length(which(colSums(Adj)==0)) && length(which(rowSums(Adj)==0)) )
            {hasnz = TRUE} else  {hasnz = FALSE}
    while(hasnz==TRUE)
      {
        coldel = which(colSums(Adj)==0)
        Adj = Adj[-coldel,-coldel]
        rowdel = which(rowSums(Adj)==0)
        Adj = Adj[-rowdel,-rowdel]
        if ( length(which(colSums(Adj)==0)) && length(which(rowSums(Adj)==0)) )
            {hasnz = TRUE} else {hasnz = FALSE}
      }

      return(Adj)
}







## ----------------------------------------------------------------------------
## Spectral Embedding

graph.spectral.embedding <- function( A,
                                      d = 4,
                                      flip.edges = FALSE,
                                      augmented = TRUE,
                                      frac  = 0,
                                      cvec  = degree(graph)/(vcount(graph) - 1) )
{
    graph <- graph_from_adjacency_matrix(A, mode = "directed", weighted=TRUE)
    cvec  = degree(graph)/(vcount(graph) - 1)

    # Add noise to A by flipping
    if (flip.edges)
    {
        A_temp <- A
        f <- floor(frac*nnzero(A))
        A[sample(which(A_temp==0),f)]=rep(1,f)
        A[sample(which(A_temp!=0),f)]=rep(0,f)
        cat("\nFraction of edges flipped: ", frac)
        cat("\nRandomly flip ",f," edges out of a total of ",nnzero(A))
    }
    B <- A+Diagonal(x=cvec)
    #cat('\nValue of k is:',d,'\n')
    if (augmented == TRUE){x <- svds(B,k=d)} else{x <- svds(A,k=d)}

    X <- x$u[,1:d]
    Y <- x$v[,1:d]

    X <- scale(X,center=FALSE,scale=1/sqrt(x$d[1:d]))
    Y <- scale(Y,center=FALSE,scale=1/sqrt(x$d[1:d]))

    data <- cbind(X,Y)
    # cat("\nDimensions of data is:",nrow(data),"x",ncol(data),"\n")

    out = list(singular_values=x, data=data)
    return(data)
}















## ----------------------------------------------------------------------------
## Running the Simulations

runSimulator <- function(
  n=19902,
  groups=c(15768,4000,1000,3000,2000,2500,2500,2000),
  num.cores= 1,
  num_runs=1,
  num.of.graphs=50,
  graph_start_offset=0,
  Grange=4:12, #(2*nrow(P)),
  modelNames='VVV', # 'VEV' , 'VVV' , etc.
  d=d_ind,
  r=Inf,
  rp=Inf,
  directed=TRUE,
  augmented=TRUE,
  sample.ratio=4, # size of sample set used for initialization
  maxP=0.2,
  scale=TRUE,
  #seed=seednum,
  thresh_d=NA, # the value at which d cuts off. But we are generally going to pick d it seems.
  compute.Pprime=FALSE, # Pprime is the proportion of verticies connected in newly assigned groups after clustering
  #Q is what we previous called P-hat. It's the centroids of the clusters, dot producted together to obtain a best estimate of the probability matrix.
  flip.edges=FALSE,
  fraction.flipped = 0)
{

    cat('\nIs the graph directed:', directed)

    cat('\nAugment adjacency matrix with weighted diagonal:', augmented)
    cat("\nValue of d is:",d)


    empty_list = vector(mode="list",length=num.of.graphs)
    all_BIC = empty_list
    all_bic=NULL; all_k=NULL
    all_classification=empty_list

    # i=1
    # while (i <= num_runs)
    # {
        # cat("\n\n\n\n#################################################################")
        # cat("\n## Run Number ##  ",i,"\t< seed",seedlist[i],">")
        # cat("\n#################################################################\n")
        # set.seed(seedlist[i])            # pick seednumber from a predefined list

        end.trial = FALSE
        mcl    = list();
        mcl.index = 0; sdIndex=0

        registerDoParallel(cores=num.cores)
        while (mcl.index < num.of.graphs)
        {
            mcl <- foreach(j=1:num.cores) %dopar%
            {

                    mcl.index = mcl.index + j
                    sd = seedlist[mcl.index]
		                graph_number = mcl.index + graph_start_offset

                    cat("\n\n\n\n#################################################################")
                    cat("\n## Graph Number ##  ",graph_number,"\t< seed",sd,">")
                    cat("\n#################################################################\n")

                    #-- pick seednumber from a predefined list
                    set.seed(sd)

                    ##-- Adjacency Matrix
                    A.name = paste("/home/kmehta7/exp/P05/binary_adj_matrix/A50/A50_",as.character(graph_number),".rds",sep="")
                    bin_A = readRDS(A.name)
                    bin_A = unname(bin_A)

                    data <- graph.spectral.embedding(A=bin_A,d=d,augmented = augmented,flip.edges=flip.edges,frac=fraction.flipped)

                    cat("\nDimensions of data of graph",graph_number,"is:",nrow(data),"x",ncol(data),"\n")
                    df=data.frame(data)


                    mclust.options('subset'=NROW(df))
                    mc <- Mclust(df[,1:(2*d)],
                                             modelNames=modelNames,
                                             G=Grange,
                                             #initialization = list(hcPairs = randomPairs(df,seed=sd),subset = 1:NROW(df) ),
                                             # initialization=list( subset=c(1:NROW(df)) ),
                                             verbose=FALSE)
            } #j foreach
            # sdIndex = sdIndex + 1

            gc()
            for (j in 1:num.cores){
                if(!is.null(mcl[[j]])){
                    mcl.index = mcl.index + 1
                    cat("\n-- Mclust successful --\n")
                    print(summary(mcl[[j]]))
                    all_BIC[[mcl.index]] = mcl[[j]]$BIC
                    all_bic[mcl.index] = mcl[[j]]$bic
                    all_k[mcl.index] = mcl[[j]]$G
                    all_classification[[mcl.index]]  = mcl[[j]]$classification
                    } else {cat("\n-- Mclust FAILED --\n");print(summary(mcl[[j]]));}
            } #j

            # save as .rData
            save( all_bic, all_BIC, all_k, all_classification,
                 file = paste(outdir,"/",filename,".rData",sep=""))

        } #mcl.index while

        # bic_classification[[i]]  = mcl[[ which.max(all_bic[[i]]) ]]$classification
        # i=i+1
    # } #i

} #func




## ----------------------------------------------------------------------------
## MAIN

# Intialize variables
tic(" Total Time")

# create output directory
outdir= "/scratch/kmehta7"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

# set filename
filename = "P05_MBHAC_4to100_d15_g100_x1"

# save as .txt file
# sink(paste(outdir,"/",filename,".txt",sep=""))

num_runs            = 1
num.of.graphs       = 100
graph_start_offset  = 0
gr                  = 4:100
d_ind               = 15
num.cores           = 20


seedlist <- readRDS("seed_list.rds")


# call main function
mainrun <- runSimulator(Grange=gr, num.of.graphs=num.of.graphs, graph_start_offset=graph_start_offset, d=d_ind, num.cores=num.cores)





cat("\n\n*****************************************************************")
cat("\n** Done! Program terminated successfully **")

toc()
sink() # end writing .txt file

cat("\n\n*****************************************************************")
cat("\n** Done! Program terminated successfully **")
