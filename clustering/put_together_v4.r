## ----------------------------------------------------------------------------
## Preamble
rm(list=ls())   # clean workspace
cat("\014")     # clear screen


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




################################################################################################
################################################################################################
##user input
tau = 0.90
file.descrip = "P05_MBHAC_4to100_d15_g101_g200"


## auto generating filenames
pi.filename = paste("/media/WDHDD/ckt/IVC/IVC_data/IVC_",file.descrip,".rData",sep="")
CM.filename = paste("/media/WDHDD/ckt/IVC/CM_data/CM_",file.descrip,".rds",sep="")
output.filename = paste("/media/WDHDD/ckt/IVC/cl_results/cl_",file.descrip,".rds",sep="")

# load data
load(pi.filename); rm(pi)
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
