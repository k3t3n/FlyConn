## ----------------------------------------------------------------------------
## Packages
library(mclust)
library(caret)

cat("\014")     # clear screen
rm(list=ls())   # clean workspace




## ----------------------------------------------------------------------------
## confusion matrix from 'caret' library
conf.Mat <- function(x_hat,x,
                     print.confm=TRUE)
{
    cat("\nThe confusion matrix is:\n")
    reference <-factor(x)
    predicted <-factor(x_hat)
    u         <- union(predicted, reference)
    t         <- table(factor(predicted, u), factor(reference, u))
    confm     <- confusionMatrix(t)

    ##--sort by cluster number
    ind = sort(as.numeric(rownames(confm$table)),index.return=TRUE)$ix

    if(print.confm)
        print(confm$table[ind,ind])

    return(confm)
} # end function






## ----------------------------------------------------------------------------
## Group label sorting when ground truth is unknown
sort.original.label <- function(origData,verbose=FALSE)
{
    mc.classification = origData
    sorted.index.orig = names( sort((table(mc.classification)),decreasing=TRUE) )
    sorted.index.orig = strtoi(sorted.index.orig)

#    sorted.index= names( sort((table(mc.classification)),decreasing=TRUE) )
#    sorted.index = strtoi(sorted.index)


    #sort in reverse order
    for (j in 1:length(table(mc.classification,exclude=0)) )
    {

        sorted.index= names( sort((table(mc.classification,exclude=0)),decreasing=TRUE) )
        sorted.index = strtoi(sorted.index)

        if (verbose) print(sorted.index)
        i = sorted.index[j]

        # switch i and j in mc.classification
        for(k in 1:length(mc.classification)) #move the i's out of the way so we can later change them to j's
        {
          if(mc.classification[k]==i){
            mc.classification[k]<--1}
        }
        #Now we change the j's to i's. Remember to do this for all of mc.classification, not just those
        #in the set of interest. We're applying a permutation to all of mc.classification.
        for(k in 1:length(mc.classification))
        {
          if(mc.classification[k]==j){
            mc.classification[k]<-i}
        }
        # Change the original i's back to j.
        for(k in 1:length(mc.classification))
        {
          if(mc.classification[k]==-1){
            mc.classification[k]<-j}
        }

    }#j

    return(mc.classification)
} # end function






## ----------------------------------------------------------------------------
## Group label reordering

group.reorder <- function(original,mc.classification)
{
    gps = unname(table(original))
    # print(table(original));cat("------------------------------------------------------\n\n\n")


    groupsizeordering = as.numeric(row.names(table(original)))
    t_max = length(table(mc.classification))
    pass_num = 1

    # print(table(mc.classification));cat("------------------------------------------------------\n\n\n")


    for(i in groupsizeordering) # we are going through from most represented vertices to least in the generated graph clustered
    {
      vertices<-which(original==i) # the original vertices wtih the number i.
      #print("The classification of these same vertices are as follows:")
      #print(mc.classification[vertices])
      #find the most commonly occurring number j >= i on those vertices in mc.classification
      # There isn't a function for mode directly so you make a table that sorts but the values that occur and counts them, and then you take the max in the table
      y<-table(mc.classification[vertices])
      # cat("The table y is:"); print(y)
      j<-as.numeric(names(y)[y==max(y)])[1]
      # cat(j, "is the most common vertex type\n")
      # browser()
      if (length(which(groupsizeordering==j))==0 )
      {
          mc.classification[which(mc.classification==j)]=i
          # cat("\nPass #",pass_num); pass_num=pass_num+1;
          # print(table(mc.classification))
      }else{
          if(  (j>length(gps))  || (which(groupsizeordering==j)>which(groupsizeordering==i))  )
          {
            # if the group with label j occurs earlier (has more vertices) than the group with label i,
            #Now we don't need to worry about j<i, but rather that we have already looked the number associated to j.
            # so that means that j occurs before i in the vector groupsizeordering

            # switch i and j in mc.classification
            for(k in 1:length(mc.classification)) #move the i's out of the way so we can later change them to j's
            {
              if(mc.classification[k]==i){
                mc.classification[k]<--1}
            }
            #Now we change the j's to i's. Remember to do this for all of mc.classification, not just those
            #in the set of interest. We're applying a permutation to all of mc.classification.
            for(k in 1:length(mc.classification))
            {
              if(mc.classification[k]==j){
                mc.classification[k]<-i}
            }
            # Change the original i's back to j.
            for(k in 1:length(mc.classification))
            {
              if(mc.classification[k]==-1){
                mc.classification[k]<-j}
            }
            # cat("\nPass #",pass_num); pass_num=pass_num+1;
            # print(table(mc.classification))
          }


          if (length(y)>1)
          {
            y_sort = sort(table(mc.classification[vertices]),decreasing=TRUE)
            y_names = as.numeric( names(y_sort) )

                for ( nn in 2:length(y) )
                {
                      #y = table(mc.classification[vertices],exclude = i)
                      #j<-as.numeric(names(y)[y==max(y)])

                      if ( unname(y_sort)[nn] > unname(table(mc.classification))[t_max] )
                      {
                        j = y_names[nn]

                        # switch i and j in mc.classification
                        for(k in 1:length(mc.classification)) #move the i's out of the way so we can later change them to j's
                        {
                          if(mc.classification[k]==t_max){
                            mc.classification[k] <- -1}
                        }
                        #Now we change the j's to i's. Remember to do this for all of mc.classification, not just those
                        #in the set of interest. We're applying a permutation to all of mc.classification.
                        for(k in 1:length(mc.classification))
                        {
                          if(mc.classification[k]==j){
                            mc.classification[k] <- t_max}
                        }
                        # Change the original i's back to j.
                        for(k in 1:length(mc.classification))
                        {
                          if(mc.classification[k]==-1){
                            mc.classification[k] <- j}
                        }

                        # cat("\nPass #",pass_num); pass_num=pass_num+1;
                        # print(table(mc.classification))
                        t_max=t_max-1
                        # cat("t_max: ",t_max,"\n")
                  }

                } #nn
          }

      }

# cat("------------------------------------------------------\n\n\n")

    } # for groupsizeordering

return(mc.classification)
} # end funtion







###############################################################################
## ----------------------------------------------------------------------------
##-User Input
file.descrip.1 = "P05_MBHAC_4to100_d11_g1_g100_mc75"
file.descrip.2 = "P05_MBHAC_4to100_d11_g1_g400_mc75"
index.range = 1:2 # no. of clusterings being compared

empty_list = vector(mode="list",length=1)
fc_cleaned = empty_list

# fc_cleaned[[1]] = readRDS('/media/WDHDD/clustering/exp/ARGO/P05/IVC/recheck/cl_95.rds')
load(paste("/media/WDHDD/ckt/IVC/fcl_results/fcl_",file.descrip.1,".rData",sep=""))
fc_cleaned[[1]] = fcl[[95]]
load(paste("/media/WDHDD/ckt/IVC/fcl_results/fcl_",file.descrip.2,".rData",sep=""))
fc_cleaned[[2]] = fcl[[381]]



##--------------------
# Aconn = readRDS('/media/WDHDD/clustering/exp/data/A_conn.rds')
# neuron_name = attributes(Aconn)$Dimnames[[1]]
# ntwk = readRDS("/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/ensemble_clustering/ntwk_d11.rds")
# # zero.vert = which(ntwk==0)
# # ntwk = ntwk[-zero.vert]
# fc_cleaned[[2]] = ntwk

# fc_cleaned[[1]] = sort.original.label(fc_cleaned[[1]])
# fc_cleaned[[1]][which(fc_cleaned[[1]]>=55)]=0





fc_clabeled = vector(mode="list",length=max(index.range))
for(k in 1:(length(index.range)-1))
{
    i =  index.range[k]
    j = index.range[k+1]
    if ( length(table(fc_cleaned[[i]])) > length(table(fc_cleaned[[j]])) )
    {
        fc_clabeled[[j]] = group.reorder(fc_cleaned[[i]],fc_cleaned[[j]])
        fc_clabeled[[i]]   = fc_cleaned[[i]] #fc_clabeled[[i]] = group.reorder(fc_clabeled[[j]],fc_cleaned[[i]])
    } else
    {
        fc_clabeled[[i]] = group.reorder(fc_cleaned[[j]],fc_cleaned[[i]])
        fc_clabeled[[j]]   = fc_cleaned[[j]] #fc_clabeled[[j]] = group.reorder(fc_clabeled[[i]],fc_cleaned[[j]])
    }

    cat('\n\n\n\n\n\n-----------------------------------------')
    cat('\nClustering:',i)
    cat('\n-----------------------------------------')
    print(table(fc_clabeled[[i]]))
    cat('\n\n-----------------------------------------')
    cat('\nClustering:',j)
    cat('\n-----------------------------------------')
    print(table(fc_clabeled[[j]]))
    cat('\n\n-----------------------------------------')
    cat('\nARI:',adjustedRandIndex(fc_clabeled[[i]],fc_clabeled[[j]]))
    cat('\n-----------------------------------------')
    cat('\n')
    confm = conf.Mat(fc_clabeled[[i]], fc_clabeled[[j]], print.confm=TRUE)
}


# for(i in index.range){ print(table(fc_cleaned[[i]]))  }
