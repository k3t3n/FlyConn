## ----------------------------------------------------------------------------
## Packages
library(mclust)
library(caret)
library(aricode)

cat("\014")     # clear screen
rm(list=ls())   # clean workspace



## ----------------------------------------------------------------------------
# INPUT:  table = The table of pairwise crossclassifications between two partitions.
#         The number of rows equals the number of clusters in the first partition.
#         The number of columns equals the number of clusters in the second partition.
#
# OUTPUT: ari = the adjusted Rand index of the table
#         vari = the variance of the adjusted Rand index under the hypergeometric distribution
#         vari_multi = the variance of the adjusted Rand index under the multinomial distribution
# ----------------------------------------------------------------------------
ari_prog <- function(P,Q,clip_zeros)
{
  if(clip_zeros)
  {
    P[which(Q==0)]=0
    Q[which(P==0)]=0
    P=P[-which(P==0)]
    Q=Q[-which(Q==0)]
  }

  table = table(P,Q)
  perc_retained = (length(P)/19902)*100

  n <- sum(table)
  rm <- colSums(table)
  cm <- rowSums(table)
  tsq <- table*table
  rmsq <- rm*rm
  cmsq <- cm*cm
  nsq <- n*n
  a <- (sum(sum(tsq)) - n)/2
  b <- (sum(rmsq) - sum(sum(tsq)))/2
  c <- (sum(cmsq) - sum(sum(tsq)))/2
  d <- (sum(sum(tsq)) + nsq - sum(rmsq) - sum(cmsq))/2
  nc2 <- choose(n,2)
  ari <- (nc2*(a + d) - ((a + b)*(a + c) + (c + d)*(b + d)))/(nc2*nc2 - ((a + b)*(a + c) + (c + d)*(b + d)))
  tcub <- table*table*table
  rmcub <- rm*rm*rm
  cmcub <- cm*cm*cm
  e <- 2*sum(rmsq) - n*(n + 1)
  f <- 2*sum(cmsq) - n*(n + 1)
  g <- 4*sum(rmcub) - 4*(n + 1)*sum(rmsq) + n*(n+1)^2
  h <- n*(n - 1)
  i <- 4*sum(cmcub) - 4*(n + 1)*sum(cmsq) + n*(n+1)^2
  vad <- 1/16*(2*n*(n-1) - ((e*f)/(n*(n-1)))^2 + (4*(g-h)*(i-h))/(n*(n-1)*(n-2)) + (e^2 - 4*g+2*h)*(f^2 - 4*i+2*h)/(n*(n-1)*(n-2)*(n-3)))
  vari <- (nc2*nc2*vad)/(nc2*nc2 - ((a + b)*(a + c) + (c + d)*(b + d)))^2
  r1 <- dim(table)[1]; c1 <- dim(table)[2]
  t <- table
  t1 <- matrix(0,r1,c1)
  t2 <- matrix(0,r1,c1)
  for (i in 1:r1){
    for (j in 1:c1){
      t1[i,j] <- t[i,j]*(2*t[i,j] - (rm[j] + cm[i]))^2
      t2[i,j] <- t[i,j]*(2*t[i,j] - (rm[j] + cm[i]))
    }
  }
  vargam <- ((2/n)^4)*(sum(t1) - (1/n)*(sum(t2))^2)
  varrand <- (1/16)*((n*(n-1))^2)*vargam
  vari_multi <- (nc2*nc2*varrand)/(nc2*nc2 - ((a + b)*(a + c) + (c + d)*(b + d)))^2

  out_list = list(ari=ari, vari=vari, vari_multi=vari_multi, perc_retained=perc_retained)
  return(out_list)
}




## ----------------------------------------------------------------------------
##--Confidence Interval
confi_int = function(ari, vari)
{
  lb = ari-1.96*sqrt(vari)
  ub = ari+1.96*sqrt(vari)
  c(lb,ub)
}




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





## ----------------------------------------------------------------------------
## Create final clustering by setting minimum cluster size
adjust.cluster.size <- function(file.descrip, thresh, remove.below)
{
  ## auto generating filenames
  # pi.filename     = paste("/media/WDHDD/ckt/IVC/IVC_data/IVC_",file.descrip,".rData",sep="")
  CM.filename     = paste("/media/WDHDD/ckt/IVC/CM_data/CM_",file.descrip,".rds",sep="")
  input.filename  = paste("/media/WDHDD/ckt/IVC/cl_results/cl_",file.descrip,".rds",sep="")

  ##--Load Data
  cl <- readRDS(input.filename)
  CM <- readRDS(CM.filename)

  num.vertices = length(cl[[100]])
  # cat('\n\nNumber of vertices:',num.vertices)

  ## --------------------------------------------------------------------------
  cl_sorted = sort.original.label(cl[[thresh]])
  lost_index = min(which(unname(table(cl_sorted))<remove.below))
  fcl = cl_sorted
  fcl[which(cl_sorted>=lost_index)] = 0

  return(fcl)

} #end function











###############################################################################
## ----------------------------------------------------------------------------
##-User Input
d     = c(15,15,15,15,15)
nG    = c(100,100,100,100,100)
tau   = c(0.95,0.95,0.95,0.95,0.95)
csize = c(150,125,100,75,50)
do_nG = FALSE
do_tau = FALSE
do_csize = TRUE


compute.ARIMatrix = TRUE
compute.confusionMatrix = FALSE


##-auto generating
tau = round(tau*nG)
empty_list = vector(mode="list",length=1)
fc_cleaned = empty_list
file.descrip = empty_list
if(do_nG){
  for(k in 1:length(d))
  {
    if(k==1){
      gname = '_g1_g100'
    } else if(k==2){
      gname = '_g101_g300'
    } else if(k==3){
      gname = '_g301_g700'
    }
    file.descrip[[k]] = paste('P05_MBHAC_4to100_d',d[k],gname,sep='')
    fc_cleaned[[k]] = adjust.cluster.size(file.descrip[[k]], tau[[k]], csize[[k]])
  }
} else {
  for(k in 1:length(d))
  {
    if(k==1){
      gname = '_g1_g100'
    } else if(k==2){
      gname = '_g101_g200'
    } else if(k==3){
      gname = '_g201_g300'
    } else if(k==4){
      gname = '_g301_g400'
    } else if(k==5){
      gname = '_g401_g500'
    }
    file.descrip[[k]] = paste('P05_MBHAC_4to100_d',d[k],gname,sep='')
    fc_cleaned[[k]] = adjust.cluster.size(file.descrip[[k]], tau[[k]], csize[[k]])
  }
}

if(do_nG){
  rnames = nG
  cnames = nG
  file.descrip[[4]] = paste('P05_MBHAC_4to100_d',d[k],'_g201_g300',sep='')
  file.descrip[[5]] = paste('P05_MBHAC_4to100_d',d[k],'_g201_g400',sep='')
  if(d[k]==15) file.descrip[[6]] = paste('P05_MBHAC_4to100_d',d[k],'_g1_g400',sep='')
  if(d[k]==11) file.descrip[[6]] = paste('P05_MBHAC_4to100_d',d[k],'_g401_g800',sep='')
  fc_cleaned[[4]]  = adjust.cluster.size(file.descrip[[4]], tau[[1]], csize[[1]])
  fc_cleaned[[5]]  = adjust.cluster.size(file.descrip[[5]], tau[[2]], csize[[2]])
  fc_cleaned[[6]]  = adjust.cluster.size(file.descrip[[6]], tau[[3]], csize[[3]])
}
if(do_tau){
  rnames = tau
  cnames = tau
  file.descrip[[4]] = paste('P05_MBHAC_4to100_d',d[k],'_g301_g400',sep='')
  file.descrip[[5]] = paste('P05_MBHAC_4to100_d',d[k],'_g401_g500',sep='')
  fc_cleaned[[4]]  = adjust.cluster.size(file.descrip[[4]], tau[[1]], csize[[1]])
  fc_cleaned[[5]]  = adjust.cluster.size(file.descrip[[5]], tau[[2]], csize[[1]])
  fc_cleaned[[6]]  = adjust.cluster.size(file.descrip[[1]], tau[[3]], csize[[1]])
}
if(do_csize){
  rnames = csize
  cnames = csize
  fc_cleaned[[6]]  = adjust.cluster.size(file.descrip[[2]], tau[[1]], csize[[1]])
  fc_cleaned[[7]]  = adjust.cluster.size(file.descrip[[3]], tau[[1]], csize[[2]])
  fc_cleaned[[8]]  = adjust.cluster.size(file.descrip[[4]], tau[[1]], csize[[3]])
  fc_cleaned[[9]]  = adjust.cluster.size(file.descrip[[5]], tau[[1]], csize[[4]])
  fc_cleaned[[10]] = adjust.cluster.size(file.descrip[[1]], tau[[1]], csize[[5]])
}


##-----------------------------------------------------------------------------
if(compute.confusionMatrix)
{
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
}


##-----------------------------------------------------------------------------
if(compute.ARIMatrix)
{
  ARI.matrix1 = matrix(0,length(d),length(d))
  ARI.matrix2  = matrix(0,length(d),length(d))
  indx = 0
  for(i in 1:length(d))
  {
    for(j in 1:length(d))
    {

      if(i<j){
        out1 = ari_prog(fc_cleaned[[i]],fc_cleaned[[j]],clip_zeros=TRUE)
        out2 = ari_prog(fc_cleaned[[i]],fc_cleaned[[j]],clip_zeros=FALSE)
      } else if(i==j) {
        if(do_tau || do_nG){
              if(i==1){
                out1 = ari_prog(fc_cleaned[[i]],fc_cleaned[[4]],clip_zeros=TRUE)
                out2 = ari_prog(fc_cleaned[[i]],fc_cleaned[[4]],clip_zeros=FALSE)
              }
              if(i==2){
                out1 = ari_prog(fc_cleaned[[i]],fc_cleaned[[5]],clip_zeros=TRUE)
                out2 = ari_prog(fc_cleaned[[i]],fc_cleaned[[5]],clip_zeros=FALSE)
              }
              if(i==3){
                out1 = ari_prog(fc_cleaned[[i]],fc_cleaned[[6]],clip_zeros=TRUE)
                out2 = ari_prog(fc_cleaned[[i]],fc_cleaned[[6]],clip_zeros=FALSE)
              }
        }
        if(do_csize){
              if(i==1){
                out1 = ari_prog(fc_cleaned[[i]],fc_cleaned[[6]],clip_zeros=TRUE)
                out2 = ari_prog(fc_cleaned[[i]],fc_cleaned[[6]],clip_zeros=FALSE)
              }
              if(i==2){
                out1 = ari_prog(fc_cleaned[[i]],fc_cleaned[[7]],clip_zeros=TRUE)
                out2 = ari_prog(fc_cleaned[[i]],fc_cleaned[[7]],clip_zeros=FALSE)
              }
              if(i==3){
                out1 = ari_prog(fc_cleaned[[i]],fc_cleaned[[8]],clip_zeros=TRUE)
                out2 = ari_prog(fc_cleaned[[i]],fc_cleaned[[8]],clip_zeros=FALSE)
              }
              if(i==4){
                out1 = ari_prog(fc_cleaned[[i]],fc_cleaned[[9]],clip_zeros=TRUE)
                out2 = ari_prog(fc_cleaned[[i]],fc_cleaned[[9]],clip_zeros=FALSE)
              }
              if(i==5){
                out1 = ari_prog(fc_cleaned[[i]],fc_cleaned[[10]],clip_zeros=TRUE)
                out2 = ari_prog(fc_cleaned[[i]],fc_cleaned[[10]],clip_zeros=FALSE)
              }
        }

      }

        if(i<=j){
          ARI.matrix1[i,j] = paste(round(out1$ari,4)," \u00B1 ",round(1.96*sqrt(out1$vari_multi),4)," (",round(out1$perc_retained,2),"%)",sep="")
          ARI.matrix2[i,j] = paste(round(out2$ari,4),"\u00B1",round(1.96*sqrt(out2$vari_multi),4))
        }


    }
  }
  rownames(ARI.matrix1)=rnames; colnames(ARI.matrix1)=cnames
  rownames(ARI.matrix2)=rnames; colnames(ARI.matrix2)=cnames
  write.csv(ARI.matrix1, 'ARI.matrix1.csv')
  write.csv(ARI.matrix2, 'ARI.matrix2.csv')
}
