################################################################################
##--calculate graph theoretic labels

gt_labels <- function(gt.name, gt.vector,
                      check.n.highest=TRUE, n=13, check.n.lowest=TRUE,
                      perc=0.01,
                      check.percentage.highest=FALSE,
                      check.thresholds=FALSE, top.thresh=200, bottom.thresh=5)
{
  # br_name = names(Py)

  gt_labels = vector(mode='list',length=length(gt.vector))

  col_max = max(gt.vector)
  col_min = min(gt.vector)

  top.n = sort(gt.vector,decreasing=TRUE,index.return=TRUE)$x[n]
  bottom.n = sort(gt.vector,index.return=TRUE)$x[n]



  if(check.thresholds)
  {
    for(i in 1:length(gt.vector))
    {
      temp = NULL

      if( gt.vector[i] >= top.thresh )
      {
        temp = c(temp,paste("(H)",gt.name,sep=""))
        temp = gsub(" ", "", temp, fixed = TRUE)
      }
      if( gt.vector[i] <= bottom.thresh )
      {
        temp = c(temp,paste("(L)",gt.name,sep=""))
        temp = gsub(" ", "", temp, fixed = TRUE)
      }

      if(is.null(temp)){
        gt_labels[[i]] = paste(" ")
      } else {
        gt_labels[[i]] = paste(temp,collapse=", ")
      }
    }
  }



  if(check.n.highest)
  {
    for(i in 1:length(gt.vector))
    {
      temp = NULL

      if( gt.vector[i] >= ((1-perc)*top.n) )
      {
        temp = c(temp,paste("(H)",gt.name,sep=""))
        temp = gsub(" ", "", temp, fixed = TRUE)
      }
      if(check.n.lowest)
      {
        if( gt.vector[i] <= ((1+perc)*bottom.n) )
        {
          if( gt.vector[i] == 0 )
          {
            temp = c(temp,paste("(Zero)",gt.name,sep=""))
            temp = gsub(" ", "", temp, fixed = TRUE)
          } else{
            temp = c(temp,paste("(L)",gt.name,sep=""))
            temp = gsub(" ", "", temp, fixed = TRUE)
          }
        }
      }

      if(is.null(temp)){
        gt_labels[[i]] = paste(" ")
      } else {
        gt_labels[[i]] = paste(temp,collapse=", ")
      }
    }#for
  }



  if(check.percentage.highest)
  {
    for(i in 1:length(gt.vector))
    {
      temp = NULL

      if( gt.vector[i] >= ((1-perc)*col_max) )
      {
        temp = c(temp,paste("(H)",gt.name,sep=""))
        temp = gsub(" ", "", temp, fixed = TRUE)
      }
      if( gt.vector[i] <= ((1+perc)*col_min) )
      {
        temp = c(temp,paste("(L)",gt.name,sep=""))
        temp = gsub(" ", "", temp, fixed = TRUE)
      }

      if(is.null(temp)){
        gt_labels[[i]] = paste(" ")
      } else {
        gt_labels[[i]] = paste(temp,collapse=", ")
      }
    }
  }


  # write.csv(unlist(gt_labels),'gt_labels.csv')
  return(gt_labels)
}

##----------------------------------------------------------------------------##




################################################################################
##-Load
# x  <- list()
# load('/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/data/results_mod_sim.rData')
# x[[1]] <- rowMeans(simplify2array(deg.clust.in))
# x[[2]] <- rowMeans(simplify2array(deg.clust.out))
# x[[3]] <- rowMeans(simplify2array(deg.clust.total))
# x[[4]] <- rowMeans(simplify2array(bet.clust))
# x[[5]] <- 1-rowMeans(simplify2array(pc.clust))
# x[[6]] <- rowMeans(simplify2array(wc.clust))
# x[[7]] <- x[[3]]-x[[6]]
#
# x.names = c('neuron_in', 'neuron_out', 'neuron_total',
#             'neuron_betCen', 'neuron_pc',
#             'neuron_wc','neuron_bc')
#
# y = NULL
# for(i in 1:length(x))
# {
#   fn.output = gt_labels(x.names[i],x[[i]])
#   y = cbind(y,unlist(fn.output))
# }
#
# write.csv(y,'gt3_labels.csv')
