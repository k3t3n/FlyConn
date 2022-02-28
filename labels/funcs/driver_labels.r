################################################################################
##--calculate bias in aprior vs conditional prob.

driver_labels <- function(Pyx,Py,Px,
                          perc=0.05,n=3)
{
  dr_name = names(Py)
  # dr_name = dr_name[-7] #remove 'unknown'

  dr_labels = vector(mode='list',length=length(Px))

  col_maxs = apply(Pyx, 2, max)
  prop_thresh = 3*Py

  top.n = rep(0,length(dr_name))
  bottom.n = rep(0,length(dr_name))
  for(j in 1:length(dr_name))
  {
    top.n[j] = sort(Pyx[,j],decreasing=TRUE,index.return=TRUE)$x[n]
    bottom.n[j] = sort(Pyx[,j],index.return=TRUE)$x[n]
  }

  for(i in 1:length(Px))
  {
    temp = NULL
    for(j in 1:length(dr_name))
    {
      if( (Pyx[i,j]>0) && (Pyx[i,j] >= ((1-perc)*top.n[j]) || Pyx[i,j]>=prop_thresh[j] || Pyx[i,j]>=0.5) )
        temp = c(temp,dr_name[j])
    }
    if(is.null(temp)){
      dr_labels[[i]] = paste(" ")
    } else {
      dr_labels[[i]] = paste(temp,collapse=", ")
    }
  }
# browser()
  return(dr_labels)
}

##----------------------------------------------------------------------------##
