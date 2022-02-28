################################################################################
##--calculate bias in aprior vs conditional prob.

neurotrans_labels <- function(Pyx,Py,Px,
                              perc=0.10,n=5)
{
  nt_name = names(Py)
  nt_name = nt_name[-7] #remove 'unknown'

  nt_labels = vector(mode='list',length=length(Px))

  col_maxs = apply(Pyx, 2, max)
  prop_thresh = 10*Py

  top.n = rep(0,length(nt_name))
  bottom.n = rep(0,length(nt_name))
  for(j in 1:length(nt_name))
  {
    top.n[j] = sort(Pyx[,j],decreasing=TRUE,index.return=TRUE)$x[n]
    bottom.n[j] = sort(Pyx[,j],index.return=TRUE)$x[n]
  }

  for(i in 1:length(Px))
  {
    temp = NULL
    for(j in 1:length(nt_name))
    {
      if( Pyx[i,j] >= ((1-perc)*top.n[j]) || Pyx[i,j]>=prop_thresh[j] || Pyx[i,j]>=0.5)
        temp = c(temp,nt_name[j])
    }
    if(is.null(temp)){
      nt_labels[[i]] = paste(" ")
    } else {
      nt_labels[[i]] = paste(temp,collapse=", ")
    }
  }
# browser()
  return(nt_labels)
}

##----------------------------------------------------------------------------##
