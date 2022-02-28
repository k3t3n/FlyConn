################################################################################
##--calculate bias in aprior vs conditional prob.

community_labels <- function(Pyx,Py,Px)
{
  nt_name = names(Py)
  # nt_name = nt_name[-7] #remove 'unknown'

  nt_labels = vector(mode='list',length=length(Px))

  col_maxs = apply(Pyx, 2, max)
  prop_thresh = 10*Py

  for(i in 1:length(Px))
  {
    temp = NULL
    for(j in 1:length(nt_name))
    {
      if( Pyx[i,j] >= 0.67 )
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
