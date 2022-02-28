################################################################################
##--calculate bias in aprior vs conditional prob.

birthtime_labels <- function(Pyx,Py,Px,
                             perc=0.10,n=3)
{
  br_name = names(Py)

  br_labels = vector(mode='list',length=length(Px))

  col_maxs = apply(Pyx, 2, max)
  col_mins = apply(Pyx, 2, min)

  top.n = rep(0,length(br_name))
  bottom.n = rep(0,length(br_name))
  for(j in 1:length(br_name))
  {
    top.n[j] = sort(Pyx[,j],decreasing=TRUE,index.return=TRUE)$x[n]
    bottom.n[j] = sort(Pyx[,j],index.return=TRUE)$x[n]
  }

  for(i in 1:length(Px))
  {
    temp = NULL
    for(j in 1:length(br_name))
    {
      if( Pyx[i,j] >= ((1-perc)*top.n[j]))
      {
        temp = c(temp,paste("(H)",br_name[j],sep=""))
        temp = gsub(" ", "", temp, fixed = TRUE)
      }
      if( Pyx[i,j] <= ((1+perc)*bottom.n[j]))
      {
        temp = c(temp,paste("(L)",br_name[j],sep=""))
        temp = gsub(" ", "", temp, fixed = TRUE)
      }
    }
    if(is.null(temp)){
      br_labels[[i]] = paste(" ")
    } else {
      br_labels[[i]] = paste(temp,collapse=", ")
    }
  }
# browser()
  return(br_labels)
}

##----------------------------------------------------------------------------##
