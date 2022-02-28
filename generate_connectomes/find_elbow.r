find_elbow <- function( singular_values, verbose=TRUE, getelb.n=3 )
{

    require(ggplot2)
    require(reshape)
    require(plotrix)
    require(grid)
    #require(extrafont)
    require(inflection)

    elbow_finder <- function(x_values, y_values) {
      # Max values to create line
      max_x_x <- max(x_values)
      max_x_y <- y_values[which.max(x_values)]
      max_y_y <- max(y_values)
      max_y_x <- x_values[which.max(y_values)]
      max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

      # Creating straight line between the max values
      fit <- lm(max_df$y ~ max_df$x)

      # Distance from point to line
      distances <- c()
      for(i in 1:length(x_values)) {
        distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
      }

      # Max distance point
      x_max_dist <- x_values[which.max(distances)]
      y_max_dist <- y_values[which.max(distances)]

      return(c(x_max_dist, y_max_dist))
    }



    maximum_curvature <- function(x_values, y_values) {
      x=x_values; y = y_values;
      delta=(max(x)-min(x)+1)/length(x)
      dy = double(length(y)); ddy=dy; curv=dy; deriv_curv = curv;
      for (t in c( 2:(length(x)-1) ) )
      {
          dy[t]   = (y[t+1]-y[t-1])/(2*delta)
          ddy[t]  = (y[t+1]-2*y[t]+y[t-1])/(delta^2)
          curv[t] = abs(ddy[t])/( (1 +  dy[t]^2)^(3/2) )
      }

      # deriv <- function(x, y) diff(y) / diff(x)
      # middle_pts <- function(x) x[-1] - diff(x) / 2
      # first_d <- deriv(x_values, y_values)
      # second_d <- deriv(middle_pts(x_values), first_d)
      # curv = abs(second_d) / ( (1 +  first_d^2)^(3/2) )

      max_x = which.max(curv)
      max_y = y[which.max(curv)]
      return(c(max_x,max_y))
    }

    closest_orig <- function(x_values, y_values) which.min(sqrt(x_values^2 + y_values^2))


    ## Maximum Profile Likelihood ## Zhu, Mu and Ghodsi, Ali (2006)
    source('/media/WDHDD/clustering/svd/elbow/getElbows.R')


    ####################################################################################################
    #singular_values=readRDS('/home/noob/WDHDD/clustering/exp/svd_values_Aconn.rds')
    #d=singular_values[1:25]
    d=singular_values
    x=c(1:length(d))
    #browser()
    e1=elbow_finder(x,d)
    e2=maximum_curvature(x,d)
    e3=ede(x,d,0)
    e4=closest_orig(x,d)
    e5=getElbows(d,n=getelb.n,plot=FALSE)
    if(verbose){
        cat("\014")
        cat("\n Maximum straight line distance: ",e1[1])
        cat("\n Inflection point: ",e3[1])
        cat("\n Maximum Curvature: ",e2[1])
        cat("\n Closest to origin: ",e4)
        cat("\n Profile Likelihood: ",e5)
        cat("\n\n")
    }

    out = list(MLD=e1[1], IP=e3[1], MC=e2[1], C2O=e4, PL=e5)
    return(out)

} #end_function
