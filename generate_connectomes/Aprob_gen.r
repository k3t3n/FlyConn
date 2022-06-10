Aprob_gen <- function(mean_conn=0.5, A_conn)
# select the desired mean connectivity
{

    ### packages
    require(Matrix)
    require(blockmodeling)
    require(caret)
    require(latex2exp)
    graphics.off()  # close plots
    cat("\014")     # clear screen


    ## load matrix
    # A_conn = readRDS('../data/A_conn.rds')


    ### remove all-zero zero rows and columns
        if ( length(which(colSums(A_conn)==0)) && length(which(rowSums(A_conn)==0)) )
        {hasnz = TRUE} else
        {hasnz = FALSE}
        while(hasnz==TRUE)
        {
            coldel = which(colSums(A_conn)==0)
            A_conn = A_conn[-coldel,-coldel]
            rowdel = which(rowSums(A_conn)==0)
            A_conn = A_conn[-rowdel,-rowdel]
            if ( length(which(colSums(A_conn)==0)) && length(which(rowSums(A_conn)==0)) )
                {hasnz = TRUE} else {hasnz = FALSE}
        }


    ### connection Probability Matrix: Binomial distribution
        i=0; p_vec = c(1:100)/100
        m=integer(length(p_vec))
        for (p in p_vec)
        {
            i=i+1
            q=(1-(1-p)^(A_conn@x))
            m[i] = mean(q)
        }

        #pick value of p which gives us the desired mean connectivity
        p_select = p_vec[ which.min(abs(m-mean_conn)) ]
        cat('\nMean connectivity desired =',mean_conn)
        cat('\nat p_conn =',p_select,'\n')

        # create matrix
        A_prob=(1-(1-p_select)^(A_conn@x))

        return(A_prob)
}
