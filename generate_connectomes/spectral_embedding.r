### Graph Spectral Embedding
spectral_embedding <- function( A,
                                d = 4,
                                flip.edges = FALSE,
                                augmented = TRUE,
                                frac  = 0,
                                return_UV = FALSE)
{
    require(igraph)
    require(RSpectra)
    require(Matrix)
    
    graph <- graph_from_adjacency_matrix(A, mode = "directed", weighted=TRUE)
    cvec  = degree(graph)/(vcount(graph) - 1)

    # Add noise to A by flipping
    if (flip.edges)
    {
        A_temp <- A
        f <- floor(frac*nnzero(A))
        A[sample(which(A_temp==0),f)]=rep(1,f)
        A[sample(which(A_temp!=0),f)]=rep(0,f)
        cat("\nFraction of edges flipped: ", frac)
        cat("\nRandomly flip ",f," edges out of a total of ",nnzero(A))
    }
    B <- A+Diagonal(x=cvec)
    #cat('\nValue of k is:',d,'\n')
    if (augmented == TRUE){x <- svds(B,k=d)} else{x <- svds(A,k=d)}

    if (return_UV)
    {
            X <- x$u[,1:d]
            Y <- x$v[,1:d]

            X <- scale(X,center=FALSE,scale=1/sqrt(x$d[1:d]))
            Y <- scale(Y,center=FALSE,scale=1/sqrt(x$d[1:d]))

            data <- cbind(X,Y)
            cat("\nDimensions of data is:",nrow(data),"x",ncol(data),"\n")

            out = list(singular_values=x$d, data=data)
        } else{
            out = x$d
        }

    return(out)

}
