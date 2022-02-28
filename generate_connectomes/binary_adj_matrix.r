### Adjacency Binary Matrix
binary_adj_matrix <- function(A_conn, A_prob)
{
    #create the binary Adjacency Matrix
    u   = runif(length(A_prob), min = 0, max = 1) #sample from unifrom distribution
    Adj_entries = (A_prob>u)*1
    Adj = A_conn
    Adj@x = Adj_entries


    # remove all-zero zero rows and columns
    if ( length(which(colSums(Adj)==0))!=0 || length(which(rowSums(Adj)==0))!=0 )
            {hasnz = TRUE} else  {hasnz = FALSE}
    while(hasnz==TRUE)
      {
        if ( length(which(colSums(Adj)==0))!=0 )
        {
            coldel = which(colSums(Adj)==0)
            Adj = Adj[-coldel,-coldel]
        }
        if ( length(which(rowSums(Adj)==0))!=0 )
        {
            rowdel = which(rowSums(Adj)==0)
            Adj = Adj[-rowdel,-rowdel]
        }
        if ( length(which(colSums(Adj)==0))!=0 || length(which(rowSums(Adj)==0))!=0 )
            {hasnz = TRUE} else {hasnz = FALSE}
      }

      return(Adj)
}
