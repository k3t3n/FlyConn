## user set parameters
mean_conn = 0.5
num.graphs = 5
gen_new_Aprob = FALSE
save_new_Aprob = FALSE
gen_new_graphs = FALSE
save_new_graphs = FALSE
d = 200
calc.svd = TRUE
calc.elbow = TRUE
##---------------------


seedlist = seq(19,19*num.graphs,by = 19)

##-- load Aconn
Aconn = readRDS('/media/WDHDD/clustering/exp/data/A_conn.rds')


##-- generate or load Aprob
cat('\n#####################################################\n')
aprob_filename = paste('/media/WDHDD/clustering/exp/data/A_prob_',as.character(mean_conn*100),'.rds',sep="")
if (gen_new_Aprob)
    {
        cat('\nGenerating probability matrix..')
        source('/media/WDHDD/clustering/exp/Aprob_gen.r')
        Aprob = Aprob_gen(mean_conn)


        if(save_new_Aprob)
        {
            cat('\nSaving probability matrix..')
            saveRDS(Aprob, aprob_filename)
        }
    } else
    {
        cat('\nLoaded existing probability matrix..')
        Aprob = readRDS(aprob_filename)
    }


##-- start
dir.name = paste('/media/WDHDD/clustering/exp/data/binary_adj_matrix/A',as.character(mean_conn*100),sep="")
for (i in 1:num.graphs)
{
    cat('\n\n#####################################################\n\n')
    graph_filename = paste(dir.name,'/A',as.character(mean_conn*100),'_',as.character(i),'.rds',sep="")
    ##-- generate random binary graphs
    if (gen_new_graphs)
    {
        source('/media/WDHDD/clustering/exp/binary_adj_matrix.r')

        cat('\nGenerating random graph#',i)
        set.seed(seedlist[i])
        binA = binary_adj_matrix(Aconn,Aprob)

        if(save_new_graphs)
        {
            cat('\nSaving random graph#',i)
            dir.create(dir.name,recursive=TRUE,showWarnings=FALSE)
            saveRDS(binA, graph_filename)
        }
    } else
    {
        cat('\nLoaded existing graph#',i)
        binA=readRDS(graph_filename)
    }


    ##-- compute singular values
    if(calc.svd)
    {
        source('/media/WDHDD/clustering/exp/spectral_embedding.r')
        singular_values = spectral_embedding(A=binA, d=d)
    }


    ##-- find elbow
    if(calc.elbow)
    {
        source('/media/WDHDD/clustering/exp/find_elbow.r')
        find_elbow(singular_values)
    }

} # for
