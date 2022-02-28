## ----------------------------------------------------------------------------
## Preamble
library(mclust)
library(Matrix)
rm(list=ls())
cat("\014")





## ----------------------------------------------------------------------------
## Combine clusterings
concatenate.data <- function(file_prefix, dirname, concat.filename,
                             num.files = 8,
                             num.graphs.per.file = 100)
{
  for(fileindex in 1:num.files)
  {
      load( paste(dirname,file_prefix,fileindex,".rData",sep="") )
      cat("\nWorking on file",fileindex)
      if(fileindex == 1)
      {
        ab = all_bic
        ak = all_k
        AB = all_BIC
        aclass = all_classification
      } else{
        ab = c(ab, all_bic)
        ak = c(ak, all_k)
        AB = c(AB, all_BIC)
        aclass = c(aclass, all_classification)
      }
  }

  all_bic = ab
  all_k =ak
  all_BIC = AB
  all_classification = aclass

  save( all_bic, all_k, all_BIC, all_classification,
        file = concat.filename)
  cat("\n\nSuccessfully saved data!!")
  cat("\n#############################\n")
}#end func





## ----------------------------------------------------------------------------
## Create a data frame of all clusterings
create.df <- function(all_classification, df.filename, num.graphs)
{
        # Aconn = readRDS("/media/WDHDD/clustering/exp/data/A_conn.rds")
        # neuron_name = attributes(Aconn)$Dimnames[[1]]

        neuron_name = readRDS("/media/WDHDD/ckt/data/neuron_names.rds")
        df=data.frame("neuron_name"=neuron_name)

        for(gr.index in 1:num.graphs)
        {
            cat("\nWorking on graph",gr.index)
            temp=rep(NA,length(neuron_name))
            gr.name = paste("/home/noob/WDHDD/ckt/data/binary_adj_matrix/A50/A50_",toString(gr.index),".rds",sep="")
            Ap = readRDS(gr.name)
            for (v.index in 1:length(neuron_name))
            {
                 p.index = which( attributes(Ap)$Dimnames[[1]] == neuron_name[v.index] )
                 if(length(p.index)!=0)  temp[v.index] = all_classification[[gr.index]][p.index]
            }
            df = cbind(df,temp)
            names(df)[(gr.index+1)] = paste("G",toString(gr.index),sep="")
        }

        saveRDS(df,file=df.filename)
        return(df)
}





## ----------------------------------------------------------------------------
## User Input Variables
concat_data = TRUE
#-- only if concat_data = TRUE
        num.files=7
        num.graphs.per.file=100
num.graphs = 700
create_dataframe = TRUE
file.descrip = "P05_MBHAC_4to100_d15"
dirname = "/home/noob/WDHDD/ckt/data/ARGO/A50/"
# remove_miss_vertices = FALSE ## CAUTION!! ## Need to create new dataframe


## auto generating filenames
file_prefix     = paste(file.descrip,"_g",num.graphs.per.file,"_x",sep="")
concat.filename = paste(dirname,file.descrip,"_g",num.graphs,".rData",sep="")
df.filename     = paste("df_data/df_",file.descrip,"_g",num.graphs,".rds",sep="")
## ----------------------------------------------------------------------------





## ----------------------------------------------------------------------------
## MAIN
if(concat_data)
{
  cat("\nConcatenating individual data files")
  cat("\n------------------------------------\n")
  concatenate.data(file_prefix, dirname, concat.filename, num.files, num.graphs.per.file)
}

cat("\nLoading concatenated data file...")
load(concat.filename)

if(create_dataframe)
{
    cat("\nCreating new dataframe")
    cat("\n------------------------\n")
    create.df(all_classification, df.filename, num.graphs)
}

# ##--
# if (remove_miss_vertices)
# {
#     cat("\n**CAUTION** Need to create new dataframe!!")
#     browser()
#     df = readRDS('/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/data/df_removed_g100_d--.rds")
#     cat("\nLoaded removed/cleaned data set!")
# } else{
#     df.filename = paste("/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/data/df_g100_d",d.used,".rds",sep="")
#     df = readRDS(df.filename)
#     cat("\nLoaded entire data set!")
# }
#
# ##--create mcl list
# mcl <- list()
# for(i in 1:num.graphs)
# {
#     mcl[[i]] = df[,(i+1)]
#     if (!remove_miss_vertices)
#         mcl[[i]][which(is.na(mcl[[i]]))]=0
# }
#
# if (!remove_miss_vertices)
#     cat("\nAll missing vertices assigned as 0 !!")

cat("\n\n*****************************************************************")
cat("\n** Done! Program terminated successfully **")
