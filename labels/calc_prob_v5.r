##----------------------------------------------------------------------------##
##-- Preamble
library(RColorBrewer)
cat("\014")     # clear screen
rm(list=ls())   # clear workspace


##----------------------------------------------------------------------------##




################################################################################
# ##--correcting name errors
# A = readRDS('/home/noob/WDHDD/clustering/exp/data/A_conn.rds')
#
# A_index = c(1:nrow(A))
# A_names = A@Dimnames[[1]][A_index]
# ## correct errors
# A_names[290]    = "5-HT1B-F-000000"
# A_names[324]    = "Cha-F-000000"
# A_names[1040]   = "Cha-F-100278"
# A_names[3328]   = "E0585-F-000000"
# A_names[3492]   = "fru-F-000000"
# A_names[6068]   = "Gad1-F-000000"
# A_names[13265]  = "npf-F-000000"
# A_names[13414]  = "Tdc2-F-000000"
# A_names[13674]  = "TH-F-000000"
# A_names[13701]  = "TH-F-000034"
# A_names[14077]  = "Trh-F-000000"
# A_names[14986]  = "VGlut-F-000000"
# A_names[16175]  = "VGlut-F-200382"
# A_names[17674]  = "VGlut-F-500042"
# A_names[18021]  = "VGlut-F-500471"
# A_names[18022]  = "VGlut-F-500472"
##----------------------------------------------------------------------------##




################################################################################
ACl.create <- function(file.descrip, bottom_index, top_index, tau, min.clust, d.used)
{


  # df.filename = paste("/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/data/df_g100_d",
  #                          d.used,".rds",sep="")
  #
  # df = readRDS(df.filename)
  m = readRDS('/home/noob/WDHDD/clustering/scraping/metadata.rds')

  ##--create dataframe with neuron names
  neuron_name = readRDS('/media/WDHDD/ckt/data/neuron_names.rds')
  ACl = data.frame("neuron_name"=neuron_name)


  ##--for the chosen class-type
  ## 1. neuron_name
  ## 2. driver
  ## 3. neurotrans
  ## 4. gender
  ## 5. age
  ## 6. birthtime


  ACls.newcol <- function(neuronName, ctype)
  {
        ##--creates a single empty class vector
        class_vec = rep(0,length(neuronName))
        ## loop through all neurons
        m_index = NULL
        for(n in 1:length(neuronName))
        {
          m_index = which(m$name==neuronName[n])
          class_vec[n] = m[[ctype]][m_index]
        }
        return(class_vec)
  }


  ##-driver
  ctype = 2
  ACl$driver = ACls.newcol(neuron_name, ctype)


  ##-birthtime
  ctype = 6
  ACl$birthtime = ACls.newcol(neuron_name, ctype)


  #--morphology
  # morpho_class_vec= readRDS('/media/WDHDD/clustering/morphology/class_vec_morph.rds')
  # morpho_classification = readRDS('/media/WDHDD/clustering/morphology/ensemble_clustering/fcleaned_mMBEM8000_t33_min100.rds')
  # # load('/home/noob/WDHDD/clustering/morphology/ARGO/pvmorph_mMBEM8000.rData'); morpho_classification = all_classification[[80]]
  # morpho = rep(0,length(neuron_name))
  # for(n in 1:length(morpho)) if(morpho_class_vec[n]) morpho[n] = morpho_classification[morpho_class_vec[n]]
  morpho = rep(0,length(neuron_name))
  ACl$morpho = morpho


  #--neurotrans
  ctype = 3
  ACl$neurotrans = ACls.newcol(neuron_name, ctype)


  #--community
  anat_comm = readRDS('/media/WDHDD/clustering/community/community_classification.rds')
  ACl$community = anat_comm



  ##--network
  # combined.filename = paste("/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/data/combined_d",
  #                          d.used,"_t",tau,"_min",min.clust,".rData",sep="")
  # tempfc = readRDS(ntwk.filename)
  ntwk.filename =paste("/media/WDHDD/ckt/IVC/fcl_results/fcl_",file.descrip,"_d",d.used,"_g",bottom_index,"_g",top_index,"_mc",min.clust,".rData",sep="")
  load(ntwk.filename)
  ACl$network = fcl[[tau]]
  rm(avg_pairsum, cl_sorted, fcl, min_pairsum, num_lost)

  #--spad
  # spad = readRDS('/media/WDHDD/clustering/exp/ARGO/P05/IVC/recheck/spad/spad_clustering_tau95_k54.rds')
  spad = rep(0,length(neuron_name))
  ACl$spad = spad #spad$classification

  Aclass.filename = paste("/media/WDHDD/ckt/labels/ACl_data/Aclass_",
                           file.descrip,"_d",d.used,"_g",bottom_index,"_g",top_index,"_mc",min.clust,"_tau",tau,".rds",sep="")

  saveRDS(ACl,file=Aclass.filename)
}
##----------------------------------------------------------------------------##




################################################################################
##--joint prob two variables

joint.prob <- function(file.descrip, bottom_index, top_index, tau, min.clust, d.used, X.ind, Y.ind)
{
  ##--
  Aclass.filename = paste("/media/WDHDD/ckt/labels/ACl_data/Aclass_",
                           file.descrip,"_d",d.used,"_g",bottom_index,"_g",top_index,"_mc",min.clust,"_tau",tau,".rds",sep="")

  ACl = readRDS(Aclass.filename)

  gps <- ACl$network
  gps = gps[which(gps!=0)]
  gps = unname(table(gps))

  X = ACl[[X.ind]]
  Y = ACl[[Y.ind]]

  remove.zeros = TRUE

  ##--remove common zeros
  if(remove.zeros)
  {
    zero.vert = which(ACl$network==0)
    X=X[-zero.vert]
    Y=Y[-zero.vert]

    if(X.ind==4 || Y.ind==4)
    {
      morpho_classification = ACl$morpho[-zero.vert]
      zero.vert = which(morpho_classification==0)
      X=X[-zero.vert]
      Y=Y[-zero.vert]
    }
  }

  uX = sort(unique(X))
  uY = sort(unique(Y))

  ##--relabel network clusters
  if(X.ind==7)
  {
    k=0
    for(i in uX)
    {
      if(i!=0)
      {
          k=k+1
          X[which(X == uX[i])]   = paste('Ntwk',toString(k),sep='')
          uX[i] = paste('Ntwk',toString(k),sep='')
      }
    }
  }

  if(Y.ind==7)
  {
    k=0
    for(i in uY)
    {
      if(i!=0)
      {
          k=k+1
          Y[which(Y == uY[i])]   = paste('Ntwk',toString(k),sep='')
          uY[i] = paste('Ntwk',toString(k),sep='')
      }
    }
  }

  if(X.ind==4)
  {
    k=0
    for(i in uX)
    {
      if(i!=0)
      {
          k=k+1
          X[which(X == uX[i])]   = paste('MorPvec',toString(k),sep='')
          uX[i] = paste('MorPvec',toString(k),sep='')
      }
    }
  }

  if(Y.ind==4)
  {
    k=0
    for(i in uY)
    {
      if(i!=0)
      {
          k=k+1
          Y[which(Y == uY[i])]   = paste('MorPvec',toString(k),sep='')
          uY[i] = paste('MorPvec',toString(k),sep='')
      }
    }
  }
  #--end relabel

  ##--intersection counts
  Cxy = array( rep(0),dim=c(length(uX),length(uY)) )
  rownames(Cxy)=uX
  colnames(Cxy)=uY

  for(i in 1:nrow(Cxy))
  {
    for(j in 1:ncol(Cxy))
    {
        cmmn = intersect(which(X==uX[i]),which(Y==uY[j]))
        Cxy[uX[i],uY[j]] = length(cmmn)
    }
  }


  ##--joint prob.
  Pxy = Cxy/length(X)

  ##--marginal prob.
  Px = rowSums(Cxy)/length(X)
  Py = colSums(Cxy)/length(Y)

  ##--conditional prob.
  Py_x = Pxy; Py_x[,]=0
  Px_y = Py_x

  for(i in 1:nrow(Pxy))
  {
    Py_x[i,] = Pxy[i,]/Px[i] ##--or Py_x = Cxy[1,]/sum(Cxy[1,])
  }
  for(i in 1:ncol(Pxy))
  {
    Px_y[,i] = Pxy[,i]/Py[i] ##--or Py_x = Cxy[1,]/sum(Cxy[1,])
  }

  ##--MI
  # natstobits(mutinformation(ACl[-zero.vert,2], ACl[-zero.vert,4], method="emp"))
  # philentropy::MI(Px,Py, Pxy, unit = "log2")
  mi=0
  for(i in 1:nrow(Pxy))
  {
    for(j in 1:ncol(Pxy))
    {
        if(Pxy[i,j]!=0)
          mi = mi + ( Pxy[i,j]*(log2(Pxy[i,j]) - log2(Px[i]*Py[j])) )
    }
  }
  mi = unname(mi)

  Hx=0
  for(i in 1:length(Px))
  {
    if(Px[i]!=0)
      Hx = Hx - Px[i]*(log2(Px[i]))
  }
  Hx = unname(Hx)

  Hy=0
  for(i in 1:length(Py))
  {
    if(Py[i]!=0)
      Hy = Hy - Py[i]*(log2(Py[i]))
  }
  Hy = unname(Hy)

  Hx_y = Hx-mi
  Hy_x = Hy-mi
  nmi = mi/(min(c(Hx,Hy)))

  nzero = length(X)
  num.neurons = length(which(ACl$network!=0))
  khat = max(ACl$network)

  cat("\n\n\n\n\n\n#############################################################################")
  cat("\n-- Summary --")
  cat("\n\n   --------------------------------------------------------------------------")
  cat("\n   Minimum cluster size:",min.clust)
  cat("\n   Ensemble clustering threshold:",tau)
  cat("\n   Dimesion of spectral embedding:",d.used)
  cat("\n   --------------------------------------------------------------------------")
  cat("\n\n   Comparing X=",names(ACl)[X.ind],"versus Y=",names(ACl)[Y.ind])
  cat("\n   Mutual information:",mi ,"bits")
  cat("\n   Normalized mutual information [0-1]:",nmi )
  cat("\n   Information in",names(ACl)[X.ind],"not explained by",names(ACl)[Y.ind],"H(X|Y):",Hx_y ,"bits")
  cat("\n   Information in",names(ACl)[Y.ind],"not explained by",names(ACl)[X.ind],"H(Y|X):",Hy_x ,"bits")
  cat("\n#############################################################################\n")

  out.ls = list(Px=Px, Py=Py, Pxy=Pxy, Px_y=Px_y, Py_x=Py_x,
                Hx=Hx, Hy=Hy, Hx_y=Hx_y, Hy_x=Hy_x, nmi=nmi, mi=mi,
                nzero=nzero, khat=khat, num.neurons=num.neurons, gps=gps)
  return(out.ls)
}
##----------------------------------------------------------------------------##







################################################################################
##--plot class dependency

plot.dep <- function(m, param.text1=NULL, param.text2=NULL, meas.text=NULL, lay=layout_on_grid)
{
  require(igraph)
  g<-graph.adjacency(m, weighted = TRUE, mode='directed')
  g<-delete.edges(g, which(is.na(E(g)$weight)) )
  g<-delete.vertices(g, degree(g)==0)

  V(g)$size <- 25
  V(g)$color <- "gold"
  V(g)$frame.color <- "gray"
  V(g)$label.color <- "black"
  V(g)$label.cex <- 1.2
  E(g)$label.cex <- 1.2
  E(g)$label=round(E(g)$weight,2)
  E(g)$arrow.size  <- (0.6+E(g)$weight)
  E(g)$width <- (0.3+E(g)$weight)*1.3
  #layout = layout_on_grid

  # cols.blu<-c(brewer.pal(8, "Blues")[6:8],"#000000")
  cols.blu<-c("#A8A8A8", "#8FD744", "#35B779", "#21908C", "#44d7cc", "#1d7fc3", "#313c8e", "#443A83", "#440154", "#a52156", "#dd3636", "#000000")
  x=colorRampPalette(cols.blu, space = "Lab", bias=1)
  lvls = as.numeric(cut(E(g)$weight,breaks = 20))
  E(g)$color <- x(20)[lvls] #sapply(lvls, function(x){cols.blu[x]})
  E(g)$label.color <- x(20)[lvls]

  plot(g, layout=lay, edge.curved=0.2)
  title(main=meas.text, cex.main=1, font.main=1)
  mtext(param.text1, cex=1, side = 1, line = 0.5)
  mtext(param.text2, cex=1, side = 1, line = 1.5)

}

##----------------------------------------------------------------------------##












################################################################################
##--user set
file.descrip = "P05_MBHAC_4to100"
bottom_index = 1
top_index = 200
tau = 180
min.clust  = 100
d.used     = 11
new.ACl    = TRUE
save.image = FALSE
calc.bias  = FALSE
calc.neurotrans.labels  = TRUE
calc.driver.labels      = TRUE
calc.birthtime.labels   = TRUE
calc.community.labels   = TRUE
##-- pick X and Y
##-- pick X as network for labels!!
## 1. neuron_name
## 2. driver
## 3. birthtime
## 4. morpho
## 5. neurotrans
## 6. community
## 7. network
## 8. spad
X.vec = c(7)
Y.vec = c(2,3,5,6,7)
##----------------------------------------------------------------------------##





################################################################################
##--main code
propH.mat = matrix(rep(NA), nrow=8, ncol=8)
nmi.mat = matrix(rep(NA), nrow=8, ncol=8)
for(X.ind in X.vec)
{
    for(Y.ind in Y.vec)
    {

      if(new.ACl)
        ACl.create(file.descrip=file.descrip, bottom_index=bottom_index, top_index=top_index, tau=tau, min.clust=min.clust, d.used=d.used)

      if(X.ind==Y.ind)
      {
        propH.mat[X.ind,Y.ind] = 0
      } else {
        if(is.na(propH.mat[X.ind,Y.ind]) || is.na(propH.mat[Y.ind,X.ind]) )
        {
          prob = joint.prob(file.descrip=file.descrip, bottom_index=bottom_index, top_index=top_index, tau=tau, min.clust=min.clust, d.used=d.used, X.ind=X.ind, Y.ind=Y.ind)

          propH.mat[X.ind,Y.ind] = prob$mi/prob$Hy #1-prob$Hx_y/prob$Hx
          propH.mat[Y.ind,X.ind] = prob$mi/prob$Hx #1-prob$Hy_x/prob$Hy

          nmi.mat[X.ind,Y.ind] = prob$nmi
          nmi.mat[X.ind,Y.ind] = nmi.mat[Y.ind,X.ind]

          param.text2 = paste("no. of classified vertices=",prob$num.neurons,",  estimated khat=",prob$khat,sep="")
        }
      }



      ################################################################################
      if(X.ind==7)
      {
        if(exists("label_frame")==FALSE)
        label_frame = data.frame(c(1:length(prob$Px))) #names(prob$Px))
        names(label_frame)[1] = ""
        label_frame$gps = prob$gps

        ##----------------------------------------------------------------------------##
        if(calc.neurotrans.labels)
        {
          if(Y.ind==5)
          {
            source('/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/neurotrans_labels.r')
            nt_labels = neurotrans_labels(prob$Py_x,prob$Py,prob$Px)
            label_frame$neurotrans_label = unlist(nt_labels)
          }
        }
        ##----------------------------------------------------------------------------##
        if(calc.driver.labels)
        {
          if(Y.ind==2)
          {
            source('/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/driver_labels.r')
            dr_labels = driver_labels(prob$Py_x,prob$Py,prob$Px)
            label_frame$driver_label = unlist(dr_labels)
          }
        }
        ##----------------------------------------------------------------------------##
        if(calc.birthtime.labels)
        {
          if(Y.ind==3)
          {
            source('/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/birthtime_labels.r')
            bt_labels = birthtime_labels(prob$Py_x,prob$Py,prob$Px)
            label_frame$birthtime_label = unlist(bt_labels)
          }
        }
        ##----------------------------------------------------------------------------##
        if(calc.community.labels)
        {
          if(Y.ind==6)
          {
            source('/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/community_labels.r')
            comm_labels = community_labels(prob$Py_x,prob$Py,prob$Px)
            label_frame$community_label = unlist(comm_labels)
          }
        }
        ##----------------------------------------------------------------------------##
        csv.filename = paste("/media/WDHDD/ckt/labels/labels_csv/labels_",
                           file.descrip,"_d",d.used,"_g",bottom_index,"_g",top_index,"_mc",min.clust,"_tau",tau,".csv",sep="")
        write.csv(label_frame, csv.filename, row.names=FALSE)
      }
      ################################################################################

    }#Y
}#X

propH.mat = propH.mat[-1,];propH.mat = propH.mat[,-1]
rownames(propH.mat)=c('driver','birthtime','morphology','neurotrans','community','network','spad')
colnames(propH.mat)=c('driver','birthtime','morphology','neurotrans','community','network','spad')

param.text1 = paste("embdedding d=",d.used,",  minimum cluster size=",min.clust,",  ensemble threshold=",tau,sep="")
meas.text = "Proportion of Information Explained"

# if(save.image) {library(Cairo); CairoPS(file="Rplot", width = 7, height = 7, font='serif', family='serif')}
# if(save.image) {setEPS(); postscript(file="Rplot.eps", family='FreeSerif', fonts='serif')}
if(save.image) cairo_ps('Rplot.eps')
plot.dep(propH.mat, param.text1, param.text2, lay=layout_in_circle)
if(save.image) dev.off()
##----------------------------------------------------------------------------##
