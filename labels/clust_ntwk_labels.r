## ----------------------------------------------------------------------------
## Preamble
library(igraph)
library(qgraph)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
theme_update(plot.title = element_text(hjust = 0.5, size=16))
theme_update(axis.text=element_text(size=15),
             axis.title=element_text(size=30) )# ,face="bold")

cat("\014")     # clear screen
rm(list=ls())   # clean workspace



## ----------------------------------------------------------------------------
P   = readRDS('/home/noob/WDHDD/clustering/exp/ARGO/P05/IVC/recheck/PMAT_95.rds')
fc  = readRDS('/home/noob/WDHDD/clustering/exp/ARGO/P05/IVC/recheck/cl_95.rds')
gps = table(fc[-which(fc==0)])
k_hat = max(fc) #length(gps)


my_cols <- c("#8c8c8c", "#9000E4", "#ffe119", "#f58231", "#00b159", "#e6beff", "#42d4f4", "#e6194B",
             "#F3FCB9", "#4686CE", "#ACEB8D", "#F9CDAD", "#99B898", "#2F9599", "#C06C84", "#f29959" )

plot.ckt=FALSE
plot.bars=TRUE
save.plot=FALSE
compute.bet.deg=TRUE
compute.AD=FALSE



## ----------------------------------------------------------------------------
g<-graph.adjacency(P, weighted = TRUE, mode='directed', diag=TRUE)

P.uw = P
P.uw[which(P.uw>0)]=1
g.uw<-graph.adjacency(P.uw, mode='directed', diag=FALSE)
set.seed(60) #5, 50




## ----------------------------------------------------------------------------
##--vertex color
cols.or<-brewer.pal(9, "Oranges")[1:7]
cols.gr<-brewer.pal(9, "Greens")


V(g)$size         <- (gps+120)/18
V(g)$frame.color  <- "#594F4F" #"#D8D8D8"
V(g)$label.color  <- "#594F4F"
V(g)$label.cex    <- (gps+300)/310
E(g)$arrow.size   <- (0.6+E(g)$weight)*1.7
E(g)$width        <- 0.4+E(g)$weight*8.2
#E(g)$label.cex    <- 1



##--edge color
cols.blu<-c(brewer.pal(9, "Blues")[4:9],"#000000")
#cols.blu<-c('#A8A7A7','#363636',"#000000")
x=colorRampPalette(cols.blu, space = "Lab", bias=1.7)
lvls = as.numeric(cut(E(g)$weight,breaks = 8))
E(g)$color <- x(8)[lvls] #sapply(lvls, function(x){cols.blu[x]})



if(plot.ckt)
{
   l <- layout_with_lgl(g)
   l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
   par(mfrow=c(1,1), mar=c(10,12,12,10))
   tkplot(g, rescale=F, layout=l*2, vertex.label.degree=pi/2, vertex.label.dist=0.1, loop.angle=pi/2)
   # qgraph(P,layout=l,diag=TRUE,directed=TRUE,vsize=gps/90+1.2, edge.color="#355C7D")
}


if(compute.bet.deg)
{
  x.names = c('clust_in', 'clust_out', 'clust_total', 'clust_bet') #, "clust_particip")

  x <- list()
  x[[1]] = strength(g, mode = "in", loops = TRUE) # ,normalized = FALSE)
  x[[2]] = strength(g, mode = "out", loops = TRUE)
  x[[3]] = strength(g, mode = "total", loops = TRUE)
  x[[4]] = betweenness(g.uw, directed = TRUE)#, weights = E(g)$weight)
  # x[[5]] = NetworkToolbox::participation(P)$overall
  # lvls = as.numeric(cut(deg.tot,breaks = 9))

  x[[1]] = x[[1]]*gps; x[[2]] = x[[2]]*gps; x[[3]] = x[[3]]*gps

  df.str = data.frame(dir=rep(c("In", "Out"), each=54),
                class_num=c(1:54),
                class_str=c(x[[1]],x[[2]]))

  df.btw = data.frame(class_num=c(1:54),
                class_btw=c(x[[4]]))

  df.combined = data.frame(class_num=c(1:54),
                    class_str=c(x[[3]]),
                    class_btw=c(x[[4]]),
                    hub=rep("",54))
  df.combined$hub[c(32,50,40,46,37,39)]="Hub"

browser()
  if(plot.bars)
  {
    # par(mar=c(5.1, 4.1, 4.1, 5.1))
    c1 <- rgb(0,158,160, max = 255, alpha = 255, names = "bluish.green")
    c2 <- rgb(204,37,41, max = 255, alpha = 222, names = "light.red")
    my_cols <- c("#8c8c8c", "#9000E4", "#ffe119", "#f58231", "#00b159", "#e6beff", "#42d4f4",  "#e6194B",
                  "#F3FCB9", "#4686CE", "#ACEB8D", "#F9CDAD", "#99B898", "#2F9599", "#C06C84" )
    # plt <- barplot(rbind(x[[1]],x[[2]]), beside=FALSE,

    if(save.plot) cairo_ps('./scatter_gt.eps',width=6.5,height=5.5)
    plt = ggplot(df.combined, aes(x=class_str, y=class_btw, fill=factor(hub))) +
          geom_point(alpha=0.4, size=1.4) +
          geom_label_repel(aes(label = class_num), size = 4.4, max.overlaps=20) + labs(fill = " ") +
          scale_fill_manual(values=c("#ffffff", "#f8766d"), breaks="Hub") +
          guides(fill = guide_legend(override.aes = aes(label = ""))) +
          xlab("Strength (Weighted Degree)") +
          ylab("Betweenness Centrality") +
          theme_light() +
          theme(axis.text=element_text(size=15), axis.title=element_text(size=15),
              legend.text=element_text(size=15), legend.position=c(0.9,0.1))
    print(plt)
    if(save.plot) dev.off()

    if(save.plot) cairo_ps('./clust_str.eps',width=5,height=3)
        plt = ggplot(data=df.str, aes(x=class_num, y=class_str, fill=dir)) +  geom_bar(stat="identity") +
        scale_fill_brewer(palette="Paired")+
        theme_minimal() +
        xlab("Class No.") +
        ylab("Strength (Weighted Degree)") +
        labs(fill = " ") +
        theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
              legend.text=element_text(size=12), legend.position=c(0.1,0.9)) +
        annotate("text", x = c(32,50,40,46,37,39,42,52), y=x[[3]][c(32,50,40,46,37,39,42,52)]+20,
              label = c(rep("H",6),rep("M",2)), color="#d55e00", size=(3.2), fontface="bold") +
        scale_x_continuous(breaks = c(10,20,30,40,50))
        print(plt)
    if(save.plot) dev.off()

    if(save.plot) cairo_ps('./clust_bet.eps',width=5,height=3)
        plt = ggplot(data=df.btw, aes(x=class_num, y=class_btw)) +  geom_bar(stat="identity", fill="mediumseagreen") +
        # scale_fill_manual(values="#009e73")+
        theme_minimal() +
        xlab("Class No.") +
        ylab("Betweenness Centrality") +
        labs(fill = " ") +
        theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
              legend.text=element_text(size=12), legend.position=c(0.1,0.9)) +
        annotate("text", x = c(32,50,40,46,37,39,42,52), y=x[[4]][c(32,50,40,46,37,39,42,52)]+1.5,
              label = c(rep("H",6),rep("M",2)), color="#d55e00", size=(3.2), fontface="bold") +
        scale_x_continuous(breaks = c(10,20,30,40,50))
        print(plt)
    if(save.plot) dev.off()
  }


browser()

  source('/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/gt_labels.r')
  y <- NULL
  for(i in 1:length(x))
  {
    fn.output = gt_labels(x.names[i],x[[i]])
    y = cbind(y,unlist(fn.output))
  }
  write.csv(y,'gt3_labels.csv')
}


if(compute.AD)
{
  load('/home/noob/WDHDD/clustering/exp/ARGO/P05/IVC/recheck/absorption/absorption_v6_sum10e4_95.rData')
  x = cbind(avg_in_absorption, avg_out_absorption, avg_in_driftiness, avg_out_driftiness)
  # load('/media/WDHDD/clustering/exp/ARGO/P05/random_walk/absorption/avg/absorption_v5.rData')
  # x = cbind(x,avg_in_driftiness,avg_out_driftiness)
  source('/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/gt_labels.r')

  y <- NULL
  for(i in 1:ncol(x))
  {
    fn.output = gt_labels(colnames(x)[i],x[,i],n=13)
    y = cbind(y,unlist(fn.output))
  }
  avg_abs_label = rep(NA,nrow(y))
  avg_drift_label = rep(NA,nrow(y))
  for(i in 1:nrow(y))
  {
    temp = paste(y[i,1:2], collapse=", ")
    temp = gsub(" , ", "", temp)
    temp = gsub(",  ", "", temp)
    avg_abs_label[i] = temp

    temp = paste(y[i,3:4], collapse=", ")
    temp = gsub(" , ", "", temp)
    temp = gsub(",  ", "", temp)
    avg_drift_label[i] = temp
  }
  abs_drift_labels = cbind(avg_abs_label, avg_drift_label)
  colnames(abs_drift_labels) = c("Avg Absorption", "Avg Driftiness")
  write.csv(abs_drift_labels,'abs_drift_labels.csv')
}



# if(plot.bars)
# {
#   dev.off()
#   if(save.plot) pdf('./gps41.pdf',width=11,height=5.5) #cairo_ps("Rplot.eps", width=7, height=7)
#
#     deg.in = degree(g, mode = "in", loops = TRUE, normalized = FALSE)
#     deg.out = degree(g, mode = "out", loops = TRUE, normalized = FALSE)
#
#     plt <- barplot(rbind(gps), beside=FALSE,
#     # main = "Participation Coefficient",
#     xlab = "Cluster No.",
#     col = my_cols[c(13)],
#     # ylab = "Day",
#     names.arg = 1:41)
#
#     legend("topleft",
#             c("number of vertices"),
#             fill = my_cols[c(13)]
#           )
#
#   if(save.plot) dev.off()
# }
