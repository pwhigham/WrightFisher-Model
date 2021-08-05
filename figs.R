#
# figs.R
#
# Figures for PNAS paper
#############################################

library("Hmisc")  # error bars
library("igraph")

#############################################
#   NETWORKS FIGURE
##############################################
plot.network <- function(graph="networks/sfp1_64/1.txt",
                         vertex.size=3.5,
                         layout=layout_with_fr,...)
{
  g = as.matrix(read.table(file=graph))
  diag(g) <- 0
  op <- par(mar=c(0,0,0,0))
  g <- graph_from_adjacency_matrix(g,mode="undirected",weighted=TRUE)
  plot(g,layout=layout,vertex.size=vertex.size,vertex.label=NA,
       vertex.color="black",...)
  print(paste("Mean degree (diag=0):",mean(degree(g)),"N=",length(V(g))))
  par(op)
  g
}

plot.lattice <- function(graph="networks/lat200/1.txt",...)
{
  g = as.matrix(read.table(file=graph))
  diag(g) <- 0
  par(mar=c(2,2,2,2))
  g <- graph_from_adjacency_matrix(g,mode="undirected",weighted=TRUE)
  ll = layout_on_grid(g,width=10,height=20)
  plot(g,layout=ll,vertex.size=3.5,vertex.label=NA,
       vertex.color="black",...)
  print(paste("Mean degree (diag=0):",mean(degree(g)),"N=",length(V(g))))
  
}


mdelta.fig <- function(mv.file="results/pan200.delta.csv",
                       graph="networks/pan200/1.txt",use.ylim=TRUE,
                       fix.ylim=c(0,0.005),
                       pan.delta="results/pan200.delta.csv.sav", # always exists
                       text.gamma.xy=c(0.073,0.0045),
                       text.gval.xy=c(0.18,0.0045)
                       )
{
  saved.mv <- paste(mv.file,".sav",sep="")
  if (file.exists(saved.mv))
  {
    mv <- read.csv(saved.mv,header=T)
  }
  else
  {
    mv <- delta.p(mv.file)
    write.csv(mv,saved.mv)
  }
  main.str <- paste("N = ",unique(mv$N),
                    " s = ",unique(mv$s))
  
  model.mdelta <- mean.delta.network(N=mv$N,s=mv$s,
                                     p=mv$p,
                                     graph=graph)
  if (use.ylim)
  {
    ylimit=fix.ylim
  }
  else
  {
    miny = min(model.mdelta,mv$mx)
    maxy = max(model.mdelta,mv$mx)
    ylimit <- c(miny,maxy)
  }
  
  plot(mv$p,mv$mx,pch=16,xlab="",ylab="",
       cex=1.5,main="",ylim=ylimit)
  lines(mv$p,mv$mx,lwd=2)
  
  g <- as.matrix(read.table(file=graph)) 
  gammas <- gamma(g)
  gmx = round(gammas[1],2)
  gm <- bquote(gamma[M])
  text(text.gamma.xy[1],text.gamma.xy[2],gm,cex=1.5)
  text(text.gval.xy[1],text.gval.xy[2],paste("=",gmx,sep=""),cex=1.5)
  
  text.cex=0.8
  mtext("p",side=1,line=0.7,cex=text.cex)
  mdeltastr = (bquote("M" [delta][p]))
  mtext(mdeltastr,side=2,line=2,cex=text.cex)
  errbar(x=mv$p,y=mv$mx,mv$mx+mv$cix, mv$mx-mv$cix, add=T, pch=1, cap=.02)
  
  points(mv$p,y=model.mdelta,pch=16,col=2,cex=1)
  lines(mv$p,y=model.mdelta,lty=1,col=2,lwd=2)
  
  # add panmictic for comparison
  pan.data <- read.csv(pan.delta)
  lines(pan.data$p,pan.data$mx,lty=1,col="grey",lwd=2)
  points(pan.data$p,pan.data$mx,pch=16,col="grey",cex=1)
  
  mv
}

vdelta.fig <- function(mv.file="results/pan200.delta.csv",
                       graph="networks/pan200/1.txt",
                       use.ylim=TRUE,
                       fix.ylim=c(0,0.06),
                       pan.delta="results/pan200.delta.csv.sav",
                       text.gamma.xy=c(0.071,0.054),
                       text.gval.xy=c(0.18,0.054)) #  
                       
{
  saved.mv <- paste(mv.file,".sav",sep="")
  if (file.exists(saved.mv))
  {
    mv <- read.csv(saved.mv,header=T)
  }
  else
  {
    mv <- delta.p(mv.file)
    write.csv(mv,saved.mv)
  }
  
  main.str <- paste("N = ",unique(mv$N),
                    " s = ",unique(mv$s))
  
  print(paste(unique(mv$N),unique(mv$s)))
  model.vdelta <- var.delta.network(N=mv$N,s=mv$s,
                                     p=mv$p,
                                     graph=graph)
  
  if (use.ylim)
  {
    ylimit=fix.ylim
  }
  else
  {
    miny = min(model.vdelta,mv$vx)
    maxy = max(model.vdelta,mv$vx)
    ylimit <- c(miny,maxy)
  }
  
  plot(mv$p,mv$vx,pch=16,xlab="",ylab="",
       cex=1.5,main="",ylim=ylimit)
  lines(mv$p,mv$vx,lwd=2)
  
  g <- as.matrix(read.table(file=graph)) 
  gammas <- gamma(g)
  
  gmx = round(gammas[2],2)
  gm <- bquote(gamma[V])
  text(text.gamma.xy[1],text.gamma.xy[2],gm,cex=1.5)
  text(text.gval.xy[1],text.gval.xy[2],paste("=",gmx,sep=""),cex=1.5)
  
  text.cex=0.8
  mtext("p",side=1,line=0.7,cex=text.cex) # was line=-1.5
  vdeltastr = (bquote("V" [delta][p]))
#  vdeltastr = expression({{V^{G}}[delta][p]})
  mtext(vdeltastr,side=2,line=2,cex=text.cex)
#  errbar(x=mv$p,y=mv$vx,mv$vx+mv$civarhigh, mv$mx-mv$civarlow, add=T, pch=1, cap=.02)
  
  points(mv$p,y=model.vdelta,pch=16,col=2,cex=1)
  lines(mv$p,y=model.vdelta,lty=1,col=2,lwd=2)
  
  
  pan.data <- read.csv(pan.delta)
  lines(pan.data$p,pan.data$vx,lty=1,col="grey",lwd=2)
  points(pan.data$p,pan.data$vx,pch=16,col="grey",cex=1)
  
  
  mv
}
fix.fig <- function(fix.file="results/pan200.fix.csv",
                    graph="networks/pan200/1.txt",
                    pan.graph="networks/pan200/1.txt",
                    pan.fix="results/pan200.fix.csv",
                    gt = c(0.49,0.1), # hack for text placement
                    et = c(0.55,0.1))
{
  saved.fix <- paste(fix.file,".sav",sep="")
  if (file.exists(saved.fix))
  {
    fix <- read.csv(saved.fix,header=T)
  }
  else
  {
    fix <- load.fixation(fix.file=fix.file)
    write.csv(fix,saved.fix)
  }
  saved.fix <- paste(pan.fix,".sav",sep="")
  
  if (file.exists(saved.fix))
  {
    pan.fix <- read.csv(saved.fix,header=T)
  }
  else
  {
    pan.fix <- load.fixation(fix.file=fix.file)
  }
  
  main.str <- paste("N = ",unique(fix$N),
                    " s = ",unique(fix$s))
  print(paste(unique(fix$N),unique(fix$s)))
  plot(fix$p,fix$fix.m,pch=16,xlab="",ylab="",
       cex=1.5,main="",ylim=c(0,1),type='l')
  lines(fix$p,fix$fix.m,lwd=2)
  #mtext(main.str,side=3,line=0.2)
  text.cex=0.8
  mtext("p",side=1,line=0.7,cex=text.cex)
  mtext("u(p,G)",side=2,line=2.2,cex=text.cex)
  
  errbar(x=fix$p,y=fix$fix.m,
         fix$fix.m+fix$ci, 
         fix$fix.m-fix$ci, add=T, pch=16, cap=.015,lwd=1.5)
  
  g <- as.matrix(read.table(file=graph)) 
  gammas <- gamma(g)
  print(gammas)
  model.fix <- gamma.fix(N=unique(fix$N),
            s=unique(fix$s),
            p=fix$p,gammas)
  points(fix$p,y=model.fix,pch=16,col=2,cex=1)
  lines(fix$p,y=model.fix,lty=1,col=2,lwd=2)
  
  lines(pan.fix$p,pan.fix$fix.m,lty=1,lwd=2,col='grey')
  points(pan.fix$p,pan.fix$fix.m,pch=16,cex=1,col='grey')

  g <- as.matrix(read.table(file=graph)) 
  gammas <- gamma(g)
  gmx = round(gammas[1]/gammas[2],2)
  gm <- bquote(gamma[G])
  
  text(x=gt[1],y=gt[2],gm,cex=1.5)
  text(x=et[1],y=et[2],paste("=",gmx,sep=""),cex=1.5)
  
  
  fix
}

#
# Plot the layout plots for different spaces, Mx,Vx,Fix,Picture
#
plot.PNAS.fig1 <- function()
{
  # Since all same size can just use normal controls
  
  jpeg("figures/fig1.jpg", units="cm",width = 11.4, height = 11.4, 
       res=600,
       pointsize=6,quality=95)
  
  
  op <- par(mfrow=c(5,4),mar=c(2,3,0.5,0.5))
  
  lab.x=0.1
  lab.y=0.9
  lab.size=2.0
  
  #pan
  graph <- "networks/pan200/1.txt"
  delta.file <- "results/pan200.delta.csv"
  fix.file <- "results/pan200.fix.csv"
  mdelta.fig(mv.file=delta.file,graph=graph)
  vdelta.fig(mv.file=delta.file,graph=graph)
  fix.fig(fix.file=fix.file,graph=graph)
  text(lab.x,lab.y,labels="A",cex=lab.size)
  plot.network(graph=graph,layout=layout_with_fr,edge.width=0.000001)
  
  # ring
  graph <- "networks/ring200/1.txt"
  delta.file <- "results/ring200.delta.csv"
  fix.file <- "results/ring200.fix.csv"
  mdelta.fig(mv.file=delta.file,graph=graph)
  vdelta.fig(mv.file=delta.file,graph=graph)
  fix.fig(fix.file=fix.file,graph=graph)
  text(lab.x,lab.y,labels="B",cex=lab.size)
  plot.network(graph=graph,layout=layout_in_circle)
  
  # degree
  graph <- "networks/rnd200/1.txt"
  delta.file <- "results/rnd200.delta.csv"
  fix.file <- "results/rnd200.fix.csv"
  mdelta.fig(mv.file=delta.file,graph=graph)
  vdelta.fig(mv.file=delta.file,graph=graph)
  fix.fig(fix.file=fix.file,graph=graph)
  text(lab.x,lab.y,labels="C",cex=lab.size)
  plot.network(graph=graph,layout=layout_with_fr)
  
  # sf
  graph <- "networks/sf200/1.txt"
  delta.file <- "results/sf200.delta.csv"
  fix.file <- "results/sf200.fix.csv"
  mdelta.fig(mv.file=delta.file,graph=graph)
  vdelta.fig(mv.file=delta.file,graph=graph)
  fix.fig(fix.file=fix.file,graph=graph)
  text(lab.x,lab.y,labels="D",cex=lab.size)
  plot.network(graph=graph,layout=layout_with_fr)
  
  # star
  par(mar=c(2,3,0.3,0.5))
  
  graph <- "networks/star64/1.txt"
  delta.file <- "results/star64.delta.csv"
  fix.file <- "results/star64.fix.csv"
  mdelta.fig(mv.file=delta.file,graph=graph)
  vdelta.fig(mv.file=delta.file,graph=graph)
  fix.fig(fix.file=fix.file,graph=graph)
  text(lab.x,lab.y,labels="E",cex=lab.size)
  plot.network(graph=graph,layout=layout_with_fr)
  
 dev.off()
  
  par(op)      
}
PNAS.fix.fig <- function()
{
  jpeg("figures/figfix.jpg", units="cm",width = 6.8, height = 10, 
       res=600,
       pointsize=6,quality=95)
  
  
  op <- par(mfrow=c(3,2),mar=c(2,3,0.5,0.5))
  
  lab.x=0.1
  lab.y=0.9
  lab.size=2.0
  
  graph <- "networks/pan200/1.txt"
  fix.file <- "results/pan200.fix.csv"
  fix.fig(fix.file=fix.file,graph=graph,et = c(0.57,0.1))
  text(lab.x,lab.y,labels="A",cex=lab.size)
  
  graph <- "networks/ring200/1.txt"
  fix.file <- "results/ring200.fix.csv"
  
  fix.fig(fix.file=fix.file,graph=graph,et = c(0.57,0.1))
  text(lab.x,lab.y,labels="B",cex=lab.size)

  graph <- "networks/lat200/1.txt"
  fix.file <- "results/lat200.fix.csv"
  
  fix.fig(fix.file=fix.file,graph=graph,et = c(0.57,0.1))
  text(lab.x,lab.y,labels="C",cex=lab.size)
  
  graph <- "networks/rnd200/1.txt"
  fix.file <- "results/rnd200.fix.csv"
  
  fix.fig(fix.file=fix.file,graph=graph,et = c(0.62,0.1))
  text(lab.x,lab.y,labels="D",cex=lab.size)
  
  graph <- "networks/sf200/1.txt"
  fix.file <- "results/sf200.fix.csv"
  fix.fig(fix.file=fix.file,graph=graph,et = c(0.62,0.1))
  text(lab.x,lab.y,labels="E",cex=lab.size)

  graph <- "networks/star200/1.txt"
  fix.file <- "results/star200.fix.csv"
  fix.fig(fix.file=fix.file,graph=graph,et = c(0.62,0.1))
  text(lab.x,lab.y,labels="F",cex=lab.size)
  
  dev.off()
  par(op)
  
}
# Plot of ttf for each space 
fig.PNAS.ttf <- function()
{
  jpeg("figures/figttf.jpg", units="cm",width = 6.8, height = 6.8, 
       res=600,
       pointsize=6,quality=95)
  
  op <- par(mar=c(2,3,0.5,0.5))
  pan.gens <- load.fixation("results/pan200.fix.csv")
  ring.gens <- load.fixation("results/ring200.fix.csv")
  lat8.gens <- load.fixation("results/lat200.fix.csv")
  deg.gens <- load.fixation("results/rnd200.fix.csv")
  sf.gens <- load.fixation("results/sf200.fix.csv")
  star.gens <- load.fixation("results/star200.fix.csv")

  plot(pan.gens$p,pan.gens$gens,type='l',ylim=c(5,5000),
       xlab="",ylab="",log="y",col=1)
  errbar(x=pan.gens$p,y=pan.gens$gens,col=1,
         pan.gens$gens+pan.gens$gens.ci, 
         pan.gens$gens-pan.gens$gens.ci, add=T, pch=1, cap=.02,errbar.col=1)
  text.cex=1
  mtext("p",side=1,line=0.5,cex=text.cex)
  
  mtext("Time to Fixation log(Gens) of allele B",side=2,line=1.9)
  
  lines(ring.gens$p,ring.gens$gens,col=2)
  errbar(x=ring.gens$p,y=ring.gens$gens,col=2,
         ring.gens$gens+ring.gens$gens.ci, 
         ring.gens$gens-ring.gens$gens.ci, add=T, pch=1, cap=.02,errbar.col=2)
  
  lines(lat8.gens$p,lat8.gens$gens,col=3)
  errbar(x=lat8.gens$p,y=lat8.gens$gens,col=3,
         lat8.gens$gens+lat8.gens$gens.ci, 
         lat8.gens$gens-lat8.gens$gens.ci, add=T, pch=1, cap=.02,errbar.col=3)
  
  lines(deg.gens$p,deg.gens$gens,col=4)
  errbar(x=ring.gens$p,y=deg.gens$gens,col=4,
         deg.gens$gens+deg.gens$gens.ci, 
         deg.gens$gens-deg.gens$gens.ci, add=T, pch=1, cap=.02,errbar.col=4)
  
  lines(sf.gens$p,sf.gens$gens,col=5)
  errbar(x=ring.gens$p,y=sf.gens$gens,col=5,
         sf.gens$gens+sf.gens$gens.ci, 
         sf.gens$gens-sf.gens$gens.ci, add=T, pch=1, cap=.02,errbar.col=5)
  
  lines(star.gens$p,star.gens$gens,col=6)
  errbar(x=ring.gens$p,y=star.gens$gens,col=6,
         star.gens$gens+star.gens$gens.ci, 
         star.gens$gens-star.gens$gens.ci, add=T, pch=1, cap=.02,errbar.col=6)
  
  text(0.73,2000,"Ring",col=2)
  text(0.55,100,"Panmictic",col=1)
  text(0.2,23,"Star",col=6)
  text(0.22,725,"Scale Free (p=1.3)",col=5)
  text(0.3,267,"N(7,2)",col=4)
  text(0.5,319,"Lattice",col=3)
  
  par(op)
  dev.off()
}

PNAS.graphs.fig <- function()
{
  jpeg("figures/figgraphs.jpg", units="cm",width = 6.8, height = 10, 
       res=600,
       pointsize=6,quality=95)
  
  lab.x=-0.9
  lab.y=0.9
  lab.size=2.5
  
  op <- par(mfrow=c(3,2))
  #,mar=c(2,3,0.3,0.5))
  
  graph <- "networks/pan200/1.txt"
  plot.network(graph=graph,layout=layout_with_fr,edge.width=0.000001)
  text(lab.x,lab.y,labels="A",cex=lab.size)
  
  lines(x=c(-1.1,1.1,1.1,-1.1,-1.1),y=c(-1.05,-1.05,1.05,1.05,-1.05))
  
  # ring
  graph <- "networks/ring200/1.txt"
  plot.network(graph=graph,layout=layout_in_circle)
  text(lab.x,lab.y,labels="B",cex=lab.size)
  
  lines(x=c(-1.1,1.1,1.1,-1.1,-1.1),y=c(-1.05,-1.05,1.05,1.05,-1.05))
  
  
  graph <- "networks/lat200/1.txt"
  plot.lattice()
  par(mar=c(0,0,0,0))
  #plot.network(graph=graph,layout=layout_on_grid)
  text(-1.15,1.146,labels="C",cex=lab.size)
  
  lines(x=c(-1.37,1.37,1.37,-1.37,-1.37),y=c(-1.3,-1.3,1.3,1.3,-1.3))
#  lines(x=c(-1.1,1.1,1.1,-1.1,-1.1),y=c(-1.05,-1.05,1.05,1.05,-1.05))
  
  
  # degree
  graph <- "networks/rnd200/1.txt"
  plot.network(graph=graph,layout=layout_with_fr)
  text(lab.x,lab.y,labels="D",cex=lab.size)
  lines(x=c(-1.1,1.1,1.1,-1.1,-1.1),y=c(-1.05,-1.05,1.05,1.05,-1.05))
  # sf
  graph <- "networks/sf200/1.txt"
  plot.network(graph=graph,layout=layout_with_fr)
  text(lab.x,lab.y,labels="E",cex=lab.size)
  lines(x=c(-1.1,1.1,1.1,-1.1,-1.1),y=c(-1.05,-1.05,1.05,1.05,-1.05))
  # star
  
  graph <- "networks/star200/1.txt"
  plot.network(graph=graph,layout=layout_with_fr)
  text(lab.x,lab.y,labels="F",cex=lab.size)
  lines(x=c(-1.1,1.1,1.1,-1.1,-1.1),y=c(-1.05,-1.05,1.05,1.05,-1.05))
  
  
  dev.off()
  par(op)
  
}
PNAS.Mv.fig <- function()
{
  # Since all same size can just use normal controls
  
  jpeg("figures/figMv.jpg", units="cm",width = 6.8, height = 10, 
       res=600,
       pointsize=6,quality=95)
  
  
  op <- par(mfrow=c(3,2),mar=c(2,3.28,0.5,0.5))
  
  lab.x=0.85
  lab.y=0.0045
  lab.size=2.1
  
  #pan
  graph <- "networks/pan200/1.txt"
  delta.file <- "results/pan200.delta.csv"
  mdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.092,0.0045))
  text(lab.x,lab.y,labels="A",cex=lab.size)

  # ring
  graph <- "networks/ring200/1.txt"
  delta.file <- "results/ring200.delta.csv"
  mdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.092,0.0045),
             text.gval.xy=c(0.23,0.0045))
  text(lab.x,lab.y,labels="B",cex=lab.size)

  graph <- "networks/lat200/1.txt"
  delta.file <- "results/lat200.delta.csv"
  mdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.092,0.0045),text.gval.xy=c(0.21,0.0045))
  text(lab.x,lab.y,labels="C",cex=lab.size)
  
  # degree
  graph <- "networks/rnd200/1.txt"
  delta.file <- "results/rnd200.delta.csv"
  mdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.092,0.0045),text.gval.xy=c(0.23,0.0045))
  text(lab.x,lab.y,labels="D",cex=lab.size)

  # sf
  graph <- "networks/sf200/1.txt"
  delta.file <- "results/sf200.delta.csv"
  mdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.092,0.0045),text.gval.xy=c(0.23,0.0045))
  text(lab.x,lab.y,labels="E",cex=lab.size)

  # star

  graph <- "networks/star200/1.txt"
  delta.file <- "results/star200.delta.csv"
  mdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.092,0.0045),text.gval.xy=c(0.21,0.0045))
  text(lab.x,lab.y,labels="F",cex=lab.size)

#  graph <- "networks/fail64/1.txt"
#  delta.file <- "results/fail64.delta.csv"
#  mdelta.fig(mv.file=delta.file,graph=graph,
#             text.gamma.xy=c(0.092,0.0045))
#  text(lab.x,lab.y,labels="F",cex=lab.size)
  
  dev.off()
  
  par(op)      
  
  
  
}

PNAS.Vx.fig <- function()
{
  # Since all same size can just use normal controls
  
  jpeg("figures/figVx.jpg", units="cm",width = 6.8, height = 10, 
       res=600,
       pointsize=6,quality=95)
  
  
  op <- par(mfrow=c(3,2),mar=c(2,3.28,0.5,0.5))
  
  lab.x=0.85
  lab.y=0.0035
  lab.size=2.1
  
  #pan
  graph <- "networks/pan200/1.txt"
  delta.file <- "results/pan200.delta.csv"
  vdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.12,0.00135),
             text.gval.xy=c(0.2,0.00135),
             use.ylim=TRUE,
             fix.ylim=c(0,0.0014))
  text(0.85,0.0013,labels="A",cex=lab.size)
  
  # ring
  graph <- "networks/ring200/1.txt"
  delta.file <- "results/ring200.delta.csv"
  vdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.12,0.00135),
             text.gval.xy=c(0.25,0.00135),
             use.ylim=T,
             fix.ylim=c(0,0.0014))
  text(0.85,0.0013,labels="B",cex=lab.size)

  graph <- "networks/lat200/1.txt"
  delta.file <- "results/lat200.delta.csv"
  vdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.12,0.00135),
             text.gval.xy=c(0.23,0.00135),
             use.ylim=T,
             fix.ylim=c(0,0.0014))
  text(0.85,0.0013,labels="C",cex=lab.size)
  

  # degree
  graph <- "networks/rnd200/1.txt"
  delta.file <- "results/rnd200.delta.csv"
  vdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.12,0.00135),
             text.gval.xy=c(0.25,0.00135),
             use.ylim=T,
             fix.ylim=c(0,0.0014))
  text(0.85,0.0013,labels="D",cex=lab.size)
  
  # sf
  graph <- "networks/sf200/1.txt"
  delta.file <- "results/sf200.delta.csv"
  vdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.12,0.0019),
             text.gval.xy=c(0.25,0.0019),
             use.ylim=T,
             fix.ylim=c(0.0,0.002))
  text(lab.x,0.0018,labels="E",cex=lab.size)
  
  # star
  
  graph <- "networks/star200/1.txt"
  delta.file <- "results/star200.delta.csv"
  vdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.11,0.0615),
             text.gval.xy = c(0.25,0.0615),
             use.ylim=T,
             fix.ylim=c(0,0.065))
  text(lab.x,0.055,labels="F",cex=lab.size)
  
  dev.off()
  
  par(op)      
}

###################### starpan figure for supplementary
##########################################################

fig.starpan <- function()
{
  jpeg("figures/Figstarpan.jpg", units="cm",width = 10, height = 8, 
       res=600,
       pointsize=6,quality=95)
  lab.size=2.5
  op <- par(mfrow=c(2,2),mar=c(2,3.28,0.5,0.5))
  
  plot.network("networks/starpan200/1.txt",layout=layout_components)
  text(-0.9,0.8,labels="A",cex=lab.size)  
  graph <- "networks/starpan200/1.txt"
  delta.file <- "results/starpan200.delta.csv"
  mdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.6,0.001),
             text.gval.xy =c(0.71,0.001))
  text(0.1,0.0045,labels="B",cex=lab.size) 
  vdelta.fig(mv.file=delta.file,graph=graph,
             text.gamma.xy=c(0.6,0.004),
             text.gval.xy = c(0.72,0.004),
             use.ylim=T,
             fix.ylim=c(0,0.02))
  text(0.1,0.018,labels="C",cex=lab.size) 
  
  fix.file <- "results/starpan200.fix.csv"
  fix.fig(fix.file=fix.file,graph=graph,et = c(0.71,0.18),
          gt=c(0.6,0.18))
  text(0.1,0.9,labels="D",cex=lab.size) 
  
  dev.off()  
  par(op)
}  