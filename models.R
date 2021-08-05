###############################################
# models.R
# Fixation and Delta models
#
# Includes basic tools for analysis
###############################################
##############################################
##############################################
# delta.p
#
# Calculate Mdeltax and Vdeltax for each p
# Return table for each p:
# p Mx Vx N s
#############################################
delta.p <- function(delta.file="results/pan64.delta.csv")
{
  df <- read.csv(delta.file,header=T)
  N <- unique(df$N)
  s <- unique(df$s)
  df$fdelta <- (df$m1-df$m0)/N
  p.vals <- sort(unique(df$p))
  
  res <- NULL
  for (prob in p.vals)
  {
    dfp <- subset(df,df$p==prob)
    
    ci=qt(0.975, nrow(dfp) - 1) * sd(dfp$fdelta) / sqrt(nrow(dfp))
    
    pvar <- var(dfp$fdelta)    
    civarlow <- (pvar*(nrow(dfp)-1))/qchisq(0.975,(nrow(dfp)-1))
    civarhigh <- (pvar*(nrow(dfp)-1))/qchisq(0.025,(nrow(dfp)-1))
    
    res <- rbind(res,c(p=prob,mx=mean(dfp$fdelta),
                         vx=pvar,
                         cix = ci,
                         civarlow = civarlow,
                         civarhigh = civarhigh,
                 N=N,
                 s=s))
  }
  print(nrow(dfp))
  as.data.frame(res)
}
###################################
# Mdelta and Vdelta models using Kimura very rough approximation
# Assumes well-mixed population
###################################
#
###############
# mean.delta
###############
mean.delta <- function(N=64,s,p=0.5)
{
  s * p * (1 - p) # model correction * ((N-1)/N)
}
##################
# var.delta
##################
# Note: takes s but is independent of selection strength.
# This is just done for consistency - other versions use s
#
var.delta <- function(N=64,s=0.0,p=0.5)
{
  (p*(1-p)/N)
}
#################################################################
############  Mdelta, Vdelta using matrices for space
#################################################################

# folder - path to folder holding networks, naming 1.txt, 2.txt, ...
# num.spaces - number of spaces to use for estimate
# For regular spaces or spaces with fixed structure num.spaces=1
#
# Generalised model of mean delta -- 
# For each node, calculate mean.delta with size correction
# using the degree for the node * (1/N)
# 
#
mean.delta.network <- function(N=64,s,p=0.5,
                               graph="networks/pan64/1.txt")
{
    graph <- as.matrix(read.table(file=graph))  
    g <- gamma(graph)
    print(paste("Mean Correction:",g[1],sep=""))
    g[1]*mean.delta(N=N,s=s,p=p)
}

##############################################
# Generalisation of variance for networks
##############################################

var.delta.network <- function(N=64,s=0.02,p=0.5,
                              graph="networks/pan64/1.txt")
{
    gm <- as.matrix(read.table(file=graph))  
    g <- gamma(gm)
    print(paste("Var Correction:",g[2],sep=""))
    g[2] * var.delta(N=N,s=s,p=p) 
}

varMgamma.delta.network <- function(N=64,s,p=0.5,
                              graph="networks/pan64/1.txt")
{
  graph <- as.matrix(read.table(file=graph))  
  g <- gamma(graph)
  print(paste("Var Correction:",g[1],sep=""))
  g[1] * var.delta(N=N,s=s,p=p) # Use gamma M for variance - 
}

################################################################
# gamma(g)
###########
# Calculate gamma_m and gamma_v for a given adjacency matrix
# Return c(gamma_m,gamma_v)
################################################################
gamma <- function(g)
{
  v <- rep(1,nrow(g))
  N <- nrow(g)
  
  mgamma <- ((g%*%v-1)/(g%*%v))
  mgamma <- v%*%mgamma
  
  D = rowSums(g)
  
# New improved measure of gamma V
  tot = 0
  for (i in 1:N)
  {
    distar = D[i] * (sum(1/D[which(g[i,]>0)]))^2
    
    tot = tot + (distar - 1)/D[i]
#    print(paste("D[i]",D[i],"DIStar:",distar,"SUM:",((sum(1/D[which(g[i,]>0)]))^2),
#          "TOT:",(distar-1)/D[i]))
  }
# OLD VERSION   
#  vgamma <- g%*%g
#  vgamma <- (vgamma%*%v)/(g%*%v)
#  vgamma <- (vgamma - 1)/(g%*%v)
#  vgamma <- v%*%vgamma
  c(mgamma/N,tot/N)
}
##################################################
# FIXATION
##################################################

#####################################################
# This is the original panmictic fixation by kimura
#####################################################

kimura.ux <- function(N,s,p)
{
  (1-exp(-2*N*s*p))/(1-exp(-2*N*s))
}

gamma.fix <- function(N,s,p,gammas)
{
  g = gammas[1]/gammas[2]
  (1-exp(-2*g*N*s*p))/(1-exp(-2*g*N*s))
}

#################################################################
# load.fixation
# Load data for fixation run
#################################################################
load.fixation <- function(fix.file="results/pan64.fix.csv")
{
  f <- read.csv(fix.file,header=T)
  N <- unique(f$N)
  s <- unique(f$s)
  
  p.vals <- sort(unique(f$p))
  res <- NULL
  
  for (prob in p.vals)
  {
    dfp <- subset(f,f$p==prob)
    # and to store gens we need to use just those where p==1
    gens.count <- subset(dfp,dfp$fix==1)
    gens.ci=qt(0.975, nrow(gens.count) - 1) * sd(gens.count$gens) / sqrt(nrow(gens.count))
    
    ci=qt(0.975, nrow(dfp) - 1) * sd(dfp$fix) / sqrt(nrow(dfp))
    res <- rbind(res,c(p=prob,
					   fix.m=mean(dfp$fix), # Note 0 or 1 for fix
                       fix.v=var(dfp$fix),
					             ci = ci,
                       N=N,
                       s=s,
					             gens=mean(gens.count$gens),
					             gens.ci = gens.ci))
  }
  as.data.frame(res)
}

