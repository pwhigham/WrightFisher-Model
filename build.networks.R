#############################################
# Create the networks using rewiring
# Save as files for use with C version
############################################

library(igraph) # Used to create ring and scale-free graph

write.matrix <- function(m,file)
{
  write.table(m,file,row.names=FALSE,col.names=FALSE)
}
################### ADJACENCY GRAPHS ############################
#
# Create graphs with 0/1 (no weighting)
# These are classic adjacency graphs - used for PNAS paper
# Note that the diag(g) = 1 since treating each node on the graph as
# a location means that each node is part of its deme.
################################################################

#
# Panmictic - fully connected
#
graph.panmictic <- function(N=64,filename="networks/pan64/1.txt")
{
  res <- matrix(nrow=N,ncol=N,data=1) # all 1's
  write.matrix(res,file=filename)
}
# 
# ring
# Ring where each node is connected to two nearest neigbours
#
graph.ring <- function(N=64,filename="networks/ring64/1.txt")
{
    g <- watts.strogatz.game(1,N,1,0.0) # see igraph - lazy way to create
    network <- as.matrix(as_adjacency_matrix(g))
    diag(network) <- 1
    write.matrix(network,file=filename)
}
#
# star
# Centre node connected to all - all other nodes just connected to centre
#
graph.star <- function(N=64,num.files=1,filename="networks/star64/1.txt")
{
    network <- matrix(nrow=N,ncol=N)
    network[,] <- 0
    network[1,] <- 1
    network[,1] <- 1
    diag(network) <- 1
    write.matrix(network,file=filename)
}

# For 10x20 lattice set dimvector=(10,20)
# Create a layout using ll = layout_on_grid(g,width=20,height=10)
# and then plot.network("networks/lat200/1.txt",layout=ll)
#

#
# lattice
# See igraph for help
graph.lattice <- function(dimvector=c(8,4),N=64,filename="networks/lattice200/1.txt")
{
    g <- make_lattice(dimvector=dimvector,
                      directed=FALSE,
                      circular=TRUE)
    network <- as.matrix(as_adjacency_matrix(g)) # Fixed for now
    diag(network) <- 1
    write.matrix(network,file=filename)
}
#
# scale-free
# See igraph for help
#
graph.sf <- function(power=1.5,N=64,filename="networks/sf64/1.txt")
{
    g <- sample_pa(n=N,power=power,directed=FALSE)
	# Ensure graph is connected
    while(!is_connected(g)) g <- sample_pa(N,power=power,directed=FALSE) 
    network <- as.matrix(as_adjacency_matrix(g)) # Fixed for now
    diag(network) <- 1
    write.matrix(network,file=filename)
}
#
# line
#
graph.line <- function(N=64,filename="networks/line_64/1.txt")
{
  res <- matrix(nrow=N,ncol=N)
  res[,] <- 0
  diag(res) <- 1    
  res[1,2] <- 1 # Just connected to one 
  for (i in 2:(N-1))  # for each row...
  {
    res[i,(i-1)] <- 1  # Connect to neighbours
    res[i,(i+1)] <- 1
  }
  res[N,(N-1)] <- 1 # and the other end
  write.matrix(res,file=filename)
}
#######################################################
#
# graph.degrees
#
# Given a function that generates a degree vector, creates
# networks by randomly wiring with given degree vectors.
# The ... are additional arguments (apart from N) for the
# degree function degree.func.
# The example would call using rnorm might be:
# 			graph.degrees(N=64,degree.func=rnorm,5,1)
# would produce a graph with mean node degree 5 and sd 1
#######################################################
graph.degrees <- function(N=64,
                          degree.func=rnorm,
                          filename="networks/deg64/1.txt",
                                  ...)
{
  fn <- 0
  while (fn < 1)
  {
    network <- matrix(nrow=N,ncol=N,data=0)
    nodevals <- as.integer(degree.func(N,...))
    while(length(which(nodevals<=0)) > 0)
      nodevals <- as.integer(degree.func(N,...))
    diag(network) <- 1
    fails <- 0
    while( (length(which(nodevals > 0)) > 1) && 
           (fails < N))
    {
      nodes <- which(nodevals > 0) # candidate nodes
      picknodes <- sample(nodes,2) # pick 2 at random to join
      
      if (network[picknodes[1],picknodes[2]] == 0)  # Not done before
      {
        network[picknodes[1],picknodes[2]] <- 1
        network[picknodes[2],picknodes[1]] <- 1
        
        nodevals[picknodes[1]] <- nodevals[picknodes[1]] - 1
        nodevals[picknodes[2]] <- nodevals[picknodes[2]] - 1
        fails <- 0
      }
      else
      {
        fails <- fails + 1
      }
    }
    if (fails < N)
    {
      # Confirm connected network
      g <- graph_from_adjacency_matrix(network,mode="directed")
      if (is_connected(g))
      {
        fn <- fn + 1
        write.matrix(network,file=filename)        
      }
    }
  }  
}

panline <- function(N=64,filename="networks/fail64/1.txt")
{
  network <- matrix(nrow=64,ncol=64)
  
  # 30,25 individuals per pan connected by a line...
  
  network[1:30,1:30] <- 1
  network[39:64,39:64] <- 1
  
  # and connect using a line
  for (i in 30:39)
  {
    network[i,i+1] = 1
    network[i+1,i] = 1
  }
  diag(network) <- 1
  g <- graph_from_adjacency_matrix(network,mode="directed")
  print(paste ("Connected?",is_connected(g))) 
  write.matrix(network,file=filename)
}
#
# Hardwired - star at one end, panmictic at other, linked by line
# N=200
#
starpan <- function(filename="networks/starpan200/1.txt")
{
  N = 200
  network <- matrix(nrow=N,ncol=N,data=0)
  # Make a star connected to node 1 with 30 nodes
  #
  # Make panmictic at one end, star at the other
  network[1:90,1:90] <- 1
  
   for (i in 90:101)
   {
     network[i,i-1] <- 1
     network[i-1,i] <- 1
   }
   
   network[100,101:200] <- 1
   network[101:200,100] <- 1
#   
  diag(network) <- 1
  g <- graph_from_adjacency_matrix(network,mode="directed")
  print(paste ("Connected?",is_connected(g))) 
  write.matrix(network,file=filename)
}
##########################
# Star -> LINE -> Star
#
# Hardwired for N=200
###########################
sls <- function(filename="networks/sls200/1.txt")
{
  N = 200
  network <- matrix(nrow=N,ncol=N,data=0)
  # Make a star connected to node 1 with 30 nodes
  #
  # Make panmictic at one end, star at the other
  
  network[1,1:90] <- 1
  network[2:90,1] <- 1
  
  for (i in 90:101)
  {
    network[i,i-1] <- 1
    network[i-1,i] <- 1
  }
  
  network[100,101:200] <- 1
  network[101:200,100] <- 1
  #   
  diag(network) <- 1
  g <- graph_from_adjacency_matrix(network,mode="directed")
  print(paste ("Connected?",is_connected(g))) 
  write.matrix(network,file=filename)
}


#######################################################
# Create islands of island.size and connect them
# with random migrant links between nearest islands
# Assumes island size is chosen to be multiple of N
#######################################################
islands <- function(N=200,filename="networks/island200/1.txt",
                    island.size=20,migrants=5)
{
  network <- matrix(nrow=N,ncol=N,data=0)
  diag(network) <- 1
  
  num.islands <- N/island.size
  island.size <- island.size-1
  first = 1
  for (i in 1:(num.islands))
  {
    network[first:(first+island.size),first:(first+island.size)] <- 1
    first <- first + island.size+1
  }
  # now randomly connect consecutive islands with migrants
  first = 1
  for (i in 1:(num.islands-1))
  {
    firstspots <- sample(first:(first+island.size),migrants,replace=FALSE)
    secspots <- sample((first+island.size+1):(first+(2*island.size)+1),migrants)
    
    network[firstspots,secspots] <- 1
    network[secspots,firstspots] <- 1
    first <- first + island.size+1
  }
  
  write.matrix(network,file=filename)
}

############################## WEIGHTED NETWORKS #####################
######################################################################
#
# We don't have any theory for weighting locations, so initial experiments
# will just set each outward edge to some value w (w <= 1.0)
# We will create the network(s) by loading in a previous network will all 1's and setting
# all values = 1 to w, then setting diag(g) = 1 since each location should 
# definitely be used in a deme.
######################################################################

weight.network <- function(graph="networks/pan200/1.txt", wgraph="Wnetworks/pan200/1.txt",
               w = 0.5)
{
  n <- as.matrix(read.table(graph))
  # Since n is adjacency matrix 0,1's, we just need to multiply by w
  n <- n * w
  diag(n) <- 1
  write.matrix(n,file=wgraph)
}
