bridges <- function (x)
{
  # Find the bridges in an undirected graph.

  # IMPLEMENTATION OF THE ALGORITHM FOR BRIDGES PROPOSED BY SCHMIDT (2013)
  # A simple test on 2-vertex- and 2-edge-connectivity. Information processing letters 113, 241-244.
  mode ="undirected"
  x <- getWmat(x)
  
  # 1. compute the depth-first search tree
  # 1.a depth first search
  igx <- graph.adjacency(x, mode = mode)
  if(!igraph::is.connected(igx)) stop ("The graph is disconnected, please enter a connected graph")
  dfs <-graph.dfs(igx, root = 1, neimode = "all", father = TRUE)
  # 1.b buid the tree
  # edges in the tree are toward the father
  tree.dir <- getWmat(data.frame(cbind(1:ncol(x), dfs$father)[-1,]))
  tree.und <- symmetrize(tree.dir, rule ="weak")
  
  # 2. orientation and network recomposition
  # 2.a reorient the edges outside the tree in x toward the dfs childs
  ntree.und <- x - tree.und
  ntree.dir <- ntree.und [dfs$order, dfs$order]
  ntree.dir[lower.tri(ntree.dir)] <- 0
  ntree.dir <- ntree.dir[match(1:ncol(x), dfs$order), match(1:ncol(x), dfs$order)]
  # 2.b recompose the network as directed
  net <- ntree.dir + tree.dir
  
  # 3. chain decomposition
  visited <- rep(FALSE, ncol(x))
  chains <- x
  
  for(i in dfs$order)
  {
    # - for any nodes i, in order of dfs
    
    nei <- (1:ncol(x))[as.logical(ntree.dir[i,])]
    
    if(length(nei) >= 1)
    {
      visited[i] <- TRUE
      for (j in nei)
      {
        # - for every backedge (i.e., edges not in the tree) that starts at i
        
        # We traverse the chain that starts with the backedge and stop at the ﬁrst
        #   vertex that is already marked as visited. During such a traversal,
        #   every traversed vertex is marked as visited.
        already <- FALSE
        chain <- matrix(c(i, j), ncol = 2, nrow = 1)
        visited[j] <- TRUE
        
        # travel the chain until a visited node is found
        from <- j
        while (already == FALSE)
        {
          to <- (1:ncol(x))[as.logical(tree.dir[from,])][1]
          chain <- rbind(chain, c(from, to))
          if(visited[to] == TRUE) already <- TRUE
          visited[to] <- TRUE
          from <- to
        }
        # in chains, delete all the edges in in the chain, so that only the edges outside any chain will remain
        chain <- (-1) * as.matrix(1*sparseMatrix(i=chain[,1], j=chain[,2], dims = c(nrow(x), ncol(x))))
        chain <- symmetrize(chain, rule="weak")
        chains<- chains + chain
      }
    }
  }
  
  # 5. output
  # the edges not included in any chain are bridges
  bridges <- chains * (chains == 1)
  bridges
}

negativedges<-function(ggm, negpro=.5)
{
  diag(ggm) <- 0
  if(negpro == 0) return (ggm)
  edges <- sum(ggm)/2
  N <- ncol(ggm)
  pos <- round(negpro*edges)
  neg <- edges-pos
  sgn <- sample(c(rep(1,pos), rep(-1, neg)), size = edges, replace = FALSE)
  ggm[lower.tri(ggm)][ggm[lower.tri(ggm)] == 1] <- sgn
  ggm <- symmetrize(ggm, rule = "lower")
  ggm
}


# convert the adjacency matrix ggm into a partial correlation matrix pcm
ggm2pcm<-function(ggm, minpcor = .2,  maxpcor = .9, maxiter = 1000, verbose = FALSE)
{
  # Starting from an adjacency matrix (ggm), the function generates a partial correlation
  # matrix (pcm) by replacing the ones in ggm with random numbers from minpcor to maxpcor.
  # from the pcm, a correlation matrix (cm) is generated
  # Finally it checks for positive-definiteness: if a positive definite cm cannot be
  # reached after maxiter iterations, the function returns -1
  # Parameters:
  # -ggm (binary matrix): adjacency matrix, the output of function Topology
  # -minpcor, maxpcor: the minimum and the maximum values of the partial correlations 
  # -maxiter: maximum number of iterations allowed in searching for a positive-definite matrix.

  N <- ncol(ggm)
  iter <- 0
  check <- FALSE
  
  while(check == FALSE)
  {
    if(iter > maxiter & check == FALSE)
    {
      if(verbose) warning(paste("Positive-definiteness was not reached after ", iter-1, " iterations"))
      return("pcm" = NA)
    }
    
    iter <- iter + 1
    if(verbose) print(paste("ggm2pcm iter =", iter))
    pcr = runif(N*(N-1)/2, minpcor, maxpcor) # generate uniform numbers in the allowed range
    ld.ggm <- ggm[lower.tri(ggm)]
    pcm <- matrix(0,ncol = N,nrow = N) # empty partial correlation matrix
    pcm[lower.tri(pcm)] <- ld.ggm*pcr # put the uniform numbers in the partial cor matrix
    pcm[upper.tri(pcm)]<-t(pcm)[upper.tri(pcm)] # this is the partial correlation matrix
    diag(pcm) <- 1
    
    if(!is.positive.definite(pcm)) next
    
    # check that also the correlation matrix is positive definite
    suppressWarnings(cm <- pcor2cor(pcm)) 
    if(any(is.na(cm))) 
    {next} else check <- is.positive.definite(cm) # TRUE # check for pd
  }
  pcm
}


ERtopology<-function(N, m, force.connected = TRUE, plot.top = FALSE, maxiter = 1000)
{
  # Generate a random Erdos-Renyi network with N nodes and m edges. if force.connected = TRUE,
  # a network is redrawn until a connected network emerges or until maxiter networks have been
  # drawn without one being connected.
  # consider that a network with m < N-1 cannot be connected, in this case a disconnected
  # network is returned and if force.connected == TRUE, also a warning.
  
  
  if(m > (N*(N-1))/2)
  {
    warning (paste("A network of N = ", N, " nodes can have at most", (N*(N-1))/2, " edges. A network with", (N*(N-1))/2, "edges is returned instead of m = ", m, " edges."))
    m <- (N*(N-1))/2
  }  
  if (force.connected == TRUE & m < (N-1))
  {
    warning(paste("A network cannot be connected with N = ", N, " nodes and m = ", m, " edges. Please select at least m = ", N-1,  "edges. A disconnected network is returned."))
    force.connected <- FALSE
  }
  
  # generate a random network with specified number of edges
  # generate the network
  ggm<-erdos.renyi.game(n = N, p.or.m = m, type = "gnm", directed = FALSE, loops = FALSE)
  
  if (force.connected == TRUE)
  {
    # if force.connected == TRUE and the network is not connected, draw another network
    iter <- 0
    while(!igraph::is.connected(ggm) & iter < maxiter)
    {
      iter <- iter + 1
      ggm <- erdos.renyi.game(n = N, p.or.m = m, type = "gnm", directed = FALSE, loops = FALSE)
    } 
    if (!igraph::is.connected(ggm))
    {
      warning(paste("A connected network was not retrieved after iter =", iter, " iterations. A disconnected network is therefore returned."))
      force.connected <- FALSE
    }
  }
  
  am <- getWmat(ggm)
  if(plot.top) qgraph(am)
  am
  
#   # TEST OF FUNCTION EStopology, with 50 nodes
#   # case A) m < N-1
#   m = 48
#   # A1) force.connected = TRUE: a warning is expected with disconnected networks with m edges
#   system.time(
#     nets<-lapply(seq(from=.01, to=1, by=.01), function(x) ERtopology(N=50, m=m, force.connected=TRUE, plot.top=FALSE))
#               )
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected) == FALSE) # all disconnected
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x)) == m) # all with m edges
#   
#   # A2) force.connected = FALSE: no warnings expected, all disconnected networks with m edges
#   system.time(
#     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) ERtopology(N=50, m=m, force.connected=FALSE, plot.top=FALSE))
#   )
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected) == FALSE)
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x)) == m)
#    
#   # case B) m >= N-1
#   m = 55 # a few edges, this makes a disconnected network really likely.
#   # B1) force.connected = TRUE. In many cases warnings are returned with disconnected networks because a connected one is not find after maxiter = 1000 attempts
#   # m=55 are few for a graph with N=50 to be connected by chance.
#   system.time(
#     nets<-lapply(seq(from=.01, to=1, by=.01), function(x) ERtopology(N=50, m=m, force.connected=TRUE, plot.top=FALSE))
#   )
#   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected))
#   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
#   
#   # B2) force.connected = FALSE: no warnings expected, some disconnected and perhaps few connected (by chance) networks with m edges
#   system.time(
#     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) ERtopology(N=50, m=m, force.connected=FALSE, plot.top=FALSE))
#   )
#   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected))
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x)) == m)
#     
#   # case C) m > N(N-1)/2. A network with m = (N-1)/2 edges is returned with a warning
#   m = (50*49/2)+1
#   system.time(
#     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) ERtopology(N=50, m=m, force.connected=FALSE, plot.top=FALSE))
#   )
#   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
}

WStopology <- function(N, m, p_ws = .1, exact.m = TRUE, force.connected = TRUE, plot.top = FALSE, maxiter=1000)
{
  # generate a Watts&Strogatz network with specified number of nodes (N),
  # probability of random rewiring (p_ws) and the total number of edges (m).
  
  # force.connected: the process of random rewiring in WS model does not ensure that the
  # resulting network is connected. If force.connected == TRUE the network is resampled
  # until a connected network is found.
  # consider that a network with m < N-1 cannot be connected, in this case a disconnected
  # network is returned and if force.connected == TRUE, also a warning.
  
  # in the WS model, the neighborhood of each node, nei (2 times the nodes' degree)
  # is equal for all nodes in the initial regular lattice, therefore the actual number of
  # edges will be q = N*nei. If m is chosen such as m/N is integer (i.e., m is a multiple of N)
  # then m edges will be present, independent of parameter exact.m.
  
  # exact.m:  If m is chosen such as m/N is not integer, the function picks
  # nei = ceiling (m/N), therefore the true number of edges q is higher than m.
  # if exact.m == FALSE, a number of edges q > m is kept.
  # if exact.m == TRUE and force.connected = FALSE, q-m edges are removed randomly
  #   the network can result disconnected as the end of this process.
  # if exact.m == TRUE and force.connected = TRUE, q-m edges are removed randomly, with
  #   the constraint that the network remains connected. function bridges
  #   is used at each step to ensure that no bridge is removed (a bridge is an edge
  #   that, if removed, leaves the network disconnected).
  
  # maxiter: if a connected network is not found after maxiter attempts, a 
  # disconnected network is returned with a warning

  nei <- ceiling(m/N)
  if (force.connected == TRUE & m < (N-1) & exact.m == TRUE)
  {
    warning(paste("A network cannot be connected with N = ", N, " nodes and m = ", m, " edges. Please select at least m = ", N-1,  "edges. A disconnected network is returned."))
    force.connected <- FALSE
  }
  if(m > (N*(N-1))/2)
  {
    warning (paste("A network of N = ", N, " nodes can have at most", (N*(N-1))/2, " edges. A network with", (N*(N-1))/2, "edges is returned instead of m = ", m, " edges."))
    m <- (N*(N-1))/2
  }
  
  # generate the graph
  ggm <- watts.strogatz.game(dim = 1, size = N, nei = nei, p = p_ws)
  if (force.connected == TRUE)
  {
    iter <- 0
    while(!igraph::is.connected(ggm) & iter < maxiter)
    {
      iter <- iter + 1
      ggm <- watts.strogatz.game(dim = 1, size = N, nei = nei, p = p_ws)
    } 
    if (!igraph::is.connected(ggm))
    {
      warning(paste("A connected network was not retrieved after iter =", iter, " iterations. A disconnected network is therefore returned."))
      force.connected <- FALSE
    }
  }
  
  am <- getWmat(ggm) # convert ggm into an adjacency matrix
  q <- ecount(ggm)
  
  if(exact.m == TRUE & m != q & force.connected == FALSE)
  {
    # a vector of all the edges
    vec <- mat2vec(am)
    # remove m-q edges randomly
    vec[vec == 1][sample(x = 1:q, size = (q-m), replace = FALSE)] <- 0
    # rebuild the adjacency matrix from vec
    am2 <- matrix(0, ncol = N, nrow = N)
    am2[upper.tri(am2)] <- vec
    am2[lower.tri(am2)] <- t(am2)[lower.tri(am2)]
    am <- am2
  } else if (exact.m == TRUE & m != q & force.connected == TRUE)
  {
      am2 <- am
      for(i in 1:(q-m))
      {
        # at each step remove an edge that is not a bridge
        bri <- bridges(am2)
        vec <- cbind(mat2vec(am2), mat2vec(bri))
        vec[vec[,1] == 1 & vec[,2] == 0, 1][sample(x = 1:nrow(vec[vec[,1] == 1 & vec[,2] == 0,]), size = 1)] <- 0
        
        # and rebuild am2
        am2 <- matrix(0, ncol = N, nrow = N)
        am2[upper.tri(am2)] <- vec[,1]
        am2[lower.tri(am2)] <- t(am2)[lower.tri(am2)]
        am <- am2
      }
    }
  if(plot.top) qgraph(am)
  am
  
  
  #   # TEST OF FUNCTION WStopology, with 50 nodes
  #   # case A) m < N-1
  #   m = 48
  #   # A1) exact.m = TRUE & force.connected = TRUE: a warning is expected with disconnected networks and m edges
  #   system.time(
  #     nets<-lapply(seq(from=.01, to=1, by=.01), function(x) WStopology(N=50, m=m, p_ws=x, exact.m=TRUE, force.connected=TRUE, plot.top=FALSE))
  #               )
  #   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected) == FALSE)
  #   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x)) == m)
  #   
  #   # A2) exact.m = TRUE & force.connected = FALSE: no warnings expected, all disconnected networks with m edges
  #   system.time(
  #     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) WStopology(N=50, m=m, p_ws=x, exact.m=TRUE, force.connected=FALSE, plot.top=FALSE))
  #   )
  #   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected) == FALSE)
  #   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x)) == m)
  #   
  #   # A3) exact.m = FALSE & force.connected = TRUE: no warnings expected, mostly connected network with q > m and some disconnected
  #   # network with a warning if a connected network is not recovered after maxiter attempts.
  #   system.time(
  #     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) WStopology(N=50, m=m, p_ws=x, exact.m=FALSE, force.connected=TRUE, plot.top=FALSE))
  #   )
  #   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected)) # how many are connected?
  #   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
  #   
  #   # A4) exact.m = FALSE & force.connected = FALSE: no warnings expected, since the true number of edges can be now
  #   # q > m, the network can be (but does not have to be) connected.
  #   system.time(
  #     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) WStopology(N=50, m=m, p_ws=x, exact.m=FALSE, force.connected=FALSE, plot.top=FALSE))
  #   )
  #   sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected)
  #   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
  #   
  #   # case B) m >= N-1 and is not a multiple of N
  #   m = 98 #m is not a multiple of N, so exact.m will act here
  #   # B1) exact.m = TRUE & force.connected = TRUE: connected networks with exactly m edges. In some cases warnings can be returned with disconnected networks
  #   system.time(
  #     nets<-lapply(seq(from=.01, to=1, by=.01), function(x) WStopology(N=50, m=m, p_ws=x, exact.m=TRUE, force.connected=TRUE, plot.top=FALSE))
  #   )
  #   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected))
  #   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
  #   
  #   # B2) exact.m = TRUE & force.connected = FALSE: no warnings expected, some disconnected and some connected (by chance) networks with m edges
  #   system.time(
  #     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) WStopology(N=50, m=m, p_ws=x, exact.m=TRUE, force.connected=FALSE, plot.top=FALSE))
  #   )
  #   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected))
  #   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x)) == m)
  #   
  #   # B3) exact.m = FALSE & force.connected = TRUE: no warnings expected, connected networks, all with q > m
  #   system.time(
  #     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) WStopology(N=50, m=m, p_ws=x, exact.m=FALSE, force.connected=TRUE, plot.top=FALSE))
  #   )
  #   sum(sapply(lapply(nets, graph.adjacency), igraph::is.connected))
  #   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
  #   
  #   # B4) exact.m = FALSE & force.connected = FALSE: no warnings expected, since the true number of edges can be now
  #   # q > m, the network can be (but does not have to be) connected.
  #   system.time(
  #     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) WStopology(N=50, m=m, p_ws=x, exact.m=FALSE, force.connected=FALSE, plot.top=FALSE))
  #   )
  #   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected))
  #   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
  #   
  #   # case C) m > N(N-1)/2. A network with m = (N-1)/2 edges (full network) is returned with a warning
  #   m = (50*49/2)+1
  #   system.time(
  #     nets <- lapply(seq(from=.01, to=1, by=.01), function(x) WStopology(N=50, m=m, p_ws=x, exact.m=FALSE, force.connected=FALSE, plot.top=FALSE))
  #   )
  #   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
}
  
BAtopology <- function(N, m, pw_ba = 1, exact.m = TRUE, force.connected = TRUE, plot.top = FALSE)
{
  # Generate a Barabasi-Albert network with specified number of nodes (N) and
  # power law with exponent pw_ba.
  # force.connected: different from WS model, the BA model is always connected if all edges
  # are retained.
  
  # The final number of edges in the barabasi.game is determined by
  # the fact that k edges are added at each step, however this is not true
  # for the first k steps: at step 1 no edge is attached, at step 2 only
  # one edge can be attached, 2 edges at step 3 etc., at kth step and subsequent k-1 edges are added.
  # Therefore the true number of edges q is determined by q = sum(1:(k-1))+k*(N-k)
  # with some algebra, this can be expressed for k as a function of q:
  # k = ((2*N-1)- sqrt((2*N-1)^2-8*q))/2.
  
  # If m is chosen such as k is integer, then the desired number of edges is always
  # present and the network will be connected. This can be ensured by selecting 
  # an integer k and then computing m as m = sum(1:(k-1))+k*(N-k).
  
  # exact.m: if m is chosen such as k is not integer, the function picks
  # k <- ceiling(((2*N-1)- sqrt((2*N-1)^2-8*m))/2), therefore the true number of edges
  # q is higher than m.
  # if exact.m == FALSE, a number of edges q > m is kept.
  # if exact.m == TRUE and force.connected = FALSE, q-m edges are removed randomly
  #   the network can result disconnected as the end of this process.
  # if exact.m == TRUE and force.connected = TRUE, q-m edges are removed randomly, with
  #   the condition that the network remains connected.
  
  if (force.connected == TRUE & exact.m == TRUE & m < (N-1))
  {
    warning(paste("A network cannot be connected with N = ", N, " nodes and m = ", m, " edges. Please select at least m = ", N-1,  "edges. A disconnected network is returned."))
    force.connected <- FALSE
  }
  if(m > (N*(N-1))/2)
  {
    warning (paste("A network of N = ", N, " nodes can have at most", (N*(N-1))/2, " edges. A network with", (N*(N-1))/2, "edges is returned instead of m = ", m, " edges."))
    m <- (N*(N-1))/2
  }  
  
  k <- ceiling(((2*N-1)- sqrt((2*N-1)^2-8*m))/2)  # number of edges added at each step after the kth step
  # generate the network
  ggm <- barabasi.game(n = N, power = pw_ba, directed = FALSE, m = k)
  am <- getWmat(ggm) # convert ggm into an adjacency matrix
  q <- ecount(ggm)
  
  if(exact.m == TRUE & m != q & force.connected == FALSE)
  {
    # a vector of all the edges
    vec <- mat2vec(am)
    # remove m-q edges randomly
    vec[vec == 1][sample(x = 1:q, size = (q-m), replace = FALSE)] <- 0
    # rebuild the adjacency matrix from vec
    am2 <- matrix(0, ncol = N, nrow = N)
    am2[upper.tri(am2)] <- vec
    am2[lower.tri(am2)] <- t(am2)[lower.tri(am2)]
    am <- am2
  } else if (exact.m == TRUE & m != q & force.connected == TRUE)
  {
    am2 <- am
    for(i in 1:(q-m))
    {
      # at each step remove an edge that is not a bridge
      bri <- bridges(am2)
      vec <- cbind(mat2vec(am2), mat2vec(bri))
      vec[vec[,1] == 1 & vec[,2] == 0, 1][sample(x = 1:nrow(vec[vec[,1] == 1 & vec[,2] == 0,]), size = 1)] <- 0
      
      # and rebuild am2
      am2 <- matrix(0, ncol = N, nrow = N)
      am2[upper.tri(am2)] <- vec[,1]
      am2[lower.tri(am2)] <- t(am2)[lower.tri(am2)]
      am <- am2
    }
  }
  if(plot.top) qgraph(am)
  am
#   
#   # TEST OF FUNCTION BAtopology, with 50 nodes
#   # case A) m < N-1
#   m = 48
#   # A1) exact.m = TRUE & force.connected = TRUE: a warning is expected with disconnected networks and m edges
#   system.time(
#     nets<-lapply(seq(from=.5, to=3, by=.1), function(x) BAtopology(N=50, m=m, pw_ba = x, exact.m = TRUE, force.connected = TRUE))
#   )
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected) == FALSE)
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x)) == m)
#   
#   # A2) exact.m = TRUE & force.connected = FALSE: no warnings expected, all disconnected networks with m edges
#   system.time(
#     nets<-lapply(seq(from=.5, to=3, by=.1), function(x) BAtopology(N=50, m=m, pw_ba = x, exact.m = TRUE, force.connected = FALSE))
#   )
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected) == FALSE)
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x)) == m)
#   
#   # A3) exact.m = FALSE & force.connected = TRUE: no warnings expected, all connected network with q > m
#   system.time(
#     nets<-lapply(seq(from=.5, to=3, by=.1), function(x) BAtopology(N=50, m=m, pw_ba = x, exact.m = FALSE, force.connected = TRUE))
#   )
#   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected)) # how many are connected?
#   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
#   
#   # A4) exact.m = FALSE & force.connected = FALSE: no warnings expected, all connected.
#   system.time(
#     nets<-lapply(seq(from=.5, to=3, by=.1), function(x) BAtopology(N=50, m=m, pw_ba = x, exact.m = FALSE, force.connected = FALSE))
#   )
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected) == TRUE)
#   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
#   
#   # case B) m >= N-1 and is not a multiple of N
#   m = 95 #  m =95 determines a fractional k, so exact.m will act here
#   # B1) exact.m = TRUE & force.connected = TRUE: connected networks with exactly m edges. In some cases warnings can be returned with disconnected networks
#   system.time(
#     nets<-lapply(seq(from=.5, to=3, by=.1), function(x) BAtopology(N=50, m=m, pw_ba = x, exact.m = TRUE, force.connected = TRUE))
#   )
#   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected))
#   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
#   
#   # B2) exact.m = TRUE & force.connected = FALSE: no warnings expected, some disconnected (by chance) and some connected networks with m edges
#   system.time(
#     nets<-lapply(seq(from=.5, to=3, by=.1), function(x) BAtopology(N=50, m=m, pw_ba = x, exact.m = TRUE, force.connected = FALSE))
#   )
#   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected))
#   all(sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x)) == m)
#   
#   # B3) exact.m = FALSE & force.connected = TRUE: no warnings expected, connected networks, all with q > m
#   system.time(
#     nets<-lapply(seq(from=.5, to=3, by=.1), function(x) BAtopology(N=50, m=m, pw_ba = x, exact.m = FALSE, force.connected = TRUE))
#   )
#   sum(sapply(lapply(nets, graph.adjacency), igraph::is.connected))
#   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
#   
#   # B4) exact.m = FALSE & force.connected = FALSE: no warnings expected, since the true number of edges can be now
#   # q > m, the network will always be connected.
#   system.time(
#     nets<-lapply(seq(from=.5, to=3, by=.1), function(x) BAtopology(N=50, m=m, pw_ba = x, exact.m = FALSE, force.connected = FALSE))
#   )
#   sum(sapply(lapply(nets, graph.adjacency, mode="undirected"), igraph::is.connected))
#   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
#   
#   # case C) m > N(N-1)/2. A network with m = (N-1)/2 edges (full network) is returned with a warning
#   m = (50*49/2)+1
#   system.time(
#     nets<-lapply(seq(from=.5, to=3, by=.1), function(x) BAtopology(N=50, m=m, pw_ba = x, exact.m = FALSE, force.connected = FALSE))
#   )
#   sapply(lapply(nets, graph.adjacency, mode="undirected"), function(x) ecount(x))
}  


# 
# # example of usage.
# # caution: some combination of values seem to lead to systematic failures.
# net<-Topology(N=50, m=50, topology="BA", exact.m=TRUE, force.connected=TRUE, negpro=0, not=NULL, pw_ba=.1, p_ws=.1, minpcor=.1, maxpcor=.6)
# par(mfrow=c(1,2))
# lay<-qgraph(net$ggm)
# qgraph(net$pcm, layout = lay)
# 
Topology<-function(N, m, topology = c("ER", "BA", "WS"), exact.m = TRUE, force.connected = TRUE, negpro = 0, not = NULL, pw_ba = 1, p_ws = .1,  minpcor = .2,  maxpcor = .9, maxiter = c(50, 1000, 50), verbose = TRUE)
{
  # This function creates a network with the desired topology, allowing to regulate directly the number of edges
  #   (which is not always possible using igraph functions).
  # It also tests for graph isomorphisms: if the obtained graph is isomorphic to a graph in the "not" list
  #   the function iterate until a non-isomorphic graph is found (or until maxiter calls are done)

  # N (scalar): number of nodes
  # m (scalar): number of edges
  # topology (string): ER=random, BA=scale-free, WS=small world.
  # pw_ba: the power of the preferential attachment in the barabasi-albert power law (see igraph::barabasi.game).
  #   It is relevant only if topology = "BA"
  # p_ws: the probability of random rewiring in the WS game. It is relevant only if topology = "WS"
  
  # minpcor = minimum value for the weights.
  # maxpcor = maximum value for the weights.  
  
  # exact.m (logical): WS and BA games do not allow any possible value of m given a value of N. In WS game each node has the
  #   same number of edges in the original lattice; In BA game each new node is added with the same number
  #   of new connections to other nodes (with some exception, see BAtopology function for details).
  #   if exact.m = FALSE, a number of edges q >= m will be present. If exact.m = TRUE, the edges will be deleted
  #   randomly until the number of edges is exactly m.
  # force.connected: if TRUE the network will be connected, if FALSE the network can be disconnected.
  #   
  #   In ER game, a network is redrawn until a conncetd network emerges.
  #   In BA game, the initial network is always connected, therefore if exact.m = FALSE, force.connected
  #     will have no effect. However, if exact.m = TRUE, the network can become disconnected after the
  #     removal of the exceeding edges. force.connected = TRUE ensures that at every edge deletion, no
  #     brige is removed. This verification can be however time-consuming.
  #   In WS game, the initial network can be disconnected as an effect of random rewiring and it can also
  #     become disconnected as an effect of edge deletion if exact.m = TRUE. force.conneced ensures that
  #     the network remains connected at each step.
  
  # negpro (integer in [0,1]): proportion of negative edges.
  # not (a list of graph matrices). If the network extracted is equal to one of these graphs, a new network
  #   is extracted. This has been implemented to reduce the risk of having multiple identical networks
  #   in simulation studies.

  # maxiter (vector of integers): maximum number of iterations for different parts of the function.
  #   maxiter[1]: max attempts to search for a network that is not isomorphic to other networks in "not" and that 
  #     is also a positive definite correlation matrix.
  #   maxiter[2]: attempts to find a connected network when topology = "ER" or "WS"
  #   maxiter[3]: attemtps to find a positive definite correlation matrix. I reccommend a low value for this
  #               because in case of repeated failure (maxiter[3] failures), a new network is drawn (up to maxiter[1] times).

  if(length(topology) != 1 | !topology %in% c("ER", "BA", "WS")) stop("Wrong topology specification")
  
  check <- FALSE
  iter <- 0
  while(check == FALSE & iter < maxiter[1])
  {
    iter <- iter + 1
    print(iter)
    if(topology=="ER") ggm <- ERtopology(N = N, m = m, force.connected = force.connected, plot.top = FALSE, maxiter = maxiter[2])
    if(topology=="BA") ggm <- BAtopology(N = N, m = m, pw_ba = pw_ba, exact.m = exact.m, force.connected = force.connected, plot.top = FALSE) 
    if(topology=="WS") ggm <- WStopology(N = N, m = m, p_ws = p_ws, exact.m = exact.m, force.connected = force.connected, plot.top = FALSE, maxiter = maxiter[2]) 
    if(verbose) cat("  topology done...")
    
    ##########################
    # verify nonisomorphism ##
    if(!is.null(not))
    {
      isomorphic <- sapply(lapply(not, graph.adjacency, mode="undirected", diag=F), graph.isomorphic, graph.adjacency(ggm,  mode="undirected", diag=F))
      if(any(isomorphic)) {
        next
        if(verbose) cat(" nonisomorphic failed.")
        } else
        {
          check <- TRUE
          if(verbose) cat("  nonisomorphic done...")
        }
    } else(check <- TRUE)
    

    ############################################
    # verify positive - definiteness of adjmat #
    if (check == TRUE)
    {
      ggm <- negativedges(ggm, negpro=negpro) # put negative edges in the matrix
      pcm <- ggm2pcm(ggm, minpcor = minpcor, maxpcor = maxpcor, maxiter = maxiter[3], verbose = verbose)
      if(verbose) cat("  partial correlation matrix done")
      if(any(is.na(pcm)))
      {
        check <- FALSE
        if(verbose) cat("  positive-definiteness failed.")
        next
      }
      if(verbose) cat("  partial correlation matrix ok!")
    }
  }

  if(iter >= maxiter[1])
  {
    if(!is.null(not))
    {
      if(any(isomorphic)) warning(paste("A new network was not identified in", iter, "iterations"))
    }
    if(all(is.na(pcm))) warning("Positive-definiteness was not reached")
  }
  
  list("ggm"=ggm, "pcm"=pcm)
}















###############################
# OLD STUFF TO FIX EVENTUALLY #
# ###############################
# Topology_vec<-function(vec, plot.top=F, not=NULL)
# {
#   Topology(N=vec$N, m=vec$m, topology=vec$topology, negpro=vec$negpro, pw_ba=vec$pw_ba, minpcor=vec$minpcor, maxpcor=vec$maxpcor, plot.top=plot.top, not=not)
# }
# 
# 
# test.topology<-function(ggm)
# {
#   degs<-sna::degree(ggm)
#   xmin<-ifelse(min(degs)<1, 1, min(degs))
#   powerlaw<-unlist(power.law.fit(degs, implementation="plfit", xmin=xmin))
#   smallw<-smallworldness(ggm)
#   c(powerlaw, smallw)
# }
# 
# 
# test.cm<-function(cm)
# {
#   # summary of correlations
#   veccm<-vectorizeMatrix(cm)
#   vectom<-vectorizeMatrix(TOMsimilarity(cm, TOMType="signed", TOMDenom="min", verbose=F))
# 
#   output<-c("cor_mean"=mean(veccm), "cor_sd"=sd(veccm), "cor_max"=max(veccm),
#             "tom_mean"=mean(vectom), "tom_sd"=sd(vectom), "tom_max"=max(vectom))
#   output
# }
# 
# 
# 
# test.clon<-function(cm, assi)
# {
#   library(WGCNA)
#   library(Matrix)
#   tom<-TOMsimilarity(adjMat=cm, TOMType="signed", TOMDenom="min", verbose=F)
#   logic<-as.matrix(bdiag(lapply(assi, function(x) matrix(T, ncol=x, nrow=x))))
#   intra<-logic
#   diag(intra)<-F
#   inter<-!logic
#   
#   out<-c("mean_cor_intra"=mean(cm[intra]),
#          "max_cor_intra"=max(cm[intra]),
#          "mean_cor_inter"=mean(cm[inter]),
#          "max_cor_inter"=max(cm[inter]),
#          "mean_tom_intra"=mean(tom[intra]),
#          "max_tom_intra"=max(tom[intra]),
#          "mean_tom_inter"=mean(tom[inter]),
#          "max_tom_inter"=max(tom[inter]))
#  out
# }
# 
# 
# 
# 
# smallworldness<-function(ggm, B=50)
# {
#   #compute the small worldness of Humphries & Gurney (2008) (sigma)
#   #and the index of Telesford et al. (2011) (omega)
#   # transitivity of ggm
#   ggm<-graph.adjacency(ggm, mode="undirected", diag=F)
#   N<-vcount(ggm)
#   m<-ecount(ggm)
#   
#   clusttrg<-transitivity(ggm, type="global", isolates="zero")
#   lengthtrg<-average.path.length(graph=ggm, directed=F, unconnected=F)
#    
#   #generate B rnd networks with the same degree distribution of ggm
#   deg.dist<-igraph::degree(ggm, mode="all", loops=F)
#   rndggm<-lapply(1:B, function(x)degree.sequence.game(deg.dist, method="simple.no.multiple"))
#   # compute the average (global) clustering coefficient over the B random networks
#   clustrnd<-mean(sapply(rndggm, transitivity, type="global", isolates="zero"))
#   
#   # compute the average shortest path length in random networks, the shortest path
#   # length among unconnected nodes is computed as N, i.e., 1 plus the max possible path length
#   lengthrnd<-mean(sapply(rndggm, average.path.length, directed=F, unconnected=F))
#   
#   # generate a regular lattice with the same number of nodes and edges
#   # this uses a different procedure then Telesford et al. (2011)
#   nei<-ceiling(m/N)
#   latt<-watts.strogatz.game(dim=1, size=N, nei=nei, p=0)
#   mdiff<-ecount(latt)-m
#   latt<-as.matrix(get.adjacency(latt, type="both"))
#   lattv<-vectorizeMatrix(latt)
#   # remove the edges in excess
#   lattv[lattv==1][sample(x=1:sum(lattv),size=mdiff, replace=F)]<-0
#   # rebuild the adjacency matrix
#   latt_mat<-matrix(ncol=N, nrow=N)
#   latt_mat[upper.tri(latt_mat)]<-lattv
#   latt_mat<-sna::symmetrize(latt_mat, rule="upper")
#   diag(latt_mat)<-0
#   # compute clustering coefficient and shortest path length for lattice graph
#   latt_final<-graph.adjacency(latt_mat, mode="undirected", diag=F)
#   clustlatt<-transitivity(latt_final, type="global", isolates="zero")
#   
#   if(clustrnd==0) return(c("sigma"=0, "omega"=0))
#   
#   # compute humphries&gourney(2008) smallworld-ness
#   sigma<-(clusttrg/clustrnd)/(lengthtrg/lengthrnd)
#   
#   # compute Telesford et al's (2011) smallworld-ness
#   omega<-(lengthrnd/lengthtrg)-(clusttrg/clustrnd)
#   
#   c("sigma"=sigma, "omega"=omega)
# }
# 
# 
# 
# 
# # convert the adjacency matrix ggm into a partial correlation matrix pcm
# ggm2cm<-function(ggm, minpcor=.2,  maxpcor=.9, maxiter=1000)
# {
#   # Starting from an adjacency matrix (ggm), the function generates a partial correlation
#   # matrix (pcm) by replacing the ones in ggm with random numbers from minpcor to maxpcor.
#   # from the pcm, a correlation matrix (cm) is generated
#   # Finally it checks for positive-definiteness: if a positive definite cm cannot be
#   # reached after maxiter iterations, the function returns -1
#   # Parameters:
#   # -ggm (binary matrix): adjacency matrix, the output of function Topology
#   # -minpcor, maxpcor: the minimum and the maximum values of the partial correlations 
#   # -maxiter: maximum number of iterations allowed in searching for a positive-definite matrix.
# 
#   N<-ncol(ggm)
#   iter<-0
#   check<-F
#   
#   while(check==F)
#   {
#     pcr=runif(N*(N-1)/2,minpcor,1)
#     ld.ggm <- ggm[lower.tri(ggm)]
#     pcm <- matrix(0,ncol=N,nrow=N) 
#     pcm[lower.tri(pcm)] <- ld.ggm*pcr # put the uniform numbers in the partial cor matrix
#     pcm <- (pcm+t(pcm))/2 + diag(1,N) # this is the partial correlation matrix
#     suppressWarnings(cm<-pcor2cor(pcm)) # this is the correlation matrix
#     if(any(is.na(cm))) 
#     {next} else check<-is.positive.definite(cm) # TRUE # check for pd
#     iter<-iter+1
#     if(iter>=maxiter & check==F)
#     {
#       warning(paste("Positive-definiteness was not reached after ", maxiter, " iterations"))
#       return(NA)
#     }
#   }
#   return(list("cm"=cm, "pcm"=pcm))
# }
