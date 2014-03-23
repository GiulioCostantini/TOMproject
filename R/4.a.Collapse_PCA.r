# scree test and parallel analysis in one command
collapse_PCA<-function(data, N = NULL, assignments, around = 5, rep = 1000, quantile = .95, rotate = "oblimin", stoppingrule = c("easystop", "parallel", "optimal"))
{  
  module <- list()
  if("easystop" %in% stoppingrule)
  {
    pca <- principal(data, nfactors = N, rotate = rotate)  
    load <- pca$loadings
    # assign each node to the module according to its highest loading
    mod <- (abs(load) == apply(abs(load), 1, max))
    
    # solve (the rare possibility of) ties
    ties <- apply(mod, 1, sum)
    if (any(ties > 1))
    {
      select <- sapply(ties[ties > 1], function(x) sample(1:x, 1))
      for(i in 1:sum(ties > 1))
      {
        mod[ties > 1, ][i,][mod[ties > 1, ][i,] == TRUE][-select[i]] <- FALSE
      }
    }
    module[[length(module) + 1]] <- matrix(1:nrow(t(mod)), ncol = ncol(t(mod)), nrow = nrow(t(mod)), byrow = FALSE)[t(mod)]
    names(module)[length(module)]<-"PCA_easystop"
  } 
  
  if("parallel"%in%stoppingrule)
  {
    ev <- eigen(cor(data, use="pairwise.complete.obs")) # get eigenvalues
    pa <- nFactors::parallel(subject=nrow(data),var=ncol(data), rep=rep, quantile=quantile, model="components")
    nS <- nScree(x=ev$values, aparallel=pa$eigen$qevpea)
    nfac<-nS$Components$nparallel
    pca<-principal(data, nfactors = nfac, rotate=rotate)  
    load<-pca$loadings
    mod<-(abs(load)==apply(abs(load), 1, max))
    
    # solve (the rare possibility of) ties
    ties <- apply(mod, 1, sum)
    if (any(ties > 1))
    {
      select <- sapply(ties[ties > 1], function(x) sample(1:x, 1))
      for(i in 1:sum(ties > 1))
      {
        mod[ties > 1, ][i,][mod[ties > 1, ][i,] == TRUE][-select[i]] <- FALSE
      }
    }
    
    module[[length(module)+1]]<-matrix(1:nrow(t(mod)), ncol=ncol(t(mod)), nrow=nrow(t(mod)), byrow=F)[t(mod)]
    names(module)[length(module)]<-"PCA_parallel"
  }
  
  if("optimal" %in% stoppingrule)
  {
    optimod <- NA
    # first check whether a perfect partition has already been recovered
    if("easystop" %in% stoppingrule | "parallel" %in% stoppingrule)
    {
      for(i in 1:length(module))
      {
        if(adjustedRandIndex(assignments, module[[i]]) == 1)
        {
          optimod <- module[[i]]
          break
        }
      }
    }
    if(any(is.na(optimod)))
    {
      optimod <- list()
      # nf is a vector of the numbers of factors considered
      nf <- (N-around):(N+around)
      nf <- nf[nf > 1 & nf < ncol(data)]
      if("easystop" %in% stoppingrule) nf <- nf[nf != N]

      for(i in nf)
      {
        pca <- principal(data, nfactors = i, rotate = rotate)  
        load <- pca$loadings
        # assign each node to the module according to its highest loading
        mod <- (abs(load) == apply(abs(load), 1, max))
        
        # solve (the rare possibility of) ties
        ties <- apply(mod, 1, sum)
        if (any(ties > 1))
        {
          select <- sapply(ties[ties > 1], function(x) sample(1:x, 1))
          for(i in 1:sum(ties > 1))
          {
            optimod[[i]][ties > 1, ][i,][optimod[[i]][ties > 1, ][i,] == TRUE][-select[i]] <- FALSE
          }
        }
        optimod[[length(optimod) + 1]] <- matrix(1:nrow(t(mod)), ncol = ncol(t(mod)), nrow = nrow(t(mod)), byrow = FALSE)[t(mod)]
        names(optimod)[length(optimod)] <- paste0("PCA_optimod_", i)
      }
      # pick the best module assignment (the first best module, in case of ties)
      adjrand <- lapply(optimod, adjustedRandIndex, assignments)
      optimod <- optimod[[match(max(unlist(adjrand)), unlist(adjrand))]]
    }
       
    module[[length(module) + 1]] <- optimod
    names(module)[length(module)]<-"PCA_optimal"
  }
  
  modulematrix<-data.frame(matrix(ncol=length(module), nrow=ncol(data)))
  for(i in 1:length(module)) modulematrix[,i]<-module[[i]]
  names(modulematrix)<-names(module)
  modulematrix
}

