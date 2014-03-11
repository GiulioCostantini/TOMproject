# scree test and parallel analysis in one command
collapse_PCA<-function(data, N = NULL, rep = 1000, quantile = .95, rotate = "oblimin", stoppingrule_pca = c("easystop", "parallel"))
{  
  require(nFactors)
  module <- list()
  
  if("easystop" %in% stoppingrule_pca)
  {
    pca <- principal(data, nfac = N, rotate = rotate)  
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
  
  if("parallel"%in%stoppingrule_pca)
  {
    ev <- eigen(cor(data, use="pairwise.complete.obs")) # get eigenvalues
    pa <- parallel(subject=nrow(data),var=ncol(data), rep=rep, quantile=quantile, model="components")
    nS <- nScree(x=ev$values, aparallel=pa$eigen$qevpea)
    nfac<-nS$Components$nparallel
    pca<-principal(data, nfac=nfac, rotate=rotate)  
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
  
  modulematrix<-data.frame(matrix(ncol=length(module), nrow=ncol(data)))
  for(i in 1:length(module)) modulematrix[,i]<-module[[i]]
  names(modulematrix)<-names(module)
  modulematrix
}
