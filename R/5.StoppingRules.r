StoppingRules <- function(dist, data, method = "average")
{
  require(dynamicTreeCut)
  
  tree <- hclust(d=as.dist(dist), method=method)
  clusterings <- data.frame(matrix(nrow = ncol(data), ncol =0))
  values <- data.frame(matrix(nrow = 12, ncol=0, dimnames=list(c(
    "index", "value", "Hybrid", "pamStage", "pamRespectsDendro", "useMedioids",
    "respectSmallClusters", "maxPAMdist", "x", "deepSplit", "dmax", "gmin"), c())))
  
  NBC <- NbClust_mod(data=t(data), diss = as.dist(dist), min.nc=1,
                     max.nc=ncol(data)-1, index="alllong", method = method)
  
  for(i in 1:ncol(NBC$Best.nc))
  {
    nc <- NBC$Best.nc[1,i]
    index <- NBC$Best.nc[2,i]
    nm <- colnames(NBC$Best.nc)[i]
    
    clusterings[,nm] <- cutree(tree, k=nc)
    values["index", nm] <- nm
    values["value", nm] <- index
  }
  
  # cutreeDynamicTree
  clusterings[,"cutTree_1"] <- dynamicTreeCut::cutreeDynamicTree(dendro=tree, deepSplit=TRUE, minModuleSize=2)
  values["index","cutTree_1"] <- "cutTree_1"
  values["deepSplit","cutTree_1"] <- TRUE
  values["Hybrid", nm] <- FALSE
  clusterings[,"cutTree_2"] <- dynamicTreeCut::cutreeDynamicTree(dendro=tree, deepSplit=FALSE, minModuleSize=2)
  values["index","cutTree_2"] <- "cutTree_2"
  values["deepSplit","cutTree_2"] <- FALSE
  values["Hybrid", nm] <- FALSE
  # cutreeHybrid
  x <- c(.01, .10, .19, .28, .37, .46, .55, .64, .73, .82, .91, .99)
  hmin <- min(tree$height)
  hmax <- hmin + .99*(max(tree$height)-hmin)
  dmax <- hmin+x*(hmax-hmin)
  gmin <- (1-x)*3/4*(hmax - hmin)
  counter <- 1
  
  for(pamStage in c(TRUE, FALSE))
  {
    if(pamStage == FALSE)
    {
      for(i in 1:length(x))
      {
        nm <- paste0("cutTreeHyb_", counter)
        clusterings[, nm] <-  cutreeHybrid(dendro = tree, distM = dist, minClusterSize = 1,
                                           maxCoreScatter = dmax[i], minGap = gmin[i], 
                                           pamStage = FALSE, verbose = FALSE)$labels
        values["index", nm] <- nm
        values["pamStage", nm] <- FALSE
        values["dmax", nm] <- dmax[i]
        values["gmin", nm] <- gmin[i]
        values["Hybrid", nm] <- TRUE
        
        if(i %in% 8:11)
          values["deepSplit", nm] <- i-8
        counter <- counter+1
        
      }
    }
    
    if(pamStage == TRUE)
    {
      for(pamRespectsDendro in c(TRUE, FALSE))
      {
        for(respectSmallClusters in c(TRUE, FALSE))
        {
          for(useMedoids in c(TRUE, FALSE))
          {
            for(maxPamDist in c(TRUE, FALSE))
            {
              for(i in 1:length(x))
              {
                
                nm <- paste0("cutTreeHyb_", counter)
                clusterings[, nm] <- cutreeHybrid(dendro=tree, distM=dist, minClusterSize=1,
                                                  maxCoreScatter=dmax[i],minGap=gmin[i], 
                                                  pamStage= TRUE, pamRespectsDendro=pamRespectsDendro,
                                                  useMedoids=useMedoids, respectSmallClusters=respectSmallClusters,
                                                  maxPamDist=maxPamDist)$labels
                values["Hybrid", nm] <- TRUE
                values["index", nm] <- nm
                values["pamStage", nm] <- TRUE
                values["dmax", nm] <- dmax[i]
                values["gmin", nm] <- gmin[i]
                values["pamRespectsDendro", nm] <- pamRespectsDendro
                values["useMedoids", nm] <- useMedoids
                values["respectSmallClusters", nm] <- respectSmallClusters
                values["maxPamDist", nm] <- maxPamDist
                if(i %in% 8:11)
                  values["deepSplit", nm] <- i-8
                counter <- counter+1
              }
            }
          }
        }
      }
    }
  }
  list("clusterings" = clusterings, "values" = values)  
}


