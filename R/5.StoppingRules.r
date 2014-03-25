StoppingRules <- function(dist, data, method = "average", rule = c("static", "dynamic", "optimal"))
{
  tree <- hclust(d=as.dist(dist), method=method)
  clusterings <- data.frame(matrix(nrow = ncol(data), ncol =0))
  values <- data.frame(matrix(nrow = 12, ncol=0, dimnames=list(c(
    "index", "value", "Hybrid", "pamStage", "pamRespectsDendro", "useMedioids",
    "respectSmallClusters", "maxPAMdist", "x", "deepSplit", "dmax", "gmin"), c())))
  
  if("static" %in% rule)
  {
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
  }
  if("dynamic" %in% rule)
  {
    # cutreeDynamicTree
    clusterings[,"cutTree_1"] <- dynamicTreeCut::cutreeDynamicTree(dendro=tree, deepSplit=TRUE, minModuleSize=2)
    values["index","cutTree_1"] <- "cutTree_1"
    values["deepSplit","cutTree_1"] <- TRUE
    values["Hybrid", "cutTree_1"] <- FALSE
    clusterings[,"cutTree_2"] <- dynamicTreeCut::cutreeDynamicTree(dendro=tree, deepSplit=FALSE, minModuleSize=2)
    values["index","cutTree_2"] <- "cutTree_2"
    values["deepSplit","cutTree_2"] <- FALSE
    values["Hybrid", "cutTree_2"] <- FALSE
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
              for(maxPamDist in c("default", "zero"))
              {
                for(i in 1:length(x))
                {
                  
                  nm <- paste0("cutTreeHyb_", counter)
                  if(maxPamDist)
                  {
                    # leave the default value for maxPamDist
                    clusterings[, nm] <- cutreeHybrid(dendro=tree, distM=dist, minClusterSize=1,
                                                      maxCoreScatter=dmax[i],minGap=gmin[i], 
                                                      pamStage= TRUE, pamRespectsDendro=pamRespectsDendro,
                                                      useMedoids=useMedoids, respectSmallClusters=respectSmallClusters,
                                                      verbose = 0)$labels
                  } else {
                    # give value zero to maxPamDist
                    clusterings[, nm] <- cutreeHybrid(dendro=tree, distM=dist, minClusterSize=1,
                                                      maxCoreScatter=dmax[i],minGap=gmin[i], 
                                                      pamStage= TRUE, pamRespectsDendro=pamRespectsDendro,
                                                      useMedoids=useMedoids, respectSmallClusters=respectSmallClusters,
                                                      maxPamDist = 0, verbose = 0)$labels
                  }
                  
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
  }
  if("optimal" %in% rule)
  {
    # easystop
    values["index","easystop"] <- "easystop"
    clusterings[,"easystop"] <- cutree(tree, k = length(assi))
    
    # optimal
    values["index", "optimal"] <- "optimal"
    optclust <- rep(1, sum(assi))
    rand <- adjustedRandIndex(assignments, optclust)
    
    assignments <- c()
    for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))
    
    for(i in 2:sum(assi))
    {
      clust <- cutree(tree, k = i)
      rand2 <- adjustedRandIndex(assignments, clust)
      if(rand2 > rand)
      {
        optclust <- clust
        rand <- rand2
      }
    }
    clusterings[,"optimal"] <-optclust
    
  }
  
  
  # when an object is not assigned by dynamic tree cut, it is turned into zero
  # we define instead a single-object cluster
  fixclusterings <- function(x)
  {
    if(sum(x == 0) >= 1)
    {
      frm <- max(x)+1
      to <- frm + sum(x == 0) -1
      x[x == 0] <- seq(from = frm, to = to, by = 1)
    }
    x
  }
  
  clusterings <- apply(clusterings, 2, fixclusterings)
  
  list("clusterings" = clusterings, "values" = values)  
}

