Simulation1 <- function(pv)
{
  # generate the network (ggm and pcm)
  net <- Topology(N = pv$N, m = pv$N, topology = pv$topology, 
                  exact.m = pv$exact.m, force.connected = pv$force.connected,
                  negpro = pv$negpro, pw_ba = pv$pw_ba, p_ws = pv$p_ws,
                  minpcor = pv$minpcor, maxpcor = pv$maxpcor)
  pcm <- net$pcm
  
  # assign the clones
  assi <- Assign(N = pv$N, N2 = pv$N2, method_assign = pv$method_assign,
                 pwr = pv$pwr, shape1 = pv$shape1, shape2 = pv$shape2)
  assignments <- c()
  for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))
  
  # clone the network
  cloned <- Clone(pcm = pcm, assi = assi, n = pv$n,
                  furthersharedvar = pv$furthersharedvar, rndvar = pv$rndvar,
                  p = pv$p, keeporig = pv$keeporig)
  
  # collapse the network with PCA
  resultpca <- PCA(data = cloned$cloned, N = pv$N, assignments = assignments,
      around = pv$around, rep = 1000, quantile = .95, rotate = "oblimin",
      stoppingrule = c("easystop", "parallel", "optimal"))
  
  # evaluate network collapsed with PCA
  randpca <- apply(resultpca, 2, adjustedRandIndex, assignments)
  jaccpca <- apply(resultpca, 2, cluster_similarity, assignments, similarity = "jaccard")
  
  # collapse the network with all the variations of TOM
  a <- list()
  vls <- data.frame(matrix(ncol=5, nrow=0, dimnames=list(c(), c("tt", "td", "sq", "tpc", "mneg"))))
  mineigen <- c(0, .7, 1)
  i <- 0
  for(tt in c("unsigned", "signed"))
  {
    for(td in c("min", "mean", "squared"))
    {
      for(squared in c(TRUE, FALSE))
      {
       for(TOMPCA in c(TRUE, FALSE)) 
       {
         if(TOMPCA == TRUE)
         {
           for(mg in 1:length(mineigen))
           {
             i <- i+1
             a[[i]] <- TOM(data = cloned$cloned, TOMType = tt, TOMDenom = td,
                           squared = squared, TOMPCA = TOMPCA,
                           mineigen = mineigen[mg], verbose = FALSE)
             vls[i,] <- c(tt, td, squared, TOMPCA, mineigen[mg])
           }
         } else if(TOMPCA == FALSE)
         {
           i <- i+1
           a[[i]] <- TOM(data = cloned$cloned, TOMType = tt, TOMDenom = td,
                         squared = squared, TOMPCA = TOMPCA, verbose = FALSE)
           vls[i,] <- c(tt, td, squared, TOMPCA, NA)
         }
       }
      }
    }
  }
  
  # unlist the nested lists
  for(i in 1:length(a))
  {
    if(typeof(a[[i]])=="list")
      if(length(a[[i]])==1)
        a[[i]] <- a[[i]][[1]]
  }
  
  # stopping rules
  res<-list()
  for(i in 1:length(a))
  {
    res[[i]] <- StoppingRules(dist=a[[i]], assi = assi, data=cloned$cloned, method=pv$method,rule=c("dynamic", "optimal"))
  }
  
  # score all TOM x stopping rules with jaccard index
  randtom <- sapply(res, function(x) apply(x[[1]], 2, adjustedRandIndex, assignments))
  jacctom <- sapply(res, function(x) apply(x[[1]], 2, cluster_similarity, assignments, similarity = "jaccard"))
  
  # only optimal stopping rule
  randtom_opti <- randtom["optimal",]
  jacctom_opti <- jacctom["optimal",]
   
  # only best dynamic stopping rule (the jaccard index is ancillary to adjrand)
  randtom_bestdy <- apply(randtom[!rownames(randtom) %in% c("easystop", "optimal"),], 2, max)
  whichrnd <- apply(randtom[!rownames(randtom) %in% c("easystop", "optimal"),], 2, which.max)
  jacctom_bestdy <- jacctom[cbind(whichrnd, 1:48)] 

  names(randtom_opti) <- names(jacctom_opti) <- names(randtom_bestdy) <-
    names(jacctom_bestdy) <- paste0("TOM", 1:length(randtom_opti))
  
  # prepare the output of rand and jaccard for both TOM and PCA
  eval <- data.frame(matrix(nrow = (length(randtom_opti)+length(randpca)), ncol = 4,
                 dimnames = list(c(names(randtom_opti), names(randpca)), 
                                 c("randopti", "randbest", "jaccopti",
                                   "jaccbest"))))
  eval$randopti <- c(randtom_opti, randpca)
  eval$randbest  <- c(randtom_bestdy, randpca)
  eval$jaccopti <- c(jacctom_opti, jaccpca)
  eval$jaccbest <- c(jacctom_bestdy, jaccpca)
  
  # prepare the overall output
  output <- list("results" = eval, "TOMs" = vls, "pcm" = net$pcm, "ggm" = net$ggm, "assi" = assi)
}