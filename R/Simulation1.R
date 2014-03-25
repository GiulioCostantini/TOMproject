# library(devtools)
# load_all("C:/Users/giulio.costantini/GitHub/dynamicTreeCut_Gmod/")
# 
# param <- read.csv2("ParametersSim1.csv", as.is = TRUE)
# param[,c("exact.m", "force.connected", "keeporig")] <- apply(param[,c("exact.m", "force.connected", "keeporig")], 2, as.logical)
# pv <- param[1,]
# 
# 
# Simulation1 <- function(pv)
# {
#   # generate the network (ggm and pcm)
#   net <- Topology(N = pv$N, m = pv$N, topology = pv$topology, 
#                   exact.m = pv$exact.m, force.connected = pv$force.connected,
#                   negpro = pv$negpro, pw_ba = pv$pw_ba, p_ws = pv$p_ws,
#                   minpcor = pv$minpcor, maxpcor = pv$maxpcor)
#   pcm <- net$pcm
#   
#   # assign the clones
#   assi <- Assign(N = pv$N, N2 = pv$N2, method_assign = pv$method_assign,
#                  pwr = pv$pwr, shape1 = pv$shape1, shape2 = pv$shape2)
#   assignments <- c()
#   for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))
#   
#   # clone the network
#   cloned <- Clone(pcm = pcm, assi = assi, n = pv$n,
#                   furthersharedvar = pv$furthersharedvar, rndvar = pv$rndvar,
#                   p = pv$p, keeporig = pv$keeporig)
#   
#   # collapse the network with PCA
#   resultpca <- PCA(data = cloned$cloned, N = pv$N, assignments = assignments,
#       around = pv$around, rep = 1000, quantile = .95, rotate = "oblimin",
#       stoppingrule = c("easystop", "parallel", "optimal"))
#   
#   # collapse the network with all the variations of TOM
#   a <- list()
#   vls <- data.frame(matrix(ncol=5, nrow=0, dimnames=list(c(), c("tt", "td", "sq", "tpc", "mneg"))))
#   mineigen <- c(0, .7, 1)
#   i <- 0
#   for(tt in c("unsigned", "signed"))
#   {
#     for(td in c("min", "mean", "squared"))
#     {
#       for(squared in c(TRUE, FALSE))
#       {
#        for(TOMPCA in c(TRUE, FALSE)) 
#        {
#          if(TOMPCA == TRUE)
#          {
#            for(mg in 1:length(mineigen))
#            {
#              i <- i+1
#              a[[i]] <- TOM(data = cloned$cloned, TOMType = tt, TOMDenom = td,
#                            squared = squared, TOMPCA = TOMPCA,
#                            mineigen = mineigen[mg])
#              vls[i,] <- c(tt, td, squared, TOMPCA, mineigen[mg])
#            }
#          } else if(TOMPCA == FALSE)
#          {
#            i <- i+1
#            a[[i]] <- TOM(data = cloned$cloned, TOMType = tt, TOMDenom = td,
#                          squared = squared, TOMPCA = TOMPCA)
#            vls[i,] <- c(tt, td, squared, TOMPCA, NA)
#          }
#        }
#       }
#     }
#   }
#   
#   bk <- a
#   
#   # unlist the nested lists
#   for(i in 1:length(a))
#   {
#     if(typeof(a[[i]])=="list")
#       if(length(a[[i]])==1)
#         a[[i]] <- a[[i]][[1]]
#   }
#   
#   # stopping rules
#   res<-list()
#   
#   for(i in 1:length(a))
#   {
#     res[[i]] <- StoppingRules(dist=a[[i]], data=cloned$cloned, method=pv$method,rule=pv[c("rule1", "rule2", "rule3")])
#   }
#   
#   
#   resulttom <- lapply(a, StoppingRules, data=cloned$cloned, method=pv$method,rule=pv[c("rule1", "rule2", "rule3")])
#   
#   # evaluation
#   Rand <- apply
#   
#   
#   # independent variables
#   
# }