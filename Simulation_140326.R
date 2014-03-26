# load the required packages, including the modded dynamicTreeCut
library(devtools)
#install_github("TOMproject", "GiulioCostantini")
install.packages("dynamicTreeCut_1.60-1.zip", repos = NULL)
library(TOMproject)
library(dynamicTreeCut)
library(parallel)

# load the parameters
param <- read.csv2("Parameters_140326.csv", as.is = TRUE)
param[,c("exact.m", "force.connected", "keeporig")] <- apply(param[,c("exact.m", "force.connected", "keeporig")], 2, as.logical)

rep <- 5
counter <- 1
myoutput1 <- list()
for(i in 1:nrow(param))
{
  for(j in 1:rep)
  {
    myoutput1[[counter]] <-tryCatch(Simulation1(param[i,]), error = function(q) NA)
    counter <- counter + 1
  }
}


# 1. find the best TOM #
#########################
result1 <- array(dim = c(sum(!is.na(myoutput1)), dim(myoutput1[[1]]$results)),
                 dimnames = list(c(), rownames(myoutput1[[1]]$results),
                          colnames(myoutput1[[1]]$results)))

counter <- 1
for(i in 1:length(myoutput1))
{
  if(!is.na(myoutput1[[i]])[1])
  {
    result1[counter,,] <- as.matrix(myoutput1[[i]]$results)
  }
  counter <- counter + 1
}

# check the general trend
randopti <- result1[,, "randopti"]
randbest <- result1[,, "randbest"] 

boxplot(randopti)
boxplot(randbest)

# check specific factors in TOM
library(stringr)
library(plotrix)
randoptiTOM <- randopti[,str_detect(colnames(randopti), "TOM")]
TOMIVs <- myoutput1[[1]]$TOMs
TOMIVs$mneg[is.na(TOMIVs$mneg)] <- -1

# TOMType
boxplot(randoptiTOM[,order(TOMIVs$tt)], xaxt="n", col=c(rep("red", 24), rep("green", 24)))
axis(side=1, at=1:ncol(randoptiTOM), labels=sort(TOMIVs$tt))

# TOMDenom
boxplot(randoptiTOM[,order(TOMIVs$td)], xaxt="n", col=c(rep("red", 16), rep("green", 16), rep("blue", 16)))
axis(side=1, at=1:ncol(randoptiTOM), labels=sort(TOMIVs$td))

# TOMDenom
boxplot(randoptiTOM[,order(TOMIVs$sq)], xaxt="n", col=c(rep("red", 24), rep("green", 24)))
axis(side=1, at=1:ncol(randoptiTOM), labels=sort(TOMIVs$sq))

# TOMpcas (x4)
boxplot(randoptiTOM[,order(TOMIVs$mneg)], xaxt="n", col=c(rep("red", 12), rep("green", 12), rep("blue", 12), rep("white", 12)))
axis(side=1, at=1:ncol(randoptiTOM), labels=sort(TOMIVs$mneg))


# which is the best TOM?
bestTOM <- which.max(apply(randoptiTOM, 2, mean))
boxplot(cbind(randoptiTOM[,bestTOM], randopti[,"PCA_optimal"]))



# parallelized version of the simulation: IT DOES NOT WORK!
sim1 <- function(x) tryCatch(Simulation1(x), error = function(q) NA)

parlist <- list()
for(i in 1:nrow(param)) parlist[[i]] <- param[i,]
cl <- makeCluster(getOption("cl.cores", 2))
myoutput <- parallel::parLapply(cl = cl, X=parlist, fun=sim1)



