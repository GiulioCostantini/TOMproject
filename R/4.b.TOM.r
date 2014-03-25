TOM <- function(data, TOMType="signed", TOMDenom = "mean", squared = FALSE, TOMPCA = FALSE, rotate="varimax", mineigen=c(0, .7, 1) )
{
  
  if(!TOMPCA)
  {
    out <- TOMnopca(adjMat=cor(data), TOMType=TOMType, TOMDenom = TOMDenom, squared = squared, paironly = FALSE)
  } else if(TOMPCA)
  {
    out <- TOMpca(data = data, rotate=rotate, mineigen=mineigen, TOMType=TOMType, TOMDenom = TOMDenom, squared = squared, verbose = TRUE)
  }
out


# 
# # systematic test of conditions that are common to WGCNA
# library(WGCNA)
# discord <- c()
# 
# for(tt in c("unsigned", "signed"))
# {
#   for(td in c("min", "mean"))
#   {
#     data <-matrix(rnorm(10000), ncol=100, nrow = 100)
#     am<- cor(data)
#     tm <- am
#     if(tt == "unsigned") tm <- abs(am)
#     a <- TOM(data = data, TOMType=tt, TOMDenom=td, squared=FALSE, TOMPCA=FALSE)
#     b <- TOMdist(adjMat=tm, TOMType=tt, TOMDenom=td, verbose = FALSE)
#     discord <- c(discord, max(abs(a-b)))
#     
#   }
# }
# 
# # test all the options: how symilar are the results
# 
# data <-matrix(rnorm(400), ncol=20, nrow = 20)
# i <- 0
# a <- list()
# vls <- data.frame(matrix(ncol=5, nrow=0, dimnames=list(c(), c("tt", "td", "sq", "tpc", "mneg"))))
# mineigen = c(0, .7, 1)
# 
# for(tt in c("unsigned", "signed"))
# {
#   for(td in c("min", "mean", "squared"))
#   {
#     for(squared in c(TRUE, FALSE))
#     {
#      for(TOMPCA in c(TRUE, FALSE)) 
#      {
#        if(TOMPCA == FALSE)
#        {
#          i <- i+1
#          a[[i]] <- TOM(data = data, TOMType=tt, TOMDenom=td, squared=squared, TOMPCA=TOMPCA)
#          vls[i,] <- c(tt, td, squared, TOMPCA, NA)
#          
#        } else if(TOMPCA == TRUE) {
#          res <- TOM(data = data, TOMType=tt, TOMDenom=td, squared=squared, TOMPCA=TOMPCA)
#          for(qq in 1:length(res))
#          {
#            i <- i+1
#            a[[i]] <- res[[qq]]
#            vls[i,] <- c(tt, td, squared, TOMPCA, mineigen[qq])
#          }
#        }
#        
#      }
#     }
#   }
# }
#
# library(ade4)
# mantest <- array(dim=c(length(a), length(a), 2))
# 
# for(i in 1:(length(a)-1))
# {
#   for(j in (i+1):length(a))
#   {
#     mt <- mantel.randtest(as.dist(a[[i]]), as.dist(a[[j]]), nrepet = 999)
#     mantest[j, i, ] <- c(mt$obs, mt$pvalue)
#   }
# }
# 
# colMeans(mantest[,,1], na.rm=T)
# mantest[,,2]
# results suggest that all matrices differ from each other, but that
# the differences are not soo marked that would make us suspicious of a
# mistake somewhere

}



TOMpca<-function(data, rotate="varimax", mineigen=c(0, .7, 1), TOMType = "signed", TOMDenom="mean", squared = FALSE, verbose = TRUE)
{
  N2 <- ncol(data)
  TOMeigmat <- list()
  for(i in 1:length(mineigen))
  {
    TOMeigmat[[i]] <- matrix(0, ncol=N2, nrow = N2)
  }
  
  for(j in 1:(N2-1))
  {
    for(i in (j+1):N2)
    {
      pcadata <- psych::principal(r = data[, -c(i,j)], nfactors = (N2-2), rotate = rotate, scores = TRUE)
      newdata <- list()
      for(k in 1:length(mineigen))
      {
        newdata[[k]]<-cbind(data[,c(i,j)], pcadata$scores[,pcadata$values > mineigen[k]])
        TOMeigmat[[k]][i, j] <- TOMnopca(adjMat = cor(newdata[[k]]), paironly = TRUE, a = 1, b = 2, TOMType = TOMType, TOMDenom = TOMDenom, squared = squared)  
      }
      if(verbose) print(paste(i,j))
    }
  }
  TOM <- lapply(TOMeigmat, sna::symmetrize, rule = "lower")
  names(TOM) <- paste0("Eigenvalues_", mineigen)
  TOM
}

TOMnopca <- function(adjMat, TOMType="signed", TOMDenom = "mean", squared = TRUE, paironly = FALSE, a = NULL, b = NULL)
{
  if(TOMType == "unsigned") adjMat <- abs(adjMat)
  
  if(!paironly)
  {
    #diag(adjMat) <- 0
    N2 <- ncol(adjMat)
    degs <- sna::degree(abs(adjMat), gmode="graph")
    num <- den <- matrix(ncol=N2, nrow = N2)
    if(TOMDenom == "mean") Tfun <- mean else if(TOMDenom == "min") Tfun <- min
    if(squared) nummat <- adjMat^2*sign(adjMat) else nummat <- adjMat
    
    for(j in 1:(N2-1))
    {
      for(i in (j+1):N2)
      {
        num[i,j] <- abs(sum(adjMat[i, -c(i,j)]*adjMat[j, -c(i,j)]) + nummat[i,j])
        if(TOMDenom %in% c("mean", "min"))
        {
          den[i,j] <- Tfun(degs[c(i,j)]) + 1 - abs(adjMat[i,j])
        } else if(TOMDenom == "squared") {
          den[i,j] <- sum(apply(abs(adjMat[-c(i,j),c(i,j)]), 1, max)^2)+1
        }

      }
    }
    
    num <- sna::symmetrize(num, rule = "lower")
    den <- sna::symmetrize(den, rule = "lower")
    diag(den) <- 1
    diag(num) <- 1
    TOMs <- num/den
    TOMd <- 1-TOMs
    
    } else if(paironly) {
      if(a == b) return(0)
    # compute TOMsimilarity for only two variables
    if(!squared)
    {
      num <- sum(apply(adjMat[c(a,b),-c(a, b)],2,prod)) + adjMat[a, b]
      num <- abs(num)
    } else if(squared){
      num <- sum(apply(adjMat[c(a,b),-c(a, b)],2,prod)) + adjMat[a, b]^2 * sign(adjMat[a, b])
      num <- abs(num)
    }

    if(TOMDenom == "min")
    {
      den <- min(c(sum(abs(adjMat[a, -a])), sum(abs(adjMat[b, -b])))) + 1 - abs(adjMat[a, b])
    }
    if(TOMDenom == "mean")
    {
      den <- mean(c(sum(abs(adjMat[a, -a])),sum(abs(adjMat[b, -b])))) + 1 - abs(adjMat[a, b])
    }
    if(TOMDenom == "squared")
    {
      den <- sum(apply(abs(adjMat[-c(a,b),c(a,b)]),1,max)^2)+1
    }
    TOMd <- 1 - num/den
  }
  
  TOMd
  
  # 
  # 
  # # systematic test of conditions that are common to WGCNA
  # library(WGCNA)
  # discord <- c()
  # 
  # for(tt in c("unsigned", "signed"))
  # {
  #   for(td in c("min", "mean"))
  #   {
  #     am <- cor(matrix(rnorm(10000), ncol=100, nrow = 100))
  #     tm <- am
  #     if(tt == "unsigned") tm <- abs(am)
  #     a <- TOM(adjMat=am, TOMType=tt, TOMDenom=td, squared=FALSE, paironly=FALSE)
  #     b <- TOMdist(adjMat=tm, TOMType=tt, TOMDenom=td, verbose = FALSE)
  #     discord <- c(discord, max(abs(a-b)))
  #     
  #   }
  # }
  # 
  # # systematic test of paironly
  # am <- cor(matrix(rnorm(10000), ncol=100, nrow = 100))
  # 
  # discord_di <- discord_od <- c()
  # a <- list()
  # b <- list()
  # counter <- 0
  # for(tt in c("unsigned", "signed"))
  # {
  #   for(td in c("min", "mean", "squared"))
  #   {
  #     for(sq in c(FALSE, TRUE))
  #     {
  #       counter <- counter+1
  #       a[[counter]] <- TOM(adjMat=am, TOMType=tt, TOMDenom=td, squared=sq, paironly=FALSE)
  #       
  #       # do the same pairwise
  #       b[[counter]] <- a[[counter]]
  #       b[[counter]][,] <- 0
  #       for(i in 1:ncol(b[[counter]]))
  #       {
  #         for(j in 1:ncol(b[[counter]]))
  #         {
  #           b[[counter]][i,j] <- TOM(adjMat=am, TOMType=tt, TOMDenom=td, squared=sq, paironly=TRUE,a=i, b = j)
  #         }
  #       }
  #       diff <- abs(a[[counter]]-b[[counter]])
  #       discord_od <- c(discord_od, max(diff - diag(diag(diff))))  
  #       discord_di <- c(discord_di, max(diag(diff)))
  #     }
  #   }
  # }
  # discord_od
  # discord_di
}
