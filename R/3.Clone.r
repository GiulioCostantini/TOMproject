##########################################
# FUNCTIONS TO DEFINE GENERATOR NETWORK ##
##########################################
pcm2data <- function(pcm, n=10000)
{
  N <- ncol(pcm)
  # cholesky decomposition of the correlation matrix.
  cm.chol <- chol(pcor2cor(pcm))
  #normal with mu = 0, sd = sd
  z <- matrix(rnorm(n*N),ncol = n)
  # Generate the data by multiplying the cholesky decomposed matrix
  # to the random matrix z
  y <- t(t(cm.chol) %*% z)
  scale(y)
}


Clone <- function(pcm, assi, method_clo = c("copy", "clone", "pclone"), n = 10000, fitvar = .5, furthersharedvar = .1, uniqueness = .2, p = .5)
{
  N <- ncol(pcm)
  if(length(assi)!= N) stop("The assignment entered is not compatible with the size of the correlation matrix")
  # generate an n*N data matrix encoding the same correlation structure
  y <- pcm2data(pcm, n)
  # generate the copied nodes
  if(method_clo == "copy") output <- copy(y = y, assi = assi, fitvar = fitvar, uniqueness = uniqueness)
  # generate the cloned nodes
  if(method_clo == "clone") output <- fclone(y = y, assi = assi, fitvar = fitvar, furthersharedvar = furthersharedvar, uniqueness = uniqueness)
  # generate the pclone nodes
  if(method_clo == "pclone") output <- pclone(y = y, assi = assi, fitvar = fitvar, furthersharedvar = furthersharedvar, uniqueness = uniqueness, p = p)
  
  output
}

copy <- function(y, assi, fitvar = .5, uniqueness = .5)
{
  # the clones are identical to the original node, plus some random noise
  # the portion of random variance added to each clone is defined by uniqueness.
  
  N <- ncol(y)
  N2 <- sum(assi)
  n <- nrow(y)
  # cloned/copied nodes are the original nodes plus random noise 
  y.cloned <- matrix(ncol = N2, nrow = n)
  
  for(i in 1:N)
  {
    shared <- as.vector(scale(y[,i])*fitvar)
    unique <- scale(matrix(rnorm(assi[i]*n), ncol = assi[i]))*uniqueness
    from <- sum(assi[0:(i-1)])+1
    to <- sum(assi[1:i])
    y.cloned[,from:to] <- shared + unique
  }
  scale(y.cloned)
}

fclone <- function(y, assi, fitvar = .5, furthersharedvar = .1, uniqueness = .3, p = .8)
{
  # Clone the columns of the data matrix in a way that preserves (completely or partially) their relations
  #     with all the other nodes in the network.
  
  # INPUT: 
  # y = original data matrix.
  # assi = assignment of the clones to each variable
  # fitvar = portion of the variance shared by all the clones of a node,
  #     which is related to other nodes (and clones of other nodes) in the network.
  # furthersharedvar = further portion of variance shared by the clones of a node
  #     which is not related to other nodes (and clones of other nodes) in the network,
  #     but the original node.
  # uniqueness = portion of variance unique to each clone (i.e., nonshared with
  #      other nodes, clones of other nodes and clones of the same node).
  # Ideally fitvar + furthersharedvar + uniqueness sum up to one. If this is not
  #     the case, their are rescaled such as their sum is one.
  # Ideally uniqueness should be higher than zero, to ensure positive-definiteness.
   
  # PROCEDURE:
  # A) The variance of the original node which is shared with other nodes in the
  #       network is saved in a variable, "fitval"
  # 
  # B) Each of the original nodes is sequentially cloned. The clones of a node are
  #     built by re-assembling the variance in fitvar plus other variance
  #     according to parameters fitvar, furthersharedvar, and uniqueness.
  #     each clone is assembled as fitvar*fitval + furthersharedvar * furthersharedval +
  #     uniqueness * rndval.
  
  N <- ncol(y)
  N2 <- sum(assi)
  n <- nrow(y)
  assignments <- c()
  for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))
  y.cloned <- matrix(ncol = sum(assi), nrow = n)
    
  for(i in 1:N)
  {
    # save the variance of the original node shared with other nodes (fitval)
    fit <- lm(y[,i]~cbind(y[,(1:ncol(y))[-c(1:i)]], y.cloned[,assignments<i]))
    fitval <- fit$fitted.values # this is the shared variance of the clones

    # define the random variable "furthersharedval" and clean it from
    # traces of correlations with fitval and with all the other nodes
    fsv <- rnorm(n)
    furthersharedval <- scale(lm(fsv~cbind(y[,(1:ncol(y))], fitval, y.cloned[,assignments<i]))$residuals)
    
    # define some random variance for each clone, cleaned from traces of correlations with
    # fitval, furthersharedval and with all the other nodes.
    rndval <- matrix(rnorm(n*assi[i]), ncol = assi[i])
    for(j in 1:ncol(rndval)) rndval[,j] <- scale(lm(rndval[,j]~cbind(y[,(1:ncol(y))], y.cloned[,assignments<i], furthersharedval))$residuals)
    
    # Assemble the clones of the ith node
    y.cloned[,assignments == i] <- fitvar*fitval+as.vector(furthersharedval)*furthersharedvar+rndval*uniqueness
  }
  scale(y.cloned)
}


pclone <- function(y, assi, fitvar = .5, furthersharedvar = .1, uniqueness = .3, p = .8)
{
  # Clone the columns of the data matrix in a way that preserves (completely or partially) their relations
  #     with all the other nodes in the network.
  
  # INPUT: 
  # y = original data matrix.
  # assi = assignment of the clones to each variable
  # fitvar = portion of the variance shared by all the clones of a node,
  #     which is related to other nodes (and clones of other nodes) in the network.
  # furthersharedvar = further portion of variance shared by the clones of a node
  #     which is not related to other nodes (and clones of other nodes) in the network,
  #     but the original node.
  # uniqueness = portion of variance unique to each clone (i.e., nonshared with
  #      other nodes, clones of other nodes and clones of the same node).
  # Ideally fitvar + furthersharedvar + uniqueness sum up to one. If this is not
  #     the case, their are rescaled such as their sum is one.
  # Ideally uniqueness should be higher than zero, to ensure positive-definiteness.
  
  # PROCEDURE:
  # A) The variance of the original node which is shared with other nodes in the
  #       network is saved in a variable, "fitval"
  # 
  # B) Each of the original nodes is sequentially cloned. The clones of a node are
  #     built by re-assembling the variance in fitvar plus other variance
  #     according to parameters fitvar, furthersharedvar, and uniqueness.
  #     each clone is assembled as fitvar*fitval + furthersharedvar * furthersharedval +
  #     uniqueness * rndval.
  N <- ncol(y)
  N2 <- sum(assi)
  n <- nrow(y)
  assignments <- c()
  for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))
  y.cloned <- matrix(ncol = sum(assi), nrow = n)
  
  for(i in 1:N)
  {
    fitval <- matrix(nrow = n, ncol = assi[i])
    for(j in 1:assi[i])
    {
      # select a subset of predictors in the ggm
      # each predictor is selected with probability p
      set <- ncol(y) - i + sum(assignments < i)
      subset <- sample(c(TRUE, FALSE), size = set, replace = TRUE, prob= c(p, 1-p))
      
      # save the variance of the original node shared with other nodes (fitval)
      fit <- lm(y[,i]~cbind(y[,(1:ncol(y))[-c(1:i)]], y.cloned[,assignments<i])[,subset])
      fitval[,j] <- fit$fitted.values # this is the shared variance of the clones
    }
    
    # define the random variable "furthersharedval" and clean it from
    # traces of correlations with fitval and with all the other nodes
    fsv <- rnorm(n)
    furthersharedval <- scale(lm(fsv~cbind(y[,(1:ncol(y))], fitval, y.cloned[,assignments<i]))$residuals)
    
    # define some random variance for each clone, cleaned from traces of correlations with
    # fitval, furthersharedval and with all the other nodes.
    rndval <- matrix(rnorm(n*assi[i]), ncol = assi[i])
    for(j in 1:ncol(rndval)) rndval[,j] <- scale(lm(rndval[,j]~cbind(y[,(1:ncol(y))], y.cloned[,assignments<i], fitval, furthersharedval))$residuals)
    
    # Assemble the clones of the ith node
    y.cloned[,assignments == i] <- fitvar*fitval+as.vector(furthersharedval)*furthersharedvar+rndval*uniqueness
    
  }
  
  scale(y.cloned)
  
  #  # check that on average the reproduced cor matrix is very very similar to the original one
  #     cm <- cor(y)
  #     cm_cloned <- cor(y.cloned)
  #     cm_reprod <- matrix(ncol = N, nrow = N)
  #     for(i in 1:N)
  #     {
  #       for(j in 1:N)
  #       {
  #         submat <- cm_cloned[assignments == i, assignments == j]
  #         if(i == j) submat <- submat[lower.tri(submat)]
  #         cm_reprod[i,j] <- mean(submat)
  #       }
  #     }
  #     round(cm_reprod,2)
  #     round(cm,2)
  #     qgraph(cm, layout = "spring")
  #     require(adegenet)
  #     qgraph(cm_cloned, layout = "spring", labels = assignments, colors = fac2col(assignments))  
  
}
