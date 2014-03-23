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

clean <- function(x, y = NULL, itself = TRUE)
{
  # clean a set of variables (x) from the variance shared with 
  # another set of variables (y). If itself = TRUE, the x variables
  # are also cleaned from each other.
  # if y = NULL, cleans only the x variables from one another.
  
  if(is.null(dim(x)[2])) x <- matrix(x, ncol=1)
 # clean from y
  if(!is.null(y))
  {
    if(is.null(dim(y)[2])) y <- matrix(y, ncol=1)
    y <- y[,!apply(apply(y, 2, is.na),2,any)] #remove NA columns from Y
    fit <- lm(x~y)
    x_pure <- fit$residuals
  } else x_pure <- x
  # clean from themselves, using a recursive call
  if(itself & ncol(x) > 1)
  {
    for(i in 1:ncol(x))
    {
      x_pure[,i] <- clean(x = x_pure[,i], y = x_pure[,-i], itself = FALSE)
    }
  }
  x_pure
}


Clone <- function(pcm, assi, n = 1000, furthersharedvar = 1, rndvar = 0, p = 1, keeporig = FALSE)
{
  N <- ncol(pcm)
  if(length(assi)!= N) stop("The assignment entered is not compatible with the size of the correlation matrix")
  # generate an n*N data matrix encoding the same correlation structure
  y <- pcm2data(pcm, n)
  output <- list()
  output$orig <- y
  N2 <- sum(assi)
  origvar <- 1-rndvar
  assignments <- c()
  for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))
  y.cloned <- matrix(ncol = sum(assi), nrow = n)
  if(keeporig)
  {
    orig <- 1
    for(i in 1:(length(assi)-1)) orig[i+1] <- sum(assi[1:i])+1
  }
  
  for(i in 1:N)
  {
    # decompose the variance of the original node in two parts, the fitted
    # values of variance R^2 and the residuals, of variance 1-R^2
    set <- 1:N
    fit <- lm(y[,i]~cbind(y.cloned[,assignments%in%set[set<i]], y[,set>i]))
    fitval <- fit$fitted.values
    residualval <- fit$residuals
    R2 <- var(fitval)
    
    if(p == 1)
    {
      fitmat <- matrix(rep(fitval, assi[i]), ncol = assi[i], nrow = n, byrow = FALSE)
    } else
    {
      fitmat <- matrix(ncol = assi[i], nrow = n)
      for(j in 1:assi[i])
      {
        # identify the subset of predictors, given p
        set <- (1:N)[-i]
        select <- sample(c(TRUE, FALSE), size = length(set), replace = TRUE, prob= c(p, 1-p))
        subset <- set[select]
        # ensure that at least one predictor is preserved
        if(length(subset) == 0) subset <- sample(set, 1)
        
        # decompose the variance of the original node in two parts, the fitted
        # values of variance R^2 and the residuals, of variance 1-R^2
        fit <- lm(y[,i]~cbind(y[,subset[subset>i]], y.cloned[,assignments%in%subset[subset<i]]))
        fitmat[,j] <- fit$fitted.values
      }
      # the maximum theoretical variance in fitmat is R2. To make all the
      # clones have the same variance, the variance of all the columns of
      # fitmat is transformed into R2, by adding a random component to each,
      # of variance = "extravariance".
      
      tol = 1e-15
      extravariance <- R2 - apply(fitmat,2,var)
      tocorrect <- (1:length(extravariance))[extravariance > tol]
      if(length(tocorrect) > 0)
      {
        extravariance <- extravariance[tocorrect]
        extraSDs <- sqrt(extravariance)
        rndmat <- matrix(rnorm(n*length(tocorrect)), ncol=length(tocorrect), nrow = n)
        rndmat <- clean(rndmat, cbind(y[,set>=i], y.cloned[,assignments < i], fitmat))
        rndmat <- scale(rndmat)[,]
        rndmat <- t(t(rndmat)*extraSDs)
        fitmat[,tocorrect] <- fitmat[,tocorrect] + rndmat
      }
    }
    
    # Further decompose the residual variance
    # A portion of the residual variance of furthersharedvar
    # is kept intact in the clones of the same node
    intactresiduals <- residualval*sqrt(furthersharedvar)
    
    # While another portion is replaced. The random variable "newresiduals"
    # is obtained by replacing the remaining residual variance with a new 
    # component of variance, shared by all nodes
    newresiduals <- matrix(rnorm(n*assi[i]), ncol = assi[i])
    set <- 1:N
    newresiduals <- scale(clean(newresiduals, cbind(y[,set>=i], y.cloned[,assignments < i], fitmat),
                                itself = TRUE))
    newresiduals <- newresiduals*sqrt(1-furthersharedvar)*sqrt(1-R2)
    
    # Assemble the clones of the ith node
    y.cloned[,assignments == i] <- (fitmat + intactresiduals + newresiduals)
    
    # if the original node should be kepts, include it
    if(keeporig) y.cloned[,orig[i]] <- y[,i]
  }
  
  # rndvar: completely random noise, different for each node
  unique <- scale(matrix(rnorm(n*N2), ncol = N2))
  y.cloned <- scale(y.cloned)*sqrt(origvar) + unique*sqrt(rndvar)
  
  output$cloned <- scale(y.cloned)
  output
}