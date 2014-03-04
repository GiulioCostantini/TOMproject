#######################################################################
# FUNCTIONS FOR ASSIGNING A NUMBER OF CLONES TO EACH GENERATING NODE #
#######################################################################
Assign<-function(N, N2, method_assign = c("unif", "pwr", "beta"), pwr = 1, shape1 = .1, shape2 = 1)
{
  # Assign a number of clones to each node
  # N: number of original nodes
  # N2: total number of clones
  # method_assign: the method for assignnig clones to the original nodes
  #       unif: same number of clones to each node
  #       pwr: assign nodes according to a power function, bounded from .1 to .9
  #             with exponent pwr
  #       beta: assign nodes accordnig to a beta distribution, with parameters
  #             shape1 (alpha) and shape2 (beta)
  #...: parameters for the specific assignment
  #     pwr (method_assign="pwr"): scalar, the exponent of the power distribution
  #     shape1 (method_assign="beta"): scalar, the alpha parameter of the beta distribution
  #     shape2 (method_assign="beta"): scalar, the beta parameter of the beta distribution
  
  if(length(method_assign) !=1 | !method_assign %in% c("unif", "pwr", "beta")) stop("Wrong type of assignment specification")
  if(method_assign == "unif") assi <- assign_unif(N, N2)
  if(method_assign == "pwr") assi <- assign_pwr(N, N2, pwr = pwr)
  if(method_assign == "beta") assi <- assign_beta(N, N2, shape1 = shape1, shape2 = shape2)
  
#   assignments <- c()
#   for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))
#   list("assi" = assi, "assignments" = assignments)
assi
}


assign_unif<-function(N, N2)
{
  # - N -> the number of variables
  # - N2-> the number of clones
  x <- rep(floor(N2/N), N)
  if(N2 %% N>0)
  {
    rem <- N2 %% N
    x[1:rem] <- x[1:rem]+1    
  }
  x[order(rnorm(length(x)))] 
}

assign_pwr <- function(N, N2, pwr = 0)
{
  # - N -> the number of variables
  # - N2-> the number of clones
  # - pwr-> the exponent
  #   pwr=0 -> uniform assignment
  #   pwr>>0 -> high asimmetry
  
  x <- seq(from = .1, to = .9, length.out = N)
  distr <- (x^pwr/sum(x^pwr))*(N2-N)
  distr.r <- round(distr)
  # if the number of clones is not right because of rounding errors,
  # add/delete a clone to the most cloned node
  if(sum(distr.r) != (N2-N))
  {
    missing <- (N2-N) - sum(distr.r)
    for(i in 1:abs(missing))
    {
      nd <- sample((1:N)[distr.r >= 1], 1)
      distr.r[nd] <- distr.r[nd] + sign(missing)
    }
   }
  distr.r <- distr.r+1 # to ensure that each node has at least one clone
  distr.r[order(rnorm(length(distr.r)))]  
}

assign_beta <- function(N, N2, shape1 = .1, shape2 = 1)
{
  # function to assign to each of the n variables a
  # certain number of clones.
  # - N -> the number of variables
  # - N2-> the number of clones
  # - shape1 and shape2 are the parameters of a beta distributions
  #   that is used to weigh the probabilities that a variable is
  #   cloned. For shape1>>1 and shape2<<1 you have a uniform
  #   distribution of clones to nodes, for shape1<<1 and shape2>>1
  #   you have that most nodes do not have clones at all and a few have
  #   all the clones.
  
  distr <- rbeta(N, shape1 = shape1, shape2 = shape2)
  distr <- distr/sum(distr)
  distr.r <- round(distr*(N2-N))
  if(sum(distr.r) != (N2-N))
  {
    missing <- (N2-N) - sum(distr.r)
    for(i in 1:abs(missing))
    {
      nd <- sample((1:N)[distr.r >= 1],1)
      distr.r[nd]<-distr.r[nd] + sign(missing)
    }
  }
  distr.r <- distr.r+1 # to ensure that each node has at least one clone
  distr.r[order(rnorm(length(distr.r)))]  
}

