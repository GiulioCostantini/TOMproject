\name{Clone}
\alias{Clone}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clone a network.
}
\description{
This function starts from an original network, represented by a partial correlation matrix (\code{pcm}), and from a vector (\code{assi}. Assi is the output of function \code{\link{Assign}}) and specifies, for each node in the pcm, how many clones will correspond to it. The result is a new network in which each node is a clone of a node in the original network. The output is not in the form of a network, but of a data matrix.
}
\usage{
Clone(pcm, assi, n = 1000, furthersharedvar = 1, rndvar = 0, p = 1, keeporig = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcm}{Matrix of numerics. Partial correlation matrix.
}
  \item{assi}{Vector of integers, same size as a dimension of \code{pcm}. The i-th element specifies the number of clones of the i-th node of the \code{pcm}.
}
  \item{n}{Integer. Desired number of observations.
}
  \item{furthersharedvar}{Numeric in [0,1]. Portion of the residual variance of the original node that is kept intact in the clones. See details.
}
  \item{rndvar}{Numeric in [0,1]. Portion of the variance of the clones that will be composed by random noise. See details.
}
  \item{p}{Numeric in [0,1]. Probability that a neighbor is considered in the process of clonation.
  }
  \item{keeporig}{Logical. Should the original nodes be preserved in the cloned network.
}
}
\details{
The function takes a partial correlation matrix in input and returns a data matrix and a matrix of clones.
The backbone of the function (i.e., if parameter \code{p}==1 and parameter \code{rndvar}==0) proceedes in this way:

1) Generate a data matrix from the pcm (\code{orig}), with M=0 and SD=1, using Cholesky decomposition (This matrix is returned in output).

2) Pick the i-th variable (starting from the first) in \code{orig} and decompose it in in two parts: The "systematic" one is the variance shared with other nodes in the network. The "residual" one is its unique variance. This split is performed using multiple regression and saving the fitted values ("systematic" variance) and the residuals ("residual" variance).

3) Compose \code{assi[i]} clones of the i-th variable, by assembling their variance in this way:

 (a) The systematic component, identical to the "systematic" variance of the original node (i.e., R^2 in the multiple regression of the original node).
 
 (b) The remaining variance, equal to err=1-R^2, is made by two parts.
 
   (b1) A shared residual component: this amounts to err*\code{furthershared} and is equal to the residual component of the original node, multiplied by sqrt((1-R^2)\code{*furthershared}) to get the desired variance. This furhter component is shared by all the clones of the same node.
   
   (b2) The remaining variance, err*(1-\code{furthershared}), is selected to be perfectly unique for each clone, i.e., completely linearly independent from each of the other nodes and clones in the network and from the clones of the same node. No linear dipendence is allowed, not even by chance.
   
4) Repeat the steps 2-3 for the next node, until all the nodes in \code{pcm} have been cloned. In decomposing the variance of the node, consider the original ndoes if they have not been cloned yet, their clones otherwise.

As a result a matrix of clones (\code{cloned}, see output) is returned, with \code{n} rows and \code{sum(assi)} columns.

This function is very generic, allowing different kinds of clones to be returned: we describe some important kinds of clones.

Duplicates) If parameter \code{p} = 1, parameter \code{rndvar} = 0, and parameter \code{furthersharedvar} = 1, the clones of a node are exact replicas of the node and therefore the off-diagonal elements of the resulting correlation network are identical to the corresponding ones in the original network. The parameter \code{furthersharedvar} = 1 determines that all the original variance of the node is kept in the clones.

Same Matrix) If parameter \code{p} = 1, parameter \code{rndvar} = 0, and parameter \code{furthersharedvar} = 0, the correlations among the clones have the two following properties: the correlation among two clones of two different nodes i and j is exactly identical to that of the original i and j. The correlation among two clones of the same nodes is the minimum possible such as the first property is obtained and is equal to the R^2 of the multiple regression at step 2.

Using parameter \code{furthersharedvar} it is possible to create a bridge between the ideal cases "Duplicates" and "Same Matrix". For instance, if \code{furthersharevar}=.5, then the 50\% of the error variance of the original nodes is preserved in their clones, therefore the off-diagonal elements of the resulting correlation matrix are equal to those of the original network, but the node share some variance that is irrelevant to strictly obtain this aim.

Parameter \code{p} allows to modify the multiple regression at Step 2, such as instead of considering all the predictors (i.e., neighbors of the node in the correlation matrix), each predictor is considered with probability \code{p}. If a predictor is selected, but it has been already cloned, its clones are considered instead. Notice that however at least one predictor is always guaranteed, even if zero predictors are extracted. Parameter \code{p} has been inspired by the partial duplication process described in Chung et al. (2003); however while in Chung's network each new node (corresponding to our clones) inherits only a subset of the neighbors of the original node, in our procedure each clone "is influenced" only by a subset of the neighbors of the original node. If the network is dense enough however, it is very likely that it will correlate also with the other neigbors (by spurious correlations). If parameter \code{p} is < 1, then the component (a) of the variance of each clone, as detemrmined at step 3, is less than R^2. The remaining variance is "filled" with random variance, perfecly unique for each clone.

The parameter \code{rndvar} allows to add further random variance to the clones. For instance, if parameter \code{rndvar} = .4, the 60\% of the variance of the clones will be the one defined in step 4, while the remaining 40\% is random. This random variance however is not cleaned from random covariance with other nodes in the network. It is meant to make the cloned matrices less "perfect" and closer to the real networks and also to make the task to recover the original network structure progressively higher.

Finally, if parameter \code{keeporig} = TRUE, the nodes in the original network are kept also in the cloned network. This step is done at the very end, by replacing one clone by node with the original variable in \code{orig}.
}

\value{A list of two matrices.
\code{orig} is a data matrix generated from \code{pcm} using Cholesky decomposition. It has \code{n} observations and as many variables as there are in \code{pcm}. A partial correlaiton matrix obtained by these data should reproduce \code{pcm}, even though with some minor discrepancies.
\code{cloned} is a matrix of cloned nodes, \code{n} observations by \code{sum(assi)} variables. Each i-th node in \code{pcm} has \code{assi[i]} clones in \code{cloned}.
}
\references{
Chung, F., Linyuan, L., Dewey, G., Galas, D. J. (2003), Duplication models for biological networks. Journal of Computational Biology, 10(5), p. 677-687}
\author{Giulio Costantini
}

\examples{
\dontrun{
require(adegenet)

# generate a network with 10 (15 edges in the original adjacency matrix)
# and then generate several cloned networks with 30 clones
N <- 10
N2 <- 30
net <- Topology(N=N, m=15, topology="BA", exact.m=F, force.connected=T, negpro=0, pw_ba=1, minpcor=.1, maxpcor=.3)
pcm <- net$pcm
assi <- Assign(N=N, N2=N2, method_assign="pwr", shape1= .0001, shape2 = 1000)
assignments <- c()
for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))

par(mfrow = c(4,2))
qgraph(pcm, layout = "spring")
title(sub = "original pcm")

# The extreme case "Duplicates"
duplicates <- Clone(pcm=pcm, assi=assi, n=1000, rndvar=0, furthersharedvar=1, p=1)
 # verify that the correlation of the clones with the original nodes are all = 1
  round(cor(duplicates$cloned, duplicates$orig)^2,2)
qgraph(cor(duplicates$cloned), layout = "spring", colors = adegenet::fac2col(assignments), labels = assignments)
title(sub = "duplicates")

# "Duplicates", with .50 noise
noisy50dup <- Clone(pcm=pcm, assi=assi, n=1000, rndvar=.5, furthersharedvar=1, p=1)
qgraph(cor(noisy50dup$cloned), layout = "spring", colors = adegenet::fac2col(assignments), labels = assignments)
title(sub = "duplicates with .50 noise")

# "Duplicates", with .90 noise
noisy90dup <- Clone(pcm=pcm, assi=assi, n=1000, rndvar=.9, furthersharedvar=1, p=1)
qgraph(cor(noisy90dup$cloned), layout = "spring", colors = adegenet::fac2col(assignments), labels = assignments)
title(sub = "duplicates with .90 noise")

# The extreme case "SameMatrix".
samemat <- Clone(pcm=pcm, assi=assi, n=1000, rndvar=0, furthersharedvar=0, p=1)
  # verify that the resulting collapsed correlation matrix is identical, besides the diagonal elements, which are the R^2.
  cm <- cor(samemat$orig)
  cm_cloned <- cor(samemat$cloned)
  cm_reprod <- matrix(ncol = N, nrow = N)
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      submat <- cm_cloned[assignments == i, assignments == j]
      if(i == j & assi[i] > 1) submat <- submat[lower.tri(submat)]
      cm_reprod[i,j] <- mean(submat)
    }
  }
  round(cm_reprod, 2)
  round(cm_reprod-cm, 4)

qgraph(cor(samemat$cloned), layout = "spring", colors = adegenet::fac2col(assignments), labels = assignments)
title(sub = "same matrix")

# a Mix of cases Duplicates and SameMatrix obtained by varying parameter furthersharedvar
dup50same50 <- Clone(pcm=pcm, assi=assi, n=1000, rndvar=0, furthersharedvar=.5, p=1)
  # verify that the resulting collapsed correlation matrix is identical, besides the diagonal elements, which are the higher than R^2.
  cm <- cor(dup50same50$orig)
  cm_cloned <- cor(dup50same50$cloned)
  cm_reprod <- matrix(ncol = N, nrow = N)
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      submat <- cm_cloned[assignments == i, assignments == j]
      if(i == j & assi[i] > 1) submat <- submat[lower.tri(submat)]
      cm_reprod[i,j] <- mean(submat)
    }
  }
  round(cm_reprod, 2)
  round(cm_reprod-cm, 4)

qgraph(cor(dup50same50$cloned), layout = "spring", colors = adegenet::fac2col(assignments), labels = assignments)
title(sub = ".5 duplicate .5 samemat")

# Use parameter p to consider only the 50\% of the predictors for each clone
net50 <- Clone(pcm=pcm, assi=assi, n=1000, rndvar=0, furthersharedvar=0, p=.5)
qgraph(cor(net50$cloned), layout = "spring", colors = adegenet::fac2col(assignments), labels = assignments)
title(sub = "half predictors by clone")
}
}