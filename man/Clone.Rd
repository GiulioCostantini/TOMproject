\name{Clone}
\alias{Clone}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clone a network using one of three different methods.
}
\description{
This function starts from an original network, represented by a partial correlation matrix (\code{pcm}), and from a vector (\code{assi}; es. the output of function \code{\link{Assign}}) specifying for each node in the original network, how many clones will correspond to it. The result is a new network in which each node is a clone of a node in the original network. The output is not in the form of a network, but of a data matrix.
}
\usage{
Clone(pcm, assi, method_clo = c("copy", "clone", "pclone"), n = 10000, fitvar = 0.5, furthersharedvar = 0.1, uniqueness = 0.2, p = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcm}{Matrix of numerics. Partial correlation matrix.
}
  \item{assi}{Vector of integers. The i-th element specifies the number of clones of the i-th node.
}
  \item{method_clo}{Character string. The procedure to clone the network, see details.
  }
  \item{n}{Integer. Number of observations in the data matrix.
}
  \item{fitvar}{Numeric. Portion of systematic variance.
}
  \item{furthersharedvar}{Numeric. Portion of random variance shared among the clones of the same node. Only relevant for methods \code{fclone} and \code{pclone}
}
  \item{uniqueness}{Numeric. Portion of unique variance for each clone.
}
  \item{p}{Numeric in [0,1]. Probability of including a predictor. Relevan only for method \code{pclone} (see details).
}
}
\details{
The parameter method_clo regulates how the variance of each clone node is composed.

- metod \code{copy}: the variance of each clone in the output network is equal to the same variance of the original node, plus random noise. The amounts of variance are proportional to the parameters \code{fitvar} and \code{uniqueness}.

- method \code{fclone}: full clonation. Each clone in the original network is equal to the weighted sum of three components. 1) The variance of the original node connected to all the other nodes in the network (the amount is determined by parameter \code{fitvar}). This component is obtained by running a regression of the original node on all the other nodes in the network (plus the clones created so far) and keeping only the fitted values. 2) A further portion of variance shared with all the other clones of the same node, but not with other nodes, even if connected to the original node (the amount is determined by parameter  \code{furthersharedvar}), 3) random noise (determined by parameter \code{uniqueness}).
  
- method  \code{pclone}: partial clonation. The same as method  \code{fclone}, but for the first component of variance, which instead of being determined by regressing the original node on all the other nodes in the network, is determined for each clone by regressing it only on a subset of such nodes. Each predictor is selected with probability \code{p}.
  
The methods  \code{fclone} and  \code{pclone} are based on the idea of Chung et al. (2003) of full and partial duplication processes in a network respectively, the most important differences being that we use weighted networks and start from a network with N nodes.

Ideally the parameters \code{fitvar}, \code{furthershared} and \code{uniqueness} should sum up to one, thus specifying the proportion of each component of variance, parameter  \code{furthershared} being zero for method \code{copy}. If their sum is not one, they are rescaled such as they sum up to one.

}
\value{A matrix of numeric. A data matrix \code{n*sum(assi)}, rows are observations and columns are (standardized) variables. Each variable is a clone of a variable in the original network.
}
\references{
Chung, F., Linyuan, L., Dewey, G., Galas, D. J. (2003), Duplication models for biological networks. Journal of Computational Biology, 10(5), p. 677-687}
\author{Giulio Costantini
}

\examples{
\dontrun{
N <- 10
N2 <- 30
m <- 9
# define the original network
net <- Topology(N=N, m=m, topology="BA", exact.m=F, force.connected=T, negpro=0, pw_ba=1, minpcor=.1, maxpcor=.6)
# extract the partial correlation matrix
pcm <- net$pcm
# clone the network with different methods
assi <- Assign(N=N, N2=N2, method_assign="pwr", shape1= .0001, shape2 = 1000)
netcop <- Clone(pcm=pcm, assi=assi, method_clo="copy", n=1000, fitvar=.5, uniqueness = .5)
netclo <- Clone(pcm=pcm, assi=assi, method_clo="clone", n=1000, fitvar=.3, uniqueness = .4, furthershared = .3)
netpclo <- Clone(pcm=pcm, assi=assi, method_clo="pclone", n=1000, fitvar=1, uniqueness = 0, furthershared = 0, p = .5)

## plot the original network and the three cloned networks
assignments <- c()
for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))
par(mfrow = c(2,2))
library(adegenet)
qgraph(pcm, layout = "spring", color = fac2col(unique(assignments)))
title(main = "original")
qgraph(cor(netcop), layout="spring", color = fac2col(assignments))
title(main = "method copy")
qgraph(cor(netclo), layout="spring", color = fac2col(assignments))
title(main = "method fclone")
qgraph(cor(netpclo), layout="spring", color = fac2col(assignments))
title(main = "method pclone")
}
}