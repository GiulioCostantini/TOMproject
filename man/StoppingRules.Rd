\name{StoppingRules}
\alias{StoppingRules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stopping rules for hierarchical cluster analysis
}
\description{
Examines several stopping rules included in the package \code{\link[NbClust]{NbClust}} and \code{\link{dynamicTreeCut}} and produce a matrix of partitions.
}
\usage{
StoppingRules(dist, data, method = "average")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{
A matrix of distances, given as a matrix, not as a dist object.
}
  \item{data}{
The original data matrix, observations by variables (clusters of variables are created).
}
  \item{method}{
Method for hierarchical cluster analysis, for the moment only method "average" has been tested.
}
}
\details{
The rule "easystop" means that the correct number of clusters is considered, the rule "optimal" means that the number of clusters with the highest value of adjusted rand index is considered.
}

\value{
\item{clusterings}{A matrix of clusterings, an element by row, a stopping rule by column}
\item{values}{A matrix of information on each stopping rule, its value or its parameters. A row by parameter, a column by stopping rule. Most values are \code{NA}s}

}
\references{
Langfelder, P., Zhang, B., & Horvath, S. (2008). Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. Bioinformatics, 24(5), 719-720.
Milligan, G. W., & Cooper, M. C. (1985). An examination of procedures for determining the number of clusters in a data set. Psychometri, 50(2), 159-179.
}
\author{
Giulio Costantini
}
\note{
A large part of the code has been taken from package \code{\link[NbClust]{NbClust}} and has been modified to fit our necessities
}

\examples{
\dontrun{
library(TOMproject)
N <- 10
N2 <- 30
m <- 15

# define a network
net <- Topology(N=N, m=m, topology="BA", exact.m=F, force.connected=T, negpro=0, pw_ba=1, minpcor=.1, maxpcor=.3)
pcm <- net$pcm
# assign clones
assi <- Assign(N=N, N2=N2, method_assign="pwr", shape1= .0001, shape2 = 1000)
assignments <- c()
for(i in 1:length(assi)) assignments <- c(assignments, rep(i, assi[i]))

# clone the network
cloned <- Clone(pcm=pcm, assi=assi, n=1000, furthersharedvar=.5, rndvar=.5, p=.9, keeporig=TRUE)

## TOMZhang and stopping rules
library(WGCNA)
dist <- TOMdist(adjMat=cor(cloned$cloned), TOMType="signed", TOMDenom="mean")
res <- StoppingRules(dist = dist, data = cloned$cloned, method = "average")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
