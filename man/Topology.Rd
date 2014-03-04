\name{Topology}
\alias{Topology}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate an adjacency matrix.
}
\description{
Generate an adjacency matrix and a partial correlation matrix encoding the same conditional independencies. It allows to specify the desired number of nodes, the number of edges and the topology.
}
\usage{
Topology(N, m, topology = c("ER", "BA", "WS"), exact.m = TRUE, force.connected = TRUE, negpro = 0, not = NULL, pw_ba = 1, p_ws = 0.1, minpcor = 0.2, maxpcor = 0.9, maxiter = c(50, 1000, 50), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Integer. Number of nodes
}
  \item{m}{Integer. Number of edges.
  }
  \item{topology}{Character. The topology on which the network is based. One of
  "ER" = Erdos-Renyi random network.
  "BA" = Barabasi & Albert scale-free network.
  "WS" = Watts & Strogatz random networks.
}
  \item{exact.m}{Logical. Specify whether the network should have exactly \code{m} edges (\code{TRUE}) or at least \code{m} edges. If \code{TRUE} the computation can be very slow in some cases.
}
  \item{force.connected}{Logical. Specify whether the network should be connected.
  }
  \item{negpro}{Integer in [0,1]. The desired proprtion of negative edges. Suggested value is 0 for the simulation, otherwise positive-definiteness problems are more likely to occur.
  }
  \item{not}{List of matrices. A list of adjacency matrices to prevent extracting networks isomorphic to networks extraced at earlier stages.}
  \item{pw_ba}{Numeric, specify the power of the preferential attachment. Only relevant if \code{topology = "BA"}.
}
  \item{p_ws}{Numeric in [0,1]. Probability of random rewiring. Only relevant if \code{topology = "WS"}
}
  \item{minpcor}{Numeric. Minimum allowed value for partial correlation.
}
  \item{maxpcor}{Numeric. Maximum allowed value for partial correlation.
}
\item{maxiter}{Vector of three integers, specifying the maximum allowed iterations for three stages of computation.

\code{maxiter[1]}. Maximum number of attempts to search for a network that is not isomorphic to other networks in "not" and that results also into a positive definite partial correlation matrix. This is the most global level of iteration.

\code{maxiter[2]}. Maximum number of attempts to find a connected network, when \code{topology = "ER"} or \code{topology = "WS"} are specified.

\code{maxiter[3]}. Maximum number of attemtps to find a positive definite correlation matrix in \code{\link{ggm2pcm}}. I reccommend a relatively low value for this because in case of repeated failure (maxiter[3] failures), a new network is drawn (up to maxiter[1] times) and the procedure is repeated: insisting too much on a topological configuration that does not appear to work may not be the best solution.
}
\item{verbose}{Logical. Do you want the function to tell you how it is going and where it is struggling?}
}
\details{
This function creates a network with the desired topology, allowing to regulate directly the number of edges (which is not always possible using igraph functions). The number of edges is ensured for \code{topology = "ER"} or \code{topology = "WS"} by specifiying the nearest number of edge larger than m and then removing edges randomly. If \code{force.connected = TRUE}, then whenever an edge is removed, a check is done to preven removing a bridge, therefore disconnecting the network (see \code{\link{bridges}}).
The function also tests for graph isomorphisms: if the obtained graph is isomorphic to a graph in the "not" list, the function keeps searching until a non-isomorphic graph is found (or until \code{maxiter[1]} calls are done).
}
\value{A list of two elements
\item{ggm}{Adjacency matrix}
\item{pcm}{Partial correlation matrix}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Giulio Costantini
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# this should make it struggle a bit, but then come up with a positive definite matrix
net <- Topology(N=40, m=60, topology="WS", exact.m=TRUE, force.connected=TRUE, negpro=0, p_ws=.2, minpcor=.1, maxpcor=.5)
par(mfrow = c(1,2))
lay <- qgraph(net$ggm)
qgraph(net$pcm, layout = lay$layout)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
