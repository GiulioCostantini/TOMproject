\name{bridges}
\alias{bridges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bridges.
}
\description{
Find bridges in an undirected graph.
}
\usage{
bridges(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An undirected graph. Can be an igraph object, an adjacency matrix, an edgelist.
}
}
\details{
Find bridges using the algorithm by Schmidt (2013)}
\value{
A binary matrix of dimensions n = m = number of nodes. Each value one in position i,j indicates that the edge connecting nodes i and j is a bridge.
}
\references{
Schmidt (2013), A simple test on 2-vertex- and 2-edge-connectivity. Information Processing Letters 113, 241-244.}
\author{
Giulio Costantini
}
\examples{
x <- igraph::barabasi.game(n=30)
bridges(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bridges }
\keyword{ graphs }% __ONLY ONE__ keyword per line
