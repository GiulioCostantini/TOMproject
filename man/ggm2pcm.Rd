\name{ggm2pcm}
\alias{ggm2pcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert an adjacency matrix to a partial correlation matrix.
}
\description{
An adjacency matrix is converted to a weighted partial correlation matrix, with positive-definiteness as a condition.
}
\usage{
ggm2pcm(ggm, minpcor = 0.2, maxpcor = 0.9, maxiter = 1000, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ggm}{Adjacecny matrix, elements in {-1, 0, 1}.
}
  \item{minpcor}{Numeric. Minimum allowed value for partial correlation.
}
  \item{maxpcor}{Numeric. Maximum allowed value for partial correlation.
}
  \item{maxiter}{Integer. Maximum number of attempts to get a positive definite matrix.
}
  \item{verbose}{Logical. Specify whether the function inform you of each iteration.}
}
\details{
At each iteration, a candidate partial correlation matrix is created by drawing numbers from an uniform distribution from \code{minpcor} to \code{maxpcor}. If the matrix is positive definite and the correlation matrix implied is also positive definite, the matrix is returned in output. If positive-definiteness cannot be reaced after \code{maxiter} iterations, a non positive-definite matrix is returned with a warning.
}
\value{
\code{pcm}, a partial correlation matrix with the same conditional independencies specified in the adjacency matrix \code{ggm} given as input.
}
\author{
Giulio Costantini, Lourens Waldorp}
\examples{
ggm <- matrix(
c(0, 1, 1, 0,
  1, 0, 1, 1,
  1, 1, 0, 1,
  0, 1, 1, 0), ncol = 4)
pcm <- ggm2pcm(ggm, minpcor = 0.2, maxpcor = 0.9, maxiter = 1000, verbose = FALSE)

# plot the two matrices
par(mfrow = c(1,2))
net1 <- qgraph(ggm)
title(main = "adjacency matrix")
qgraph(pcm, layout = net1$layout)
title(main = "partial correlation matrix")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ partial correlation }
\keyword{ positive definite }% __ONLY ONE__ keyword per line
