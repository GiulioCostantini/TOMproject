\name{pcm2data}
\alias{pcm2data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial correlation matrix to data.}
\description{
Convert a partial correlation matrix to a data matrix using Cholesky decomposition.}
\usage{
pcm2data(pcm, n = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pcm}{A partial correlation matrix.
}
  \item{n}{Integer. The  number of observations in the output data.
}
}
\details{Use the Cholesky decomposition to generate data from a partial correlation matrix.
}
\value{
A data matrix, of dimensions \code{n} and \code{ncol(pcm)}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
ggm <- matrix(
c(0, 1, 1, 0,
  1, 0, 1, 1,
  1, 1, 0, 1,
  0, 1, 1, 0), ncol = 4)
pcm <- ggm2pcm(ggm, minpcor = 0.2, maxpcor = 0.9, maxiter = 1000, verbose = FALSE)
dat <- pcm2data(pcm, n = 1000)

# plot
par(mfrow = c(3,1))
net1 <- qgraph(ggm)
title(main = "adjacency matrix")
qgraph(pcm, layout = net1$layout)
title(main = "original pcm")
qgraph(cor2pcor(cor(dat)), layout = net1$layout)
title(main = "pcm reproduced by data")
}