\name{collapse_PCA}
\alias{collapse_PCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Automatic module detection with PCA
}
\description{
Perform a principal component analysis to recover the module structure of a cloned network.
}
\usage{
collapse_PCA(data, N = NULL, assignments, around = 5, rep = 1000, quantile = .95, rotate = "oblimin", stoppingrule = c("easystop", "parallel", "optimal"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Matrix of numerics. The data file.
}
  \item{N}{Integer. The dimension of the original network. Only relevant if \code{"easystop"} or \code{"optimal"} are selected as a stopping rule.
}
  \item{assignments}{Vector of integers. The correct clustering to recover, the output of function \code{\link{Assign}}, used to compute the adjusted Rand Index. Only relevant if \code{"optimal"} is selected as a stoppingrule. See details.
}
\item{around}{Integer. If \code{"optimal"} is selected as a stopping rule, severel numbers of components are considered, from \code{N-around} to \code{N+around}. See details.}
\item{rep}{Integer. Iterations for parallel analysis.
}
  \item{quantile}{Numeric in [0,1]. Quantile for the parallel analysis. Conventional values are .95 or .99.}
  \item{rotate}{Character string. The kind of rotation. An "oblimin" rotation is suggested.}
  \item{stoppingrule}{Vector of strings. Specify the stopping rule to determine the number of components. Can include one or more methods among \code{"easystop"}, \code{"parallel"}, and \code{"optimal"}.
}
}
\details{
The function performs a principal component analysis. Once the PCA is performed, each node is assigned to a module according to its highest component loading.

The number of components can be determined in three ways. If \code{"easystop"} is included in \code{stoppinngrule}, the correct number of factors (as specified by \code{N}) are kept. If \code{"parallel"} is included in \code{stoppingrule}, the number of components is determined using parallel analysis. If mehod \code{"optimal"} is included in \code{stoppingrule}, several numbers of factors are considered in the surroundings of \code{N}, from \code{N-around} to \code{N+around}. The number of factors that results in the best adjusted Rand index is retained (see \code{\link[mclust]{adjustedRandIndex}}).
}

\value{Matrix of integers. A row for each cloned node (= a row for each column of the input matrix \code{data}), a column for each stopping rule. An integer indicates the belonging of each node to a module. Consider that the choice of the numbers is arbitrary, they are on a nominal scale. The only important thing is whether two nodes are in the same module or not.
}
\author{
Giulio Costantini
}

\seealso{
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Principal Component Analysis }
\keyword{ parallel analysis }% __ONLY ONE__ keyword per line
