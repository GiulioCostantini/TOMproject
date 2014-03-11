\name{collapse_PCA}
\alias{collapse_PCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Automatic Module Detection with principal component analysis
}
\description{
Perform a principal component analysis to recover the module structure of a cloned network.
}
\usage{
collapse_PCA(data, N = NULL, rep = 1000, quantile = 0.95, rotate = "oblimin", stoppingrule_pca = c("easystop", "parallel"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Matrix of numerics. The data file.
}
  \item{N}{Integer. The dimension of the original network. Only relevant if "easystop" is selected as a stopping rule.
}
  \item{rep}{Integer. Iterations for parallel analysis.
}
  \item{quantile}{Numeric in [0,1]. Quantile for the parallel analysis. Conventional values are .95 or .99.}
  \item{rotate}{Character string. The kind of rotation. An "oblimin" rotation is suggested.}
  \item{stoppingrule_pca}{Character vector. Specify the stopping rule to determine the number of components. Can be one or both of "easystop" and "parallel".
}
}
\details{
The function performs a principal component analysis. The number of components can be determined in two ways. If \code{"easystop"} is included in \code{stoppinngrule}, the correct number of factors (as specified by \code{N}) are kept. If \code{"parallel"} is included in \code{stoppingrule}, the number of components is determined using parallel analysis. Both methods can be specified for the stoppingrule.
Once the PCA is performed, each node is assigned to the component on which it has the highest loading.
}
\value{Matrix of integers. A row for each cloned node (= a row for each column of the input matrix \code{data}), a column for each stopping rule. An integer indicates the belonging of each node to a module. Consider that the choice of the numbers is arbitrary, they are on a nominal scale. The only important thing is whether two nodes are in the same module or not.
}
\author{
Giulio Costantini
}

\seealso{See also \code{\link{collapse_TOM}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Principal Component Analysis }
\keyword{ parallel analysis }% __ONLY ONE__ keyword per line
