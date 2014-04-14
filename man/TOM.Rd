\name{TOM}
\alias{TOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variations on Topological dissimilarity Matrix (opposite of the Topological Overlap Matrix).
}
\description{
Compute several variations of the (opposite of the) Topological Overlap Matrix (Zhang & Horvath, 2005).
}
\usage{
TOM(data, TOMType="signed", TOMDenom = "mean", squared = FALSE, TOMPCA = FALSE, rotate="varimax", use = "pairwise.complete.obs", mineigen=c(0, .7, 1) )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A dataframe, observations by variables.
  }
  \item{TOMType}{String. One of "signed" and "unsigned". If signed, the signs in the correlation matrix are considered (see Langfelder, 2013).
  }
  \item{TOMDenom}{String. One of \code{"min"}, \code{"mean"} and \code{"squared"}. The options \code{"min"} and \code{"mean"} are described in Langfelder (2012). The option squared is a possible variation that I am still considering. It might not work properly.
  }
  \item{squared}{Logical. If \code{TRUE} the squared direct connection is considered in the numerator.
  }
  \item{TOMPCA}{Logical. If \code{TRUE} a technique that combines both TOM and PCA is used. It can be particularly effective in presence of clones that can inflate the overlap among a pair of nodes.
  }
  \item{rotate}{String. Only relevant if \code{TOMPCA = TRUE}. A rotation technique for function \link[psych]{principal}.
  }
  \item{mineigen}{Vector of numerics. Only relevant if \code{TOMPCA = TRUE}. Only the component scores with eigenvalue larger than \code{mineigen} are considered in computing TOM. If the length of \code{mineigen} is >1, then a TOM matrix will be given in output for each value in \code{mineigen}.
  }
}
  
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A matrix, each element i,j including the topological dissimilarity between a pair of elements.
}

\references{
Langfelder, P. (2012). Topological overlap and the case of disappearing small modules. Technical report.
Langfelder, P. (2013). Signed vs. unsigned topological overlap matrix. Technical Report.
Zhang, B., & Horvath, S. (2005). A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology, 4(1).
}

\author{
Giulio Costantini
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[WGCNA]{TOMdist}}
}

\examples{
\dontrun{
data <-matrix(rnorm(100), ncol=10, nrow = 10)

TOM(data)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
