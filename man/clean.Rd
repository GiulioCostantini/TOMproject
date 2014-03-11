\name{clean}
\alias{clean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cleans variables from the variance shared with other variables.
}
\description{The function takes in input one or two data matrices, \code{x} and \code{y}, and uses multiple regression to remove from the variables in \code{x} all the variance shared with variables in \code{y}. A parameter allows to remove also the variance shared among the variables in \code{x}.
}
\usage{
clean(x, y=NULL, itself = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data matrix, observations by variables.
}
  \item{y}{A data matrix, observations by variables. The number of observations must be the same as in x.
}
 \item{itself}{Logical. Should be the variables in x cleaned from their own shared variance? 
 }
}
\details{
There are three main ways to use this function:

The first is to give a set of variables \code{x} in input, with parameter \code{itself} = TRUE. The variables included in the output data matrix are the same as \code{x}, but with all their shared variance removed, so their correlation is perfectly zero.

The second is to specify a matrix \code{y} in input and parameter \code{itself} = FALSE. In this case, the \code{x} variables are cleaned from the variance that they share with those in \code{y}, but not with the variance that they share with other variables in \code{x}.

The third is to specify a matrix \code{y} in input and parameter \code{itself} = TRUE. In this case, the variables are cleaned both from the shared variance with the variables in \code{y} and from that shared with other variables in \code{x}.
}
\value{A data matrix, same dimensions as x.
}
\author{Giulio Costantini
}
\examples{
# generate x, a matrix of correlated variables
pcm <- matrix(c(1, .2, .3,
         .2, 1, .1,
         .3, .1, 1), ncol = 3)
x <- pcm2data(pcm, 1000)
round(cor(x), 3)

# clean the variables in x from their shared variance
# and save them in variable x1.
x1 <- clean(x, itself = TRUE)
round(cor(x1), 3)

# generate another matrix, y, of variables that correlate with x
y <- matrix(rnorm(3000), ncol = 3)
round(cor(cbind(x, y)), 3)

# clean the variables in x only from the variance that they share with variables in y
x2 <- clean(x, y, itself = FALSE)
round(cor(cbind(x2, y)),3)

# clean the variables in x both from the variance shared with y and from the variance shared with other variables in x
x3 <- clean(x, y, itself = TRUE)
round(cor(cbind(x3, y)),3)
}

\keyword{ residuals }
\keyword{ multiple regression }