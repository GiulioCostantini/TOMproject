\name{Assign}
\alias{Assign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Assign a certain amount of clones to each node.
}
\description{Assign clones to each node according to different distributions, uniform (same number of clones for each node), pwr (power function) and beta (beta distribution).
}
\usage{Assign(N, N2, method_assign = c("unif", "pwr", "beta"), pwr = 1, shape1 = .1, shape2 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Integer. Number of nodes in the original network.
}
  \item{N2}{Integer. Number of nodes in the cloned network. N2 > N.
}
  \item{method_assign}{Character string. One of "unif", "pwr" or "beta".
}
  \item{pwr}{Numeric. The exponent. Only relevant if method_assign = "pwr".
}
  \item{shape1}{Numeric. The first parameter of the beta distribution. Only relevant if method_assign = "beta".
}
  \item{shape2}{Numeric. The second parameter of the beta distribution. Only relevant if method_assign = "beta".
}
}
\details{Independent of the method chosen, each node has always at least one clone.
Method "unif" just assigns the same number of clones to each node (+-1, if N2 is not multiple of N).
Method "pwr" draws N equally spaced numbers from .1 to .9, takes their pwr-th power, reshuffles them and distributes the clones to the N nodes proportionally to the resulting vector. Method "pwr" with pwr = 0 is equivalent to method "unif".
Method "beta" draws N random numbers from a beta distribution with parameters shape1 and shape2 and distributes the clones to the N nodes, proportionally to the resulting vector. For shape1>>1 and shape2<<1 you get a uniform distribution of clones to nodes, for shape1<<1 and shape2>>1 you get that most nodes have only one clone and a few nodes have most clones.
For both methods "pwr" and "beta" rounding issues may cause the presence of few more/less clones than desired. In this case, missing clones are added randomly or clones in excess are removed randomly.
}
\value{
A vector of integers, indicating how many clones are assigned to each node.
}
\author{Giulio Costantini
}

\examples{
# uniform assignment
Assign(N = 10, N2 = 100, method_assign = "unif")
# very asymmetric power assignment
Assign(N = 10, N2 = 100, method_assign = "pwr", pwr = 10)
# very symmetric power assignment
Assign(N = 10, N2 = 100, method_assign = "pwr", pwr = 0.2)
# very asymmetric beta assignment
Assign(N = 10, N2 = 100, method_assign = "beta", shape1 = .1, shape2 = 1)
# very symmetric beta assignment
Assign(N = 10, N2 = 100, method_assign = "beta", shape1 = 100, shape2 = .001)
}
