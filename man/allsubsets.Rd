\name{allsubsets}
\alias{allsubsets}
\title{All Subsets}
\description{
The function generates all subsets of the set \eqn{\{1,\dots,n\}}.
}
\usage{
allsubsets(n)
}
\arguments{
\item{n}{the number of elements of the set whose subsets are to be returned.}
}
\value{
Returns a \eqn{2^n} by \eqn{n} matrix of logical values.
Each row corresponds to one subset.
}
\author{David Kraus, \url{http://www.davidkraus.net/}}
\note{This functions is used in the all subsets variant of the
data-driven smooth tests of the proportional hazards assumption,
see \code{\link{smoothproptest}}.}
\examples{
allsubsets(4)
}
\keyword{survival}

