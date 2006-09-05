\name{mini.survival}
\alias{mini.survival}
\alias{mini.agreg.fit}
\alias{mini.coxph.detail}
\alias{mini.survfit.coxph}
\title{Minimalist Versions of Functions from the Package survival}
\description{
Minimalist (lightweight) versions of some functions from
the package \pkg{survival} used in \code{\link{smoothproptest}}
and \code{\link{scoreproptest}}. For internal use only.
}
\usage{
mini.agreg.fit(x, y, init, control, sort.start, sort.end)

mini.coxph.detail(object)

mini.survfit.coxph(object)
}
\arguments{
  \item{x}{covariate matrix.}
  \item{y}{response (matrix with three columns: start time, stop time,
        event indicator).}
  \item{init}{initial estimates.}
  \item{control}{control options (an object of class
        \code{\link{coxph.control}}).}
  \item{sort.start, sort.end}{order.}
  \item{object}{a Cox model fit.}
}
\details{
These functions are minimalist versions of the functions
\code{\link{agreg.fit}}, \code{\link{coxph.detail}} and
\code{\link{survfit.coxph}} extracted from the package \pkg{survival}.
They do not perform error checks and do not work with data
frames, model frames etc., which makes them much faster.
They are used in \code{\link{smoothproptest}}
and \code{\link{scoreproptest}}.
}
\author{David Kraus, \url{http://www.davidkraus.net/}}
\section{Warning }{
These functions are for internal use only.
}
\keyword{survival}
