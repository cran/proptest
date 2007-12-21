\name{smoothproptest}
\alias{smoothproptest}
\alias{summary.smoothproptest}
\title{Data-driven Smooth Test of the Proportional Hazards Assumption}
\description{
The function performs the Neyman type smooth test of the proportional hazards
assumption for an \emph{individual} covariate or \emph{globally} for all of the
covariates in the Cox model for right censored survival data. Both a fixed and
data-driven choice of the alternative model is possible.
}
\usage{
smoothproptest(fit, covariate = 1, global = FALSE, dims = 4,
               basis = "legendre", time.transf = "F",
               data.driven = TRUE, nsim = 1000)

\method{summary}{smoothproptest}(object, ...)
}
\arguments{
  \item{fit}{a Cox model fit (an output of \code{\link{coxph}}).}
  \item{covariate}{integer determining which covariate is to be tested
        for proportionality. Only used if \code{global} is \code{FALSE}.}
  \item{global}{logical. Should the global test be performed?}
  \item{dims}{a vector or a single value. \code{dims} gives dimensions for
        smooth modelling of the effects of respective covariates. If \code{dims}
        is a single value and there is more than one covariate, the value is
        replicated.}
  \item{basis}{a character string. Which basis of smooth functions is to be
        used? Possible values are \code{"legendre"} and \code{"cos"} (or
        \code{"cosine"}).}
  \item{time.transf}{a character string. The basis functions are evaluated at
        transformed times. With \code{time.transf="F"}, the transformation
        is \eqn{F_0(t)/F_0(\tau)}{F_0(t)/F_0(tau)} (\eqn{F_0} is the
        distribution function corresponding to the baseline hazard).
        For \code{time.transf="L"}, the transformation is
        \eqn{\Lambda_0(t)/\Lambda_0(\tau)}{Lambda_0(t)/Lambda_0(tau)}
        (\eqn{\Lambda_0}{Lambda_0} is the cumulative baseline hazard).
        \eqn{F_0} and \eqn{\Lambda_0}{Lambda_0} are estimated.}
  \item{data.driven}{logical. Should the BIC be used to select the alternative?}
  \item{nsim}{the number of simulations to be carried out to compute the
        two-term approximation of the \eqn{p}-value. This is only used when
        the test is global.}
  \item{object}{an object of class \code{"smoothproptest"}, as returned
        by the function \code{smoothproptest}.}
  \item{...}{further parameters for printing.}
}
\details{
The Neyman type smooth test of proportionality (time-constancy of the coefficient)
of an \emph{individual} covariate against the alternative of the time-varying
coefficient consists of expressing the coefficient of the tested covariate as a
linear combination of basis functions and testing significance of the new
artificial time-dependent covariates using the partial likelihood score test.
The potentially time-varying effects of the covariates that are not tested
should be modelled as combinations of basis functions too (\code{dims} gives
dimensions for smooth modelling of the effects of the covariates which are not
tested, and of the tested covariate). This makes it possible to perform individual
covariate tests. Not modelling the other covariates would be dangerous: the test
generally would not be able to distinguish which covariate is proportional and
which not. See Kraus (2008).

In the \emph{global} test the coefficients of all of the covariates are expressed
as linear combinations of basis functions and tested for significance.
The vector \code{dims} gives the number of the basis functions for each covariate.

Both \emph{individual} and \emph{global} tests can be data-driven which means
that the alternative is selected by a BIC-like rule. For individual
(covariate-specific) tests the number of basis functions for the tested
covariate is selected (see Kraus, 2007a). For global tests, the order of
approximation of possibly time-varying effects is selected for each of
the covariates (see Kraus, 2007b, Chapter 6).

The distribution of data-driven test statistics is approximated by the
two-term \eqn{H}-approximation (Kraus, 2007a,b) because the one-term
asymptotic approximation is inaccurate.

If the test is data-driven, the \code{summary} method prints details
on the selection procedure (statistics and penalised statistics for
each alternative model). This is equivalent to \code{print(x, detail=TRUE, ...)}.
}
\value{
A list (an object of class \code{"smoothproptest"}) containing some of input values,
the test statistic(s) and \eqn{p}-value(s). The most important components are:
  \item{stat}{the test statistic (\code{stat.bic} if \code{data.driven=TRUE},
        \code{stat.d} otherwise).}
  \item{p}{the \eqn{p}-value corresponding to \code{stat} (one of \code{p.d}
        or \code{p.bic.h} below).}
  \item{stat.d, stat.bic}{the test statistic of the fixed dimension and data-driven
        test (the latter present only if \code{data.driven} is \code{TRUE}).}
  \item{p.d}{the \eqn{p}-value of the fixed dimension test based on the asymptotic
        chi-square with \eqn{d} df (\code{p.d}). Here \eqn{d} is either
        \code{dims[covariate]} for the individual covariate test or \code{sum(dims)}
        for the global test.}
  \item{p.bic.asympt, p.bic.h}{\eqn{p}-values for the data-driven test. \code{p.bic.asympt}
        is based on the asymptotic distribution (chi-square with 1 df for covariate-specific
        tests, max of chi-square with 1 df for global tests) (inaccurate).
        \code{p.bic.h} is based on the two-term approximation.}
  \item{scorestats, scorestats.penal}{statistics and penalised statistics
        for all alternative models (only for data-driven tests).}
  \item{alt}{all alternative models (only for data-driven tests).}
  \item{S}{the index of the selected alternative (only for data-driven tests).
        \code{alt[S,]} is the selected alternative.}
}
\references{
Kraus, D. (2007a) Data-driven smooth tests of the proportional
hazards assumption. \emph{Lifetime Data Anal.} \bold{13}, 1--16.

Kraus, D. (2007b) Neyman's smooth tests in survival analysis.
PhD thesis. Charles University in Prague, Dept. of Statistics.

Kraus, D. (2008) Identifying nonproportional covariates in the
Cox model. \emph{Comm. Statist. Theory Methods} \bold{37}.
To appear.
}
\author{David Kraus, \url{http://www.davidkraus.net/}}
\seealso{\code{\link{scoreproptest}}, \code{\link{coxph}}}
\examples{
## Case 4 of Kvaloy & Neef (2004, Lifetime Data Anal.):
## data generated from the distribution with hazard rate
## \lambda(t)=\exp(0.5tZ_1+Z_2-8)
## (Z_1,Z_2) jointly normal with E=4, var=1, cor=rho
## censoring times uniform(0,5)

n = 200
rho = .3
z = matrix(rnorm(n*2),ncol=2) \%*\% chol(matrix(c(1,rho,rho,1),2)) + 4
a = .5
tim = 1/(a*z[,1]) * log(1-a*z[,1]*exp(-z[,2]+8)*log(runif(n)))
ct = 5*runif(n)
nc = tim<=ct
tim = pmin(tim,ct)
fit = coxph(Surv(tim,nc)~z)

## Covariate-specific tests

test1 = smoothproptest(fit,covariate=1)   # testing Z_1 (nonproportional)
summary(test1)   # print details on the alternative models

test2 = smoothproptest(fit,covariate=2)   # testing Z_2 (proportional)
summary(test2)   # print details on the alternative models

## Global test

test = smoothproptest(fit,global=TRUE)   # global test
summary(test)   # print details on the alternative models
}
\keyword{survival}

