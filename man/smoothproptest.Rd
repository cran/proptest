\name{smoothproptest}
\alias{smoothproptest}
\alias{print.smoothproptest}
\title{Data-driven Smooth Test of the Proportional Hazards Assumption}
\description{
The function performs the Neyman type smooth test of the proportional hazards
assumption for an \emph{individual} covariate in the Cox model for right censored
survival data. Both a fixed and data-driven choice of the alternative model is
possible.
}
\usage{
smoothproptest(fit, covariate = 1, dims = 4, basis = "legendre",
               time.transf = "F", data.driven = TRUE,
               all.subsets = FALSE, h.approx = TRUE, sim = FALSE,
               nsim = 1000)
}
\arguments{
  \item{fit}{a Cox model fit (an output of \code{\link{coxph}}).}
  \item{covariate}{integer determining which covariate is to be tested
        for proportionality.}
  \item{dims}{a vector or a single value. \code{dims} gives dimensions
        for smooth modeling of the effects of the covariates that are not
        tested, and of the tested covariate. If \code{dims} is a single value
        and there is more than one covariate, the value is replicated.}
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
  \item{data.driven}{logical. Should the BIC be used?}
  \item{all.subsets}{logical. If \code{TRUE} then the BIC selects out of all
        the nonempty subsets, otherwise only the nested subsets are used.}
  \item{h.approx}{logical. Should we compute the \eqn{p}-value using the
        \eqn{H}-approximation? This is only valid for nested subsets
        (\code{all.subsets=FALSE}).}
  \item{sim}{logical. Should we compute the \eqn{p}-value using the
        simulation approximation?}
  \item{nsim}{the number of simulations to be carried out to approximate
        the \eqn{p}-value.}
}
\details{
The Neyman type smooth test of proportionality (time-constancy of the coefficient)
of a covariate against the alternative of the time-varying coefficient consists of
expressing the coefficient of the tested covariate as a linear combination of
basis functions and testing significancy of the new artificial time-dependent
covariates using the partial likelihood score test.

In the data-driven version, the alternative is selected by a BIC-like rule. The
distribution of the test statistic then may be approximated by the
\eqn{H}-approximation if the nested alternatives are used, or by simulations if
either the nested subsets or all the subsets are used. See Kraus (2005).

The potentially time-varying effects of the covariates that are not tested
should be modeled as combinations of basis functions too. This makes it possible
to perform individual covariate tests. Not doing so would be dangerous: the test
generally would not be capable to distinguish which covariate is proportional
and which not. See Kraus (2006).

The vector \code{dims} gives the number of the basis functions for each covariate.
}
\value{
A list (an object of class \code{"smoothproptest"}) containing some of input values,
test statistics and \eqn{p}-values computed by various methods (some of them are
\code{NULL} if not computed). The most important components are:
  \item{stat}{the test statistic (\code{stat.bic} if \code{data.driven=TRUE},
        \code{stat.d} otherwise).}
  \item{p}{the \eqn{p}-value corresponding to \code{stat} (one of the \eqn{p}-values
        below).}
  \item{stat.d}{the test statistic of the fixed dimension test.}
  \item{p.d.chisqd, p.d.sim}{the \eqn{p}-value of the fixed dimension test based on
        the asymptotic \eqn{\chi^2_d} distribution and on simulations.}
  \item{stat.bic}{the test statistic of the data-driven test.}
  \item{p.bic.h, p.bic.sim, p.bic.chisq1}{the \eqn{p}-value of the data-driven test
        based on the \eqn{H}-approximation (only for nested alternatives), simulations,
        and asymptotic \eqn{\chi^2_1} (only for nested alternatives; theoretically
        works, practically not).}
}
\references{
Kraus, D. (2005).
Data-driven smooth tests of the proportional hazards assumption.
Research Report 2134, Institute of Information Theory and Automation,
Prague.

Kraus, D. (2006).
Identifying nonproportional covariates in the Cox model.
Research Report 2170, Institute of Information Theory and Automation,
Prague.
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

## Tests using nested subsets; p-values computed by H-approximation

test1 = smoothproptest(fit,covariate=1)   # testing Z_1 (nonproportional)
print(test1,print.alt=TRUE)   # print details on the alternative models

test2 = smoothproptest(fit,covariate=2)   # testing Z_2 (proportional)
print(test2,print.alt=TRUE)   # print details on the alternative models

## Tests using all subsets; now simulations must be used for p-values

test1 = smoothproptest(fit,covariate=1,all.subsets=TRUE,sim=TRUE)   # Z_1
print(test1,print.alt=TRUE)   # print details on the alternative models

test2 = smoothproptest(fit,covariate=2,all.subsets=TRUE,sim=TRUE)   # Z_2
print(test2,print.alt=TRUE)   # print details on the alternative models
}
\keyword{survival}

