\name{scoreproptest}
\alias{scoreproptest}
\alias{print.scoreproptest}
\title{Test of the Proportional Hazards Assumption Based on the Score Process}
\description{
The function performs tests of the proportional hazards assumption for an
\emph{individual} covariate in the Cox model for right censored survival data.
The tests of the Kolmogorov--Smirnov, Cramer--von Mises and Anderson--Darling
type based on the component of the score process are computed and simulated
\eqn{p}-values are returned.
}
\usage{
scoreproptest(fit, covariate = 1, dims = 4, basis = "legendre",
              time.transf = "F", nsim = 1000, nsim.plot = 50,
              weight = "unit")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{a Cox model fit (an output of \code{\link{coxph}}).}
  \item{covariate}{integer determining which covariate is to be tested
        for proportionality.}
  \item{dims}{a vector or a single value. \code{dims} gives dimensions
        for smooth modeling of the effects of the covariates that are not
        tested. If \code{dims} is a single value and there is more than one
        covariate, the value is replicated.}
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
        \eqn{F_0} and \eqn{\Lambda_0}{Lambda_0} are estimated from the
        input model \code{fit}.}
  \item{nsim}{the number of simulations to be carried out to approximate
        the \eqn{p}-value.}
  \item{nsim.plot}{the number of simulated paths of the score process to be
        returned (intended for plotting with \code{\link{plot.scoreproptest}}).}
  \item{weight}{a character string. The weighted process can be used. Possible values
        are \code{"unit"} (unweighted, default), \code{"gehan"}
        (\var{nrisk/nevent}) and \code{"mm"}
        (\eqn{(\var{nrisk}-\var{nrisk}[\var{nevent}])^2/(\var{nevent}^2)}{(\var{nrisk}-\var{nrisk}[\var{nevent}])^2/(\var{nevent}^2)})
        (Marzec and Marzec, 1997).}
}
\details{
The score process is used for assessment of the proportional hazards assumption
(Lin, Wei and Ying, 1993). Each component of the score process reflects
departures from proportionality (time-constancy of the effect) of the corresponding
covariate.

However, tests based on individual components of the process are generally
not capable to distinguish which covariate is proportional and which not. The method
is only valid if the other covariates are proportional.

Therefore, the potentially time-varying effects of the covariates that are not tested
are modeled as combinations of basis functions. The test is then based on the score
process from this large model with artificial time-dependent covariates. This makes
it possible to perform individual covariate tests. See Kraus (2006).

The vector \code{dims} gives the number of the basis functions for each covariate.
}
\value{
A list (an object of class \code{"scoreproptest"}). The most important components are:
  \item{score.process}{the component of the score process corresponding to the
        tested covariate.}
  \item{time}{the vector of the event times.}
  \item{stat.ks}{the Kolmogorov--Smirnov test statistic.}
  \item{p.ks}{the simulated \eqn{p}-value for the Kolmogorov--Smirnov test.}
  \item{stat.cm}{the Cramer--von Mises test statistic.}
  \item{p.cm}{the simulated \eqn{p}-value for the Cramer--von Mises test.}
  \item{stat.ad}{the Anderson--Darling test statistic.}
  \item{p.ad}{the simulated \eqn{p}-value for the Anderson--Darling test.}
  \item{score.process.sim}{a matrix with \code{nsim.plot} columns containing
        simulated paths of the score process. These may be plotted with
        \code{\link{plot.scoreproptest}}.}
}
\references{
Kraus, D. (2006).
Identifying nonproportional covariates in the Cox model.
Research Report 2170, Institute of Information Theory and Automation,
Prague.

Lin, D.Y., Wei, L.J. and Ying, Z. (1993).
Checking the Cox model with cumulative sums of martingale-based
residuals. \emph{Biometrika}, \bold{80}, 557--572.

Marzec, L. and Marzec, P. (1997).
Generalized martingale-residual processes for goodness-of-fit
inference in Cox's type regression models. \emph{Ann. Statist.},
\bold{25}, 683--714.
}
\author{David Kraus, \url{http://www.davidkraus.net/}}
\seealso{\code{\link{plot.scoreproptest}}, \code{\link{smoothproptest}},
\code{\link{coxph}}}
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

par(mfrow=c(2,1))

test1 = scoreproptest(fit,covariate=1)   # testing Z_1 (nonproportional)
print(test1)
plot(test1,main="Score process for z1")

test2 = scoreproptest(fit,covariate=2)   # testing Z_2 (proportional)
print(test2)
plot(test2,main="Score process for z2")

par(mfrow=c(1,1))
}
\keyword{survival}% at least one, from doc/KEYWORDS
