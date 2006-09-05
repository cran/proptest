\name{proptest-package}
\alias{proptest-package}
\alias{proptest}
\docType{package}
\title{Tests of the Proportional Hazards Assumption}
\description{
Tests of the proportional hazards assumption in the Cox
model: data-driven Neyman type smooth tests and score process based
tests for identifying nonproportional covariates.
}
\details{
\tabular{ll}{
Package: \tab proptest\cr
Type: \tab Package\cr
Version: \tab 0.1-0\cr
Date: \tab 2006-09-05\cr
License: \tab GPL\cr
}

The package provides two functions for testing proportional hazards:
\code{\link{smoothproptest}} for smooth tests, and
\code{\link{scoreproptest}} for tests based on the score process.

Smooth tests consist of expressing the coefficient of the tested
covariate as a linear combination of some basis functions, i.e.,
they consist of introducing artificial time-dependent covariates
and testing their significancy. A data-driven choice of these
covariates is implemented.

Score process tests are tests based on functionals (Kolmogorov--Smirnov,
Cramer--von Mises, Anderson--Darling) of components of the score process.
Both numerical and visual assessment based on simulations is possible.

Both of the tests may be used for testing \emph{individual} covariates
because the covariates that are not tested may be modeled by smooth
functions.
}
\author{David Kraus, \url{http://www.davidkraus.net/}}
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
\keyword{package}
\keyword{survival}
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

## Data-driven smooth tests

test1 = smoothproptest(fit,covariate=1)   # testing Z_1 (nonproportional)
print(test1,print.alt=TRUE)   # print details on the alternative models

test2 = smoothproptest(fit,covariate=2)   # testing Z_2 (proportional)
print(test2,print.alt=TRUE)   # print details on the alternative models

## Tests based on the score process

par(mfrow=c(2,1))

test1 = scoreproptest(fit,covariate=1)   # testing Z_1 (nonproportional)
print(test1)
plot(test1,main="Score process for z1")

test2 = scoreproptest(fit,covariate=2)   # testing Z_2 (proportional)
print(test2)
plot(test2,main="Score process for z2")

par(mfrow=c(1,1))
}
