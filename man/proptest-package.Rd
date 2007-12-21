\name{proptest-package}
\alias{proptest-package}
\alias{proptest}
\docType{package}
\title{Tests of the Proportional Hazards Assumption}
\description{
Tests of the proportional hazards assumption in the Cox
model: data-driven Neyman type smooth tests and score process based
tests for identifying nonproportional covariates and for global checks.
}
\details{
The package provides two functions for testing proportional hazards:
\code{\link{smoothproptest}} for smooth tests, and
\code{\link{scoreproptest}} for tests based on the score process.

Smooth tests consist of expressing the coefficient of the tested
covariate or coefficients of all of the covariates as linear
combinations of some basis functions, i.e., they consist of
introducing artificial time-dependent covariates and testing their
significancy. A data-driven choice of these covariates is implemented.

Score process tests are tests based on functionals (Kolmogorov--Smirnov,
Cramer--von Mises, Anderson--Darling) of components of the score process
or on the supremum of a test process formed from the whole score process.
Both numerical and visual assessment based on simulations is possible.

Both types of tests may be used for \emph{global} verification as well as
for testing \emph{individual} covariates because the covariates that are
not tested may be modeled by smooth functions.
}
\author{David Kraus, \url{http://www.davidkraus.net/}}
\references{
Kraus, D. (2007) Data-driven smooth tests of the proportional
hazards assumption. \emph{Lifetime Data Anal.} \bold{13}, 1--16.

Kraus, D. (2008) Identifying nonproportional covariates in the
Cox model. \emph{Comm. Statist. Theory Methods} \bold{37}.
To appear.

Lin, D.Y., Wei, L.J. and Ying, Z. (1993)
Checking the Cox model with cumulative sums of martingale-based
residuals. \emph{Biometrika}, \bold{80}, 557--572.
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
summary(test1)   # print details on the alternative models

test2 = smoothproptest(fit,covariate=2)   # testing Z_2 (proportional)
summary(test2)   # print details on the alternative models

test = smoothproptest(fit,global=TRUE)   # global test
summary(test)   # print details on the alternative models

## Tests based on the score process

par(mfrow=c(3,1))

test1 = scoreproptest(fit,covariate=1)   # testing Z_1 (nonproportional)
print(test1)
plot(test1,main="Score process for z1")

test2 = scoreproptest(fit,covariate=2)   # testing Z_2 (proportional)
print(test2)
plot(test2,main="Score process for z2")

test = scoreproptest(fit,global=TRUE)   # global test of proportional hazards
print(test)
plot(test,main="Global test process")

par(mfrow=c(1,1))
}
