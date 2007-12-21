\name{plot.scoreproptest}
\alias{plot.scoreproptest}
\title{Plotting the Observed Test Process and its Simulations}
\description{
The function plots the observed test process and a number of
its realisations simulated under the hypothesis of proportional
hazards.
}
\usage{
\method{plot}{scoreproptest}(x, nsim.plot = x$nsim.plot, ...)
}
\arguments{
  \item{x}{an object of class \code{"scoreproptest"} (output of
   \code{\link{scoreproptest}}).}
  \item{nsim.plot}{the number of simulated paths of the test process
   to be plotted. It must not be greater than \code{x$nsim.plot}.}
  \item{...}{further plotting parameters.}
}
\details{
By plotting the observed path of the test process along with
its simulations, one can visually assess the time-constancy
of the effect of the corresponding covariate or the global
proportional hazards assumption.

The function plots \code{x$score.process} or \code{x$test.process}
and the first \code{nsim.plot} realisations contained in
\code{x$score.process.sim} or \code{x$test.process.sim}.
}
\references{
Lin, D.Y., Wei, L.J. and Ying, Z. (1993)
Checking the Cox model with cumulative sums of martingale-based
residuals. \emph{Biometrika}, \bold{80}, 557--572.
}
\author{David Kraus, \url{http://www.davidkraus.net/}}
\seealso{\code{\link{scoreproptest}}}
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
\keyword{survival}

