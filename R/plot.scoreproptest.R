"plot.scoreproptest" <-
function(x,nsim.plot=x$nsim.plot,...)
{
	time = x$time
	if (x$global) {
		test.process = x$test.process
		test.process.sim = matrix(x$test.process.sim,nrow=length(test.process))
		if (nsim.plot>x$nsim.plot) stop("can't plot ",nsim.plot," realisations, object contains only ",x$nsim.plot)
		ylim = range(c(test.process,test.process.sim))
		lwd=5
		plot(time,test.process,type="l",lwd=lwd,ylim=ylim,...)
		for (i in (1:nsim.plot)) lines(time,test.process.sim[,i],col="grey")
		lines(time,test.process,type="l",lwd=lwd)
	} else {
		score.process = x$score.process
		score.process.sim = matrix(x$score.process.sim,nrow=length(score.process))
		if (nsim.plot>x$nsim.plot) stop("can't plot ",nsim.plot," realisations, object contains only ",x$nsim.plot)
		ylim = range(c(score.process,score.process.sim))
		lwd=5
		plot(time,score.process,type="l",lwd=lwd,ylim=ylim,...)
		for (i in (1:nsim.plot)) lines(time,score.process.sim[,i],col="grey")
		lines(time,score.process,type="l",lwd=lwd)
	}
	invisible(x)
}

