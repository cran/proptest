"mini.agreg.fit" <-
function(x, y, init, control, sort.start, sort.end)
## returns object (list) that can be handled directly by (mini.)coxph.detail (goal: to avoid coxph with its model frames)
{
	n = nrow(y)
	nvar = ncol(x)
	start = y[,1]
	stopp = y[,2]
	event = y[,3]
	if(all(event==0)) stop("Can't fit a Cox model with zero failures")

	method = "efron"

	# Sort the data (or rather, get a list of sorted indices)
	#  For both stop and start times, the indices go from last to first
	# indices may be provided on input; goal: to avoid repeated sorting of the same thing
	if (missing(sort.start) || missing(sort.end)) {
		sort.end  = order(-stopp, event)
		sort.start= order(-start)
	}
	if (missing(control)) control = coxph.control()
	if (missing(init)) init = rep(0,nvar)
	maxiter = control$iter.max
	if (!is.null(init)) {
		if (length(init) != nvar) stop("Wrong length for inital values")
	} else
		init = rep(0,nvar)

	agfit = .C("agfit3", iter= as.integer(maxiter),
		as.integer(n),
		as.integer(nvar), 
		as.double(start), 
		as.double(stopp),
		as.integer(event),
		as.double(x),
		as.double(rep(0, n)), # offset
		as.double(rep(1, n)), # weights
		as.integer(1), # nstrat
		as.integer(n), # newstrat
		as.integer(sort.end-1),
		as.integer(sort.start-1),
		means = double(nvar),
		coef= as.double(init),
		u = double(nvar),
		imat= double(nvar*nvar), loglik=double(2),
		flag=integer(1),
		double(2*nvar*nvar +nvar*3 + n),
		as.double(control$eps),
		as.double(control$toler.chol),
		sctest=as.double(method=='efron'),PACKAGE="survival" )

	var = matrix(agfit$imat,nvar,nvar)
	coef = agfit$coef
	if (agfit$flag < nvar) which.sing = diag(var)==0
	else which.sing = rep(FALSE,nvar)

	infs = abs(agfit$u %*% var)
	if (maxiter >1) {
		if (agfit$flag == 1000)
			warning("Ran out of iterations and did not converge")
		else {
			infs = ((infs > control$eps) & infs > control$toler.inf*abs(coef))
			if (any(infs))
				warning(paste("Loglik converged before variable ",
					paste((1:nvar)[infs],collapse=","),
					"; beta may be infinite. "))
		}
	}
	lp  = x %*% coef - sum(coef *agfit$means)
	score = as.double(exp(lp))

	coef[which.sing] = NA

	list(coefficients  = coef,
	var = var,
	loglik = agfit$loglik,
	score  = agfit$sctest,
	iter   = agfit$iter,
	linear.predictors = as.vector(lp),
	means = agfit$means,
	method= 'coxph',
	x=x,
	y=y)
}

