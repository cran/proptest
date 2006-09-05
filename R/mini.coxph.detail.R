"mini.coxph.detail" <-
function(object) {
## object is a coxph-like object such as the output of mini.agreg.fit
## object must contain x, y
	#method = object$method
	method = "efron"
	x = object$x
	y = object$y

	nvar = ncol(x)
	if (ncol(y)==2) {
		mintime = min(y[,1])
		if (mintime < 0) y = cbind( 2*mintime -1, y)
		else 	y = cbind(-1,y)
	}

	ord = order(y[,2], -y[,3])

	# sort the data
	x = x[ord,]
	y = y[ord,]
	storage.mode(y) = 'double'
	if (!is.null(object$lp)) { # object is the output of (mini.)agreg.fit
		n = length(object$lp)
		score = exp(object$lp)[ord]
	} else { # object is the output of coxph
		n = length(object$linear.predictor)
		score = exp(object$linear.predictor)[ord]
	}

	ndeath = sum(y[,3])
	ff = .C("coxdetail", as.integer(n),
		as.integer(nvar),
		ndeath= as.integer(ndeath),
		y = y,
		as.double(x),
		as.integer(rep(0,n)), # newstrat
		index =as.double(score),
		weights = as.double(rep(1,n)), # weights
		means= c(method=='efron', double(ndeath*nvar)),
		u = double(ndeath*nvar),
		i = double(ndeath*nvar*nvar),
		double(nvar*(3 + 2*nvar)),
		PACKAGE="survival")
	keep = 1:ff$ndeath
	time = y[ff$index[keep],2]
	means= (matrix(ff$means,ndeath, nvar))[keep,]
	score=  matrix(ff$u, ndeath, nvar)[keep,]
	var = array(ff$i, c(nvar, nvar, ndeath))[,,keep]
	temp = list(time = time, means=means, nevent=ff$y[keep,1],
	 nrisk = ff$y[keep,2], hazard= ff$y[keep,3], score= score,  imat=var,
	 varhaz=ff$weights[keep], y=y, x=x)
	temp
}

