"mini.survfit.coxph" <-
function(object) {
## object is a coxph-like object such as the output of mini.agreg.fit
## object must contain x, y

	se.fit=FALSE
	if (!is.null(object$lp)) { # object is the output of (mini.)agreg.fit
		n = length(object$lp)
	} else { # object is the output of coxph
		n = length(object$linear.predictor)
	}
	nvar = length(object$coef)
	score = exp(object$linear.predictor)
	
	temp = c('aalen', 'kalbfleisch-prentice', 'efron',
			   'tsiatis', 'breslow', 'kaplan-meier', 'fleming-harringon',
			   'greenwood', 'exact')
	temp2 = c(2,1,3,2,2,1,3,1,1)
	type="efron"
	vartype = type
	method = temp2[match(match.arg(type, temp), temp)]
	vartype = temp2[match(match.arg(vartype, temp), temp)]
	conf.type = 'none'
	y = object$y
	ny = ncol(y)

	type = attr(y, 'type')
	if (type=='counting') {
		if (method=='kaplan-meier')
			stop ("KM method not valid for counting type data")
		ord = order(y[,2], -y[,3]) 
	} else
		if (type=='right') {
			ord = order(y[,1], -y[,2]) 
			miny = min(y[,1])
			if (miny < 0) y = cbind(2*miny -1, y)
			else y = cbind(-1, y)
		} else
			stop("Cannot handle \"", type, "\" type survival data")

	if (!is.null(object$x)) x = object$x[ord,]
	else x = 0
	weights = rep(1,n)
	newstrat = as.integer(rep(0,n))
	newstrat[n] = 1

	stype = 2
	offset2 = 0   #offset variable for the new data set
	x2 = matrix(object$means, nrow=1)
	n2 = nrow(x2)

	coef = ifelse(is.na(object$coef), 0, object$coef)
	newrisk = exp(c(x2 %*% coef) + offset2 - sum(coef*object$means))

	dimnames(y) = NULL   #I only use part of Y, so names become invalid
	storage.mode(y) = 'double'
	ndead = sum(y[,3])
	surv = .C('agsurv2', as.integer(n),
				  as.integer(nvar* se.fit),
				  y = y[ord,],
				  as.double(score[ord]),
				  strata = as.integer(newstrat),
                          wt = as.double(weights),
				  surv = double(ndead*n2),
				  varhaz = double(ndead*n2),
				  as.double(x),
				  as.double(object$var),
				  nsurv = as.integer(c(method, vartype)),
				  double(3*nvar),
				  as.integer(n2),
				  as.double(x2),
				  as.double(newrisk), PACKAGE="proptest")
	nsurv = surv$nsurv[1]
	ntime = 1:nsurv
	if (n2>1) {
		tsurv = matrix(surv$surv[1:(nsurv*n2)], ncol=n2)
		tvar  = matrix(surv$varhaz[1:(nsurv*n2)], ncol=n2)
		dimnames(tsurv) = list(NULL, dimnames(x2)[[1]])
	} else {
		tsurv = surv$surv[ntime]
		tvar  = surv$varhaz[ntime]
	}
	if (surv$strata[1] <=1)
		temp = list(n=n,time=surv$y[ntime,1],
			 n.risk=surv$y[ntime,2],
			 n.event=surv$y[ntime,3],
			 surv=tsurv,
			type=type)
	else {
		temp = surv$strata[1:(1+surv$strata[1])]
		tstrat = diff(c(0, temp[-1])) #n in each strata
		names(tstrat) = levels(object$strata)
		temp = list(n=n, time=surv$y[ntime,1],
			 n.risk=surv$y[ntime,2],
			 n.event=surv$y[ntime,3],
			 surv=tsurv,
			 strata= tstrat, ntimes.strata=tstrat,
			strata.all=NULL,
			type=type)
	}

	class(temp) = c('survfit.cox', 'survfit')
	temp
}

