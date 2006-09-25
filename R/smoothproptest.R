"smoothproptest" <-
function(fit,covariate=1,dims=4,basis="legendre",time.transf="F",data.driven=TRUE,all.subsets=FALSE,h.approx=TRUE,sim=FALSE,nsim=1000)
# Neyman type smooth test of the proportional hazards assumptions
# other (than tested) covariates can (should!) be modeled by smooth functions too (dimesnions in dims)
{
	n = fit$n
	logn = log(n)
	x = fit$x
	y = fit$y
	if (is.null(y)  ||  is.null(x)) {
		temp = getFromNamespace("coxph.getdata","survival")(fit, y=TRUE, x=TRUE)
		fit$y = temp$y
		fit$x = temp$x
		x = fit$x
		y = fit$y
	}
	nvar = ncol(x)

	if (attributes(y)$type!="right") stop("data type must be 'right' censored")
##	if (all.subsets & (!sim)) stop("when all subsets used, simulations must be done (use sim=TRUE)")
	if (data.driven && (!h.approx) && (!sim)) stop("at least one of H-approximation (h.approx=TRUE) and simulation (sim=TRUE) must be used")
	if (!data.driven) h.approx=FALSE
	if (length(dims)==1) {
		dims=rep(dims,nvar)
	} else {
		if (length(dims)!=nvar) stop("length of d must equal number of covariates or 1")
	}
	d=dims[covariate]

	# test statistics, pvalues,...; all of them set to NA, then some of them computed
	stat.d = stat.bic = NULL
	p.d = p.d.sim = p.bic.h = p.bic.sim = p.bic.w = p.bic.chisq1 = NULL
	wts = S = scorestats = scorestats.penal = NULL


	other.timevar = any(dims[-covariate]>0) # is there any covariate to be modeled smoothly?  (among those not to be tested)

	control.iter0 = coxph.control(iter.max=0)

	# building the basis based on the original model (all covariates assumed proportional):
	sfit0 = mini.survfit.coxph(fit)
	surv0 = sfit0$surv
	nevent = length(surv0)
	if (time.transf=="F") {
		psitimes = (1-surv0)/(1-surv0[nevent])
	} else {
		if (time.transf=="L") {
			psitimes = log(surv0)/log(surv0[nevent])
		} else {
			stop("unknown time transformation: ",time.transf)
		}
	}

	if (basis=="cos") basis = "cosine"
	if (basis=="legendre") {
		psi = poly(psitimes,degree=max(dims))
	} else {
		if (basis=="cosine") {
			psi = cos(matrix(psitimes,nrow=nevent,ncol=max(dims))*matrix(1:max(dims),nrow=nevent,ncol=max(dims),byr=TRUE)*pi)
		} else {
			stop("unknown basis of smooth functions: ",basis)
# 			if (basis=="spline") {
# 				psi = bs(psitimes,df=max(dims),degree=1,Boundary.knots=c(min(0,psitimes[1]),psitimes[length(psitimes)]))
# 			}
		}
	}

	endtimes = y[,ncol(y)-1]
	r = sort(endtimes,index.return=TRUE)$ix
	endtimes = endtimes[r]
	endtimes0 = c(0,endtimes)
	jumptimes = sfit0$time ## lam0[,2]
	psi2 = matrix(0,n,max(dims))
	for (i in 1:n) psi2[i,] = psi[nevent-which.max(endtimes[i]>=jumptimes[nevent:1])+1,]  #defining psi-basis even at censoring timepoints
	psi2 = cbind(1,psi2)
	x = matrix(x[r,],n,nvar)
	event = y[r,ncol(y)]
	ymax = matrix(0,n*(n+1)/2,3)
	xmax.other = matrix(0,n*(n+1)/2,sum(dims)-d+nvar-1)
	xmax.covariate = matrix(0,n*(n+1)/2,d+1)

	froms = cumsum(1+c(0,dims[-covariate],dims[covariate])) # beginnings of covariates (corresp.onding columns in xmax.other);
	                                                          # beginning from the 0-th smoth fnc, i.e., original covar.; tested covar. last
	other=(1:nvar)[-covariate]

	# building time-dependent covariates (for the other covariates, not the tested one):
	v = 0
	for (i in 1:n) {
		u = 1:i+v
		if (nvar>1) {
			for (j in 1:(nvar-1)) {
				xmax.other[u,froms[j]:(froms[j+1]-1)]=x[i,other[j]]*psi2[1:i,1:(dims[other[j]]+1)]
			}
		}
		# xmax.covariate[u,]=x[i,covariate]*psi2[1:i,1:(d+1)]	
		xmax.covariate[u,1]=x[i,covariate]
		ymax[u,1] = endtimes0[1:i]
		ymax[u,2] = endtimes[1:i]
		ymax[u,3] = c(rep(0,i-1),event[i])
		v=v+i
	}
	wh = which(ymax[,1]!=ymax[,2])
	ymax = ymax[wh,]
	ymax = Surv(ymax[,1],ymax[,2],ymax[,3])
	xmax.other = xmax.other[wh,]

	if (nvar==1) {
		fitmax = fit
	} else {
		fitmax = mini.agreg.fit(x=cbind(xmax.other,xmax.covariate[wh,1]),y=ymax)
	}
	betamax = fitmax$coef 

	# if more than 1 covariate and some is to be modeled smoothly, the time transformation
	# for the tested covariate is based on the larger smooth model:
	if ((nvar>1) && other.timevar) {
		# mini.survfit.coxph can handle output of (mini.)agreg.fit, survfit.coxph can't
		# sfitmax = survfit.coxph(fitmax)
		sfitmax = mini.survfit.coxph(fitmax)
		survmax = sfitmax$surv
		if (time.transf=="F") {
			psitimes = (1-survmax)/(1-survmax[nevent])
		} else {
			if (time.transf=="L") {
				psitimes = log(survmax)/log(survmax[nevent])
			} else {
				stop("unknown time transformation: ",time.transf)
			}
		}
		if (basis=="cos") basis = "cosine"
		if (basis=="legendre") {
			psi = poly(psitimes,degree=d)
		} else {
			if (basis=="cosine") {
				psi = cos(matrix(psitimes,nrow=nevent,ncol=d)*matrix(1:d,nrow=nevent,ncol=d,byr=TRUE)*pi)
			} else {
# 				if (basis=="spline") {
# 					psi = bs(psitimes,df=d,degree=1,Boundary.knots=c(min(0,psitimes[1]),psitimes[length(psitimes)]))
# 				}
			}
		}
		psi2 = matrix(0,n,d)
		for (i in 1:n) psi2[i,] = psi[nevent-which.max(endtimes[i]>=jumptimes[nevent:1])+1,]
		psi2 = cbind(1,psi2)
	}

	# now the time transf. is ready whatever nvar is; let's compute the smooth functions for the tested covar.
	v = 0
	for (i in 1:n) {
		u = 1:i+v
		xmax.covariate[u,]=x[i,covariate]*psi2[1:i,]
		v = v+i
	}
	xmax.covariate = xmax.covariate[wh,]

	if (data.driven) {
		if (all.subsets) {
			nalt = 2^d-1 # number of alternative models
			alt = allsubsets(d)[1:nalt,] # alternative models
		} else {
			nalt = d
			alt = lower.tri(matrix(0,d,d),diag=TRUE)[d:1,]
		}
	} else {
		nalt = 1
		alt = matrix(TRUE,1,d)
	}

	dimalt = apply(alt,1,sum)
	scorestats = rep(0,nalt)
	V.inv = list() # list of inverse V matrices for all alternatives (needed in simulations)
	sort.end = order(-ymax[,2], ymax[,3]) # to avoid repeated sorting in each run of mini.agreg.fit
	sort.start = order(-ymax[,1])		 # ...
	for (k in nalt:1) { # from nalt downto 1, because alt[1,] is the largest
		betax.init = c(betamax,rep(0,dimalt[k]))
		# coxph was slow:
		# fitx = coxph(ymax~cbind(xmax.other,xmax.covariate[,c(T,alt[k,])]),init=betax.init,iter.max=0,x=TRUE,y=TRUE)
		# instead: direct call to agreg.fit is 3 times faster:
		fitx = mini.agreg.fit(x=cbind(xmax.other,xmax.covariate[,c(TRUE,alt[k,])]), y=ymax, init=betax.init, control=control.iter0, sort.start=sort.start, sort.end=sort.end)
		# what follows was about 15% slower than mini.agreg.fit:
		# fitx = agreg.fit(x=cbind(xmax.other,xmax.covariate[,c(T,alt[k,])]), y=ymax, strata=NULL, init=betax.init,
		#				   control=control, method="efron", rownames="")
		scorestats[k] = fitx$score
		V.inv[[k]] = as.matrix(fitx$var[(froms[nvar]+1):(froms[nvar]+dimalt[k]),(froms[nvar]+1):(froms[nvar]+dimalt[k])])
	}

	scorestats.penal = scorestats - dimalt*logn
	S = which.max(scorestats.penal)
	stat.bic = scorestats[S]
	stat.d = scorestats[1]   # fixed dimension test statistic T_d
	p.d.chisqd = 1-pchisq(stat.d,d)

	if (sim) { # LWY simulation approximation
		suppressWarnings(fitx.d.detail <- mini.coxph.detail(fitx)) # mini.coxph.detail is 15% faster than coxph.detail
		#D = as.matrix(apply(array(fitx.d.detail$imat[(froms[nvar]+1):(froms[nvar]+dimalt[k]),1:froms[nvar],],c(d,froms[nvar],nevent)),c(1,2),sum))

		imat = solve(fitx$var) # imat in the largest alternative
		D = imat[(froms[nvar]+1):(froms[nvar]+dimalt[k]),1:froms[nvar]] # (2,1) submatrix of imat
		J.inv = as.matrix(fitmax$var) # is equivalent to inv of submatrix of imat: solve(imat[1,1]))
		D.J.inv = D %*% J.inv
		nvartotal = sum(dims)+nvar
		dU=fitx.d.detail$score
		temp=.C("lwy_smoothproptest",
				 as.double(dU),
				 as.integer(nevent),
				 as.integer(nvartotal),
				 as.double(D.J.inv),
				 as.double(unlist(V.inv)),
				 as.integer(alt),
				 as.integer(nalt),
				 as.integer(nsim),
				 as.integer(froms[nvar]),
				 as.double(stat.d),		# T_d
				 as.double(stat.bic),		# T_S
				 as.double(logn),
				 p.d.sim=double(1),		# simul. pval. for T_d
				 p.bic.sim=double(1),		# simul. pval. for T_S
				 wts=double(d),		# simul. weights for mixture below
				 p.bic.w=double(1),		# chisq-mixture pval. for T_S
				 PACKAGE="proptest")
		p.d.sim=temp$p.d.sim
		p.bic.sim=temp$p.bic.sim
		wts=temp$wts
		p.bic.w=temp$p.bic.w
	} # end of LWY simulation

	if (h.approx) { # two-term approximation H (can be used only for nested alternatives, for all subsets doesn't work)
		if (!all.subsets) {
			h1 = function(x,n) {
				(2*pnorm(sqrt(x))-1)*(2*pnorm(sqrt(log(n)))-1)
			}
			h2 = function(x,n) {
				(2*pnorm(sqrt(x))-1)*(2*pnorm(sqrt(log(n)))-1) + 2*(1-pnorm(sqrt(log(n))))
			}
			nu = 2 ## 2
			if (stat.bic<=logn) {
				p.bic.h = 1-h1(stat.bic,n)
			} else {
				if (stat.bic>=nu*logn) {
					p.bic.h = 1-h2(stat.bic,n)
				} else {
					p.bic.h = 1 - ( (nu*logn-stat.bic)/((nu-1)*logn)*h1(logn,n) + (stat.bic-logn)/((nu-1)*logn)*h2(nu*logn,n) )
				}
			}
		} else { # all subsets
			V.inv.1 = V.inv[dimalt==1]
			V.inv.2 = V.inv[dimalt==2]
			chol.V = chol(solve(V.inv[[1]]))
			temp=.C("h_approx_sim",
				as.double(stat.bic),
				as.double(unlist(V.inv.1)),
				as.double(unlist(V.inv.2)),
				as.double(chol.V),
				as.integer(d),
				as.double(logn),
				as.integer(nsim),
				h=double(1),
				PACKAGE="proptest")
			p.bic.h = 1-temp$h
		}
	} # end of H approximation

	if (data.driven && (!all.subsets)) {
		p.bic.chisq1 = 1-pchisq(stat.bic,1)
	}

	out = list(stat.d = stat.d,
		stat.bic = stat.bic,
		p.d.chisqd = p.d.chisqd,
		p.d.sim = p.d.sim,
		p.bic.h = p.bic.h,
		p.bic.sim = p.bic.sim,
		p.bic.w = p.bic.w,
		p.bic.chisq1 = p.bic.chisq1,
		wts=wts)
	if (data.driven) {
		out$stat = stat.bic
		out$p = ifelse(h.approx, p.bic.h, p.bic.sim)
	} else {
		out$stat = stat.d
		out$p = p.d.chisqd
	}

	out$data.driven = data.driven
	out$covariate = covariate
	out$dims = dims
	out$d = d
	out$all.subsets=all.subsets
	out$basis = basis
	out$time.transf=time.transf
	out$n = n
	out$nvar = nvar

	out$h.approx = h.approx
	out$sim = sim
	if (sim) out$nsim=nsim
	out$S = S
	out$scorestats = scorestats
	out$scorestats.penal = scorestats.penal
	out$nalt = nalt
	out$alt = alt

	class(out) = "smoothproptest"
	out
}

