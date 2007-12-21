"smoothproptest" <-
function(fit,covariate=1,global=FALSE,dims=4,basis="legendre",time.transf="F",data.driven=TRUE,nsim=1000)
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
	if (length(dims)==1) {
		dims=rep(dims,nvar)
	} else {
		if (length(dims)!=nvar) stop("length of d must equal number of covariates or 1")
	}
	
	if (global) {
		d = sum(dims)
	} else {
		d = dims[covariate]
	}

	# test statistics, pvalues,...; all of them set to NA, then some of them computed
	stat.d = stat.bic = NULL
	p.d = p.bic.h = p.bic.asympt = NULL
	S = scorestats = scorestats.penal = NULL


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
	
	if (global) { # to keep it similar with individual covariate test, we imagine that covar. 1 is tested and put it into xmax.covariate
		xmax.original = matrix(0,n*(n+1)/2,nvar)
		xmax.smooth = matrix(0,n*(n+1)/2,sum(dims))
# 		froms = nvar+cumsum(c(0,dims)) # beginnings of smooth functions covariates 
# 		other=(1:nvar)[-covariate]
		v = 0
		for (i in 1:n) {
			u = 1:i+v
			from = 1
			for (j in 1:nvar) {
				xmax.original[u,j] = x[i,j]
				if (dims[j]>0) {
					xmax.smooth[u,from:(from+dims[j]-1)]=x[i,j]*psi2[1:i,2:(dims[j]+1)]
					from = from+dims[j]
				}
			}
			# xmax.covariate[u,]=x[i,covariate]*psi2[1:i,1:(d+1)]	
			ymax[u,1] = endtimes0[1:i]
			ymax[u,2] = endtimes[1:i]
			ymax[u,3] = c(rep(0,i-1),event[i])
			v=v+i
		}

		# let's make ymax look like a 'Surv' object of type "counting"
		class(ymax) = "Surv"
		attr(ymax, "type") = "counting"
		# remove intervals of zero length (tied obs.)
		wh = (ymax[,1]==ymax[,2]) # tied observations, to be removed
		# but we must avoid removing events if they are tied with censorings
		v = 0
		for (i in 1:n) {
			j = i
			while ((wh[v+j])&&(j>0)) { # for each subject go through ties at the end
				j = j-1
			}
			# now j is the last nonzero interval for the ith subject
			ymax[v+j,3] = event[i]
			v = v+i
		}
		wh = !wh # intervals of non-zero length
		
		ymax = ymax[wh,]
# 		ymax = Surv(ymax[,1],ymax[,2],ymax[,3])
		xmax.original = xmax.original[wh,]
		xmax.smooth = xmax.smooth[wh,]

		fitmax = fit
		betamax = fitmax$coef
	} else {
		xmax.other = matrix(0,n*(n+1)/2,sum(dims)-d+nvar-1)
		xmax.covariate = matrix(0,n*(n+1)/2,d+1)
		froms = cumsum(1+c(0,dims[-covariate],dims[covariate])) # beginnings of covariates (corresp.onding columns in xmax.other);
	                                                          # beginning from the 0-th smooth fnc, i.e., original covar.; tested covar. last
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
		
		# let's make ymax look like a 'Surv' object of type "counting"
		class(ymax) = "Surv"
		attr(ymax, "type") = "counting"
		# remove intervals of zero length (tied obs.)
		wh = (ymax[,1]==ymax[,2]) # tied observations, to be removed
		# but we must avoid removing events if they are tied with censorings
		v = 0
		for (i in 1:n) {
			j = i
			while ((wh[v+j])&&(j>0)) { # for each subject go through ties at the end
				j = j-1
			}
			# now j is the last nonzero interval for the ith subject
			ymax[v+j,3] = event[i]
			v = v+i
		}
		wh = !wh # intervals of non-zero length

		ymax = ymax[wh,]
#		ymax = Surv(ymax[,1],ymax[,2],ymax[,3])
		xmax.other = xmax.other[wh,]

		if (nvar==1) {
			fitmax = fit
		} else {
			fitmax = mini.agreg.fit(x=cbind(xmax.other,xmax.covariate[wh,1]),y=ymax)
		}
		betamax = fitmax$coef
		
		
		# if more than 1 covariate and some is to be modeled smoothly and test is not global,
		# the time transformation for the tested covariate is based on the larger smooth model:
		if ((nvar>1) && other.timevar) {
			# mini.survfit.coxph can handle output of (mini.)agreg.fit, survfit.coxph can't
			# sfitmax = survfit(fitmax)
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
		v = 0
		for (i in 1:n) {
			u = 1:i+v
			xmax.covariate[u,]=x[i,covariate]*psi2[1:i,]
			v = v+i
		}
		xmax.covariate = xmax.covariate[wh,]
	}

	if (data.driven) {
 		if (global) {
# 		if (all.subsets) {
# 			nalt = 2^d-1 # number of alternative models
# 			alt = allsubsets(d)[1:nalt,] # alternative models
			# nested-like subsets
			nalt = prod(dims+1)-1 # number of alternative models
			alt = globnestedsubsets(dims)[1:nalt,] # alternative models
		} else {
			nalt = d
			alt = matrix(lower.tri(matrix(0,d,d),diag=TRUE)[d:1,],d,d)
		}
	} else {
		nalt = 1
		alt = matrix(TRUE,1,d)
	}
	
	dimalt = rowSums(alt)
	scorestats = rep(0,nalt)
	V.inv = list() # list of inverse V matrices for all alternatives (needed in simulations)
	sort.end = order(-ymax[,2], ymax[,3]) # to avoid repeated sorting in each run of mini.agreg.fit
	sort.start = order(-ymax[,1])		 # ...
	if (global) { # for the global test, we compute the score test statistic with matrix of xmax.original and xmax.smooth
		for (k in nalt:1) { # from nalt downto 1, because alt[1,] is the largest
			betax.init = c(betamax,rep(0,dimalt[k]))
			fitx = mini.agreg.fit(x=cbind(xmax.original,xmax.smooth[,alt[k,]]), y=ymax, init=betax.init, control=control.iter0, sort.start=sort.start, sort.end=sort.end)
			scorestats[k] = fitx$score
			V.inv[[k]] = as.matrix(fitx$var[(nvar+1):(nvar+dimalt[k]),(nvar+1):(nvar+dimalt[k])])
		}
	} else { # for individual test, the matrix consists of xmax.other and 
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
	}


	scorestats.penal = scorestats - dimalt*logn
	S = which.max(scorestats.penal)
	stat.bic = scorestats[S]
	stat.d = scorestats[1]   # fixed dimension test statistic T_d
	p.d = 1-pchisq(stat.d,d)
	
	
	if (data.driven) {
		if (!global) { # not global, thus nested subsets, thus two-term approximation H
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
			p.bic.asympt = 1-pchisq(stat.bic,1) # asymptotic chisq with 1 d.f.
		} else { # global, thus nested-like subsets, thus two-term approx. via simulations
			# nested-like subsets for global tests
			V.inv.1 = V.inv[dimalt==1]
			V.inv.2 = V.inv[dimalt==2]
			s2 = rep(FALSE, d) ## ss will have TRUE on the first two (at most) positions for each covar.
			for (k in 1:nalt) if (dimalt[k]<=2) s2 = s2 | alt[k,]
			V2 = solve(V.inv[[1]])[s2,s2]
			chol.V2 = chol(V2)
			d2 = sum(s2)
			alt2 = alt[(dimalt==1)|(dimalt==2),s2]
			dimalt2 = rowSums(alt2)
			nalt2 = length(dimalt2)
			diag.V2 = diag(V2)
			p.bic.asympt = 0
			p.bic.h = 0
			for (i in 1:nsim) {
				u2.sim = t(chol.V2)%*%rnorm(d2)
				m1 = max(u2.sim*u2.sim/diag.V2)
				p.bic.asympt = p.bic.asympt + (m1<=stat.bic)
				m2 = 0
				k = 0
				for (j in 1:nalt2)
					if (dimalt2[j]==2) {
						k = k+1
						m2 = max(m2, u2.sim[alt2[j,]]%*%V.inv.2[[k]]%*%u2.sim[alt2[j,]])
					}
				p.bic.h = p.bic.h + ((m1<=stat.bic)&&(m2-m1<=logn)) + ((m2<=stat.bic)&&(m2-m1>logn));
			}
			p.bic.asympt = 1 - p.bic.asympt/nsim # asymptotic p-val. (max of dep. chisq1)
			p.bic.h = 1 - p.bic.h/nsim # two-term approx.
		}
	} # end of the two-term approximation for the data-driven test

	out = list(stat.d = stat.d,
		p.d = p.d)
	if (data.driven) {
		out$stat.bic = stat.bic
		out$stat = stat.bic
		out$p = p.bic.h
		out$p.bic.h = p.bic.h
		out$p.bic.asympt = p.bic.asympt
	} else {
		out$stat = stat.d
		out$p = p.d
	}
	
	out$data.driven = data.driven
	out$global = global
	if (!global) {
		out$covariate = covariate
	}
	out$dims = dims
	out$d = d
	out$basis = basis
	out$time.transf=time.transf
	out$n = n
	out$nvar = nvar

	out$nsim=nsim
	out$S = S
	out$scorestats = scorestats
	out$scorestats.penal = scorestats.penal
	out$nalt = nalt
	out$alt = alt

	class(out) = "smoothproptest"
	out
}

