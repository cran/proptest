"scoreproptest" <-
function(fit,covariate=1,dims=4,basis="legendre",time.transf="F",nsim=1000,nsim.plot=50,weight="unit")
# score process based test of the proportional hazards assumptions
# other (than tested) covariates can (should!) be modeled by smooth functions (dimesnions in dims)
# dims = dimesnions of smooth time-varying covariates other than the tested one (hence, dims[covariate] can be arbitrary, is ignored)
{
	n = fit$n
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
	dims[covariate]=0

	if (nsim.plot>nsim) stop("number of simulations: nsim (=",nsim,") can't be less than nsim.plot (=",nsim.plot,")")

	# I will always use simulations
# 	if ((nvar>1)&(!sim)) {
# 		warning("number of covariates > 1, simulations should be used (use sim=TRUE)")
# 	}

	# suppressWarnings(fit.detail <- coxph.detail(fit));
	suppressWarnings(fit.detail <- mini.coxph.detail(fit));
	nevent = length(fit.detail$time)
	other.timevar = any(dims[-covariate]>0)

	if (other.timevar) { # preparing large model:
		# survfit.coxph = getFromNamespace("survfit.coxph","survival")
		# sfit0 = survfit.coxph(fit)
		sfit0 = mini.survfit.coxph(fit)
		surv0 = sfit0$surv
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
# 				if (basis=="spline") {
# 					psi = bs(psitimes,df=max(dims),degree=1,Boundary.knots=c(min(0,psitimes[1]),psitimes[length(psitimes)]))
# 				}
			}
		}

		endtimes = y[,ncol(y)-1]
		r = sort(endtimes,index.return=TRUE)$ix
		endtimes = endtimes[r]
		endtimes0 = c(0,endtimes)
		jumptimes = sfit0$time
		psi2 = matrix(0,n,max(dims))
		for (i in 1:n) psi2[i,] = psi[nevent-which.max(endtimes[i]>=jumptimes[nevent:1])+1,]  #defining psi-basis even at censoring timepoints
		psi2 = cbind(1,psi2) # v prvnim sloupci 1 (pro puvodni kovariatu)
		x = matrix(x[r,],n,nvar)
		event = y[r,ncol(y)]
		ymax = matrix(0,n*(n+1)/2,3)
		xmax.other = matrix(0,n*(n+1)/2,sum(dims)+nvar-1)
		xmax.covariate = rep(0,n*(n+1)/2)

		froms = cumsum(1+c(0,dims[-covariate],dims[covariate])) # zacatky kovariat (odpovidajicih sloupcu v xmax.other); zacina od 0-te hladke fce, tj. od puvod. kovar.; testovana kov. na konci
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
			xmax.covariate[u]=x[i,covariate]
			ymax[u,1] = endtimes0[1:i]
			ymax[u,2] = endtimes[1:i]
			ymax[u,3] = c(rep(0,i-1),event[i])
			v=v+i
		}
		wh = which(ymax[,1]!=ymax[,2])
		ymax = ymax[wh,]
		ymax = Surv(ymax[,1],ymax[,2],ymax[,3])
		xmax.other = xmax.other[wh,]
		xmax.covariate=xmax.covariate[wh]

		# coxph with x=T,y=T makes coxph.detail 3 times faster
		# fitmax = coxph(ymax~cbind(xmax.other,xmax.covariate),x=TRUE,y=TRUE)
		# but mini.agreg.fit is 3x faster than coxph
		fitmax = mini.agreg.fit(y=ymax,x=cbind(xmax.other,xmax.covariate))
		covariatemax = froms[nvar] # covariatemax is the index of the column with tested covariate which is the last
		nvarmax = covariatemax
		# suppressWarnings(fitmax.detail <- coxph.detail(fitmax))
		suppressWarnings(fitmax.detail <- mini.coxph.detail(fitmax))

		# the large model is now ready

	} else { # other.timevar is FALSE, no large model
		fitmax = fit
		covariatemax = covariate
		nvarmax = nvar
		fitmax.detail = fit.detail
	}

	# L = weight for weighted tests (Marzec & Marzec, 1997)
	if (weight=="gehan") {
		L = fit.detail$nrisk/nevent
	} else {
		if (weight=="mm") {
			L = (fit.detail$nrisk-fit.detail$nrisk[nevent])^2/(nevent^2)
		} else {
			L = rep(1,nevent)
		}
	}

	#print(array(fit.detail$imat,c(nvar,nvar,n)))
	##sigma2 = cumsum(array(fit.detail$imat,c(nvar,nvar,nev))[covariate,covariate,])
	imat.end.inv = fitmax$var ## sigma11 = sigma[1:nvar,1:nvar] # OK
	imat = array(fitmax.detail$imat,c(nvarmax,nvarmax,nevent))
	sigma2 = cumsum((L^2)*imat[covariatemax,covariatemax,])
	sigma2.end = sigma2[nevent]
	sigma.end = sqrt(sigma2.end)
	qqq = sigma2/sigma2.end
	dqqq = diff(c(0,qqq))
	first.last = range(which(dqqq!=0)) #the first and last point with nonzero increment of variance (=> no division by zero in AD)
	#imat = apply(array(fit.detail$imat,c(nvar,nvar,n)),1:2,sum)
	
	dU=as.matrix(fitmax.detail$score)
	#U = apply(,2,cumsum)
	score.process = apply(L*dU,2,cumsum)[,covariatemax]
	stat.ks = max(abs(score.process))/sigma.end
	edge.cut = first.last[1]:(first.last[2]-1) # edge points with zero increments of variance removed (because of division in AD)
	stat.cm = sum(score.process[edge.cut]^2*dqqq[edge.cut])/(sigma2.end)
	stat.ad = sum(score.process[edge.cut]^2*dqqq[edge.cut]/qqq[edge.cut]/(1-qqq[edge.cut]))/(sigma2.end)
	
	# simulation approximation (LWY, 1993)
	for (t in 1:nevent) imat[,,t] = imat[,,t]*L[t]
	imat = apply(imat,c(1,2),cumsum) # attention: now time is the 1st index (caused by apply)
	imat.imat.end.inv = array(0,c(nvarmax,nvarmax,nevent)) # = imat x imat.end.inv; here time will be the 3rd index
	for (t in 1:nevent) imat.imat.end.inv[,,t] = imat[t,,] %*% imat.end.inv
		temp = .C("lwy_scoreproptest",
				as.double(L*dU),
				as.integer(nevent),
				as.integer(nvarmax),
				#as.double(L),
				as.double(t(imat.imat.end.inv[covariatemax,,])),
				as.double(dqqq),
				as.double(qqq),
				as.integer(c(first.last[1],first.last[2]-1)),
				as.integer(nsim),
				as.integer(covariatemax),
				as.double(stat.ks*sigma.end),
				as.double(stat.cm*sigma2.end),
				as.double(stat.ad*sigma2.end),
				pks=double(1),
				pcm=double(1),
				pad=double(1),
				as.integer(nsim.plot),
				score.process.sim=matrix(0,nevent,nsim.plot),
				PACKAGE="proptest")
		p.ks = temp$pks
		p.cm = temp$pcm
		p.ad = temp$pad
		score.process.sim = temp$score.process.sim

# Original R version of LWY simulations; now in C it's 20 times faster
#		 for (b in 1:nsim) {
#			 U.G = dU*rnorm(nevent)
#			 #sc.G.end = apply(U.G,2,sum)
#			 #U.G =apply(L * U.G,2,cumsum)   # apply is terribly slow
#			 for (j in 1:nvarmax) U.G[,j]=cumsum(L*U.G[,j])	# for is 40 % faster
#			 sc.G.end=U.G[nevent,]
#			 #U.G = U.G - t( apply( imat.imat.end.inv, 3, function(x) x%*%sc.G.end ) )   # this is sooo sloooow
# #			for (i in 1:nevent) U.G[i,] = U.G[i,] - imat.imat.end.inv[,,i]%*%sc.G.end	  # this spares over 33% of total time; 4x faster
#			 for (i in 1:nevent) score.process.G[i] = U.G[i,covariatemax] - sum(imat.imat.end.inv[covariatemax,,i]*sc.G.end)	  # this spares a bit more
#			 
# #			score.process.G = U.G[,covariatemax]
#			 stat.G.ks[b] = max(abs(score.process.G))/sigma.end
#			 stat.G.cm[b] = sum(score.process.G[edge.cut]^2*dqqq[edge.cut])/(sigma2.end)
#			 stat.G.ad[b] = sum(score.process.G[edge.cut]^2*dqqq[edge.cut]/qqq[edge.cut]/(1-qqq[edge.cut]))/(sigma2.end)
#		 }
#		 p.simul.ks = 1-rank(c(stat.ks,stat.G.ks))[1]/(nsim+1)
#		 p.simul.cm = 1-rank(c(stat.cm,stat.G.cm))[1]/(nsim+1)
#		 p.simul.ad = 1-rank(c(stat.ad,stat.G.ad))[1]/(nsim+1)

	time = fitmax.detail$time
	out = list(score.process = score.process,
		time = time,
		covariate = covariate,
		stat.ks = stat.ks,
		p.ks = p.ks,
		stat.cm = stat.cm,
		p.cm = p.cm,
		stat.ad = stat.ad,
		p.ad = p.ad,
		nsim = nsim,
		score.process.sim = score.process.sim,
		nsim.plot = nsim.plot,
		dims = dims,
		basis = basis,
		time.transf = time.transf,
		weight = weight)

	class(out) = "scoreproptest"
	out
}

