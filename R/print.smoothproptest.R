"print.smoothproptest" <-
function(x,detail=FALSE,...)
{
# detail: print alternative models and their scores?
	if (!inherits(x,"smoothproptest"))
		stop("must be an object of class 'smoothproptest'")
	if (x$data.driven) {
		cat("\nData-driven smooth test of proportional hazards\n")

		cat("\nBasis of smooth functions: ",x$basis,sep="")
		cat("\nDimensions for covariates:")
		dims = matrix(c(1:length(x$dims),x$dims),nrow=2,byrow=TRUE)
		colnames(dims) = rep("",length(x$dims))
		rownames(dims)=c("   Covariate   ","   Dimension")
		print(dims,quote=FALSE)
		if (x$global) {
			cat("Global test")
		} else {
			cat("Tested covariate: ",x$covariate,",   max. dimension: ",x$d,sep="")
		}
# 		cat(ifelse(x$all.subsets,"\nUsed BIC for all subsets","\nUsed BIC for nested subsets"))
		if (detail) {
			cat("\nAll",x$nalt,"alternative models:\n")
			selected = rep(NA,(x$nalt))
			selected[x$S] = "<="
			altlabels = matrix("",x$nalt,x$d)
			if (x$global) {
				colnames(altlabels) = rep("",sum(x$dims))
				k=1
				for (i in 1:x$nvar) for (j in 1:x$dims[i]) {colnames(altlabels)[k] = paste(i,",",j,sep=""); k = k+1}
			} else {
				colnames(altlabels) = paste("phi_",1:x$d,sep="")
			}
			for (i in 1:x$nalt) for(j in 1:x$d) altlabels[i,j] = ifelse(x$alt[i,j],"*","")
			m = cbind(altlabels,round(cbind(scorestat=x$scorestats,scorestat.penalised=x$scorestats.penal),digits=4),selected)
			rownames(m) = rep("",x$nalt)
			print(m,na.print="",quote=FALSE,width=6)
			if (x$global) {
		 		cat("Covariates in selected alternative:",paste(rep("phi_{",sum(x$dims)),colnames(altlabels),rep("}",sum(x$dims)),sep="")[x$alt[x$S,]],sep=" ")
			} else {
		 		cat("Covariates in selected alternative:",colnames(altlabels)[x$alt[x$S,]],sep=" ")
			}
		} else {
			cat("\n")
		}
 		cat("\nScore test statistic: ",x$stat,sep="")
 		cat("\np-value: ",format.pval(x$p),"   (based on two-term approx.",ifelse(x$global,paste(" with",x$nsim,"simulations"),""),")",sep="")

	} else {
		cat("\nSmooth test of proportional hazards\n")
		cat("\nBasis of smooth functions: ",x$basis,sep="")
		cat("\nDimensions for covariates:")
		dims = matrix(c(1:length(x$dims),x$dims),nrow=2,byrow=TRUE)
		colnames(dims) = rep("",length(x$dims))
		rownames(dims)=c("   Covariate   ","   Dimension")
		print(dims,quote=FALSE)
		if (x$global)
			cat("Global test")
		else {
			cat("Tested covariate: ",x$covariate,",   dimension: ",x$d,sep="")
		}
		cat("\nTest statistic: ",x$stat,sep="") #" (",x$d," df)",sep="")
		cat("\np-value: ",format.pval(x$p),"   (based on chi^2 with ",x$d," df)",sep="")
	}
	cat("\n\n")
	invisible()
}

summary.smoothproptest = function(object, ...)
{
	if (!inherits(object,"smoothproptest"))
		stop("must be an object of class 'smoothproptest'")
	
	print.smoothproptest(object,detail=TRUE,...)

	invisible()
}
