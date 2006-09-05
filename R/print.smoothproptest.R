"print.smoothproptest" <-
function(x,print.alt=FALSE,...)
{
# alt: print alternative models and their scores?
	if (x$data.driven) {
		cat("\nData-driven smooth test of proportional hazards\n")

		cat("\nBasis of smooth functions: ",x$basis,sep="")
		cat("\nDimensions for covariates:")
		dims = matrix(c(1:length(x$dims),x$dims),nrow=2,byrow=TRUE)
		colnames(dims) = rep("",length(x$dims))
		rownames(dims)=c("   Covariate   ","   Dimension")
		print(dims,quote=FALSE)
		cat("Tested covariate: ",x$covariate,",   max. dimension: ",x$d,sep="")
		cat(ifelse(x$all.subsets,"\nUsed BIC for all subsets","\nUsed BIC for nested subsets"))
		if (print.alt) {
			cat("\nAll",x$nalt,"alternative models:\n")
			selected = rep(NA,(x$nalt))
			selected[x$S] = "<="
			altlabels = matrix("",x$nalt,x$d)
			colnames(altlabels) = paste("phi_",1:x$d,sep="")
			for (i in 1:x$nalt) for(j in 1:x$d) altlabels[i,j] = ifelse(x$alt[i,j],"*","")
			m = cbind(altlabels,round(cbind(scorestat=x$scorestats,scorestat.penalized=x$scorestats.penal),digits=4),selected)
			rownames(m) = rep("",x$nalt)
			print(m,na.print="",quote=FALSE,width=6)
		} else {
			cat("\n")
		}
 		cat("Covariates in selected alternative:",paste(rep("phi_",x$d),1:x$d,sep="")[x$alt[x$S,]],sep=" ")
 		cat("\nScore test statistic: ",x$stat,sep="")
 		cat("\np-value: ",format.pval(x$p),"   (based on ",ifelse(x$sim,paste(x$nsim,"simulations"),"H-approximation"),")",sep="")

	} else {
		cat("\nSmooth test of proportional hazards\n")
		cat("\nBasis of smooth functions: ",x$basis,sep="")
		cat("\nDimensions for covariates:")
		dims = matrix(c(1:length(x$dims),x$dims),nrow=2,byrow=TRUE)
		colnames(dims) = rep("",length(x$dims))
		rownames(dims)=c("   Covariate   ","   Dimension")
		print(dims,quote=FALSE)
		cat("Tested covariate: ",x$covariate,",   dimension: ",x$d,sep="")
		cat("\nScore test statistic: ",x$stat,sep="") #" (",x$d," df)",sep="")
		cat("\np-value: ",format.pval(x$p),"   (based on ",ifelse(x$sim,paste(x$nsim,"simulations"),paste("chi^2 with",x$d,"df")),")",sep="")
	}
	cat("\n\n")
	invisible(x)
}

