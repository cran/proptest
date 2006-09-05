"print.scoreproptest" <-
function(x,...)
{
	cat("\nScore process test of proportional hazards\n")
	cat("\nTested covariate: ",x$covariate,"\n",sep="")
	if (length(x$dims)>1) {
		cat("Basis of smooth functions for other covariates: ",x$basis,sep="")
		cat("\nDimensions for other covariates:")
		dims = as.matrix(matrix(c((1:length(x$dims)),x$dims),nrow=2,byrow=TRUE)[,-x$covariate])
		colnames(dims) = rep("",length(x$dims)-1)
		rownames(dims)=c("   Covariate   ","   Dimension")
		print(dims,quote=FALSE)
	}
	cat("Kolmogorov--Smirnov statistic: ",x$stat.ks,",   p-value: ",format.pval(x$p.ks),sep="")
	cat("\nCramer--von Mises statistic: ",x$stat.cm,",   p-value: ",format.pval(x$p.cm),sep="")
	cat("\nAnderson--Darling statistic: ",x$stat.ad,",   p-value: ",format.pval(x$p.ad),sep="")
	cat("\np-values based on",x$nsim,"simulations")
	cat("\n\n")
	invisible(x)
}

