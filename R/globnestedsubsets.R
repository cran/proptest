"globnestedsubsets" <-
function(d)
# d = vector (of length p, say)
# the number of returned sets (nrow of matrix) is prod(d+1)
# class of subsets containing all possible combinations of nested subsets for
# each of p covariates
# used in the "nested subset" variant of the _global_ data-driven PH test
{
	if (length(d)==1)
		a = lower.tri(matrix(0,d+1,d+1),diag=FALSE)[(d+1):1,1:d]
	else {
		a = matrix(NA,0,sum(d))
		for (k in d[1]:0)
 			a = rbind(a, cbind( cbind(matrix(TRUE,prod(d[-1]+1),k),matrix(FALSE,prod(d[-1]+1),d[1]-k)), Recall(d[-1])) )
	}
	a
}

