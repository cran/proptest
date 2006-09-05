"allsubsets" <-
function(n)
# all subsets of the set of n elements; output: logical 2^n by n matrix
# used in the all subset variant of the data-driven PH test
{
	if(n > 0) rbind(cbind(TRUE,Recall(n-1)), cbind(FALSE,Recall(n-1)))
}

