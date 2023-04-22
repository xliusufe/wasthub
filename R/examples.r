exams <- function(method = "wast", B = 1000, K = 1000){
	data(simulatedData_gaussian)
	pvals   = pvalhuber(data = data_gaussian, method = method, B = B, K = K)

	data(simulatedData_quantile)
	pvals   = pvalhuber(data = data_quantile, method = method, B = B, K = K)
	return(pvals)
}