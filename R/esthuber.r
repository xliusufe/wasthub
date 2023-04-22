esthub <- function(x, y, method = "adaptive", maxIter = 100, tol = 0.00001, M = 100, tauk = NULL) {
	n 	= length(y)
	p 	= ifelse(is.null(ncol(x)) , 1, ncol(x))

	scal_x 	= 1/sqrt(colSums(x^2))
	x 		= x*matrix(rep(scal_x,each=n), nrow=n, ncol=p)

	tau0 	= 1.345
	if(method == "adaptive"){
		qq 	= -1
	}
	else if(method == "CV"){
		if(is.null(tauk)){
			tauk 	= seq(0.01,10,length.out = M)
		}
		y1 		= y - mean(y)

		fit 	= lm(y1~x-1)
		halpha 	= fit$coefficients
		hsigma2 = sqrt(sum( (y1 - x%*%halpha)^2 )/(n-p))

		loss 	= rep(1, length(tauk))
		for(k in 1:length(tauk)){
			qq 	= hsigma2*tauk[k]*sqrt(n/log(n))
			fit <- .Call("_EST_HUBER_2step",
						as.numeric(y),
						as.numeric(x),
						as.integer(c(n, p, maxIter)),
						as.numeric(c(qq, tau0, tol))
					)
			loss[k] = sum( (y1 - x%*%fit$coef)^2 )
		}
		ko 	= which.min(loss)
		qq 	= hsigma2*tauk[ko]*sqrt(n/log(n))
	}
	else{
		qq 	= tau0*median( abs(y - median(y)) )/qnorm(0.75)
	}

	fit <- .Call("_EST_HUBER_2step",
				as.numeric(y),
				as.numeric(x),
				as.integer(c(n, p, maxIter)),
				as.numeric(c(qq, tau0, tol))
			)

	halpha 	= fit$coef
	tau 	= fit$tau
	hsigma2 = sum( (y-x%*%halpha)^2 )/(n-p)
	halpha 	= halpha*scal_x

	return(list(alpha	= halpha,
				tau 	= tau,
				sigma2 	= hsigma2
				)
			)
}