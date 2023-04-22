estolscp <- function(data, smooth = 'sigmoid', isBoot = FALSE, isWB = FALSE, h = NULL, maxIter = 100, tol = 0.00001, B = 1000) {

	if(!(smooth %in% c('sigmoid', 'pnorm','mixnorm'))){
		stop("Smooth function must be one of {'sigmoid', 'pnorm','mixnorm'} !")
	}
	y 	= data$Y
	n 	= length(y)
	tx 	= data$X
	x 	= data$Z
	z 	= data$U
	p1 	= ifelse(is.null(ncol(tx)), 1, ncol(tx))
	p2 	= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 	= ifelse(is.null(ncol(z)) , 1, ncol(z))

	if(is.null(h)){
		h 	= n^(1/2)/log(n)
	}

	type = switch(smooth,
					'sigmoid' 	= 1,
					'pnorm' 	= 2,
					'mixnorm' 	= 3
					)

	theta0 	= rep(1, p3)
	dims 	= c(n, p1, p2, p3-1, maxIter, type)
	params 	= c(tol, h)
	fit	= .Call("_EST_LinearCP",
				as.numeric(y),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(theta0),
				as.integer(dims),
				as.numeric(params)
			)

	htheta 	= c(1, fit$theta)
	tbeta 	= fit$beta
	halpha 	= tbeta[1:p1]
	hbeta 	= tbeta[-c(1:p1)]
	hdelta 	= z%*%htheta>0
	muhat 	= tx%*%halpha + (x%*%hbeta)*hdelta
	resids 	= y - muhat
	hsigma2 = sum( (y-muhat)^2 )/(n-(p1+p2+p3))

	if(isBoot){
		alphaB = rep(0, B)
		for(k in 1:B){
			if(isWB){
				yb 	= muhat + resids*rnorm(n)
			}
			else{
				yb 	= muhat + resids[sample.int(n, n, replace = TRUE)]
			}
			fit = .Call("_EST_LinearCP",
						as.numeric(yb),
						as.numeric(tx),
						as.numeric(x),
						as.numeric(z),
						as.integer(dims),
						as.numeric(params)
					)
			alphaB[k] = fit$beta[1]
		}
		hsigma2 = var(alphaB)
	}

	return(list(alpha	= halpha,
				beta 	= hbeta,
				theta 	= htheta,
				delta 	= hdelta,
				sigma2 	= hsigma2
				)
			)
}