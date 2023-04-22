estolsmcp <- function(data, ng = 2, smooth = 'sigmoid', isBoot = FALSE, isWB = FALSE, h = NULL, maxIter = 100, tol = 0.00001, B = 1000) {

	if(!(smooth %in% c('sigmoid', 'pnorm','mixnorm'))){
		stop("Smooth function must be one of {'sigmoid', 'pnorm','mixnorm'} !")
	}
	if(ng<1){
		warning("The number of change-planes must be equal or greater than one.  The default ng = 2 is used!")
		ng = 2
	}
	else if(ng==1){
		warning("The number of multiple change-planes is equal or greater than two. If ng = 1, the single change-plane model is used, which is same as function 'estglmBoot'!")
		fit = estolscp(data, smooth = smooth, isBoot = isBoot, isWB = isWB, h = h, maxIter = maxIter, tol = tol, B = B)
		return(fit)
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

	abeta0 	= c(rep(0, p1+ng*p2), rep(1,p3-1),c(-0.5, 0.5))
	dims 	= c(n, p1, p2, p3-1, maxIter, type, ng)
	params 	= c(tol, h)
	fit	= .Call("_EST_LINEAR_MULTI",
				as.numeric(abeta0),
				as.numeric(y),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.integer(dims),
				as.numeric(params)
			)

	htheta 	= c(1, fit$theta)
	tbeta 	= fit$beta
	halpha 	= tbeta[1:p1]
	ha 		= fit$ha
	hbeta 	= rep(0, p2*ng)
	xb 		= rep(0, n)
	hdelta 	= matrix(0, nrow = n,  ncol = ng)
	for(k in 1:ng){
		zt 		= (z%*%htheta-ha[k])*h
		if(smooth == 'sigmoid'){
			sh 	= 1.0/(1.0+exp(-zt))
		}
		else if(smooth == 'pnorm'){
			sh 	= pnorm(zt)
		}
		else{
			sh 	= pnorm(zt) + zt*dnorm(zt)
		}
		hbeta[((k-1)*p2+1):(k*p2)] = tbeta[(p1+(k-1)*p2+1):(p1+k*p2)]
		hdelta[,k] 	= z%*%htheta>ha[k]
		xb 			= xb + (x%*%hbeta[((k-1)*p2+1):(k*p2)])*sh
	}
	muhat 	= tx%*%halpha + xb
	resids 	= y - muhat
	hsigma2 = sum( (y-muhat)^2 )/(n-(p1+p2))


	if(isBoot){
		alphaB = rep(0, B)
		for(k in 1:B){
			if(isWB){
				yb 	= muhat + resids*rnorm(n)
			}
			else{
				yb 	= muhat + resids[sample.int(n, n, replace = TRUE)]
			}
			fit = .Call("_EST_LINEAR_MULTI",
						as.numeric(abeta0),
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
				ha		= ha,
				sigma2 	= hsigma2
				)
			)
}