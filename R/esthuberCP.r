esthubcp <- function(data, method = "adaptive", smooth = 'sigmoid', h = NULL, maxIter = 100, tol = 0.00001, M = 100, tauk = NULL) {

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

	tau0 	= 1.345
	if(method == "adaptive"){
		qq 	= -1
	}
	else if(method == "CV"){
		if(is.null(tauk)){
			tauk 	= seq(0.01,10,length.out = M)
		}
		y1	= y - mean(y)

		theta0 	= rep(1, p3)
		dims 	= c(n, p1, p2, p3-1, maxIter, type)
		params 	= c(tol, h)
		fit	= .Call("_EST_LinearCP",
					as.numeric(y1),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(theta0),
					as.integer(dims),
					as.numeric(params)
				)

		htheta 	= c(1,fit$theta)
		zt 		= (z%*%htheta)*h
		if(smooth == 'sigmoid'){
			sh 	= 1.0/(1.0+exp(-zt))
		}
		else if(smooth == 'pnorm'){
			sh 	= pnorm(zt)
		}
		else{
			sh 	= pnorm(zt) + zt*dnorm(zt)
		}
		w1 		= cbind(tx, x*as.numeric(sh))

		resids 	= y1 - w1%*%fit$beta
		hsigma2 = sqrt( sum(resids^2)/(n-(p1+p2+p3)) )

		beta0 	= fit$beta
		theta0 	= htheta[-1]
		dims 	= c(n, p1, p2, p3-1, maxIter, type)
		loss 	= rep(1, length(tauk))
		for(k in 1:length(tauk)){
			qq 	= hsigma2*tauk[k]*sqrt(n/log(n))
			params 	= c(tol, h, qq, tau0)
			fit <- .Call("_EST_HuberCP_2step",
						as.numeric(y),
						as.numeric(tx),
						as.numeric(x),
						as.numeric(z),
						as.numeric(beta0),
						as.numeric(theta0),
						as.integer(dims),
						as.numeric(params)
					)
			loss[k] = sum( (y1 - w1%*%fit$beta)^2 )
		}
		ko 	= which.min(loss)
		qq 	= hsigma2*tauk[ko]*sqrt(n)/log(n)
	}
	else{
		qq 	= tau0*median( abs(y - median(y)) )/qnorm(0.75)
	}

	theta0 	= rep(1, p3)
	dims 	= c(n, p1, p2, p3-1, maxIter, type)
	params 	= c(tol, h, qq, tau0)
	fit 	= .Call("_EST_HuberCP",
				as.numeric(y),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(theta0),
				as.integer(dims),
				as.numeric(params)
			)

	htheta 	= c(1,fit$theta)
	tbeta 	= fit$beta
	halpha 	= tbeta[1:p1]
	hbeta 	= tbeta[-c(1:p1)]
	hdelta 	= z%*%htheta>0
	tau 	= fit$tau
	meany 	= mean(y)
	hsigma2 = sum( (y -meany -tx%*%halpha - (x%*%hbeta)*hdelta)^2 )/(n-(p1+p2+p3))

	return(list(alpha	= halpha,
				beta 	= hbeta,
				theta 	= htheta,
				delta 	= hdelta,
				tau 	= tau,
				sigma2 	= hsigma2
				)
			)
}