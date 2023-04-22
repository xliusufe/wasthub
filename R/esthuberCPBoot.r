esthubcpBoot <- function(data, method = "adaptive", smooth = 'sigmoid', weights = 'exponential', h = NULL, maxIter = 100, tol = 0.00001, B = 1000) {

	if(!(smooth %in% c('sigmoid', 'pnorm','mixnorm'))){
		stop("Smooth function must be one of {'sigmoid', 'pnorm','mixnorm'} !")
	}
	if(!(weights %in% c('exponential', 'norm','bernoulli'))){
		stop("Weights must be one of {'exponential', 'norm','bernoulli'} !")
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
	hsigma2 = sum( (y -w1%*%tbeta)^2 )/(n-(p1+p2+p3))


	#--------- bootstrap ---------------------------------------------
	halphaB	= matrix(0, nrow = p1, ncol = B)
	hbetaB	= matrix(0, nrow = p2, ncol = B)
	hthetaB	= matrix(0, nrow = p3, ncol = B)

	dims 	= c(n, p1, p2, p3-1, maxIter, type)
	params 	= c(tol, h, qq, tau0)
	for(b in 1:B){
		G = switch(weights,
				'exponential'	= rexp(n),
				'norm'			= 1 + rnorm(n),
				'bernoulli'		= 2*rbinom(n,1,prob=0.5)
			)

		fit = .Call("_EST_HuberCP_Weight",
				as.numeric(y),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(G),
				as.numeric(theta0),
				as.integer(dims),
				as.numeric(params)
			)

		tbeta 	= fit$beta
		halphaB[, b]	= tbeta[1:p1]
		hbetaB[, b]		= tbeta[-c(1:p1)]
		hthetaB[, b]	= c(1,fit$theta)
	}
	hstd = sqrt(n)*apply(rbind(halphaB, hbetaB, hthetaB), 1, sd)

	return(list(alpha	= halpha,
				beta 	= hbeta,
				theta 	= htheta,
				delta 	= hdelta,
				tau 	= tau,
				sigma2 	= hsigma2,
				std 	= hstd,
				alphaB	= halphaB,
				betaB	= hbetaB,
				thetaB	= hthetaB
				)
			)
}