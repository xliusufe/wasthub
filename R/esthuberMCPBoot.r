esthubmcpBoot <- function(data, ng = 2, method = "adaptive", smooth = 'sigmoid', weights = 'exponential', h = NULL, maxIter = 100, tol = 0.00001, B = 1000) {

	if(!(smooth %in% c('sigmoid', 'pnorm','mixnorm'))){
		stop("Smooth function must be one of {'sigmoid', 'pnorm','mixnorm'} !")
	}
	if(ng<1){
		warning("The number of change-planes must be equal or greater than one.  The default ng = 2 is used!")
		ng = 2
	}
	else if(ng==1){
		warning("The number of multiple change-planes is equal or greater than two. If ng = 1, the single change-plane model is used, which is same as function 'estglmBoot'!")
		fit = esthubcpBoot(data, method = method, smooth = smooth, weights = weights, h = h, maxIter = maxIter, tol = tol, B = B)
		return(fit)
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


	abeta0 	= c(rep(0, p1+ng*p2), rep(1,p3-1),c(-0.5, 0.5))
	dims 	= c(n, p1, p2, p3-1, maxIter, type, ng)
	params 	= c(tol, h, qq, tau0)
	fit 	= .Call("_EST_HUBER_MULTI",
				as.numeric(abeta0),
				as.numeric(y),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.integer(dims),
				as.numeric(params)
			)

	htheta 	= c(1,fit$theta)
	tbeta 	= fit$beta
	halpha 	= tbeta[1:p1]
	ha 		= fit$ha
	tau 	= fit$tau
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
	hsigma2 = sum( (y -tx%*%halpha - xb)^2 )/(n-(p1+p2))

	#--------- bootstrap ---------------------------------------------
	halphaB	= matrix(0, nrow = p1, ncol = B)
	hbetaB	= matrix(0, nrow = p2*ng, ncol = B)
	hthetaB	= matrix(0, nrow = p3, ncol = B)
	haB		= matrix(0, nrow = ng, ncol = B)

	dims 	= c(n, p1, p2, p3-1, maxIter, type, ng)
	params 	= c(tol, h, qq, tau0)
	for(b in 1:B){
		G = switch(weights,
				'exponential'	= rexp(n),
				'norm'			= 1 + rnorm(n),
				'bernoulli'		= 2*rbinom(n,1,prob=0.5)
			)

		fit = .Call("_EST_HUBER_MULTI_WEIGHT",
				as.numeric(abeta0),
				as.numeric(y),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(G),
				as.integer(dims),
				as.numeric(params)
			)

		tbeta 	= fit$beta
		halphaB[, b]	= tbeta[1:p1]
		hbetaB[, b]		= tbeta[-c(1:p1)]
		hthetaB[, b]	= c(1,fit$theta)
		haB[, b]		= fit$ha
	}
	hstd = sqrt(n)*apply(rbind(halphaB, hbetaB, hthetaB, haB), 1, sd)

	return(list(alpha	= halpha,
				beta 	= hbeta,
				theta 	= htheta,
				delta 	= hdelta,
				tau 	= tau,
				ha		= ha,
				sigma2 	= hsigma2,
				std 	= hstd,
				alphaB	= halphaB,
				betaB	= hbetaB,
				thetaB	= hthetaB,
				haB		= haB
				)
			)
}