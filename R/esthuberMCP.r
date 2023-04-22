esthubmcp <- function(data, ng = 2, method = "adaptive", smooth = 'sigmoid', h = NULL, maxIter = 100, tol = 0.00001, M = 100, tauk = NULL) {

	if(!(smooth %in% c('sigmoid', 'pnorm','mixnorm'))){
		stop("Smooth function must be one of {'sigmoid', 'pnorm','mixnorm'} !")
	}
	if(ng<1){
		warning("The number of change-planes must be equal or greater than one.  The default ng = 2 is used!")
		ng = 2
	}
	else if(ng==1){
		warning("The number of multiple change-planes is equal or greater than two. If ng = 1, the single change-plane model is used, which is same as function 'estglmBoot'!")
		fit = esthubcp(data, method = method, smooth = smooth, h = h, maxIter = maxIter, tol = tol, M = M, tauk = tauk)
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

	tau0 	= 1.345
	if(method == "adaptive"){
		qq 	= -1

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
	}
	else if(method == "CV"){
		if(is.null(tauk)){
			tauk 	= seq(0.01,10,length.out = M)
		}
		y1 		= y - mean(y)

		abeta0 	= c(rep(0, p1+ng*p2), rep(1,p3-1),c(-0.5, 0.5))
		dims 	= c(n, p1, p2, p3-1, maxIter, type, ng)
		params 	= c(tol, h)
		fit	= .Call("_EST_LINEAR_MULTI",
					as.numeric(abeta0),
					as.numeric(y1),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.integer(dims),
					as.numeric(params)
				)

		tbeta 	= fit$beta
		halpha 	= tbeta[1:p1]
		ha 		= fit$ha
		htheta 	= c(1,fit$theta)
		xb 		= rep(0, n)
		w1 		= tx
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
			hbeta0 = tbeta[(p1+(k-1)*p2+1):(p1+k*p2)]
			xb 			= xb + (x%*%hbeta0)*sh
			w1 		= cbind(w1, x*as.numeric(sh))
		}
		hsigma2 = sum( (y1 -tx%*%halpha - xb)^2 )/(n-(p1+p2))

		abeta0 	= c(tbeta, htheta[-1], ha)
		loss 	= rep(1, length(tauk))
		for(k in 1:length(tauk)){
			qq 	= hsigma2*tauk[k]*sqrt(n/log(n))
			params 	= c(tol, h, qq, tau0)
			fit <- .Call("_EST_HUBER_MULTI_2step",
						as.numeric(abeta0),
						as.numeric(y),
						as.numeric(tx),
						as.numeric(x),
						as.numeric(z),
						as.integer(dims),
						as.numeric(params)
					)
			loss[k] = sum( (y1 - w1%*%fit$beta)^2 )
		}
		ko 	= which.min(loss)
		qq 	= hsigma2*tauk[ko]*sqrt(n)/log(n)

		fit <- .Call("_EST_HUBER_MULTI_2step",
					as.numeric(abeta0),
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.integer(dims),
					as.numeric(params)
				)
	}
	else{
		qq 	= tau0*median( abs(y - median(y)) )/qnorm(0.75)

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
	}


	tau 	= fit$tau
	htheta 	= c(1,fit$theta)
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
	hsigma2 = sum( (y -tx%*%halpha - xb)^2 )/(n-(p1+p2))

	return(list(alpha	= halpha,
				beta 	= hbeta,
				theta 	= htheta,
				delta 	= hdelta,
				tau 	= tau,
				ha		= ha,
				sigma2 	= hsigma2
				)
			)
}