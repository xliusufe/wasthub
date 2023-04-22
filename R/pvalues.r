gam.init = function(n.initials,q,Z,lb.quantile,ub.quantile,ss=1){
	if(q==1){
		gamma.initials = matrix(1,n.initials,q+1)
		gamma.initials[,1] = -quantile(Z,seq(lb.quantile,ub.quantile,length=n.initials))
	}else{
		n.initials = n.initials/ss
		out = matrix(rnorm(n.initials*q), nrow = n.initials)
		gamma.initials	= out/sqrt(apply(out^2,1,sum))
		Z.gamma.initials = Z %*% t(gamma.initials)
		ll=round(n.initials/ss)
		qtile = sample(seq(lb.quantile,ub.quantile,length=n.initials),n.initials)
		gamma.initials.1 = sapply(1:n.initials,function(x)return(
							-quantile(Z.gamma.initials[,x-floor((x-0.1)/ll)*ll],qtile[x])
						))

		gamma.initials.1=ifelse(gamma.initials.1==(-1)*apply(Z.gamma.initials,2,min),gamma.initials.1-0.001,gamma.initials.1)
		gamma.initials.1=ifelse(gamma.initials.1==(-1)*apply(Z.gamma.initials,2,max),gamma.initials.1+0.001,gamma.initials.1)
		gamma.initials.aug=do.call("rbind", rep(list(gamma.initials), ss))
		gamma.initials = cbind(gamma.initials.1,gamma.initials.aug)
	}
	return(gamma.initials)
}

EstTn_Huber_sst <- function(data, isAda = TRUE, K = 1000, M = 1000) {
	y 		= data$Y
	n 		= length(y)
	tx 		= data$X
	x 		= data$Z
	z 		= data$U
	p1 		= ifelse(is.null(ncol(tx)), 1, ncol(tx))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))
	maxIter = 50
	tol 	= 0.00001
	tau0 	= 1.345
	if(isAda){
		tau1 	= -1
	}
	else{
		tau1 	= tau0*median( abs(y - median(y)) )/qnorm(0.75)
	}

	fit <- .Call("_EST_HUBER_2step",
				as.numeric(y),
				as.numeric(tx),
				as.integer(c(n,p1,maxIter)),
				as.numeric(c(tau1, tau0, tol))
			)

	tau 	= fit$tau
	alphahat = fit$coef
	muhat 	= tx%*%alphahat
	resids 	= y - muhat

	if(p3==1){
		z1 		= quantile(z, probs = c(0.1, 0.9))
		rtheta 	= z[(z>z1[1]) & (z<z1[2])]
		K 		= length(rtheta)
	}
	else{
		rtheta = gam.init(K, p3-1, z[,-1], lb.quantile=.1, ub.quantile=0.9, ss=1)
		rtheta = t(rtheta)
	}
	# G  		= matrix(rnorm(M*n), nrow = n, ncol = M)
	G = matrix(0, n, M)
	for(k in 1:M){
		G[,k] 	= resids[sample.int(n, n, replace = TRUE)]
	}
	dims 	= c(n, p1, p2, p3, K, M, maxIter)
	params 	= c(tau, tol)
	fit 	<- .Call("_HUBER_SST",
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(resids),
					as.numeric(rtheta),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params)
				)
	pvals 	= fit$pvals
	theta 	= fit$theta
	return(pvals)
}

EstTn_Huber_wast <- function(data, isAda = TRUE, isWB = FALSE, isBeta = 0, shape1 = 1, shape2 = 1, K = 1000, M = 1000) {
	y 		= data$Y
	n 		= length(y)
	tx 		= data$X
	x 		= data$Z
	z 		= data$U
	p1 		= ifelse(is.null(ncol(tx)), 1, ncol(tx))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))

	maxIter = 100
	tol 	= 0.0001
	tau0 	= 1.345
	if(isAda){
		tau1 	= -1
	}
	else{
		tau1 	= tau0*median( abs(y - median(y)) )/qnorm(0.75)
	}
	dims 	= c(n, p1, p2, p3, M, isBeta, maxIter)


	fit <- .Call("_EST_HUBER_2step",
				as.numeric(y),
				as.numeric(tx),
				as.integer(c(n,p1,maxIter)),
				as.numeric(c(tau1, tau0, tol))
			)

	tau 	= fit$tau
	alphahat = fit$coef
	muhat 	= tx%*%alphahat
	resids 	= y - muhat

	params 	= c(tau, shape1, shape2, tol)
	yb = matrix(0, n, M)
	if(isWB){
		for(k in 1:M){
			yb[,k] 	= muhat + resids*rnorm(n);
		}
	}
	else{
		for(k in 1:M){
			yb[,k] 	= muhat + resids[sample.int(n, n, replace = TRUE)]
		}
	}
	yb 	<- cbind(y, yb)
	fit <- .Call("_HUBER_WAST",
				as.numeric(yb),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(resids),
				as.integer(dims),
				as.numeric(params)
			)

	teststat	= fit$Tn0
	teststat_p	= fit$Tn
	pvals   	= mean(teststat_p >= teststat)

	return(pvals)
}

EstTn_Huber_wastCV <- function(data, isWB = FALSE, isBeta = 0, shape1 = 1, shape2 = 1, K = 1000, M = 1000) {
	y1 		= data$Y
	n 		= length(y1)
	tx 		= data$X
	x 		= data$Z
	z 		= data$U
	p1 		= ifelse(is.null(ncol(tx)), 1, ncol(tx))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))

	maxIter = 100
	tol 	= 0.0001
	tauk 	= seq(0.01,10,length.out = K)
	tau0 	= 1.345
	y 		= y1 - mean(y1)

	fit 	= lm(y~tx-1)
	halpha 	= fit$coefficients
	hsigma2 = sqrt(sum( (y-tx%*%halpha)^2 )/(n-p1))

	loss 	= rep(0, length(tauk))
	for(k in 1:length(tauk)){
		qq 	= hsigma2*tauk[k]*sqrt(n/log(n))
		fit <- .Call("_EST_HUBER_2step",
					as.numeric(y1),
					as.numeric(tx),
					as.integer(c(n,p1,maxIter)),
					as.numeric(c(qq, tau0, tol))
				)
		loss[k] = sum( (y - tx%*%fit$coef)^2 )
	}
	ko 	= which.min(loss)
	qq 	= hsigma2*tauk[ko]*sqrt(n/log(n))

	fit <- .Call("_EST_HUBER_2step",
				as.numeric(y1),
				as.numeric(tx),
				as.integer(c(n,p1,maxIter)),
				as.numeric(c(qq, tau0, tol))
			)

	tau 	= qq
	muhat 	= tx%*%fit$coef
	resids 	= y1 - muhat

	dims 	= c(n, p1, p2, p3, M, isBeta, maxIter)
	params 	= c(tau, shape1, shape2, tol)
	yb = matrix(0, n, M)
	if(isWB){
		for(k in 1:M){
			yb[,k] 	= muhat + resids*rnorm(n);
		}
	}
	else{
		for(k in 1:M){
			yb[,k] 	= muhat + resids[sample.int(n, n, replace = TRUE)]
		}
	}
	yb 	<- cbind(y, yb)
	fit <- .Call("_HUBER_WAST",
				as.numeric(yb),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(resids),
				as.integer(dims),
				as.numeric(params)
			)

	teststat	= fit$Tn0
	teststat_p	= fit$Tn
	pvals   	= mean(teststat_p >= teststat)

	return(pvals)
}

EstTn_Huber_slr <- function(data, isAda = TRUE, K = 1000, M = 1000) {
	y 		= data$Y
	n 		= length(y)
	tx 		= data$X
	x 		= data$Z
	z 		= data$U
	p1 		= ifelse(is.null(ncol(tx)), 1, ncol(tx))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))
	maxIter = 50
	tol 	= 0.00001

	tau0 	= 1.345
	if(isAda){
		tau1 	= -1
	}
	else{
		tau1 	= tau0*median( abs(y - median(y)) )/qnorm(0.75)
	}

	fit <- .Call("_EST_HUBER_2step",
				as.numeric(y),
				as.numeric(tx),
				as.integer(c(n,p1,maxIter)),
				as.numeric(c(tau1, tau0, tol))
			)

	tau 	= fit$tau
	muhat 	= tx%*%fit$coef
	resids 	= y - muhat

	if(p3==1){
		z1 		= quantile(z, probs = c(0.1, 0.9))
		rtheta 	= z[(z>z1[1]) & (z<z1[2])]
		K 		= length(rtheta)
	}
	else{
		rtheta = gam.init(K, p3-1, z[,-1], lb.quantile=.1, ub.quantile=0.9, ss=1)
		rtheta = t(rtheta)
	}

	G = matrix(0, n, M)
	for(k in 1:M){
		G[,k] 	= resids[sample.int(n, n, replace = TRUE)]
	}

	dims 	= c(n, p1, p2, p3, K, M, maxIter)
	params 	= c(tau, tol)
	fit 	<- .Call("_HUBER_SLR",
					as.numeric(y),
					as.numeric(tx),
					as.numeric(x),
					as.numeric(z),
					as.numeric(resids),
					as.numeric(rtheta),
					as.numeric(G),
					as.integer(dims),
					as.numeric(params)
				)
	pvals 	= fit$pvals[1]
	return(pvals)
}

EstTn_Huber_approx <- function(data, isAda = TRUE, isWB = FALSE, isBeta = 0, shape1 = 1, shape2 = 1, M = 1000, N0 = 5000, MU0 = NULL, Z_K = NULL) {
	y 		= data$Y
	n 		= length(y)
	tx 		= data$X
	x 		= data$Z
	z 		= data$U
	p1 		= ifelse(is.null(ncol(tx)), 1, ncol(tx))
	p2 		= ifelse(is.null(ncol(x)) , 1, ncol(x))
	p3 		= ifelse(is.null(ncol(z)) , 1, ncol(z))

	maxIter = 100
	tol 	= 0.0001
	tau0 	= 1.345
	if(isAda){
		tau1 	= -1
	}
	else{
		tau1 	= tau0*median( abs(y - median(y)) )/qnorm(0.75)
	}
	dims 	= c(n, p1, p2, p3, M, isBeta, maxIter)


	fit <- .Call("_EST_HUBER_2step",
				as.numeric(y),
				as.numeric(tx),
				as.integer(c(n,p1,maxIter)),
				as.numeric(c(tau1, tau0, tol))
			)

	tau 	= fit$tau
	alphahat = fit$coef
	muhat 	= tx%*%alphahat
	resids 	= y - muhat

	params 	= c(tau, shape1, shape2, tol)
	yb = matrix(0, n, M)
	if(isWB){
		for(k in 1:M){
			yb[,k] 	= muhat + resids*rnorm(n);
		}
	}
	else{
		for(k in 1:M){
			yb[,k] 	= muhat + resids[sample.int(n, n, replace = TRUE)]
		}
	}
	yb 	<- cbind(y, yb)

	if(is.null(MU0)){
		MU0 	= runif(p3) - 0.5
	}
	if(is.null(Z_K)){
		Z_K		= rnorm(N0)
	}

	fit <- .Call("_HUBER_WAST_APPROX",
				as.numeric(yb),
				as.numeric(tx),
				as.numeric(x),
				as.numeric(z),
				as.numeric(resids),
				as.numeric(Z_K),
				as.numeric(MU0),
				as.integer(dims),
				as.numeric(params)
			)

	teststat	= fit$Tn0
	teststat_p	= fit$Tn
	pvals   	= mean(teststat_p >= teststat)

	return(pvals)
}

pvalhuber <- function(data, method = "wast", isAda = TRUE, isWB = FALSE, B = 1000, K = 1000, isBeta = FALSE, shape1 = 1, shape2 = 1, N0 = 5000, MU = NULL, ZK = NULL){
	isBeta = ifelse(isBeta, 1, 0)

	if(method=='wast') {
	   pvals  	= EstTn_Huber_wast(data, isAda = isAda, isWB = isWB, isBeta = isBeta, shape1 = shape1, shape2 = shape2, K = K, M = B)
	}
	else if(method=='sst'){
		pvals  	= EstTn_Huber_sst(data, isAda = isAda, K = K, M = B)
	}
	else if(method=='slrt'){
		pvals  	= EstTn_Huber_slr(data, isAda = isAda, K = K, M = B)
	}
	else if(method=="wastcv"){
		K = 100
		pvals 	= EstTn_Huber_wastCV(data, isWB = isWB, isBeta = isBeta, shape1 = shape1, shape2 = shape2, K = K, M = B)
	}
	else if(method=='wastapprox'){
		pvals  	= EstTn_Huber_approx(data, isAda = isAda, isWB = isWB, isBeta = isBeta, shape1 = shape1, shape2 = shape2, M = B, N0 = N0, MU0 = MU, Z_K = ZK)
	}
	else{
		warning("Input method is not one of {'wast', 'wastcv', 'wastapprox', 'sst', and 'slrt'}. The default method 'wast' is used!")
		pvals  	= EstTn_Huber_wast(data, isAda = isAda, isWB = isWB, isBeta = isBeta, shape1 = shape1, shape2 = shape2, K = K, M = B)
	}
	return(pvals)
}