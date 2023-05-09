#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "_WASTHUB_HEAD.h"

//-------------- Estimate robust regression --------------------------------
void Huber_d1(double *x, double *g, double tau, int n){
	int i;
	for(i=0;i<n;i++){
		if(fabs(x[i])>tau){
			g[i] = tau*SGN(x[i]);
		}
		else{
			g[i] = x[i];
		}
	}
}

double Huber_d1sum(double *x, double tau, int n){
	int i;
	double sum = 0.0;
	for(i=0;i<n;i++){
		if(fabs(x[i])>tau){
			sum += tau*SGN(x[i]);
		}
		else{
			sum += x[i];
		}
	}
	return sum/n;
}

double Huber_scale(double *resid2, double low, double up, int n, double rate, int maxstep, double eps){
	int i, step=0;
	double tmp, mid;

	while(step < maxstep){
		step++;

		mid = 0.5*(low + up);
		tmp = 0.0;
		for(i=0;i<n;i++){
			if(resid2[i] < mid){
				tmp += resid2[i]/mid;
			}
			else{
				tmp += 1.0;
			}
		}
		tmp = tmp/n - rate;
		if(tmp < 0){
			up = mid;
		}
		else{
			low = mid;
		}
		if(up-low<eps){
			break;
		}

	}
	return 0.5 * (low + up);

}

double Huber_mean(double *x, double *tau0, int n, int maxstep, double eps){
	int i, step=0;
	double meanx, stdx, tau, psiSum, psiSum0, psiSumDiff, mu, rate, low, up;
	double alpha, muDiff, a1, a2, cross;
	double *x1, *resid, *resid2;

	x1	 	= (double*)malloc(sizeof(double)*n);
	resid	= (double*)malloc(sizeof(double)*n);
	resid2	= (double*)malloc(sizeof(double)*n);

	rate = log(1.0*n)/n;
	meanx = 0.0;
	for(i=0;i<n;i++){
		meanx += x[i];
	}
	meanx /= n;

	stdx = 0.0;
	for(i=0;i<n;i++){
		x1[i] = x[i] - meanx;
		stdx += x1[i]*x1[i];
	}
	stdx = sqrt(stdx/(n-1));
	tau = stdx*sqrt(n/log(1.0*n));


	psiSum0 = Huber_d1sum(x1, tau, n);
	mu = psiSum0;
	muDiff = psiSum0;

	for(i=0;i<n;i++){
		resid[i] = x1[i] - mu;
		resid2[i] = resid[i]*resid[i];
	}
	low = minx(resid2, n);
	up 	= cumsum(resid2, n);
	tau = sqrt( Huber_scale(resid2, low, up, n, rate, maxstep, eps) );

	psiSum = Huber_d1sum(resid, tau, n);
	psiSumDiff = psiSum0 - psiSum;


	while(fabs(psiSum)> eps && step < maxstep){
		step++;

		alpha = 1.0;
		cross = muDiff*psiSumDiff;
		if(cross>0){
			a1 = cross/(psiSumDiff * psiSumDiff);
			a2 = muDiff*muDiff/cross;
			alpha = MIN( MIN(a1, a2), 100.0);
		}
		psiSum0 = psiSum;
		muDiff = alpha*psiSum;
		mu += muDiff;
		for(i=0;i<n;i++){
			resid[i] = x1[i] - mu;
			resid2[i] = resid[i]*resid[i];
		}


		low = minx(resid2, n);
		up 	= cumsum(resid2, n);
		tau = sqrt( Huber_scale(resid2, low, up, n, rate, maxstep, eps) );
		psiSum = Huber_d1sum(resid, tau, n);
		psiSumDiff = psiSum0 - psiSum;
	}

	*tau0 = tau;
	free(x1);
	free(resid);
	free(resid2);

	return mu + meanx;

}

double EstHuberR_2Step(double *x, double *y1, double *beta, double *residual, double qq, double tau0, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta0, *hess, *qy, *dpsi, *y, *med;
	double tmp, bnorm, yk, wk, tau, meany, qr;

	y 		= (double*)malloc(sizeof(double)*n);
	med 	= (double*)malloc(sizeof(double)*n);
	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	meany = 0.0;
	for(i=0;i<n;i++){
		meany += y1[i];
	}
	meany /= n;
	for(i=0;i<n;i++){
		y[i] = y1[i] - meany;
	}

	for(j=0;j<p;j++) qy[j] = 0.0;
	for(i=0;i<n;i++){
		yk = y[i];
		for(j=0;j<p;j++){
			qy[j] 	+= x[j*n+i]*yk;
		}
	}
	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}
	MatrixInvSymmetric(hess,p);
	AbyB(beta0, hess, qy, p, p, 1);
	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++) tmp += x[j*n+i]*beta0[j];
		residual[i]	= y[i] - tmp;
	}

	if(qq < 0){
		qr = SampleQuantile1(residual, n, 0.5);
		for(i=0;i<n;i++){
			med[i] = fabs(residual[i] - qr);
		}
		qr = SampleQuantile1(med, n, 0.5);
		qq = tau0 * qr;
	}

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			if(fabs(residual[i])>qq){
				wk = qq/fabs(residual[i]);
			}
			else{
				wk = 1.0;
			}
			yk = y[i];
			for(j=0;j<p;j++){
				tmp 	= x[j*n+i]*wk;
				qy[j] 	+= tmp*yk;
				dpsi[j*n+i] = tmp;
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp = beta[j] - beta0[j];
			bnorm += tmp*tmp;
		}

		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)	beta0[j] = beta[j];
			for(i=0;i<n;i++){
				tmp = 0.0;
				for(j=0;j<p;j++) tmp += x[j*n+i]*beta[j];
				residual[i] = y[i] - tmp;
			}
		}
	}
	beta[0] += meany;

	for(i=0;i<n;i++){
		residual[i] += beta[0];
	}
	beta[0] = Huber_mean(residual, &tau, n, maxstep, eps);

	qq = tau;
	for(i=0;i<n;i++){
		tmp = residual[i] - beta[0];
		if(fabs(tmp)<qq){
			residual[i] = tmp;
		}
		else{
			residual[i] = qq*SGN(tmp);
		}
	}

	free(y);
	free(med);
	free(beta0);
	free(hess);
	free(dpsi);
	free(qy);
	return tau;
}

SEXP _EST_HUBER(SEXP Y_, SEXP X_, SEXP PARA_INT_, SEXP PARA_DOUBLE_){
	// dimensions
	int *para 	= INTEGER(PARA_INT_);
	int n     	= para[0];
	int p     	= para[1];
	int maxstep	= para[2];

	double *para1 	= REAL(PARA_DOUBLE_);
	double tau   	= para1[0];
	double tau0 	= para1[1];
	double eps    	= para1[2];

	// Outcome
	SEXP rBeta, rTau, rResids, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p));
	PROTECT(rResids		= allocVector(REALSXP, 	n));
	PROTECT(rTau 		= allocVector(REALSXP, 	1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	if(p==1)
		REAL(rBeta)[0]	= EstHuberR1(REAL(Y_), REAL(rResids), tau, n);
	else
		REAL(rTau)[0] = EstHuberR(REAL(X_), REAL(Y_), REAL(rBeta), REAL(rResids), tau, tau0, n, p, maxstep, eps);

	SET_STRING_ELT(list_names, 	0,	mkChar("coef"));
	SET_STRING_ELT(list_names, 	1,  mkChar("residuals"));
	SET_STRING_ELT(list_names, 	2,	mkChar("tau"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rResids);
	SET_VECTOR_ELT(list, 		2, 	rTau);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_HUBER_2step(SEXP Y_, SEXP X_, SEXP PARA_INT_, SEXP PARA_DOUBLE_){
	// dimensions
	int *para 	= INTEGER(PARA_INT_);
	int n     	= para[0];
	int p     	= para[1];
	int maxstep	= para[2];

	double *para1 	= REAL(PARA_DOUBLE_);
	double tau   	= para1[0];
	double tau0 	= para1[1];
	double eps    	= para1[2];

	// Outcome
	SEXP rBeta, rTau, rResids, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p));
	PROTECT(rResids		= allocVector(REALSXP, 	n));
	PROTECT(rTau 		= allocVector(REALSXP, 	1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	if(p==1)
		REAL(rBeta)[0]	= Huber_mean(REAL(Y_), REAL(rTau), n, maxstep, eps);
	else
		REAL(rTau)[0] = EstHuberR_2Step(REAL(X_), REAL(Y_), REAL(rBeta), REAL(rResids), tau, tau0, n, p, maxstep, eps);

	SET_STRING_ELT(list_names, 	0,	mkChar("coef"));
	SET_STRING_ELT(list_names, 	1,  mkChar("residuals"));
	SET_STRING_ELT(list_names, 	2,	mkChar("tau"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rResids);
	SET_VECTOR_ELT(list, 		2, 	rTau);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_CATONI(SEXP Y_, SEXP X_, SEXP PARA_INT_, SEXP PARA_DOUBLE_){
	// dimensions
	int *para 	= INTEGER(PARA_INT_);
	int n     	= para[0];
	int p     	= para[1];
	int maxstep	= para[2];

	double *para1 	= REAL(PARA_DOUBLE_);
	double tau   	= para1[0];
	double eps    	= para1[1];

	// Outcome
	SEXP rBeta, rResids, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p));
	PROTECT(rResids		= allocVector(REALSXP, 	n));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	if(p==1)
		REAL(rBeta)[0]	= EstCatoniR1(REAL(Y_), REAL(rResids), tau, n);
	else
		EstCatoniR(REAL(X_), REAL(Y_), REAL(rBeta), REAL(rResids), tau, n, p, maxstep, eps);

	SET_STRING_ELT(list_names, 	0,	mkChar("coef"));
	SET_STRING_ELT(list_names, 	1,  mkChar("residuals"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rResids);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- SST of robust regression --------------------------------
double Huber_SST0(double *y, double *tx1, double *x, double *z, double *alphahat, double* thetahat, double *resid1,
		int n, int p11, int p2, int p3, double tau, int maxIter, double tol, int K, int M, double *G, double* theta){
	int i,j,k,s,t, *subg, maxk=0, count=0, sumsb = 0;
	double tmp, tmp1, Tn=0.0, Tn0=0.0, Tn_star, pvals;
	double *score1, *psi, *psin;
	double *C11, *psi_h, *K1, *Vh, *Tns, *weight;

	subg	= (int*)malloc(sizeof(int)*n);
	psin	= (double*)malloc(sizeof(double)*p2);
	psi  	= (double*)malloc(sizeof(double)*n*p2);
	psi_h  	= (double*)malloc(sizeof(double)*n*p2);
	score1 	= (double*)malloc(sizeof(double)*n*p11);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	K1		= (double*)malloc(sizeof(double)*p2*p11);
	Vh		= (double*)malloc(sizeof(double)*p2*p2);
	Tns		= (double*)malloc(sizeof(double)*M);
	weight 	= (double*)malloc(sizeof(double)*n);


	for(j=0; j<M; j++){
		for(i=0; i<n; i++){
			tmp1 = fabs(G[j*n+i]);
			if(tmp1 > tau){
				G[j*n+i] = tau*SGN(G[j*n+i]);
			}
		}
	}
	for(i=0; i<n; i++){
		tmp1 = fabs(resid1[i]);
		if(tmp1 > tau){
			resid1[i] = tau*SGN(resid1[i]);
			weight[i] = tau/fabs(tmp1);
			for(j=0; j<p11; j++){
				score1[j*n+i] = tx1[j*n+i]*weight[i];
			}
		}
		else{
			weight[i] = 1.0;
			for(j=0; j<p11; j++){
				score1[j*n+i] = tx1[j*n+i];
			}
		}
	}

	for (s = 0; s < p11; s++){
		for (t = 0; t < p11; t++){
			tmp = 0.0;
			for (i = 0; i < n; i++){
				tmp += tx1[s*n+i]*score1[t*n+i];
			}
			C11[s*p11+t] 	= tmp/n;
		}
	}

	if(p11==1){
		C11[0] = 1.0/C11[0];
	}
	else{
		MatrixInvSymmetric(C11, p11);
	}

	for(i=0; i<n; i++){
		for(j=0; j<p11; j++){
			tmp = 0.0;
			for (s = 0; s < p11; s++){
				tmp += C11[j*p11+s]*tx1[s*n+i];
			}
			score1[j*n+i] = tmp*resid1[i];
		}

		for(j=0; j<p2; j++){
			psi[j*n+i] 	= x[j*n+i]*resid1[i];
		}
	}

	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}

	for(k=0; k<K; k++){
		sumsb 	= 0;
		if(p3==1){
			for(i=0; i<n; i++){
				subg[i] = IDEX(theta[k], z[i]);
				sumsb += subg[i];
			}
		}
		else{
			for(i=0; i<n; i++){
				tmp = 0.0;
				for (j = 0; j < p3; j++){
					tmp += z[j*n+i]*theta[k*p3 + j];
				}
				subg[i] = IDEX(0.0, tmp);
				sumsb += subg[i];
			}
		}

		if (sumsb==0){
			continue;
		}
		for (s = 0; s < p2; s++){
			for (t = 0; t < p11; t++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					if(subg[i]){
						tmp += x[s*n+i]*tx1[t*n+i]*weight[i];
					}
				}
				K1[s*p11+t] = tmp/n;
			}
		}

		for (s = 0; s < p2; s++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += psi[s*n+i]*subg[i];
				tmp1 = 0.0;
				for (t = 0; t < p11; t++){
					tmp1 += K1[s*p11+t]*score1[t*n+i];
				}
				psi_h[s*n+i] = psi[s*n+i]*subg[i] - tmp1;

			}
			psin[s] = tmp;
		}

		for (s = 0; s < p2; s++){
			for (t = s; t < p2; t++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += psi_h[s*n+i]*psi_h[t*n+i];
				}
				Vh[s*p2+t] = tmp/n;
			}
		}
		for (s = 1; s < p2; s++){
			for (t = 0; t < s; t++){
				Vh[s*p2+t] = Vh[t*p2+s];
			}
		}

		if(p2 ==1){
			Vh[0] = 1.0/Vh[0];
		}
		else{
			MatrixInvSymmetric(Vh, p2);
		}
		Tn = 0.0;
		for (s = 0; s < p2; s++){
			for (t = 0; t < p2; t++){
				Tn += psin[s]*Vh[s*p2+t]*psin[t];
			}
		}

		if(Tn>Tn0){
			Tn0 = Tn;
			maxk = k;
		}
		for(j=0; j< M; j++){
			Tn_star = 0.0;
			for (s = 0; s < p2; s++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += psi_h[s*n+i]*G[j*n+i];
				}
				psin[s] = tmp;
			}


			for (s = 0; s < p2; s++){
				for (t = 0; t < p2; t++){
					Tn_star += psin[s]*Vh[s*p2+t]*psin[t];
				}
			}
			if(Tn_star>Tns[j]){
				Tns[j] = Tn_star;
			}
		}
	}

	for(j=0; j< M; j++){
		if(Tn0<=Tns[j]){
			count++;
		}
	}
	pvals = 1.0*count/M;
	for(j=0; j<p3; j++){
		thetahat[j] = theta[maxk*p3 + j];
	}


	free(subg);
	free(psi);
	free(score1);
	free(C11);
	free(Vh);
	free(psi_h);
	free(psin);
	free(K1);
	free(Tns);
	free(weight);
	return pvals;
}

double Huber_SST(double *y, double *tx1, double *x, double *z, double *alphahat, double* thetahat, double *resid1,
		int n, int p11, int p2, int p3, double tau, int maxIter, double tol, int K, int M, double *G, double* theta){
	int i,j,k,s,t, *subg, maxk=0, count=0, sumsb = 0;
	double tmp, tmp1, Tn=0.0, Tn0=0.0, Tn_star, pvals;
	double *score1, *psi, *psin;
	double *C11, *psi_h, *K1, *Vh, *Tns, *weight;

	subg	= (int*)malloc(sizeof(int)*n);
	psin	= (double*)malloc(sizeof(double)*p2);
	psi  	= (double*)malloc(sizeof(double)*n*p2);
	psi_h  	= (double*)malloc(sizeof(double)*n*p2);
	score1 	= (double*)malloc(sizeof(double)*n*p11);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	K1		= (double*)malloc(sizeof(double)*p2*p11);
	Vh		= (double*)malloc(sizeof(double)*p2*p2);
	Tns		= (double*)malloc(sizeof(double)*M);
	weight 	= (double*)malloc(sizeof(double)*n);


	for(i=0; i<n; i++){
		tmp1 = fabs(resid1[i]);
		if(tmp1 > tau){
			resid1[i] = tau*SGN(resid1[i]);
			weight[i] = 0.0;
		}
		else{
			weight[i] = 1.0;
		}
	}


	for (s = 0; s < p11; s++){
		for (t = 0; t < p11; t++){
			tmp = 0.0;
			for (i = 0; i < n; i++){
				tmp += tx1[s*n+i]*tx1[t*n+i];
			}
			C11[s*p11+t] 	= tmp/n;
		}
	}
	if(p11==1){
		C11[0] = 1.0/C11[0];
	}
	else{
		MatrixInvSymmetric(C11, p11);
	}

	for(i=0; i<n; i++){
		for(j=0; j<p11; j++){
			tmp = 0.0;
			for (s = 0; s < p11; s++){
				tmp += C11[j*p11+s]*tx1[s*n+i];
			}
			score1[j*n+i] = tmp*resid1[i];
		}

		for(j=0; j<p2; j++){
			psi[j*n+i] 	= x[j*n+i]*resid1[i];
		}
	}

	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}

	for(k=0; k<K; k++){
		sumsb 	= 0;
		if(p3==1){
			for(i=0; i<n; i++){
				subg[i] = IDEX(theta[k], z[i]);
				sumsb += subg[i];
			}
		}
		else{
			for(i=0; i<n; i++){
				tmp = 0.0;
				for (j = 0; j < p3; j++){
					tmp += z[j*n+i]*theta[k*p3 + j];
				}
				subg[i] = IDEX(0.0, tmp);
				sumsb += subg[i];
			}
		}

		if (sumsb==0){
			continue;
		}
		for (s = 0; s < p2; s++){
			for (t = 0; t < p11; t++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					if(subg[i] && fabs(resid1[i])<=tau){
						tmp += x[s*n+i]*tx1[t*n+i];
					}
				}
				K1[s*p11+t] = tmp/n;
			}
		}

		for (s = 0; s < p2; s++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += psi[s*n+i]*subg[i];
				tmp1 = 0.0;
				for (t = 0; t < p11; t++){
					tmp1 += K1[s*p11+t]*score1[t*n+i];
				}
				psi_h[s*n+i] = psi[s*n+i]*subg[i] - tmp1;

			}
			psin[s] = tmp;
		}

		for (s = 0; s < p2; s++){
			for (t = s; t < p2; t++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += psi_h[s*n+i]*psi_h[t*n+i];
				}
				Vh[s*p2+t] = tmp/n;
			}
		}
		for (s = 1; s < p2; s++){
			for (t = 0; t < s; t++){
				Vh[s*p2+t] = Vh[t*p2+s];
			}
		}

		if(p2 ==1){
			Vh[0] = 1.0/Vh[0];
		}
		else{
			MatrixInvSymmetric(Vh, p2);
		}
		Tn = 0.0;
		for (s = 0; s < p2; s++){
			for (t = 0; t < p2; t++){
				Tn += psin[s]*Vh[s*p2+t]*psin[t];
			}
		}

		if(Tn>Tn0){
			Tn0 = Tn;
			maxk = k;
		}
		for(j=0; j< M; j++){
			Tn_star = 0.0;
			for (s = 0; s < p2; s++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += psi_h[s*n+i]*G[j*n+i];
				}
				psin[s] = tmp;
			}


			for (s = 0; s < p2; s++){
				for (t = 0; t < p2; t++){
					Tn_star += psin[s]*Vh[s*p2+t]*psin[t];
				}
			}
			if(Tn_star>Tns[j]){
				Tns[j] = Tn_star;
			}
		}
	}

	for(j=0; j< M; j++){
		if(Tn0<=Tns[j]){
			count++;
		}
	}
	pvals = 1.0*count/M;
	for(j=0; j<p3; j++){
		thetahat[j] = theta[maxk*p3 + j];
	}

	free(subg);
	free(psi);
	free(score1);
	free(C11);
	free(Vh);
	free(psi_h);
	free(psin);
	free(K1);
	free(Tns);
	free(weight);
	return pvals;
}

SEXP _HUBER_SST(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP RESIDS, SEXP THETA, SEXP G, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, K, M, maxIter;
	double tol, tau;

	n 		= INTEGER(DIMs)[0];
	p1 		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	K 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];

	tau		= REAL(PARAMs)[0];
	tol 	= REAL(PARAMs)[1];

	SEXP rpvals, rAlpha, rthetahat, list, list_names;
  	PROTECT(rpvals 		= allocVector(REALSXP, 	1));
	PROTECT(rthetahat	= allocVector(REALSXP, 	p3));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));


	REAL(rpvals)[0] = Huber_SST(REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(rAlpha), REAL(rthetahat), REAL(RESIDS),
					n, p1, p2, p3, tau, maxIter, tol, K, M, REAL(G), REAL(THETA));
	SET_STRING_ELT(list_names, 	0,  mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,  mkChar("alpha"));
	SET_STRING_ELT(list_names, 	2,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rAlpha);
	SET_VECTOR_ELT(list, 		2, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

//-------------- SLR of robust regression --------------------------------
void Huber_SLR(double *pvals, double *y, double *tx, double *x, double *z, double *resids, double tau, double *theta, double* G, double *thetahat,
		int n, int p, int p2, int p3, int K, int M, int maxstep, double tol){
	int i,j,k,s,t,p11=p+p2,*subg,sumsb=0,count=0,maxk=0;
	double tmp1, Tn0;
	double *psi0, *psi1, *psi2, *psin, *C11, *alpha0, *Tns, *I0, *In, *weight;

	subg	= (int*)malloc(sizeof(int)*n);
	weight 	= (double*)malloc(sizeof(double)*n);
	alpha0 	= (double*)malloc(sizeof(double)*p11);
	Tns		= (double*)malloc(sizeof(double)*M);
	psi0 	= (double*)malloc(sizeof(double)*n*p11);
	psi1 	= (double*)malloc(sizeof(double)*n*p11);
	psi2 	= (double*)malloc(sizeof(double)*n*p11);
	psin 	= (double*)malloc(sizeof(double)*p11);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	I0		= (double*)malloc(sizeof(double)*p*p);
	In 		= (double*)malloc(sizeof(double)*p11*p11);


	for(j=0; j<M; j++){
		for(i=0; i<n; i++){
			tmp1 = fabs(G[j*n+i]);
			if(tmp1 > tau){
				G[j*n+i] = tau*SGN(G[j*n+i]);
			}
		}
	}
	for(i=0; i<n; i++){
		tmp1 = fabs(resids[i]);
		if(tmp1 > tau){
			resids[i] = tau*SGN(resids[i]);
			weight[i] = tau/fabs(tmp1);
		}
		else{
			weight[i] = 1.0;
		}
	}
	for(i=0; i<n; i++){
		weight[i] = sqrt(weight[i]);
	}

	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp1 = 0.0;
			for(i=0; i<n; i++){
				tmp1 += tx[j*n+i]*tx[k*n+i]*weight[i];
			}
			I0[j*p+k] = tmp1;
		}
	}
	MatrixInvSymmetric(I0,p);

	for(j=0; j<p; j++){
		for(i=0; i<n; i++){
			psi0[j*n+i]	= tx[j*n+i]*resids[i];
			psi1[j*n+i]	= tx[j*n+i]*weight[i];
			psi2[j*n+i]	= tx[j*n+i];
		}
	}

	for(j=0; j< M; j++){
		Tns[j]	= 0.0;
	}
	Tn0 = -1000000.0;
	for(k=0; k<K; k++){
		sumsb 	= 0;
		if(p3==1){
			for(i=0; i<n; i++){
				subg[i] = IDEX(theta[k], z[i]);
				sumsb += subg[i];
			}
		}
		else{
			for(i=0; i<n; i++){
				tmp1 = 0.0;
				for (s = 0; s < p3; s++){
					tmp1 += z[s*n+i]*theta[k*p3 + s];
				}
				subg[i] = IDEX(0.0, tmp1);
				sumsb += subg[i];
			}
		}

		if (sumsb==0){
			continue;
		}
		for(i=0; i<n; i++){
			if(subg[i]){
				for(j=0; j<p2; j++){
					psi0[(j+p)*n+i]	= x[j*n+i]*resids[i];
					psi1[(j+p)*n+i]	= x[j*n+i]*weight[i];
					psi2[(j+p)*n+i]	= x[j*n+i];
				}
			}
			else{
				for(j=0; j<p2; j++){
					psi0[(j+p)*n+i]	= 0.0;
					psi1[(j+p)*n+i]	= 0.0;
					psi2[(j+p)*n+i]	= 0.0;
				}
			}
		}

		for (s = 0; s < p11; s++){
			for (t = s; t < p11; t++){
				tmp1 = 0.0;
				for(i=0; i<n; i++){
					tmp1 += psi1[s*n+i]*psi1[t*n+i];
				}
				In[s*p11+t] = tmp1;
			}
		}
		for (s = 1; s < p11; s++){
			for (t = 0; t < s; t++){
				In[s*p11+t] = In[t*p11+s];
			}
		}

		MatrixInvSymmetric(In, p11);

		for (s = 0; s < p; s++){
			for (t = 0; t < p; t++){
				In[s*p11+t] -= I0[s*p+t];
			}
		}

		for (s = 0; s < p11; s++){
			tmp1 = 0.0;
			for(i=0; i<n; i++){
				tmp1 += psi0[s*n+i];
			}
			psin[s] = tmp1;
		}
		tmp1 = 0.0;
		for (s = 0; s < p11; s++){
			for (t = 0; t < p11; t++){
				tmp1 += psin[s]*In[s*p11+t]*psin[t];
			}
		}


		if(tmp1 > Tn0){
			Tn0 = tmp1;
			maxk = k;
		}

		for(j=0; j< M; j++){
			for (s = 0; s < p11; s++){
				tmp1 = 0.0;
				for(i=0; i<n; i++){
					tmp1 += psi2[s*n+i]*G[j*n+i];
				}
				psin[s] = tmp1;
			}
			tmp1 = 0.0;
			for (s = 0; s < p11; s++){
				for (t = 0; t < p11; t++){
					tmp1 += psin[s]*In[s*p11+t]*psin[t];
				}
			}

			if(tmp1 > Tns[j]){
				Tns[j] = tmp1;
			}
		}
	}

	for(j=0; j< M; j++){
		if(Tn0<=Tns[j]){
			count++;
		}
	}
	pvals[0] = 1.0*count/M;

	for(j=0; j<p3; j++){
		thetahat[j] = theta[maxk*p3 + j];
	}

	free(subg);
	free(alpha0);
	free(weight);
	free(Tns);
	free(psi0);
	free(psi1);
	free(psi2);
	free(psin);
	free(C11);
	free(I0);
	free(In);
}

SEXP _HUBER_SLR(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP RESIDS, SEXP THETA, SEXP G, SEXP DIMs, SEXP PARAMS){
	int n, p1, p2, p3, K, M, maxstep;
	double tol, tau;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	K 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];
	maxstep = INTEGER(DIMs)[6];

	tau 	= REAL(PARAMS)[0];
	tol		= REAL(PARAMS)[1];

	SEXP rpvals, rthetahat, list, list_names;
  	PROTECT(rpvals 		= allocVector(REALSXP, 	2));
	PROTECT(rthetahat	= allocVector(REALSXP, 	p3));
	PROTECT(list 		= allocVector(VECSXP, 	2));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));

	Huber_SLR(REAL(rpvals), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESIDS), tau, REAL(THETA), REAL(G), REAL(rthetahat),
			n, p1, p2, p3, K, M, maxstep, tol);

	SET_STRING_ELT(list_names, 	0,	mkChar("pvals"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rpvals);
	SET_VECTOR_ELT(list, 		1, 	rthetahat);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- WAST of robust regression --------------------------------
void EstHuber2(double *x, double *y, double *beta, double *residual, double qq, int n, int p, int maxstep, double eps){
	int i,j,k, step=0, flag=1;
	double *beta0, *hess, *qy, *dpsi;
	double tmp, bnorm, yk, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	hess	= (double*)malloc(sizeof(double)*p*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	for(j=0;j<p;j++) qy[j] = 0.0;
	for(i=0;i<n;i++){
		yk = y[i];
		for(j=0;j<p;j++){
			qy[j] 	+= x[j*n+i]*yk;
		}
	}
	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}
	MatrixInvSymmetric(hess,p);
	AbyB(beta0, hess, qy, p, p, 1);
	for(i=0;i<n;i++){
		tmp = y[i];
		for(j=0;j<p;j++) tmp -= x[j*n+i]*beta0[j];
		residual[i]	= tmp;
	}

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			if(fabs(residual[i])>qq){
				wk = qq/fabs(residual[i]);
			}
			else{
				wk = 1.0;
			}
			yk = y[i];
			for(j=0;j<p;j++){
				tmp 	= x[j*n+i]*wk;
				qy[j] 	+= tmp*yk;
				dpsi[j*n+i] = tmp;
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = qy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, qy, p, p, 1);
			}
		}

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp = beta[j] - beta0[j];
			bnorm += tmp*tmp;
		}

		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)	beta0[j] = beta[j];
			for(i=0;i<n;i++){
				tmp = y[i];
				for(j=0;j<p;j++) tmp -= x[j*n+i]*beta[j];
				residual[i] = tmp;
			}
		}
	}

	for(i=0;i<n;i++){
		tmp = y[i];
		for(j=0;j<p;j++) tmp -= x[j*n+i]*beta[j];
		if(fabs(tmp)<qq){
			residual[i] = tmp;
		}
		else{
			residual[i] = qq*SGN(tmp);
		}
	}

	free(beta0);
	free(dpsi);
	free(hess);
	free(qy);
}

double Huber_single(double *Tns, double *y, double *tx, double *x, double *z, double *resid0, double *alphahat, int n, int p, int p2,
		int M, double tau, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double tmp, omega=0.0, Tn, Tn0=0.0, xij;
	double *ty, *resid, *yb, *OMEGA;

	yb 		= (double*)malloc(sizeof(double)*n);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	ty  	= (double*)malloc(sizeof(double)*n);

	for(i=0; i<n; i++){
		tmp = fabs(resid0[i]);
		if(tmp>tau){
			resid0[i] = tau*SGN(resid0[i]);
		}
	}
	if(isBeta==1){
		for(i=0; i<n; i++){
			ty[i] = incbeta(z[i], shape1, shape2);
		}
	}
	else if(isBeta==0){
		shape2 *= SQRT2;
		for(i=0; i<n; i++){
			ty[i] = 0.5*erf( (z[i] - shape1)*shape2 ) +0.5;
		}
	}
	else{
		for(i=0; i<n; i++){
			ty[i] = z[i];
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			omega 	= (z[i]<z[j])?ty[i]:ty[j];
			OMEGA[i*n+j]	= omega*xij;
			Tn0 	+= OMEGA[i*n+j]*resid0[i]*resid0[j];
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[n+g*n+i];
		}
		EstHuber2(tx, yb, alphahat, resid, tau, n, p, maxIter, tol);
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(ty);
	return 2.0*Tn0/n;
}

double Huber_multiple(double *Tns, double *y, double *tx, double *x, double *z, double *resid0, double *alphahat, int n, int p, int p2, int p3,
		int M, double tau, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double tmp, omega=0.0, Tn, Tn0=0.0, rho, sd, xij;
	double *resid, *stdx, *yb, *OMEGA;

	yb 		= (double*)malloc(sizeof(double)*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);

	for(i=0; i<n; i++){
		tmp = fabs(resid0[i]);
		if(tmp>tau){
			resid0[i] = tau*SGN(resid0[i]);
		}
	}

	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}

			if(1-rho*rho < MEPS){
				omega = 0.5;
			}
			else{
				omega = 0.25 + atan(rho/sqrt(1-rho*rho))*MPI1;
			}
			OMEGA[i*n+j]	= omega*xij;
			Tn0 	+= OMEGA[i*n+j]*resid0[i]*resid0[j];
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[n+g*n+i];
		}
		EstHuber2(tx, yb, alphahat, resid, tau, n, p, maxIter, tol);
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(stdx);
	return 2.0*Tn0/n;
}

double Huber_multiple_approx(double *Tns, double *y, double *tx, double *x, double *z, double *resid0, double *alphahat, double *zk, double *mu0, int n, int p, int p2, int p3,
		int M, double tau, double shape1, double shape2, int maxIter, double tol, int N0){
	int i,j,s,g,count;
	double tmp, tmp1, tmp2, tmp3, omega=0.0, Tn, Tn0=0.0, rho, sd, xij;
	double *resid, *stdx, *yb, *OMEGA, *zmu;

	yb 		= (double*)malloc(sizeof(double)*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	zmu 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);

	for(i=0; i<n; i++){
		tmp = fabs(resid0[i]);
		if(tmp>tau){
			resid0[i] = tau*SGN(resid0[i]);
		}
	}

	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		tmp = 0.0;
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
			tmp += stdx[j*n+i]*mu0[j];
		}
		zmu[i] = tmp;
	}


	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}

			if(1-rho*rho < MEPS){
				tmp1 = zmu[j];
				count = 0;
				for (s = 0; s < N0; s++) {
					if(zk[s] < tmp1)
						count++;
				}
				omega = 1.0*count/N0;
			}
			else{
				tmp		= 0.0;
				tmp1 	= zmu[j];
				tmp2 	= rho/sqrt(1-rho*rho);
				tmp3 	= tmp1/sqrt(1-rho*rho);
				for (s = 0; s < N0; s++) {
					if(zk[s] < tmp1){
						tmp += erf(SQRT2*(tmp3 - zk[s]*tmp2)) + 1.0;
					}
				}
				omega = 0.5*tmp/N0;
			}

			OMEGA[i*n+j]	= omega*xij;
			Tn0 	+= OMEGA[i*n+j]*resid0[i]*resid0[j];
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[n+g*n+i];
		}
		EstHuber2(tx, yb, alphahat, resid, tau, n, p, maxIter, tol);
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(zmu);
	free(stdx);
	return 2.0*Tn0/n;
}

SEXP _HUBER_WAST(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP RESID0, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, isBeta, maxIter, M;
	double tau, shape1, shape2, tol;
	n 		= INTEGER(DIMs)[0];
	p1 		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	M 		= INTEGER(DIMs)[4];
	isBeta 	= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];

	tau 	= REAL(PARAMs)[0];
	shape1 	= REAL(PARAMs)[1];
	shape2 	= REAL(PARAMs)[2];
	tol 	= REAL(PARAMs)[3];

	SEXP rTn, rTn0, rAlpha, list, list_names;
	PROTECT(rTn0 		= allocVector(REALSXP, 	1));
	PROTECT(rTn 		= allocVector(REALSXP, 	M));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	if(p3==1){
		REAL(rTn0)[0] = Huber_single(REAL(rTn), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESID0), REAL(rAlpha),
						n, p1, p2, M, tau, isBeta, shape1, shape2, maxIter, tol);
	}
	else{
		REAL(rTn0)[0] = Huber_multiple(REAL(rTn), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESID0), REAL(rAlpha),
						n, p1, p2, p3, M, tau, shape1, shape2, maxIter, tol);
	}


	SET_STRING_ELT(list_names, 	0,  mkChar("Tn0"));
	SET_STRING_ELT(list_names, 	1,  mkChar("Tn"));
	SET_STRING_ELT(list_names, 	2,  mkChar("coef"));
	SET_VECTOR_ELT(list, 		0, 	rTn0);
	SET_VECTOR_ELT(list, 		1, 	rTn);
	SET_VECTOR_ELT(list, 		2, 	rAlpha);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _HUBER_WAST_APPROX(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP RESID0, SEXP Z_K, SEXP MU0, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, isBeta, maxIter, M, N0;
	double tau, shape1, shape2, tol;
	n 		= INTEGER(DIMs)[0];
	p1 		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	M 		= INTEGER(DIMs)[4];
	isBeta 	= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];
	N0 		= INTEGER(DIMs)[7];

	tau 	= REAL(PARAMs)[0];
	shape1 	= REAL(PARAMs)[1];
	shape2 	= REAL(PARAMs)[2];
	tol 	= REAL(PARAMs)[3];

	SEXP rTn, rTn0, rAlpha, list, list_names;
	PROTECT(rTn0 		= allocVector(REALSXP, 	1));
	PROTECT(rTn 		= allocVector(REALSXP, 	M));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	if(p3==1){
		REAL(rTn0)[0] = Huber_single(REAL(rTn), REAL(Y), REAL(tX), REAL(X),
						REAL(Z), REAL(RESID0), REAL(rAlpha),
						n, p1, p2, M, tau, isBeta, shape1, shape2, maxIter, tol);
	}
	else{
		REAL(rTn0)[0] = Huber_multiple_approx(REAL(rTn), REAL(Y), REAL(tX), REAL(X),
						REAL(Z), REAL(RESID0), REAL(rAlpha), REAL(Z_K), REAL(MU0),
						n, p1, p2, p3, M, tau, shape1, shape2, maxIter, tol, N0);
	}


	SET_STRING_ELT(list_names, 	0,  mkChar("Tn0"));
	SET_STRING_ELT(list_names, 	1,  mkChar("Tn"));
	SET_STRING_ELT(list_names, 	2,  mkChar("coef"));
	SET_VECTOR_ELT(list, 		0, 	rTn0);
	SET_VECTOR_ELT(list, 		1, 	rTn);
	SET_VECTOR_ELT(list, 		2, 	rAlpha);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

//-------------- Estimation of robust regression with Change Plane --------------------------------
void EstLinearCP_R(double *beta0, double *theta0, const double *tx, const double *x, const double *z, const double *y, int n, int p1, int p2, int p3, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, ologelr, nlogelr, wk, yk;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *wz, *xy, *zy;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*n*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			xy[j+p1]	+= x[j*n+i]*y[i]*tmp;
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}


	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += w[j*n+i]*w[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}

	while(step < maxstep){
		step++;
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta0, hess, xy, p, p, 1);
			}
		}

		ologelr = 0.0;
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			yk 	= y[i] - tmp;
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			wk = tmp*dsh[i];
			yk -= tmp*sh[i];

			for(j=0;j<p3;j++){
				zy[j] 	+= z[n+j*n+i]*wk*yk;
				wz[j*n+i] = z[n+j*n+i]*wk;
			}
			ologelr += yk*yk;
		}

		for(j=0; j < p3; j++){
			for(k=0; k < p3; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += wz[j*n+i]*wz[k*n+i];
				}
				hessz[j*p3+k] = tmp;
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0;j<p3;j++){
			theta[j] += theta0[j];
		}

		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				xy[j+p1]	+= x[j*n+i]*y[i]*tmp;
				w[(j+p1)*n+i] = x[j*n+i]*tmp;
			}
		}
		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta0[j];
			}
			nlogelr += tmp*tmp;
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(w);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
}

double EstHuberCP_R(double *beta0, double *theta0, const double *tx, const double *x, const double *z, const double *y1, int n, int p1, int p2, int p3, double h, int maxstep, double eps, double qq, double tau0, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, ologelr, nlogelr, meany, qr, wk, yk, tau;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *w1, *wz, *xy, *zy, *y, *med, *dpsi, *resids, *resids1;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	w1 		= (double*)malloc(sizeof(double)*n*p2);
	wz 		= (double*)malloc(sizeof(double)*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	y 		= (double*)malloc(sizeof(double)*n);
	med 	= (double*)malloc(sizeof(double)*n);
	dpsi	= (double*)malloc(sizeof(double)*p);
	resids	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);

	meany = 0.0;
	for(i=0;i<n;i++){
		meany += y1[i];
	}
	meany /= n;
	for(i=0;i<n;i++){
		y[i] = y1[i] - meany;
	}

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	ologelr = 0.0;
	for(i=0;i<n;i++){
		ologelr += y[i]*y[i];
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			xy[j+p1]	+= x[j*n+i]*y[i]*tmp;
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}

	// ############ update parameters by quadratic loss ####################
	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta[j];
			}
			nlogelr	+= tmp*tmp;
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p;j++){
				beta0[j] = beta[j];
			}
		}

		// ############ update theta ####################
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk = y[i];
			for(j=0;j<p1;j++){
				yk -= tx[j*n+i]*beta0[j];
			}
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			wk = tmp*dsh[i];
			yk -= tmp*sh[i];

			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*wk;
				zy[j] 	+= wz[j]*yk;
			}

			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0;j<p3;j++){
			theta[j] += theta0[j];
		}

		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}

			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i]	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[i] 	= tmp;
				dsh[i]	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i]	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				xy[j+p1]	+= x[j*n+i]*y[i]*tmp;
				w1[j*n+i]	= x[j*n+i]*tmp;
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p1;j++){
				tmp -= w[j*n+i]*beta0[j];
			}
			for(j=0;j<p2;j++){
				tmp -= w1[j*n+i]*beta0[p1+j];
			}
			nlogelr	+= tmp*tmp;
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n*p2;i++){
				w[n*p1+i]	= w1[i];
			}
		}
		else{
			break;
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}

	// ############ refine parameters by adaptive huber loss ####################
	for(i=0;i<n;i++){
		tmp = y[i];
		for(j=0;j<p;j++) {
			tmp -= w[j*n+i]*beta0[j];
		}
		resids[i] = tmp;
	}

	if(qq < 0){
		qr = SampleQuantile1(resids, n, 0.5);
		for(i=0;i<n;i++){
			med[i] = fabs(resids[i] - qr);
		}
		qr = SampleQuantile1(med, n, 0.5);
		qq = tau0 * qr;
	}

	ologelr = 0.0;
	for(i=0;i<n;i++){
		tmp = fabs(resids[i]);
		if(tmp>qq){
			ologelr	+= qq*(2*tmp-qq);
		}
		else{
			ologelr += tmp*tmp;
		}
	}


	step = 0;
	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0;j<p;j++) xy[j] = 0.0;
		for(j=0;j<p*p;j++) hess[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = fabs(resids[i]);
			if(tmp>qq){
				tmp1 = sqrt(qq/tmp);
				for(j=0;j<p;j++){
					dpsi[j] = w[j*n+i]*tmp1;
					xy[j] 	+= dpsi[j]*tmp1*y[i];
				}
			}
			else{
				for(j=0;j<p;j++){
					dpsi[j] = w[j*n+i];
					xy[j] 	+= dpsi[j]*y[i];
				}
			}
			for(j=0; j < p; j++){
				for(k=0; k < p; k++){
					hess[j*p+k] += dpsi[j]*dpsi[k];
				}
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta[j];
			}
			resids1[i] = tmp;
			tmp = fabs(tmp);
			if(tmp>qq){
				nlogelr	+= qq*(2*tmp-qq);
			}
			else{
				nlogelr += tmp*tmp;
			}
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p;j++){
				beta0[j] = beta[j];
			}
			for(i=0;i<n;i++){
				resids[i] = resids1[i];
			}
		}

		// ############ update theta ####################
		// ologelr = 0.0;
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			yk = resids[i];
			wk = tmp*dsh[i];
			if(fabs(yk)>qq){
				tmp1 = sqrt(qq/fabs(yk));
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk*tmp1;
					zy[j] 	+= wz[j]*tmp1*yk;
				}
			}
			else{
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk;
					zy[j] 	+= wz[j]*yk;
				}
			}

			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}
		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}

			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i]	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i]	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				w1[j*n+i] = x[j*n+i]*tmp;
			}
		}
		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p1;j++){
				tmp -= w[j*n+i]*beta0[j];
			}
			for(j=0;j<p2;j++){
				tmp -= w1[j*n+i]*beta0[p1+j];
			}
			resids1[i] = tmp;
			if(fabs(tmp)>qq){
				nlogelr += qq*(2*fabs(tmp)-qq);
			}
			else{
				nlogelr += tmp*tmp;
			}

		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*p2;i++){
				w[n*p1+i]	= w1[i];
			}
		}
		else{
			break;
		}

	}


	beta0[0] += meany;
	for(i=0;i<n;i++){
		resids[i] += beta0[0];
	}

	beta0[0] = Huber_mean(resids, &tau, n, maxstep, eps);

	free(beta);
	free(theta);
	free(w);
	free(w1);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
	free(y);
	free(med);
	free(dpsi);
	free(resids);
	free(resids1);

	return tau;
}

double EstHuberCP_2step(double *beta0, double *theta0, const double *tx, const double *x, const double *z, const double *y1, int n, int p1, int p2, int p3, double h, int maxstep, double eps, double qq, double tau0, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, bnorm, ologelr, nlogelr, meany, wk, yk, tau;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *wz, *xy, *zy, *y, *dpsi, *resids;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	y 		= (double*)malloc(sizeof(double)*n);
	dpsi	= (double*)malloc(sizeof(double)*p);
	resids	= (double*)malloc(sizeof(double)*n);

	meany = 0.0;
	for(i=0;i<n;i++){
		meany += y1[i];
	}
	meany /= n;
	for(i=0;i<n;i++){
		y[i] = y1[i] - meany;
	}

	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}
	// ############ refine parameters by adaptive huber loss ####################
	for(i=0;i<n;i++){
		tmp = y[i];
		for(j=0;j<p;j++) tmp -= w[j*n+i]*beta0[j];
		resids[i] = tmp;
	}

	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0;j<p;j++) xy[j] = 0.0;
		for(j=0;j<p*p;j++) hess[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = fabs(resids[i]);
			if(tmp>qq){
				wk = sqrt(qq/tmp);
				for(j=0;j<p;j++){
					dpsi[j]	= w[j*n+i]*wk;
					xy[j] 	+= dpsi[j]*wk*y[i];
				}
			}
			else{
				for(j=0;j<p;j++){
					dpsi[j]	= w[j*n+i];
					xy[j] 	+= dpsi[j]*y[i];
				}
			}
			for(j=0; j < p; j++){
				for(k=0; k < p; k++){
					hess[j*p+k] += dpsi[j]*dpsi[k];
				}
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp = beta[j] - beta0[j];
			bnorm += tmp*tmp;
		}

		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)	beta0[j] = beta[j];
			for(i=0;i<n;i++){
				tmp = y[i];
				for(j=0;j<p;j++) tmp -= w[j*n+i]*beta[j];
				resids[i] = tmp;
			}
		}
		// ############ update theta ####################
		ologelr = 0.0;
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			yk = resids[i];
			wk = tmp*dsh[i];
			tmp = fabs(yk);
			if(tmp > qq){
				ologelr += qq*(2*tmp-qq);
				tmp = sqrt(qq/tmp);
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk*tmp;
					zy[j] 	+= wz[j]*tmp*yk;
				}
			}
			else{
				ologelr += yk*yk;
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk;
					zy[j] 	+= wz[j]*yk;
				}
			}

			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}
		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}


		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}

			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i]	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i]	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				xy[j+p1]	+= x[j*n+i]*y[i]*tmp;
				w[(j+p1)*n+i] = x[j*n+i]*tmp;
			}
		}
		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta0[j];
			}

			resids[i] = tmp;
			tmp = fabs(tmp);
			if(tmp>qq){
				nlogelr += qq*(2*tmp-qq);
			}
			else{
				nlogelr += tmp*tmp;
			}

		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
		}
		else{
			break;
		}

	}

	beta0[0] += meany;
	for(i=0;i<n;i++){
		resids[i] += beta0[0];
	}

	beta0[0] = Huber_mean(resids, &tau, n, maxstep, eps);

	free(beta);
	free(theta);
	free(w);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
	free(y);
	free(dpsi);
	free(resids);

	return tau;
}

SEXP _EST_LinearCP(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP THETA0, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	// Outcome
	SEXP rBeta, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	for(int i=0; i<p3; i++) REAL(rTheta)[i] = REAL(THETA0)[i];


	EstLinearCP_R(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _EST_HuberCP(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP THETA0, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h, qq, tau0;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	qq 		= REAL(PARAMs)[2];
	tau0 	= REAL(PARAMs)[3];
	// Outcome
	SEXP rBeta, rTau, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rTau 		= allocVector(REALSXP, 	1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(int i=0; i<p3; i++) REAL(rTheta)[i] = REAL(THETA0)[i];

	REAL(rTau)[0] = EstHuberCP_R(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, h, maxstep, tol, qq, tau0, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("tau"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rTau);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_HuberCP_2step(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP BETA0 , SEXP THETA0, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth;
	double tol, h, qq, tau0;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	qq 		= REAL(PARAMs)[2];
	tau0 	= REAL(PARAMs)[3];
	// Outcome
	SEXP rBeta, rTau, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rTau 		= allocVector(REALSXP, 	1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(i=0; i<p1+p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(THETA0)[i];

	REAL(rTau)[0] = EstHuberCP_2step(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, h, maxstep, tol, qq, tau0, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("tau"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rTau);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}


//-------------- Estimation of weighted robust regression with Change Plane --------------------------------
double Huber_d1sum_Weight(double *x, double*g, double tau, int n){
	int i;
	double sum = 0.0;
	for(i=0;i<n;i++){
		if(fabs(x[i])>tau){
			sum += tau*SGN(x[i])*g[i];
		}
		else{
			sum += x[i]*g[i];
		}
	}
	return sum/n;
}

double Huber_mean_Weight(double *x, double *g, double *tau0, int n, int maxstep, double eps){
	int i, step=0;
	double meanx, stdx, tau, psiSum, psiSum0, psiSumDiff, mu, rate, low, up;
	double alpha, muDiff, a1, a2, cross;
	double *x1, *resid, *resid2;

	x1	 	= (double*)malloc(sizeof(double)*n);
	resid	= (double*)malloc(sizeof(double)*n);
	resid2	= (double*)malloc(sizeof(double)*n);

	rate = log(1.0*n)/n;
	meanx = 0.0;
	for(i=0;i<n;i++){
		meanx += x[i];
	}
	meanx /= n;

	stdx = 0.0;
	for(i=0;i<n;i++){
		x1[i] = x[i] - meanx;
		stdx += x1[i]*x1[i];
	}
	stdx = sqrt(stdx/(n-1));
	tau = stdx*sqrt(n/log(1.0*n));

	psiSum0 = Huber_d1sum_Weight(x1, g, tau, n);
	mu = psiSum0;
	muDiff = psiSum0;


	for(i=0;i<n;i++){
		resid[i] = x1[i] - mu;
		resid2[i] = resid[i]*resid[i]*g[i]*g[i];
	}
	low = minx(resid2, n);
	up 	= cumsum(resid2, n);
	tau = sqrt( Huber_scale(resid2, low, up, n, rate, maxstep, eps) );
	psiSum = Huber_d1sum_Weight(resid, g, tau, n);
	psiSumDiff = psiSum0 - psiSum;

	while(fabs(psiSum)> eps && step < maxstep){
		step++;

		alpha = 1.0;
		cross = muDiff*psiSumDiff;
		if(cross>0){
			a1 = cross/(psiSumDiff * psiSumDiff);
			a2 = muDiff*muDiff/cross;
			alpha = MIN( MIN(a1, a2), 100.0);
		}
		psiSum0 = psiSum;
		muDiff = alpha*psiSum;
		mu += muDiff;
		for(i=0;i<n;i++){
			resid[i] = x1[i] - mu;
			resid2[i] = resid[i]*resid[i]*g[i]*g[i];
		}


		low = minx(resid2, n);
		up 	= cumsum(resid2, n);
		tau = sqrt( Huber_scale(resid2, low, up, n, rate, maxstep, eps) );
		psiSum = Huber_d1sum_Weight(resid, g, tau, n);
		psiSumDiff = psiSum0 - psiSum;
	}

	*tau0 = tau;
	free(x1);
	free(resid);
	free(resid2);

	return mu + meanx;

}

void EstLinearCP_Weight(double *beta0, double *theta0, const double *tx, const double *x, const double *z, const double *y, double *g, int n, int p1, int p2, int p3, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, ologelr, nlogelr, wk, yk;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *wz, *xy, *zy;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*n*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i]*g[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			xy[j+p1]	+= x[j*n+i]*y[i]*tmp*g[i];
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}

	while(step < maxstep){
		step++;
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i]*g[i];
				}
				hess[j*p+k] = tmp;
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta0, hess, xy, p, p, 1);
			}
		}

		ologelr = 0.0;
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			yk 	= y[i] - tmp;
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			wk = tmp*dsh[i];
			yk -= tmp*sh[i];

			for(j=0;j<p3;j++){
				zy[j] 	+= z[n+j*n+i]*wk*yk*g[i];
				wz[j*n+i] = z[n+j*n+i]*wk;
			}
			ologelr += yk*yk*g[i];
		}

		for(j=0; j < p3; j++){
			for(k=0; k < p3; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += wz[j*n+i]*wz[k*n+i]*g[i];
				}
				hessz[j*p3+k] = tmp;
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0;j<p3;j++){
			theta[j] += theta0[j];
		}
		// printArrayDouble(beta0, p, p);
		// printArrayDouble(theta, p3, p3);

		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				xy[j+p1]	+= x[j*n+i]*y[i]*tmp*g[i];
				w[(j+p1)*n+i] = x[j*n+i]*tmp;
			}
		}
		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta0[j];
			}
			nlogelr += tmp*tmp*g[i];
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(w);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
}

double EstHuberCP_Weight(double *beta0, double *theta0, const double *tx, const double *x, const double *z, const double *y1, double *g, int n, int p1, int p2, int p3, double h, int maxstep, double eps, double qq, double tau0, int smooth){
	int i,j,k, step=0, p=p1+p2, flag=1;
	double tmp, tmp1, phix, ologelr, nlogelr, meany, qr, wk, yk, tau;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *w1, *wz, *xy, *zy, *y, *med, *dpsi, *resids, *resids1;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n);
	dsh 	= (double*)malloc(sizeof(double)*n);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	w1 		= (double*)malloc(sizeof(double)*n*p2);
	wz 		= (double*)malloc(sizeof(double)*p3);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	y 		= (double*)malloc(sizeof(double)*n);
	med 	= (double*)malloc(sizeof(double)*n);
	dpsi	= (double*)malloc(sizeof(double)*p);
	resids	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);

	meany = 0.0;
	for(i=0;i<n;i++){
		meany += y1[i];
	}
	meany /= n;
	for(i=0;i<n;i++){
		y[i] = y1[i] - meany;
	}

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i]*g[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	ologelr = 0.0;
	for(i=0;i<n;i++){
		ologelr += y[i]*y[i]*g[i];
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			xy[j+p1]	+= x[j*n+i]*y[i]*tmp*g[i];
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}

	// ############ update parameters by quadratic loss ####################
	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i]*g[i];
				}
				hess[j*p+k] = tmp;
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta[j];
			}
			nlogelr	+= tmp*tmp*g[i];
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p;j++){
				beta0[j] = beta[j];
			}
		}
		// ############ update theta ####################
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			yk 	= y[i] - tmp;
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			wk = tmp*dsh[i];
			yk -= tmp*sh[i];

			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*wk;
				zy[j] 	+= wz[j]*yk*g[i];
			}

			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*g[i];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0;j<p3;j++){
			theta[j] += theta0[j];
		}


		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}

			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i]	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[i] 	= tmp;
				dsh[i]	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i]	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				xy[j+p1]	+= x[j*n+i]*y[i]*tmp*g[i];
				w1[j*n+i]	= x[j*n+i]*tmp;
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p1;j++){
				tmp -= w[j*n+i]*beta0[j];
			}
			for(j=0;j<p2;j++){
				tmp -= w1[j*n+i]*beta0[p1+j];
			}
			nlogelr	+= tmp*tmp*g[i];
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n*p2;i++){
				w[n*p1+i]	= w1[i];
			}
		}
		else{
			break;
		}
	}

	for(i=0;i<n;i++){
		tmp = z[i];
		for(j=0;j<p3;j++){
			tmp	+= z[n+j*n+i]*theta0[j];
		}
		if(smooth==1){
			tmp 	= 1.0/(1.0+exp(-tmp*h));
			sh[i] 	= tmp;
			dsh[i] 	= h*tmp*(1-tmp);
		}
		else if(smooth==2){
			tmp1 	= tmp*h;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
			sh[i] 	= tmp;
			dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
		}
		else{
			tmp1 	= tmp*h;
			phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
			tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
			sh[i] 	= tmp;
			dsh[i] 	= h*phix*(2 - tmp1*tmp1);
		}

		for(j=0;j<p2;j++){
			w[(j+p1)*n+i] = x[j*n+i]*tmp;
		}
	}

	// ############ refine parameters by adaptive huber loss ####################
	for(i=0;i<n;i++){
		tmp = y[i];
		for(j=0;j<p;j++) {
			tmp -= w[j*n+i]*beta0[j];
		}
		resids[i] = tmp;
	}

	if(qq < 0){
		qr = SampleQuantile1(resids, n, 0.5);
		for(i=0;i<n;i++){
			med[i] = fabs(resids[i] - qr);
		}
		qr = SampleQuantile1(med, n, 0.5);
		qq = tau0 * qr;
	}

	ologelr = 0.0;
	for(i=0;i<n;i++){
		tmp = fabs(resids[i]);
		if(tmp>qq){
			ologelr	+= qq*(2*tmp-qq)*g[i];
		}
		else{
			ologelr += tmp*tmp*g[i];
		}
	}

	step = 0;
	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0;j<p;j++) xy[j] = 0.0;
		for(j=0;j<p*p;j++) hess[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = fabs(resids[i]);
			if(tmp>qq){
				wk = sqrt(qq/tmp);
				for(j=0;j<p;j++){
					dpsi[j] = w[j*n+i]*wk;
					xy[j] 	+= dpsi[j]*wk*y[i]*g[i];
				}
			}
			else{
				for(j=0;j<p;j++){
					dpsi[j] = w[j*n+i];
					xy[j] 	+= dpsi[j]*y[i]*g[i];
				}
			}
			for(j=0; j < p; j++){
				for(k=0; k < p; k++){
					hess[j*p+k] += dpsi[j]*dpsi[k]*g[i];
				}
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta[j];
			}
			resids1[i] = tmp;
			tmp = fabs(tmp);
			if(tmp>qq){
				nlogelr	+= qq*(2*tmp-qq)*g[i];
			}
			else{
				nlogelr += tmp*tmp*g[i];
			}
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p;j++){
				beta0[j] = beta[j];
			}
			for(i=0;i<n;i++){
				resids[i] = resids1[i];
			}
		}


		// ############ update theta ####################
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p2;j++){
				tmp += x[j*n+i]*beta0[j+p1];
			}
			yk = resids[i];
			wk = tmp*dsh[i];
			if(fabs(yk)>qq){
				tmp1 = sqrt(qq/fabs(yk));
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk*tmp1;
					zy[j] 	+= wz[j]*yk*tmp1*g[i];
				}
			}
			else{
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk;
					zy[j] 	+= wz[j]*yk*g[i];
				}
			}

			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*g[i];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}
		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		for(i=0;i<n;i++){
			tmp = z[i];
			for(j=0;j<p3;j++){
				tmp	+= z[n+j*n+i]*theta[j];
			}

			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[i] 	= tmp;
				dsh[i]	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[i] 	= tmp;
				dsh[i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[i] 	= tmp;
				dsh[i]	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				w1[j*n+i] = x[j*n+i]*tmp;
			}
		}
		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p1;j++){
				tmp -= w[j*n+i]*beta0[j];
			}
			for(j=0;j<p2;j++){
				tmp -= w1[j*n+i]*beta0[p1+j];
			}
			resids[i] = tmp;
			if(fabs(tmp)>qq){
				nlogelr += qq*(2*fabs(tmp)-qq)*g[i];
			}
			else{
				nlogelr += tmp*tmp*g[i];
			}

		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n*p2;i++){
				w[n*p1+i]	= w1[i];
			}
		}
		else{
			break;
		}

	}


	beta0[0] += meany;
	for(i=0;i<n;i++){
		resids[i] += beta0[0];
	}

	beta0[0] = Huber_mean_Weight(resids, g, &tau, n, maxstep, eps);

	free(beta);
	free(theta);
	free(w);
	free(w1);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
	free(y);
	free(med);
	free(dpsi);
	free(resids);
	free(resids1);

	return tau;
}

SEXP _EST_LinearCP_Weight(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP THETA0, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	// Outcome
	SEXP rBeta, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	for(int i=0; i<p3; i++) REAL(rTheta)[i] = REAL(THETA0)[i];


	EstLinearCP_Weight(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _EST_HuberCP_Weight(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP THETA0, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, maxstep, smooth;
	double tol, h, qq, tau0;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	qq 		= REAL(PARAMs)[2];
	tau0 	= REAL(PARAMs)[3];
	// Outcome
	SEXP rBeta, rTau, rTheta, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rTau 		= allocVector(REALSXP, 	1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(int i=0; i<p3; i++) REAL(rTheta)[i] = REAL(THETA0)[i];

	REAL(rTau)[0] = EstHuberCP_Weight(REAL(rBeta), REAL(rTheta), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, h, maxstep, tol, qq, tau0, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("tau"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rTau);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

//-------------- Estimation of robust regression with multiple Change Plane --------------------------------
void _EstLinearCP_Multi(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+ng*p2, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, wk, yk;
	double *beta, *theta, *sh, *dsh, *sh1, *dsh1, *hess, *hessz, *hessa, *w, *wz, *wa, *ha, *xy, *zy, *ay, *resids, *resids1, *xb;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp;
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	while(step < maxstep){
		step++;

		// calculate alpha and beta
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta0, hess, xy, p, p, 1);
			}
		}

		// calculate ak
		ologelr = 0.0;
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= y[i] - tmp;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				yk -= tmp*sh[k*n+i];
			}
			resids[i] = yk;
			ologelr += yk*yk;
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk;
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k];
				}
			}
		}

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[(j+1)*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk;
			}
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}


		// calculate theta
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*wk;
				zy[j] 	+= wz[j]*yk;
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		nlogelr = 0.0;
		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp;
					w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			nlogelr += yk*yk;
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(sh);
	free(dsh);
	free(sh1);
	free(dsh1);
	free(hess);
	free(hessz);
	free(hessa);
	free(w);
	free(wz);
	free(wa);
	free(ha);
	free(xy);
	free(zy);
	free(ay);
	free(resids);
	free(resids1);
	free(xb);
}

double _EstHuberCP_Multi(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y1, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, double qq, double tau0, int smooth){
	int i,j,k, step=0, p=p1+p2*ng, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, meany, qr, wk, yk, tau;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *w1, *wz, *xy, *zy, *y, *med, *dpsi, *resids, *resids1;
	double *sh1, *dsh1, *ha, *hessa, *xb, *ay, *wa;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	w1 		= (double*)malloc(sizeof(double)*n*p2*ng);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	y 		= (double*)malloc(sizeof(double)*n);
	med 	= (double*)malloc(sizeof(double)*n);
	dpsi	= (double*)malloc(sizeof(double)*p);
	resids	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	meany = 0.0;
	for(i=0;i<n;i++){
		meany += y1[i];
	}
	meany /= n;
	for(i=0;i<n;i++){
		y[i] = y1[i] - meany;
	}

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	ologelr = 0.0;
	for(i=0;i<n;i++){
		ologelr += y[i]*y[i];
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			for(j=0;j<p2;j++){
				xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp;
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	// ############ update parameters by quadratic loss ####################
	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta[j];
			}
			nlogelr	+= tmp*tmp;
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p;j++){
				beta0[j] = beta[j];
			}
		}

		// ############ update ak ####################
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= y[i] - tmp;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				yk -= tmp*sh[k*n+i];
			}
			resids[i] = yk;
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk;
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k];
				}
			}
		}

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[n+j*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk;
			}
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}

		// ############ update theta ####################
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*wk;
				zy[j] 	+= wz[j]*yk;
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0;j<p3;j++){
			theta[j] += theta0[j];
		}


		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp;
					w1[(k*p2+j)*n+i] = x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			nlogelr += yk*yk;
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
			for(i=0;i<n*p2*ng;i++){
				w[n*p1+i]	= w1[i];
			}
		}
		else{
			break;
		}
	}


	// for(j=p1;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			for(j=0;j<p2;j++){
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	// ############ refine parameters by adaptive huber loss ####################
	for(i=0;i<n;i++){
		tmp = y[i];
		for(j=0;j<p;j++) {
			tmp -= w[j*n+i]*beta0[j];
		}
		resids[i] = tmp;
	}

	if(qq < 0){
		qr = SampleQuantile1(resids, n, 0.5);
		for(i=0;i<n;i++){
			med[i] = fabs(resids[i] - qr);
		}
		qr = SampleQuantile1(med, n, 0.5);
		qq = tau0 * qr;
	}

	ologelr = 0.0;
	for(i=0;i<n;i++){
		tmp = fabs(resids[i]);
		if(tmp>qq){
			ologelr	+= qq*(2*tmp-qq);
		}
		else{
			ologelr += tmp*tmp;
		}
	}


	step = 0;
	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0;j<p;j++) xy[j] = 0.0;
		for(j=0;j<p*p;j++) hess[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = fabs(resids[i]);
			if(tmp>qq){
				tmp1 = sqrt(qq/tmp);
				for(j=0;j<p;j++){
					dpsi[j] = w[j*n+i]*tmp1;
					xy[j] 	+= dpsi[j]*tmp1*y[i];
				}
			}
			else{
				for(j=0;j<p;j++){
					dpsi[j] = w[j*n+i];
					xy[j] 	+= dpsi[j]*y[i];
				}
			}
			for(j=0; j < p; j++){
				for(k=0; k < p; k++){
					hess[j*p+k] += dpsi[j]*dpsi[k];
				}
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta[j];
			}
			resids1[i] = tmp;
			tmp = fabs(tmp);
			if(tmp>qq){
				nlogelr	+= qq*(2*tmp-qq);
			}
			else{
				nlogelr += tmp*tmp;
			}
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p;j++){
				beta0[j] = beta[j];
			}
			for(i=0;i<n;i++){
				resids[i] = resids1[i];
			}
		}

		// ############ update ak ####################
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= resids[i];
			if(fabs(yk)>qq){
				tmp1 = sqrt(qq/fabs(yk));
				for(k=0; k < ng; k++){
					tmp = 0.0;
					for(j=0;j<p2;j++){
						tmp += x[j*n+i]*beta0[p1+k*p2+j];
					}
					xb[n+k*n+i] = tmp;
					wa[k] = tmp*dsh[k*n+i]*tmp1;
					ay[k] += wa[k]*yk*tmp1;
				}
			}
			else{
				for(k=0; k < ng; k++){
					tmp = 0.0;
					for(j=0;j<p2;j++){
						tmp += x[j*n+i]*beta0[p1+k*p2+j];
					}
					xb[n+k*n+i] = tmp;
					wa[k] = tmp*dsh[k*n+i];
					ay[k] += wa[k]*yk;
				}
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k];
				}
			}
		}

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[n+j*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk;
			}
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}

		// ############ update theta ####################
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			if(fabs(yk)>qq){
				tmp1 = sqrt(qq/fabs(yk));
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk*tmp1;
					zy[j] 	+= wz[j]*tmp1*yk;
				}
			}
			else{
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk;
					zy[j] 	+= wz[j]*yk;
				}
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}
		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					w1[(k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			if(fabs(yk)>qq){
				nlogelr += qq*(2*fabs(yk)-qq);
			}
			else{
				nlogelr += yk*yk;
			}
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*p2*ng;i++){
				w[n*p1+i]	= w1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}

	}


	beta0[0] += meany;
	for(i=0;i<n;i++){
		resids[i] += beta0[0];
	}

	beta0[0] = Huber_mean(resids, &tau, n, maxstep, eps);

	free(beta);
	free(theta);
	free(w);
	free(w1);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
	free(y);
	free(med);
	free(dpsi);
	free(resids);
	free(resids1);
	free(sh1);
	free(dsh1);
	free(ha);
	free(hessa);
	free(xb);
	free(ay);
	free(wa);

	return tau;
}

double _EstHuberCP_Multi_2step(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y1, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, double qq, int smooth){
	int i,j,k, step=0, p=p1+p2*ng, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, meany, wk, yk, tau;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *w1, *wz, *xy, *zy, *y, *dpsi, *resids, *resids1;
	double *sh1, *dsh1, *ha, *hessa, *xb, *ay, *wa;

	y 		= (double*)malloc(sizeof(double)*n);
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	w1 		= (double*)malloc(sizeof(double)*n*p2*ng);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	dpsi	= (double*)malloc(sizeof(double)*p);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));


	meany = 0.0;
	for(i=0;i<n;i++){
		meany += y1[i];
	}
	meany /= n;
	for(i=0;i<n;i++){
		y[i] = y1[i] - meany;
	}

	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}

			for(j=0;j<p2;j++){
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	// ############ refine parameters by adaptive huber loss ####################
	for(i=0;i<n;i++){
		tmp = y[i];
		for(j=0;j<p;j++) tmp -= w[j*n+i]*beta0[j];
		resids[i] = tmp;
	}

	ologelr = 0.0;
	for(i=0;i<n;i++){
		tmp = fabs(resids[i]);
		if(tmp>qq){
			ologelr	+= qq*(2*tmp-qq);
		}
		else{
			ologelr += tmp*tmp;
		}
	}

	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0;j<p;j++) xy[j] = 0.0;
		for(j=0;j<p*p;j++) hess[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = fabs(resids[i]);
			if(tmp>qq){
				wk = sqrt(qq/tmp);
				for(j=0;j<p;j++){
					dpsi[j]	= w[j*n+i]*wk;
					xy[j] 	+= dpsi[j]*wk*y[i];
				}
			}
			else{
				for(j=0;j<p;j++){
					dpsi[j]	= w[j*n+i];
					xy[j] 	+= dpsi[j]*y[i];
				}
			}
			for(j=0; j < p; j++){
				for(k=0; k < p; k++){
					hess[j*p+k] += dpsi[j]*dpsi[k];
				}
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta[j];
			}
			resids1[i] = tmp;
			tmp = fabs(tmp);
			if(tmp>qq){
				nlogelr	+= qq*(2*tmp-qq);
			}
			else{
				nlogelr += tmp*tmp;
			}
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p;j++){
				beta0[j] = beta[j];
			}
			for(i=0;i<n;i++){
				resids[i] = resids1[i];
			}
		}

		// ############ update ak ####################
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= resids[i];
			if(fabs(yk)>qq){
				tmp1 = sqrt(qq/fabs(yk));
				for(k=0; k < ng; k++){
					tmp = 0.0;
					for(j=0;j<p2;j++){
						tmp += x[j*n+i]*beta0[p1+k*p2+j];
					}
					xb[n+k*n+i] = tmp;
					wa[k] = tmp*dsh[k*n+i]*tmp1;
					ay[k] += wa[k]*yk*tmp1;
				}
			}
			else{
				for(k=0; k < ng; k++){
					tmp = 0.0;
					for(j=0;j<p2;j++){
						tmp += x[j*n+i]*beta0[p1+k*p2+j];
					}
					xb[n+k*n+i] = tmp;
					wa[k] = tmp*dsh[k*n+i];
					ay[k] += wa[k]*yk;
				}
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k];
				}
			}			
		}

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[n+j*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk;
			}
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}

		// ############ update theta ####################
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			if(fabs(yk)>qq){
				tmp1 = sqrt(qq/fabs(yk));
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk*tmp1;
					zy[j] 	+= wz[j]*tmp1*yk;
				}
			}
			else{
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk;
					zy[j] 	+= wz[j]*yk;
				}
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}
		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}


		nlogelr = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					w1[(k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			if(fabs(yk)>qq){
				nlogelr += qq*(2*fabs(yk)-qq);
			}
			else{
				nlogelr += yk*yk;
			}
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*p2*ng;i++){
				w[n*p1+i]	= w1[i];
			}
		}
		else{
			break;
		}

	}

	beta0[0] += meany;
	for(i=0;i<n;i++){
		resids[i] += beta0[0];
	}

	beta0[0] = Huber_mean(resids, &tau, n, maxstep, eps);

	free(beta);
	free(theta);
	free(w);
	free(w1);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
	free(y);
	free(dpsi);
	free(resids);
	free(resids1);
	free(sh1);
	free(dsh1);
	free(ha);
	free(hessa);
	free(xb);
	free(ay);
	free(wa);

	return tau;
}

SEXP _EST_LINEAR_MULTI(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	// Outcome
	SEXP rBeta, rTheta, rHa, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	_EstLinearCP_Multi(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, ng, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_HUBER_MULTI(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h, qq, tau0;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	qq 		= REAL(PARAMs)[2];
	tau0 	= REAL(PARAMs)[3];
	// Outcome
	SEXP rBeta, rTheta, rHa, rTau, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rTau 		= allocVector(REALSXP, 	1));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	4));
	PROTECT(list 		= allocVector(VECSXP, 	4));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	REAL(rTau)[0] = _EstHuberCP_Multi(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, ng, h, maxstep, tol, qq, tau0, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("tau"));
	SET_STRING_ELT(list_names, 	3,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rTau);
	SET_VECTOR_ELT(list, 		3, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(6);
	return list;
}

SEXP _EST_HUBER_MULTI_2step(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h, qq;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	qq 		= REAL(PARAMs)[2];

	// Outcome
	SEXP rBeta, rTheta, rHa, rTau, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rTau 		= allocVector(REALSXP, 	1));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	4));
	PROTECT(list 		= allocVector(VECSXP, 	4));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	REAL(rTau)[0] = _EstHuberCP_Multi_2step(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), n, p1, p2, p3, ng, h, maxstep, tol, qq, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("tau"));
	SET_STRING_ELT(list_names, 	3,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rTau);
	SET_VECTOR_ELT(list, 		3, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(6);
	return list;
}

//-------------- Estimation of weighted robust regression with multiple Change Plane --------------------------
void _EstLinearCP_Multi_Weight(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y, double *g, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, int smooth){
	int i,j,k, step=0, p=p1+ng*p2, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, wk, yk;
	double *beta, *theta, *sh, *dsh, *sh1, *dsh1, *hess, *hessz, *hessa, *w, *wz, *wa, *ha, *xy, *zy, *ay, *resids, *resids1, *xb;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	resids 	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i]*g[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;	
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			

			for(j=0;j<p2;j++){
				xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp*g[i];
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	while(step < maxstep){
		step++;

		// calculate alpha and beta
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i]*g[i];
				}
				hess[j*p+k] = tmp;
			}
		}
		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta0[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta0, hess, xy, p, p, 1);
			}
		}

		// calculate ak
		ologelr = 0.0;
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= y[i] - tmp;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				yk -= tmp*sh[k*n+i];
			}
			resids[i] = yk;
			ologelr += yk*yk*g[i];
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk*g[i];
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k]*g[i];
				}
			}			
		}

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[(j+1)*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk*g[i];
			}
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}


		// calculate theta
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*wk;
				zy[j] 	+= wz[j]*yk*g[i];
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*g[i];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		nlogelr = 0.0;
		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;	
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;	
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp;
					w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			nlogelr += yk*yk*g[i];
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}
	}

	free(beta);
	free(theta);
	free(sh);
	free(dsh);
	free(sh1);
	free(dsh1);
	free(hess);
	free(hessz);
	free(hessa);
	free(w);
	free(wz);
	free(wa);
	free(ha);
	free(xy);
	free(zy);
	free(ay);
	free(resids);
	free(resids1);
	free(xb);
}

double _EstHuberCP_Multi_Weight(double *beta0, double *theta0, double *ha0, double *tx, double *x, double *z, double *y1, double *g, int n, int p1, int p2, int p3, int ng, double h, int maxstep, double eps, double qq, double tau0, int smooth){
	int i,j,k, step=0, p=p1+p2*ng, flag=1, count=0;
	double tmp, tmp1, zt0, phix, ologelr, nlogelr, meany, qr, wk, yk, tau;
	double *beta, *theta, *sh, *dsh, *hess, *hessz, *w, *w1, *wz, *xy, *zy, *y, *med, *dpsi, *resids, *resids1;
	double *sh1, *dsh1, *ha, *hessa, *xb, *ay, *wa;
	xy 		= (double*)malloc(sizeof(double)*p);
	zy 		= (double*)malloc(sizeof(double)*p3);
	ay 		= (double*)malloc(sizeof(double)*p3);
	sh 		= (double*)malloc(sizeof(double)*n*ng);
	dsh 	= (double*)malloc(sizeof(double)*n*ng);
	sh1 	= (double*)malloc(sizeof(double)*n*ng);
	dsh1 	= (double*)malloc(sizeof(double)*n*ng);
	ha 		= (double*)malloc(sizeof(double)*ng);
	beta 	= (double*)malloc(sizeof(double)*p);
	theta 	= (double*)malloc(sizeof(double)*p3);
	w 		= (double*)malloc(sizeof(double)*n*p);
	w1 		= (double*)malloc(sizeof(double)*n*p2*ng);
	wz 		= (double*)malloc(sizeof(double)*p3);
	wa 		= (double*)malloc(sizeof(double)*ng);
	hess	= (double*)malloc(sizeof(double)*p*p);
	hessz	= (double*)malloc(sizeof(double)*p3*p3);
	hessa	= (double*)malloc(sizeof(double)*ng*ng);
	y 		= (double*)malloc(sizeof(double)*n);
	med 	= (double*)malloc(sizeof(double)*n);
	dpsi	= (double*)malloc(sizeof(double)*p);
	resids	= (double*)malloc(sizeof(double)*n);
	resids1	= (double*)malloc(sizeof(double)*n);
	xb 		= (double*)malloc(sizeof(double)*n*(ng+1));

	meany = 0.0;
	for(i=0;i<n;i++){
		meany += y1[i];
	}
	meany /= n;
	for(i=0;i<n;i++){
		y[i] = y1[i] - meany;
	}

	for(j=0;j<p;j++)	beta0[j] = 0.0;
	for(j=0;j<p;j++) 	xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p1;j++){
			xy[j] 	+= tx[j*n+i]*y[i]*g[i];
		}
	}
	for(j=0;j<p1;j++){
		for(i=0;i<n;i++){
			w[j*n+i] = tx[j*n+i];
		}
	}

	ologelr = 0.0;
	for(i=0;i<n;i++){
		ologelr += y[i]*y[i]*g[i];
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			for(j=0;j<p2;j++){
				xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp*g[i];
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	// ############ update parameters by quadratic loss ####################
	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += w[j*n+i]*w[k*n+i]*g[i];
				}
				hess[j*p+k] = tmp;
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta[j];
			}
			nlogelr	+= tmp*tmp*g[i];
		}
		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p;j++){
				beta0[j] = beta[j];
			}
		}

		// ############ update ak ####################
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= y[i] - tmp;
			for(k=0; k < ng; k++){
				tmp = 0.0;
				for(j=0;j<p2;j++){
					tmp += x[j*n+i]*beta0[p1+k*p2+j];
				}
				xb[n+k*n+i] = tmp;
				wa[k] = tmp*dsh[k*n+i];
				yk -= tmp*sh[k*n+i];
			}
			resids[i] = yk;
			for(k=0; k < ng; k++){
				ay[k] += wa[k]*yk*g[i];
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k]*g[i];
				}
			}
		}

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[n+j*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk*g[i];
			}
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}

		// ############ update theta ####################
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			for(j=0;j<p3;j++){
				wz[j]	= z[n+j*n+i]*wk;
				zy[j] 	+= wz[j]*yk*g[i];
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*g[i];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}

		for(j=0;j<p3;j++){
			theta[j] += theta0[j];
		}


		for(j=p1;j<p;j++) xy[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					xy[p1+k*p2+j]		+= x[j*n+i]*y[i]*tmp*g[i];
					w1[(k*p2+j)*n+i] = x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			nlogelr += yk*yk*g[i];
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
			for(i=0;i<n*p2*ng;i++){
				w[n*p1+i]	= w1[i];
			}
		}
		else{
			break;
		}
	}


	// for(j=p1;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		zt0 = z[i];
		for(j=0;j<p3;j++){
			zt0	+= z[(j+1)*n+i]*theta0[j];
		}
		for(k=0; k < ng; k++){
			tmp = zt0 - ha0[k];
			if(smooth==1){
				tmp 	= 1.0/(1.0+exp(-tmp*h));
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*tmp*(1-tmp);
			}
			else if(smooth==2){
				tmp1 	= tmp*h;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*exp(-0.5*tmp1*tmp1)*MPI2;
			}
			else{
				tmp1 	= tmp*h;
				phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
				tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
				sh[k*n+i] 	= tmp;
				dsh[k*n+i] 	= h*phix*(2 - tmp1*tmp1);
			}
			for(j=0;j<p2;j++){
				w[(p1+k*p2+j)*n+i] 	= x[j*n+i]*tmp;
			}
		}
	}

	// ############ refine parameters by adaptive huber loss ####################
	for(i=0;i<n;i++){
		tmp = y[i];
		for(j=0;j<p;j++) {
			tmp -= w[j*n+i]*beta0[j];
		}
		resids[i] = tmp;
	}

	if(qq < 0){
		qr = SampleQuantile1(resids, n, 0.5);
		for(i=0;i<n;i++){
			med[i] = fabs(resids[i] - qr);
		}
		qr = SampleQuantile1(med, n, 0.5);
		qq = tau0 * qr;
	}

	ologelr = 0.0;
	for(i=0;i<n;i++){
		tmp = fabs(resids[i]);
		if(tmp>qq){
			ologelr	+= qq*(2*tmp-qq)*g[i];
		}
		else{
			ologelr += tmp*tmp*g[i];
		}
	}

	// printArrayDouble(beta0, p, p);
	// printArrayDouble(theta, p3, p3);

	step = 0;
	while(step < maxstep){
		step++;

		// ############ update beta ####################
		for(j=0;j<p;j++) xy[j] = 0.0;
		for(j=0;j<p*p;j++) hess[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = fabs(resids[i]);
			if(tmp>qq){
				tmp1 = sqrt(qq/tmp);
				for(j=0;j<p;j++){
					dpsi[j] = w[j*n+i]*tmp1;
					xy[j] 	+= dpsi[j]*tmp1*y[i]*g[i];
				}
			}
			else{
				for(j=0;j<p;j++){
					dpsi[j] = w[j*n+i];
					xy[j] 	+= dpsi[j]*y[i]*g[i];
				}
			}
			for(j=0; j < p; j++){
				for(k=0; k < p; k++){
					hess[j*p+k] += dpsi[j]*dpsi[k]*g[i];
				}
			}
		}

		if(p<2){
			if(hess[0]<MEPS){
				break;
			}
			else{
				beta[0] = xy[0]/hess[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hess,p);
			if(flag<0){
				break;
			}
			else{
				AbyB(beta, hess, xy, p, p, 1);
			}
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			tmp = y[i];
			for(j=0;j<p;j++){
				tmp -= w[j*n+i]*beta[j];
			}
			resids1[i] = tmp;
			tmp = fabs(tmp);
			if(tmp>qq){
				nlogelr	+= qq*(2*tmp-qq)*g[i];
			}
			else{
				nlogelr += tmp*tmp*g[i];
			}
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p;j++){
				beta0[j] = beta[j];
			}
			for(i=0;i<n;i++){
				resids[i] = resids1[i];
			}
		}

		// ############ update ak ####################
		for(j=0;j<ng;j++) ay[j] = 0.0;
		for(j=0;j<ng*ng;j++) hessa[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p1;j++){
				tmp += tx[j*n+i]*beta0[j];
			}
			xb[i] = tmp;
			yk 	= resids[i];
			if(fabs(yk)>qq){
				tmp1 = sqrt(qq/fabs(yk));
				for(k=0; k < ng; k++){
					tmp = 0.0;
					for(j=0;j<p2;j++){
						tmp += x[j*n+i]*beta0[p1+k*p2+j];
					}
					xb[n+k*n+i] = tmp;
					wa[k] = tmp*dsh[k*n+i]*tmp1;
					ay[k] += wa[k]*yk*tmp1*g[i];
				}
			}
			else{
				for(k=0; k < ng; k++){
					tmp = 0.0;
					for(j=0;j<p2;j++){
						tmp += x[j*n+i]*beta0[p1+k*p2+j];
					}
					xb[n+k*n+i] = tmp;
					wa[k] = tmp*dsh[k*n+i];
					ay[k] += wa[k]*yk*g[i];
				}
			}

			for(j=0; j < ng; j++){
				for(k=0; k < ng; k++){
					hessa[j*ng+k] += wa[j]*wa[k]*g[i];
				}
			}
		}

		if(ng<2){
			if(hessa[0]<MEPS){
				break;
			}
			else{
				ha[0] = ay[0]/hessa[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessa, ng);
			if(flag<0){
				break;
			}
			else{
				AbyB(ha, hessa, ay, ng, ng, 1);
			}
		}
		for(j=0; j < ng; j++){
			ha[j] = ha0[j] - ha[j];
		}

		count = 0;
		for(k=1; k < ng; k++){
			if(ha[k-1]>ha[k]){
				count++;
				break;
			}
		}
		if(!count){
			nlogelr = 0.0;
			for(i=0;i<n;i++){
				yk 	= y[i] - xb[i];
				zt0 = z[i];
				for(j=0;j<p3;j++){
					zt0	+= z[n+j*n+i]*theta0[j];
				}
				for(k=0; k < ng; k++){
					tmp = zt0 - ha[k];
					if(smooth==1){
						tmp 	= 1.0/(1.0+exp(-tmp*h));
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*tmp*(1-tmp);
					}
					else if(smooth==2){
						tmp1 	= tmp*h;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
					}
					else{
						tmp1 	= tmp*h;
						phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
						tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
						sh1[k*n+i] 	= tmp;
						dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
					}
					yk -= xb[n+k*n+i]*tmp;
				}
				resids1[i] = yk;
				nlogelr += yk*yk*g[i];
			}
			if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
				ologelr = nlogelr;
				for(j=0;j<ng;j++){
					ha0[j] = ha[j];
				}
				for(i=0;i<n;i++){
					resids[i]	= resids1[i];
				}
				for(i=0;i<n*ng;i++){
					sh[i] 	= sh1[i];
					dsh[i] 	= dsh1[i];
				}
			}
		}

		// ############ update theta ####################
		for(j=0;j<p3;j++) zy[j] = 0.0;
		for(j=0;j<p3*p3;j++) hessz[j] = 0.0;
		for(i=0;i<n;i++){
			yk 	= resids[i];
			wk = 0.0;
			for(k=0; k < ng; k++){
				wk += xb[n+k*n+i]*dsh[k*n+i];
			}
			if(fabs(yk)>qq){
				tmp1 = sqrt(qq/fabs(yk));
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk*tmp1;
					zy[j] 	+= wz[j]*tmp1*yk*g[i];
				}
			}
			else{
				for(j=0;j<p3;j++){
					wz[j]	= z[n+j*n+i]*wk;
					zy[j] 	+= wz[j]*yk*g[i];
				}
			}
			for(j=0; j < p3; j++){
				for(k=0; k < p3; k++){
					hessz[j*p3+k] += wz[j]*wz[k]*g[i];
				}
			}
		}

		if(p3<2){
			if(hessz[0]<MEPS){
				break;
			}
			else{
				theta[0] = zy[0]/hessz[0];
			}
		}
		else{
			flag = MatrixInvSymmetric(hessz,p3);
			if(flag<0){
				break;
			}
			else{
				AbyB(theta, hessz, zy, p3, p3, 1);
			}
		}
		for(j=0; j < p3; j++){
			theta[j] += theta0[j];
		}

		nlogelr = 0.0;
		for(i=0;i<n;i++){
			yk 	= y[i] - xb[i];
			zt0 = z[i];
			for(j=0;j<p3;j++){
				zt0	+= z[n+j*n+i]*theta[j];
			}
			for(k=0; k < ng; k++){
				tmp = zt0 - ha0[k];
				if(smooth==1){
					tmp 	= 1.0/(1.0+exp(-tmp*h));
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*tmp*(1-tmp);
				}
				else if(smooth==2){
					tmp1 	= tmp*h;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*exp(-0.5*tmp1*tmp1)*MPI2;
				}
				else{
					tmp1 	= tmp*h;
					phix 	= exp(-0.5*tmp1*tmp1)*MPI2;
					tmp 	= 0.5*erf(SQRT2*tmp1) +0.5 + tmp1*phix;
					sh1[k*n+i] 	= tmp;
					dsh1[k*n+i] = h*phix*(2 - tmp1*tmp1);
				}
				yk -= xb[n+k*n+i]*tmp;
				for(j=0;j<p2;j++){
					w1[(k*p2+j)*n+i] 	= x[j*n+i]*tmp;
				}
			}
			resids1[i] = yk;
			if(fabs(yk)>qq){
				nlogelr += qq*(2*fabs(yk)-qq)*g[i];
			}
			else{
				nlogelr += yk*yk*g[i];
			}
		}

		if((nlogelr<ologelr) && (1 - nlogelr/ologelr > eps)){
			ologelr = nlogelr;
			for(j=0;j<p3;j++){
				theta0[j] = theta[j];
			}
			for(i=0;i<n;i++){
				resids[i]	= resids1[i];
			}
			for(i=0;i<n*p2*ng;i++){
				w[n*p1+i]	= w1[i];
			}
			for(i=0;i<n*ng;i++){
				sh[i] 	= sh1[i];
				dsh[i] 	= dsh1[i];
			}
		}
		else{
			break;
		}

	}


	beta0[0] += meany;
	for(i=0;i<n;i++){
		resids[i] += beta0[0];
	}

	beta0[0] = Huber_mean_Weight(resids, g, &tau, n, maxstep, eps);

	free(beta);
	free(theta);
	free(w);
	free(w1);
	free(wz);
	free(sh);
	free(dsh);
	free(hess);
	free(hessz);
	free(xy);
	free(zy);
	free(y);
	free(med);
	free(dpsi);
	free(resids);
	free(resids1);
	free(sh1);
	free(dsh1);
	free(ha);
	free(hessa);
	free(xb);
	free(ay);
	free(wa);

	return tau;
}

SEXP _EST_LINEAR_MULTI_WEIGHT(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	// Outcome
	SEXP rBeta, rTheta, rHa, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	_EstLinearCP_Multi_Weight(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, ng, h, maxstep, tol, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_HUBER_MULTI_WEIGHT(SEXP BETA0, SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP G, SEXP DIMs, SEXP PARAMs){
	int i, n, p1, p2, p3, maxstep, smooth, ng;
	double tol, h, qq, tau0;
	n 		= INTEGER(DIMs)[0];
	p1		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	maxstep	= INTEGER(DIMs)[4];
	smooth	= INTEGER(DIMs)[5];
	ng		= INTEGER(DIMs)[6];

	tol 	= REAL(PARAMs)[0];
	h   	= REAL(PARAMs)[1];
	qq 		= REAL(PARAMs)[2];
	tau0 	= REAL(PARAMs)[3];
	// Outcome
	SEXP rBeta, rTheta, rTau, rHa, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p1+ng*p2));
	PROTECT(rTheta 		= allocVector(REALSXP, 	p3));
	PROTECT(rTau 		= allocVector(REALSXP, 	1));
	PROTECT(rHa 		= allocVector(REALSXP, 	ng));
	PROTECT(list_names 	= allocVector(STRSXP, 	4));
	PROTECT(list 		= allocVector(VECSXP, 	4));

	for(i=0; i<p1+ng*p2; i++) REAL(rBeta)[i] = REAL(BETA0)[i];
	for(i=0; i<p3; i++) REAL(rTheta)[i] = REAL(BETA0)[p1+ng*p2+i];
	for(i=0; i<ng;i++) REAL(rHa)[i] = REAL(BETA0)[p1+ng*p2+p3+i];

	REAL(rTau)[0] = _EstHuberCP_Multi_Weight(REAL(rBeta), REAL(rTheta), REAL(rHa), REAL(tX), REAL(X), REAL(Z), REAL(Y), REAL(G), n, p1, p2, p3, ng, h, maxstep, tol, qq, tau0, smooth);

	SET_STRING_ELT(list_names, 	0,	mkChar("beta"));
	SET_STRING_ELT(list_names, 	1,	mkChar("theta"));
	SET_STRING_ELT(list_names, 	2,	mkChar("tau"));
	SET_STRING_ELT(list_names, 	3,	mkChar("ha"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rTheta);
	SET_VECTOR_ELT(list, 		2, 	rTau);
	SET_VECTOR_ELT(list, 		3, 	rHa);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(6);
	return list;
}
