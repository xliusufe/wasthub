#include <stdio.h>  	// required for exit
#include <stdlib.h> 	// required for malloc(), free();
#include <string.h> 	// required for memcpy()
#include <float.h>  	// required for DBL_EPSILON
#include <math.h>
#include "_WASTHUB_HEAD.h"

void AbyB(double *outVector, double *A, double *v, int n, int p, int q){
	int i,j,k;
	double tmp;
	for (i=0;i<n;i++){
		for(k=0;k<q;k++){
			tmp = 0;
			for(j=0;j<p;j++)
				tmp += A[j*n + i]*v[k*p + j];
			outVector[k*n+i] = tmp;
		}
	}
}

void tAbyB(double *outMatrix, const double *A, const double *B, int n, int p, int q){
    int i,j,k;
    double temp;
	for (i = 0; i<p; i++){
		for (k = 0; k<q; k++){
			temp = 0.0;
			for (j = 0; j < n; j++)
				temp += A[i + j*p] * B[k + j*q];
			outMatrix[i*q + k] = temp;
		}
	}
}

double exp_approx(double x, int n){
	x = 1.0 + x/256;
	for (int i = 0; i < n; i++){
		x *= x;
	}

	return x;
}

void Standarize(double *y, double *std, int n, int p, double *x, int flag)
{
	int i, j;
	double s, s1;
	for(j=0;j<p;j++)
	{
		s = 0; s1 = 0;
		for(i=0;i<n;i++) s += x[j*n+i];
		s = s/n;
		for(i=0;i<n;i++) s1  += x[j*n+i]*x[j*n+i];
		s1 = s1/n - s*s;
		if(flag)
			std[j] = sqrt(s1);
		else
			std[j] = sqrt(n*s1/(n-1));
		for(i=0;i<n;i++) y[j*n+i] = (x[j*n+i]-s)/std[j];
	}
}

void Std(double *std, double *x, int n, int p, int flag){
	int i, j;
	double s, s1;
	for(j=0;j<p;j++)
	{
		s = 0; s1 = 0;
		for(i=0;i<n;i++) s += x[j*n+i];
		s = s/n;
		for(i=0;i<n;i++) s1  += x[j*n+i]*x[j*n+i];
		s1 = s1/n - s*s;
		if(flag)
			std[j] = sqrt(s1);
		else
			std[j] = sqrt(n*s1/(n-1));
	}
}

void sortN(int *ind0, double *x, int n, int dd){
	int i, j, MaxInd, d, *ind;
	double tmp;
	ind = (int*)malloc(sizeof(int)*n);
	for(i=0;i<n;i++) ind[i] = i;

	d = (dd==n?dd-1:dd);
	for(i=0;i<d;i++)
	{
		tmp = x[0]; MaxInd = ind[0];
		for(j=1;j<n-i;j++)
		{
			if(x[j]<tmp)
			{
				x[j-1] = x[j];
				x[j] = tmp;
				ind[j-1] = ind[j];
				ind[j] = MaxInd;
			}
			else
			{
				tmp = x[j];
				MaxInd = ind[j];
			}
		}
	}
	for(j=0;j<dd;j++) ind0[j] = ind[n-j-1];
	free(ind);
}

void SortQ(double *s, int l, int r)
{
    int i, j;
	double x;
    if (l < r)
    {
        i = l;
        j = r;
        x = s[i];
        while (i < j)
        {
            while(i < j && s[j] > x) j--;
			if(i < j) s[i++] = s[j];
            while(i < j && s[i] < x) i++;
			if(i < j) s[j--] = s[i];
        }
        s[i] = x;
        SortQ(s, l, i-1);
        SortQ(s, i+1, r);
    }
}

double minx(double *x, int n){
	int i;
	double mx;
	mx = x[0];
	for(i=1; i<n; i++){
		if(x[i]<mx){
			mx = x[i];
		}
	}
	return mx;
}

double cumsum(double *x, int n){
	int i;
	double mx;
	mx = x[0];
	for(i=1; i<n; i++){
		mx += x[i];
	}
	return mx;
}

int LowTriangularInv0(double *B, int n, double *A){
	// Input:
	// A is a lower triangular matrix
	//
	// Output:
	// B = inv(A)
	//
	int i,j,k;
	for(i=0;i<n;i++)
		if(fabs(A[i*n+i])<EPS)	return(0);
	for(i=0;i<n;i++)	B[i*n+i] = 1;
	for(j=1;j<n;j++)
		for(i=0;i<j;i++)	B[j*n+i] = 0;

	for(i=n-1;i>=0;i--)//rows
	{
		if(fabs(A[i*n+i]-1)>EPS)
			for(j=i;j<n;j++)
				B[j*n+i] = B[j*n+i]/A[i*n+i];
		if(i>0)
		{
			for(j=i;j<n;j++)// columns
				for(k=0;k<i;k++)// rows
					B[j*n+k] = B[j*n+k] - A[i*n+k]*B[j*n+i];
		}
	}
	return(1);
}

int LowTriangularInv(double *B, int n, double *A){
	// Input:
	// A is a lower triangular matrix
	//
	// Output:
	// B = inv(A)
	//
	int i,j,k;
	for(i=0;i<n;i++)
		if(fabs(A[i*n+i])<EPS)	return(0);
	for(i=0;i<n;i++)	B[i*n+i] = 1;
	for(j=1;j<n;j++)
		for(i=0;i<j;i++)	B[j*n+i] = 0;

	for(i=n-1;i>=0;i--)//rows
	{
		if(fabs(A[i*n+i]-1)>EPS)
			for(j=i;j<n;j++)
				B[j*n+i] = B[j*n+i]/A[i*n+i];
		if(i>0)
		{
			for(j=i;j<n;j++)// columns
				for(k=0;k<i;k++)// rows
					B[j*n+k] = B[j*n+k] - A[i*n+k]*B[j*n+i];
		}
	}
	return(1);
}

void QRDecompN(double *E, double *R, double *x, int n, int p){
	// Input:
	// X is a p*n matrix
	//
	// Output:
	// R is a p*p lower triangular matrix
	// E is a p*n matrix satisfying E*t(E) = I_p
	//
	double *Z, *znorm;
	double  tmp, tmp1;
	int i,j, k;

	Z = (double*)malloc(sizeof(double)*n*p);
	znorm = (double*)malloc(sizeof(double)*p);

	// calculate the first column
	tmp = 0;
	for(i=0;i<n;i++){
		Z[i] = x[i];
		tmp += Z[i]*Z[i];
	}
	znorm[0] = sqrt(tmp);
	tmp = 0;
	for(i=0;i<n;i++){
		E[i] = x[i]/znorm[0];
		tmp += E[i]*x[i];
	}
	R[0] = tmp;

	//iteration from j=1...p
	for(j=1;j<p;j++){
		for(k=0;k<j;k++){
			tmp=0;	for(i=0;i<n;i++) tmp += E[k*n+i]*x[j*n+i];
			R[j*p+k] = tmp;
		}
		tmp1 = 0;
		for(i=0;i<n;i++){
			tmp = 0; for(k=0;k<j;k++) tmp += R[j*p+k]*E[k*n+i];
			Z[j*n+i] = x[j*n+i] - tmp;
			tmp1 += pow(Z[j*n+i],2);
		}
		znorm[j] = sqrt(tmp1);
		tmp1 = 0;
		for(i=0;i<n;i++) E[j*n+i] = Z[j*n+i]/znorm[j];
		for(i=0;i<n;i++) tmp1 += E[j*n+i]*x[j*n+i];
		R[j*p+j] = tmp1;
	}
	free(Z); free(znorm);
}

void SampleQuantile(double *qr, int m, double *z, int n, double *q)
{
	double *zs=(double*)malloc(sizeof(double)*n);
	int i, ind;
	for(i=0;i<n;i++) zs[i] = z[i];
	SortQ(zs, 0, n-1);
	for(i=0;i<m;i++)
	{
		ind = floor(q[i]*n);
		if (ind!=n*q[i])
			qr[i] = zs[ind];
		else
			qr[i] = (zs[ind-1] + zs[ind])/2;
	}
	free(zs);
}

double SampleQuantile1(double *z, int n, double q){
	int i, ind;
	double qr;
	double *zs=(double*)malloc(sizeof(double)*n);
	for(i=0;i<n;i++) zs[i] = z[i];
	SortQ(zs, 0, n-1);

	ind = floor(q*n);
	if (ind!=n*q){
		qr = zs[ind];
	}
	else{
		qr = (zs[ind-1] + zs[ind])/2;
	}

	free(zs);
	return qr;
}

int MatrixInvSymmetric(double *a,int n){
	int i,j,k,m;
    double w,g,*b;
    b = (double*)malloc(n*sizeof(double));

    for (k=0; k<=n-1; k++){
        w=a[0];
        if (fabs(w)+1.0==1.0){
            free(b); return(-2);
        }
        m=n-k-1;
        for (i=1; i<=n-1; i++){
            g=a[i*n]; b[i]=g/w;
            if (i<=m) b[i]=-b[i];
            for (j=1; j<=i; j++)
                a[(i-1)*n+j-1]=a[i*n+j]+g*b[j];
        }
        a[n*n-1]=1.0/w;
        for (i=1; i<=n-1; i++)
            a[(n-1)*n+i-1]=b[i];
    }
    for (i=0; i<=n-2; i++)
        for (j=i+1; j<=n-1; j++)
            a[i*n+j]=a[j*n+i];
    free(b);
    return(2);
}

int CholeskyDecomL(double *L, int n, double *a){
	// L is lowertriangle matrix satisfying a = L*L'
	int i,j,k;
	double sum, *p;
	p = (double*)malloc(sizeof(double)*n);

	for(i=0; i<n; i++)
	{
		for(j=i; j<n;j++)
		{
			for(sum=a[j*n+i], k=i-1; k>=0; k--)
			{
				sum-=a[k*n+i]*a[k*n+j];
			}
			if (i==j)
			{
				if(sum<=0.0)
					return(-1);
				p[i]=sqrt(sum);
			}
			else
				a[i*n+j]=sum/p[i];
		}
	}
	for(i=0;i<n;i++){
		for(j=i+1; j<n;j++){
			L[i*n+j] = a[i*n+j];
		}
		L[i*n+i] = p[i];
		for(j=0; j<i;j++){
			L[i*n+j] = 0.0;
		}
	}
	free(p);
	return(1);
}

int CholeskyDecomU(double *U, int n, double *a){
	// U is uppertriangle matrix satisfying a = U'*U
	int i,j,k;
	double sum, *p;
	p = (double*)malloc(sizeof(double)*n);

	for(i=0; i<n; i++)
	{
		for(j=i; j<n;j++)
		{
			for(sum=a[j*n+i], k=i-1; k>=0; k--)
			{
				sum-=a[k*n+i]*a[k*n+j];
			}
			if (i==j)
			{
				if(sum<=0.0)
					return(-1);
				p[i]=sqrt(sum);
			}
			else
				a[i*n+j]=sum/p[i];
		}
	}
	for(i=0;i<n;i++){
		for(j=i+1; j<n;j++){
			U[j*n+i] = a[i*n+j];
		}
		U[i*n+i] = p[i];
		for(j=0; j<i;j++){
			U[j*n+i] = 0.0;
		}
	}
	free(p);
	return(1);

	// CholeskyDecomU(H, p, I0);
	// for (s = 0; s < p; s++){
	// 	for (t = 0; t < p; t++){
	// 		tmp = 0.0;
	// 		for(j=0; j< t+1; j++){
	// 			tmp += H[s*p+j]*H[t*p+j];
	// 		}
	// 		printf("%f  ", tmp);
	// 	}
	// 	printf("\n");
	// }
}

void gser(double *gamser, double a, double x, double *gln){
	int n;
	double sum,del,ap;
	*gln = lgamma(a);
	if (x <= 0.0) {
		if (x < 0.0)
			fprintf(stderr, "x less than 0 in routine gser");
		*gamser=0.0;
		return;
	}
	else {
		ap = a;
		del = sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser = sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		fprintf(stderr, "a too large, ITMAX too small in routine gser");
		return;
	}
}

void gcf(double *gammcf, double a, double x, double *gln){
	int i;
	double an,b,c,d,del,h;
	*gln = lgamma(a);
	b = x+1.0-a;
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < FPMIN)	d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN)	c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS)
			break;
	}
	if (i > ITMAX) fprintf(stderr, "a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

double incgamma(double x, double a){
	double gamser,gammcf, gln;
	if (x < 0.0 || a <= 0.0) fprintf(stderr, "Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser, a, x, &gln);
		return gamser;
	}
	else {
		gcf(&gammcf, a, x, &gln);
		return 1.0-gammcf;
	}
}

double incbeta(double x, double a, double b) {
	if (x < 0.0 || x > 1.0) return 1.0/0.0;

	/*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
	if (x > (a+1.0)/(a+b+2.0)) {
		return (1.0-incbeta(1.0-x,b,a)); /*Use the fact that beta is symmetrical.*/
	}

	/*Find the first part before the continued fraction.*/
	const double lbeta_ab = lgamma(a)+lgamma(b)-lgamma(a+b);
	const double front = exp(log(x)*a+log(1.0-x)*b-lbeta_ab) / a;

	/*Use Lentz's algorithm to evaluate the continued fraction.*/
	double f = 1.0, c = 1.0, d = 0.0;

	int i, m;
	for (i = 0; i <= ITMAX; ++i) {
		m = i/2;

		double numerator;
		if (i == 0) {
			numerator = 1.0; /*First numerator is 1.0.*/
		} else if (i % 2 == 0) {
			numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
		} else {
			numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
		}

		/*Do an iteration of Lentz's algorithm.*/
		d = 1.0 + numerator * d;
		if (fabs(d) < FPMIN) d = FPMIN;
		d = 1.0 / d;

		c = 1.0 + numerator / c;
		if (fabs(c) < FPMIN) c = FPMIN;

		const double cd = c*d;
		f *= cd;

		/*Check for stop.*/
		if (fabs(1.0-cd) < EPS) {
			return front * (f-1.0);
		}
	}

    return 1.0/0.0; /*Needed more loops, did not converge.*/
}

void Kernelh(double *x, double *kern, int n, double h0, int type){
	// type = 1: Gaussian kernel
	// type = 2; Epanecknikov kernel
	// type = 3: uniform kernel
	int i;
	double tmp, h;
	Std( &h, x, n, 1, 0);
	h /= h0;
	for(i=0; i<n; i++){
		tmp = x[i]*h;
		tmp *= tmp;
		if(type==1){
			kern[i] = exp(-0.5*tmp)*MPI2*h;
		}
		else if(type==2){
			kern[i] = 0.75*(1-tmp)*(tmp<1?1:0)*h;
		}
		else{
			kern[i] = (tmp>0.5?0:1)*(tmp<-0.5?0:1)*h;
		}
	}
}

void EstLinearR(double *beta, double *residual, const double *x, const double *y, int n, int p){
	int i,j,k;
	double tmp, *hess, *xy;
	xy 		= (double*)malloc(sizeof(double)*p);
	hess	= (double*)malloc(sizeof(double)*p*p);

	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i];
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
	AbyB(beta, hess, xy, p, p, 1);

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0; j < p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - tmp;
	}

	free(hess);
	free(xy);
}

double EstCatoniR1(double *y, double *residual, double tau, const int n){
	int i;
	double beta=0.0;
	for(i=0;i<n;i++){
		if(fabs(y[i])<tau){
			beta += y[i];
		}
		else{
			beta += tau*SGN(y[i]);
		}		
	}
	beta /= n;

	for(i=0;i<n;i++){
		residual[i] = y[i] - beta;
	}

	return beta;
}

void EstCatoniR(double *x, double *y, double *beta, double *residual, double qq, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta0, *qy, *dpsi, *hess;
	double tmp, bnorm, yk, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	yk = 0.0;
	for(i=0;i<n;i++)	yk += y[i];
	yk /= n;
	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(i=0;i<n;i++)	residual[i]	= y[i] - yk;

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = residual[i]/qq;
			if(fabs(tmp)>1.414214){
				yk = SGN(tmp)*0.942809*qq;
				wk = 0.0;
			}
			else{
				yk = qq*tmp*( 1.0 - tmp*tmp/6);
				wk = 1.0 - 0.5*tmp*tmp;
			}
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*wk;
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
			tmp = beta[j];
			bnorm += tmp*tmp;
		}

		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)	beta0[j] += beta[j];
			for(i=0;i<n;i++){
				tmp = 0.0;
				for(j=0;j<p;j++) tmp += x[j*n+i]*beta0[j];
				residual[i] = y[i] - tmp;
			}
		}
	}

	for(j=0; j < p; j++){
		beta[j] = beta0[j];
	}
	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++) tmp += x[j*n+i]*beta[j];
		tmp = (y[i] - tmp)/qq;
		if(fabs(tmp)>1.414214){
			residual[i] = SGN(tmp)*0.942809;
		}
		else{
			residual[i] = qq*tmp*( 1.0 - tmp*tmp/6);
		}
	}
	

	free(beta0);
	free(dpsi);
	free(qy);
	free(hess);
}

double EstHuber(double *x, double *y, double *beta, double *residual, double *weight, double *hess, double qq, int n, int p, int maxstep, double eps, double sigma){
	int i,j,k, step=0;
	double *beta0, *qy, *dpsi;
	double tmp, bnorm, yk, wk, loglokl=0.0;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	yk = 0.0;
	for(i=0;i<n;i++)	yk += y[i];
	yk /= n;
	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(i=0;i<n;i++)	residual[i]	= y[i] - yk;

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = residual[i]/qq;
			if(fabs(tmp)>1.414214){
				yk = SGN(tmp)*0.942809;
				wk = 0.0;
			}
			else{
				yk = qq*tmp*( 1.0 - tmp*tmp/6);
				wk = 1.0 - 0.5*tmp*tmp;
			}
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*wk;
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
			tmp = beta[j];
			bnorm += tmp*tmp;
		}

		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)	beta0[j] += beta[j];
			for(i=0;i<n;i++){
				tmp = 0.0;
				for(j=0;j<p;j++) tmp += x[j*n+i]*beta0[j];
				residual[i] = y[i] - tmp;
				loglokl += 0.5*residual[i]*residual[i]*IDEX(residual[i], qq) + (qq*fabs(residual[i])-0.5*qq*qq)*IDEX(qq, residual[i]);
			}
		}
	}

	for(j=0; j < p; j++){
		beta[j] = beta0[j];
	}
	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++) tmp += x[j*n+i]*beta[j];
		tmp = (y[i] - tmp)/qq;
		if(fabs(tmp)>1.414214){
			residual[i] = SGN(tmp)*0.942809;
			weight[i] = 0.0;
		}
		else{
			residual[i] = qq*tmp*( 1.0 - tmp*tmp/6);
			weight[i] = 1.0 - 0.5*tmp*tmp;
		}
	}

	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i]*weight[i];
			}
			hess[j*p+k] = tmp;
		}
	}
	MatrixInvSymmetric(hess,p);

	free(beta0);
	free(dpsi);
	free(qy);
	return -2.0*loglokl;
}

double EstHuberR1(double *y, double *residual, double tau, const int n){
	int i;
	double beta=0.0;
	for(i=0;i<n;i++){
		if(fabs(y[i])<tau){
			beta += y[i];
		}
		else{
			beta += tau*SGN(y[i]);
		}
	}
	beta /= n;

	for(i=0;i<n;i++){
		residual[i] = y[i] - beta;
	}

	return beta;
}

double EstHuberR(double *x, double *y1, double *beta, double *residual, double qq, double tau0, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta0, *hess, *qy, *dpsi, *y, *med;
	double tmp, bnorm, yk, wk, meany, qr;

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
		tmp = 0.0;
		for(j=0;j<p;j++) tmp += x[j*n+i]*beta[j];
		tmp = y[i] - tmp;
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
	return qq;
}

double _Omega(double *OMEGA, double *x, double *z, double *resid0, int n, int p2, int p3, int typewgt){
	int i,j,s;
	double rho, xij, omega, sd, Tn0=0.0, *stdx;

	stdx 	= (double*)malloc(sizeof(double)*n*p3);


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

	if(typewgt==1){
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
	}
	else{
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
				omega = 0.5 - acos(rho)*MPI1;
				OMEGA[i*n+j]	= omega*xij;
				Tn0 	+= OMEGA[i*n+j]*resid0[i]*resid0[j];
			}
		}
	}


	free(stdx);
	return Tn0;
}

double _OmegaSingle(double *OMEGA, double *x, double *z, double *resid0, int n, int p2, double shape1, double shape2, int isBeta){
	int i,j,s;
	double xij, omega, Tn0=0.0, *ty;

	ty  = (double*)malloc(sizeof(double)*n);

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


	free(ty);
	return Tn0;
}

double _Omega_approx(double *OMEGA, double *x, double *z, double *mu0, double *zk, double *resid0, int n, int p2, int p3, int N0){
	int i,j,s, count;
	double tmp, tmp1, tmp2, tmp3, rho, xij, omega, sd, Tn0=0.0, *stdx, *zmu;

	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	zmu 	= (double*)malloc(sizeof(double)*n);


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
			tmp = omega*xij;
			OMEGA[i*n+j]	= tmp;
			Tn0 	+= tmp*resid0[i]*resid0[j];
		}
	}


	free(stdx);
	free(zmu);
	return Tn0;
}
