#ifndef WAST_H_INCLUDED
#define WAST_H_INCLUDED

#define GAMMA_M 5
#define MPI 3.1415926
#define MPI1 0.1591549   //1.0/(2*MPI)
#define MPI2 0.3989423   //1.0/sqrt(2*MPI)
#define SQRT2 0.7071068  //1.0/sqrt(2.0)
#define LOG2PI 1.837877  //log(2pi)
#define MEPS 1e-10
#define ITMAX 100
#define EPS 1.0e-8
#define FPMIN 1.0e-30
#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define SGN(a) (((a)<0)?(-1):1)
#define IDEX(a, b) (((a)<(b))?1:0)

void AbyB(double *outVector, double *A, double *v, int n, int p, int q);

void tAbyB(double *outMatrix, const double *A, const double *B, int n, int p, int q);

double exp_approx(double x, int n);

void Standarize(double *y, double *std, int n, int p, double *x, int flag);

void Std(double *std, double *x, int n, int p, int flag);

void sortN(int *ind0, double *x, int n, int dd);

void SortQ(double *s, int l, int r);

double minx(double *x, int n);

double cumsum(double *x, int n);

int LowTriangularInv(double *B, int n, double *A);

void QRDecompN(double *E, double *R, double *x, int n, int p);

void SampleQuantile(double *qr, int m, double *z, int n, double *q);

double SampleQuantile1(double *z, int n, double q);

int MatrixInvSymmetric(double *a,int n);

int CholeskyDecomL(double *L, int n, double *a);

int CholeskyDecomU(double *U, int n, double *a);

void gser(double *gamser, double a, double x, double *gln);

void gcf(double *gammcf, double a, double x, double *gln);

double incgamma(double x, double a);

double incbeta(double x, double a, double b);

void Kernelh(double *x, double *kern, int n, double h0, int type);

void EstLinearR(double *beta, double *residual, const double *x, const double *y, int n, int p);

double EstCatoniR1(double *y, double *residual, double tau, const int n);

void EstCatoniR(double *x, double *y, double *beta, double *residual, double qq, int n, int p, int maxstep, double eps);

double EstHuber(double *x, double *y, double *beta, double *residual, double *weight, double *hess, double qq, int n, int p, int maxstep, double eps, double sigma);

double EstHuberR1(double *y, double *residual, double tau, const int n);

double EstHuberR(double *x, double *y1, double *beta, double *residual, double qq, double tau0, int n, int p, int maxstep, double eps);

#endif // WAST_H_INCLUDED
