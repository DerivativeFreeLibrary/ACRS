#include "acrs_headers.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv){
	int n, m, nc, iexit;
	int maxiter, prnlev, iter, nf;
	double tol, diff_initial, omega;
	double *lb, *ub, *cl, *cu;
	double *xbest, fbest;
	double *xtrial, *c, *cent, *xtemp;
	int *ind, *mask;
	int i, j;
	workingset *w;

	setdim(&n,&nc);

	m  = 25*n;
	// Working set allocation
	w = (workingset *)malloc(m*sizeof(workingset));
	for(i=0;i<m;i++) w[i].x = (double *)malloc(n*sizeof(double));

	// Upper and Lower bound allocation and definition
	lb = (double *)malloc(n*sizeof(double));
	ub = (double *)malloc(n*sizeof(double));

	setbounds(n,lb,ub);

	// Upper and Lower bound of gen.constraints allocation and definition
	cl = (double *)malloc(nc*sizeof(double));
	cu = (double *)malloc(nc*sizeof(double));
	for(i=0;i<nc;i++){ 
		cl[i]=-1.e+24; cu[i]=1.e+24;
	}

	// xbest allocation
	xbest = (double *)malloc(n*sizeof(double));

	// working arrays allocation
	xtrial = (double *)malloc(n*sizeof(double));
	c      = (double *)malloc(nc*sizeof(double));
	cent   = (double *)malloc(n*sizeof(double));
	xtemp  = (double *)malloc(n*sizeof(double));
	ind    = (int    *)malloc((n+1)*sizeof(int));
	mask   = (int    *)malloc(m*sizeof(int));

	iexit = 0;
	prnlev = 0;
	maxiter= 10000;
	tol = 1.e-6;

	printf("RAND_MAX = %d\n",RAND_MAX);
	printf("Now calling acrs ...\n");

	acrs(n, lb, ub, nc, cl, cu, m, w, maxiter, prnlev, tol,
		 xbest, &fbest, &iter, &nf, &iexit, &diff_initial, 
		 &omega, xtrial, c, cent, xtemp, ind, mask);

	printf("acrs exited reporting iexit=%d\n",iexit);

	for(i=0;i<m;i++) free(w[i].x);
	free(w);
	free(lb);
	free(ub);
	free(cl);
	free(cu);
	free(xbest);
	free(xtrial);
	free(c);
	free(cent);
	free(xtemp);
	free(ind);
	free(mask);

	return(0);
}

