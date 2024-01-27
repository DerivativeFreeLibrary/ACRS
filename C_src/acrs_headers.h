//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <time.h>

struct type_workingset {
	double		*x;
	double		f;
};

typedef struct type_workingset	workingset;

void setdim(int *n, int* nc);
void setbounds(int n, double* lb, double* ub);
double funct(int n, double* x);
void fconstr(int n, int nc, double* x, double* c);

void acrs(int n, double *lb, double *ub, int nc, double *cl, double *cu,
		  int m, workingset *w, int maxiter, int prnlev, double tol,
		  double *xbest, double *fbest, int *iter, int *nftot, int *iexit,
		  double *diff_initial, double *omega, double *xtrial, double *c,
		  double *cent, double *xtemp, int *ind, int *mask);

void generate_set_w(int maxtgen, int n, int m, int *nftot, int outlev, 
					double *lb, double *ub, double tol, workingset *w,
					int nc, double *c, double *cl, double *cu, int *iexit, double *xin);

void price_base_iteration(int n, int nc, int m, workingset *w, int *iexit, int *nftot,
						  int outlev, double *lb, double *ub, double *cl, double *cu,
						  double *omega, double *diff, double *diff_initial, int proj,
						  int *collis_count, int *rej_count, int *viol_count,
						  double *cent, double *xtrial, double *c, int *ind, int *mask);

void insert(workingset *w, int i, int m, int n);

int maxloc(workingset *w, int m, int *mask);

//extern void random_num_(double* r);
