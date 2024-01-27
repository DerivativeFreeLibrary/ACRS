#include <math.h>

void setdim(int *n, int* nc){
	*n  = 4;
	*nc = 0;
	return;
};

void setbounds(int n, double* lb, double* ub){
	int i;
	for(i=0;i < n; i++){
		lb[i] = -10.0;
		ub[i] =  10.0;
	}
	return;
};

double funct(int n, double* x){
	int i;
	double f, A, B, C, D;

	A=x[0]+10.0*x[1];
	B=x[2]-x[3];
	C=x[1]-2.0*x[2];
	D=x[0]-x[3];
	
	f=pow(A,2) + 5.0*pow(B,2) + pow(C,4) + 10.0*pow(D,4);	

	return(f);
}

void fconstr(int n, int nc, double* x, double* c){
	return;
}

