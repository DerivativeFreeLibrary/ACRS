//============================================================================================
//    ACRS - Derivative-Free Adaptive Controlled Random Search algorithm for 
//           bound constrained global optimization problems 
//    Copyright (C) 2011  G.Liuzzi, S.Lucidi
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//    G. Liuzzi, S. Lucidi, F. Parasiliti, M. Villani. Multiobjective optimization techniques 
//    for the design of induction motors, IEEE Transactions on Magnetics, 39(3): 1261-1264 (2003)
//    DOI: 10.1109/TMAG.2003.810193
//
//    L. Cirio, S. Lucidi, F. Parasiliti, M. Villani. A global optimization approach for 
//    the synchronous motors design by finite element analysis, International Journal of 
//    Applied Electromagnetics and Mechanics, 16(1): 13-27 (2002)
//============================================================================================

#include "acrs_headers.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void acrs(int n, double *lb, double *ub, int nc, double *cl, double *cu,
		  int m, workingset *w, int maxiter, int prnlev, double tol,
		  double *xbest, double *fbest, int *iter, int *nftot, int *iexit,
		  double *diff_initial, double *omega, double *xtrial, double *c,
		  double *cent, double *xtemp, int *ind, int *mask){
/*---------------------------------------------------------------------------------
 
   PRICE ALGORITHM
   first implemented by G. Di Pillo, S. Lucidi and
   successively modified by G. Liuzzi
   November 14, 2005
 
   Reference: P. Brachetti, M. De Felice Ciccoli, G. Di Pillo, S. Lucidi
              "A new version of the Price's algorithm for global optimization"
              Journal of Global Optimization, 10  pp.165-184, 1997.
  --------------------------------------------------------------------------------- */

	float cpu_limit, maxtgen;

/*---------------------------------------------------------------------------------
!	N			: (ON ENTRY) The number of variables. Unchanged on exit
!				  it must be >= 1
!
!	LB(N)		: (ON ENTRY) The lower bounds on the variables. Unchanged on exit
!				  Every variable must have a proper lower bound
!
!	UB(N)		: (ON ENTRY) The upper bounds on the variables. Unchanged on exit
!				  Every variable must have a proper upper bound
!
!	NC			: (ON ENTRY) The number of general constraints. Unchanged on exit
!				  The gen. constraints are thought to be in the form  cl <= g <= cu.
!				  If there are no gen. constraints then set NC <= 0.
!
!	CL(NC)		: (ON ENTRY) The lower bounds of the gen. constraints. Unchanged on exit
!
!	CU(NC)		: (ON ENTRY) The upper bounds of the gen. constraints. Unchanged on exit
!
!	M			: (ON ENTRY) Dimension of the working set. (ON EXIT) Unchanged
!				  It should be a bit more than N+1. The suggested value is MAX(50, 25*N).
!				  If WKDIMEN <= N+1 then the routines terminates with an error
!
!	W			: (ON ENTRY)
!				  if IEXIT > -2 Unspecified.
!				  if IEXIT = -2, -3 initial working set.
!				  (ON EXIT) The last working set used by ACRS
!
!	MAXITER		: (ON ENTRY) The maximum allowed number of iterations
!
!	PRNLEV		: (ON ENTRY) The printing level. (ON EXIT) Unchanged
!				  If PRNLEV < -1 then OUTLEV is automatically reset to -1
!				  If PRNLEV >  2 then OUTLEV is automatically reset to  2
!				  PRNLEV = -1 -->  NO OUTPUT
!                 PRNLEV =  0 -->  CONCISE OUTPUT
!                 PRNLEV = +1 -->  VERBOSE
!                 PRNLEV = +2 -->  VERY VERBOSE (for debug purpose only)
!
!	TOL			: (ON ENTRY) Tolerance in the stopping criterion. (ON EXIT) Unchanged
!
!	XBEST		: (ON ENTRY)
!				  if IEXIT ==  0, unspecified
!				  if IEXIT == -1, specifies the point to be added to the working set
!				  (ON EXIT) The computed solution.
!
!	FBEST		: (ON ENTRY) Unspecified. (ON EXIT) The minimum computed function value.
!
!	ITER		: (ON ENTRY) Unspecified. (ON EXIT) The number of iterations to convergence
!
!	NFTOT		: (ON ENTRY) Unspecified. (ON EXIT) The number of objective function evaluations to convergence.
!
!	IEXIT		: (ON ENTRY)
!				  0 --> Generate set S from scratch
!				 -1 --> Generate set S then add point specified by XSOL to it
!				 -2 --> Do not generate set S. Use the one provided
!				 -3 --> Do not generate set S. Continue from previous set
!				  (ON EXIT) An integer specifying the type of output.
!				  0 --> Normal termination
!				  Abnormal or early termination signals
!				  1 --> N <= 0
!				  2 --> maximum number of iterations reached
!				  3 --> exceeded CPU time limit
!				  4 --> unable to generate n+1 different points
!				  5 --> exceeded number of constraints violations
!				  6 --> exceeded max. number of f. evaluations
!				  7 --> some var has an improper lower bound
!				  8 --> some var has an improper upper bound
!				  9 --> constraints are too tight to be easily satisfied
!				 10 --> provided working set dimension too small
!--------------------------------------------------------------------------------- */
	int i, j, outlev, viol_count, collis_count, rej_count;
	int viol, collision, proj;
	double diff;
	float timeset, timesol, timetot;
	 

	maxtgen   = 0.0;
	cpu_limit = 0.0;

	//CALL CPU_TIME(time)

	if( n <=  0) {
		*iexit = 1;
		goto _900;
	}

	if( m <= n+1 ) {
		*iexit = 10;
		goto _900;
	}

	outlev = prnlev;
	if( prnlev < -1 ) outlev = -1;
	if( prnlev >  2 ) outlev =  2;

	for(i=0;i<n;i++){
		if(lb[i] <= -1.e+16){
			*iexit = 7;
			break;
		}
		if(ub[i] >= 1.e+16){
			*iexit = 8;
			break;
		}
	}

	if ( *iexit > 0 ) goto _1000;

	//srand(37598);

	proj = 1;
	if ( *iexit > -3 ) (*omega) = (double)n;

	*iter       = 0;
	*nftot      = 0;

	viol       = 0;
	collision  = 0;
	rej_count  = 0;
	viol_count = 0;

	//CALL CPU_TIME(timeset)
	//timeset = timeset - time

	if ( outlev >= 1 ) printf("++++++++++++ determine the initial set s +++++++++++++\n");

	if ( *iexit > -3 ){
		generate_set_w(maxtgen,n,m,nftot,outlev,lb,ub,tol,w,nc,c,cl,cu,iexit,xbest);

		if( *iexit > 0 ) goto _1000;
	}

	
	/*
	FILE *fid;
	fid = fopen("fort.1","r");
	//fid = fopen("cort.1","w");
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			fscanf(fid,"%lf",&w[i].x[j]);
			//fprintf(fid,"%lf\n",w[i].x[j]);
		}
		fscanf(fid,"%lf",&w[i].f);
		//fprintf(fid,"%lf\n",w[i].f);
	}
	fclose(fid);
	*/

	//CALL CPU_TIME(timeset)
	//timeset = timeset - time

	diff = w[m-1].f - w[0].f;
	if ( *iexit > -3 ) *diff_initial = diff;

	*iexit = 0;

	if ( *diff_initial > 1.e+10 ){
		*diff_initial = 1.e+10;
		if ( outlev >= 1 ) printf("modified diff_initial\n");
	}

	while (1){
	
		//---------------------
		// do some printing
		//---------------------
		if(outlev >= 0){
			//if(mod(iter,30)==0){
			if ( *iter - (*iter / 30)*30 == 0){
				printf("\n");
				printf("   Iter      NF        FMIN          FMAX        DIFF\n");
				printf("-----------------------------------------------------------\n");
			}
			printf("  %d | %d | %11.4e | %11.4e | %9.2e\n",
					  *iter,*nftot,w[0].f,w[m-1].f,diff);
		}

		if ( outlev >= 1 ){
			printf("\n");
			printf("****** iteration %d ****** \n",*iter);
		}

		/*------------------------------------------*
		 * Stopping criterion: stops either because *
		 *		tolerance    reached -- iexit = 0   *
		 *		maxiter      reached -- iexit = 2   *
		 *		max cpu time reached -- iexit = 3   *
		 *------------------------------------------*/
		
		//CALL CPU_TIME(timesol)
		//timesol = timesol - time
		
		if ( diff < tol      ) break;
		if ( *iter >= maxiter ){
			*iexit = 2;
			break;
		}
		if ( timesol >= cpu_limit ){
			//*iexit = 3;
			//break;
		}

	/*int i, outlev, viol_count, collis_count, rej_count;
	int viol, collision, proj;
	double diff, ftemp;
	float timeset, timesol, timetot, time; */

		price_base_iteration(n,nc,m,w,iexit,nftot,outlev,lb,ub,cl,cu, 
					omega,&diff,diff_initial,proj,&collis_count,&rej_count,&viol_count,
					cent, xtrial, c, ind, mask);

		if ( *iexit > 0 ) break;

		*iter = *iter + 1;

	}

/*----------------------------------------------------------------*
 * Normal exit													  *
 *----------------------------------------------------------------*/

_900: 1;	

	//CALL CPU_TIME(timesol)
	//timesol = timesol - time
	//timetot = timeset + timesol

	if ( outlev >= 0 ){

		switch (*iexit){
		case 0:
			printf(" ********************* FINAL RESULTS ********************\n");
			printf(" modified improved price algorithm \n");
			printf(" dimension of the problem (n) = %d\n", n);
			printf(" number of initial points used (m) = %d\n" , m);
			printf(" tolerance used in the (global) stopping criterion = %f\n", tol );

			printf("\n");
			printf(" Final function value = %14.9e\n", w[0].f);
			printf(" %f %f\n",w[m-1].f,w[0].f);

			printf("\n");
			printf(" # iterations : %d\n" , *iter);
			printf(" # function evaluations : %d\n" , *nftot);
			printf(" Total time = %f seconds\n" ,timetot );

			printf("\n");
			printf(" Minimum point : \n");
			for(i=0;i<n;i++)  printf("  X _%d = %13.8e\n",i,w[0].x[i]);

			printf("\n");
			printf(" ********************************************************\n");
			break;
		case 1:
			printf(" *** error:  n  must be > 0        *** \n");
			printf(" *** please, correct and resubmit. *** \n");
			break;

		case 2:
			printf(" *** warning: maximum number of iterations ***   (%d)\n",*iter);
			printf(" Total time = %f seconds\n" ,timetot );
			break;

		case 3:
			printf(" *** warning: exceeded cpu time limit ***\n");
			printf(" Total time = %f seconds\n" ,timetot );
			break;

		case 4:
			printf(" *** error: unable to generate n+1 different points ***\n");
			break;

		case 5:
			printf(" *** error: exceeded number of constraints violations ***\n");
			break;

		case 6:
			printf(" *** warning: exceeded max number f.evaluations ***\n");
			break;

		case 10:
			printf(" *** error: working set dimension too small ***\n");
			break;

		}
    }

	for(i=0;i<n;i++)  xbest[i] = w[0].x[i];
	*fbest = w[0].f;

	if ( outlev >= 0 ) printf("Setup time = %f seconds\nSolve time = %f seconds\nTotal time = %f seconds\n",timeset,timesol,timetot);

_1000: 1;

	switch (*iexit){
	case 7:
		if ( outlev >= 0 ){
			printf(" *** error: some var has an improper lower bound ***\n");
			printf(" *** please, correct and resubmit.               ***\n");
		}
		break;
	case 8:
		if ( outlev >= 0 ){
			printf(" *** error: some var has an improper upper bound ***\n");
			printf(" *** please, correct and resubmit.               ***\n");
		}
		break;
	case 9:
		if ( outlev >= 0 ){
			printf(" *** error: constraints are too tight to be easily satisfied ***\n");
			printf(" *** please, correct and resubmit.                           ***\n");
		}
		break;
	}

}

int maxloc(workingset *w, int m, int *mask){
	int i, ind;
	ind = -1;

	for(i=0;i<m;i++){
		if ( (ind == -1) && (mask[i]) ){ 
			ind = i;
		}else{
			if ( ind != -1 ){
				if ( (mask[i]) && (w[i].f>w[ind].f) ) ind = i;
			}
		}
	}
	return (ind);
}

void price_base_iteration(int n, int nc, int m, workingset *w, int *iexit, int *nftot,
						  int outlev, double *lb, double *ub, double *cl, double *cu,
						  double *omega, double *diff, double *diff_initial, int proj,
						  int *collis_count, int *rej_count, int *viol_count,
						  double *cent, double *xtrial, double *c, int *ind, int *mask){

	double p=2.0;
	double r, f_min, f_max, phi, sum, f, f_w, alpha;
	int i, j, mass_i;
	int collision, viol;

	collision = 0;
	viol      = 0;
	f_min     = w[0].f;
	f_max     = w[m-1].f;
	
	for(i=0;i<m;i++) mask[i] = 0;

/*==============================================================*
 * Step 2 : Choose at random  n+1  points  over  S and determine*
 *          the centroid                                        *
 *==============================================================*/

// Random choice of n+1 points

	if ( outlev >= 2 )printf(" random choice of %d points among %d\n",n+1,m);

	r = rand();
	r = (double)r/(double)RAND_MAX;
	//random_num_(&r);
	r = (pow(p,r)-1)/(p-1);

	ind[0]       = (int)(floor((double)m*r));
	if(ind[0] > m-1) ind[0] = m-1;
	mask[ind[0]] = 1;

	for(i=1;i<n+1;i++){

_300:	        r = rand();
		r = r/(double)RAND_MAX;
		//random_num_(&r);
		r = (pow(p,r)-1)/(p-1);

		ind[i]       = (int)(floor((double)m*r));
		if(ind[i] > m-1) ind[i] = m-1;
		mask[ind[i]] = 1;
		
		for(j=0;j<i;j++){
			if ( ind[i] == ind[j] ){
				if ( !collision ) *collis_count = 0;
				collision = 1;
				*collis_count = *collis_count + 1;
				if ( *collis_count > 1000 * n ){
					*iexit = 4;
					return;
				}
				goto _300;
			}

		}

	}

	collision = 0;


/*-----------------------------------------------------------------------------*
 * Determine the trial point according to the IMPROVED+MODIFIED Price scheme:  *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * - Determine the weighted centroid                                           */

	mass_i = maxloc(w,m,mask);

	if ( *rej_count >= 10 * m ) *omega = *omega / 2.0;

	phi = (*omega) * ( ( (*diff) * (*diff) ) /  (*diff_initial) );
	if ( (*diff_initial >= 1.e+4) && (phi < 1.e-6) ){
		*diff_initial = 10.0 * (*omega) * (*diff);
		if ( outlev >= 1 ) printf(" modified diff_initial determining the centroid \n");
	}

	sum = 0.0;
	for(j=0;j<n+1;j++) sum = sum + ( 1.0 / ( w[ind[j]].f - f_min + phi ) );

	for(i=0;i<n;i++){
		cent[i] = 0.0;
		for(j=0;j<n+1;j++) cent[i] = cent[i] + w[ind[j]].x[i] / ( w[ind[j]].f - f_min + phi );
		cent[i] = cent[i] / sum;
	}

/*                                                                             *
 * - Determine the trial point by a weighted reflection                        *
 *                                                                             */

	f_w = 0.0;
	for(j=0;j<n+1;j++) f_w = f_w + w[ind[j]].f / ( w[ind[j]].f - f_min + phi );
	f_w = f_w / sum;

	alpha = 1.0 - ( ( w[mass_i].f - f_w ) / ( f_max - f_min + phi ) );

	for(i=0;i<n;i++) xtrial[i] = ( 1.0 + alpha ) * cent[i] - alpha * w[mass_i].x[i];

	for(j=0;j<n+1;j++) mask[ind[j]] = 0;

	if ( outlev >= 2 ){
		for(i=0;i<n;i++) printf(" x_trial(%d) = %f\n",i,xtrial[i]);
	}

// Check the consistency of the new trial point with the box constraints

	if ( proj ){

		for(i=0;i<n;i++){
			if ( xtrial[i] > ub[i] ){
				xtrial[i] = ub[i];
				if ( outlev >= 2 ) printf(" trial point does not satisfy ub constraints: projected \n");
			}
			if ( xtrial[i] < lb[i] ){
				xtrial[i] = lb[i];
				if ( outlev >= 2 ) printf(" trial point does not satisfy lb constraints: projected \n");
			}
		}

	}else{
		for(i=0;i<n;i++){
			if ( (xtrial[i] > ub[i]) || (xtrial[i] < lb[i]) ){

				if ( outlev >= 2 ){
					printf(" the trial point does not satisfy box constraints\n");
					printf("%f %f %f\n",lb[i],xtrial[i],ub[i]);
				}

				if ( !viol ) *viol_count = 0;
				viol       = 1;
				*viol_count = *viol_count + 1;

				if ( *viol_count > 1000 * n ){
					*iexit = 5;
					return;
				}

				return;

			}
		}

	}

//--------------------------------------------------------------------------


	if( nc > 0 ) fconstr(n,nc,xtrial,c);
	viol = 0;

/*                                                                             *
 * - Check the general constraints assumed in the form  CONSTR(nc) <= 0        *
 *                                                                             */

	for(i=0;i<nc;i++){
		if ( (c[i] > cu[i]) || (c[i] < cl[i]) ){
			if ( outlev >= 2 ) printf(" the trial point does not satisfy general constraints \n");
			return;
		}
	}

	f = funct(n,xtrial);
	*nftot = *nftot + 1;

	if ( outlev >= 2 ) printf(" f_trial = %f\n", f);

/*======================================================================*
 * Step 3-4 : Comparison of f(X_trial) with f_max                       *
 *======================================================================*/

	if ( f >= f_max ){
		if ( outlev >= 1 ) printf(" trial point rejected \n");
		rej_count = rej_count + 1;
	}else{
		rej_count = 0;

		if ( outlev >= 1 ){
			printf(" trial point accepted \n");
			printf(" diff = %f\n", *diff);
		}

		for(i=0;i<n;i++) w[m-1].x[i] = xtrial[i];
		w[m-1].f = f;
		insert (w,m-1,m,n);
		*diff = w[m-1].f - w[0].f;
	}
	return;
}

void generate_set_w(int maxtgen, int n, int m, int *nftot, int outlev, 
		    double *lb, double *ub, double tol, workingset *w,
		    int nc, double *c, double *cl, double *cu, int *iexit, double *xin){
	int i, j, i_max, viol_count, viol_tot, it;
	workingset *wt;
	double *w1;
	float timegen;
	double r, f, f_max, f_min;
	int viol;

	w1 = (double *)malloc(n*sizeof(double));

	wt = (workingset *)malloc(m*sizeof(workingset));
	for(i=0;i<m;i++) wt[i].x = (double *)malloc(n*sizeof(double));

	viol_tot = 0;
	//call cpu_time(time)

	it = 1;
	for(i=0;i<m;i++){
		viol_count = 0;
_90:	1;
		//call cpu_time(timegen)
		//timegen = timegen - time
		//if ( (timegen > maxtgen) || (viol_tot > 100000000) ){
		if (viol_tot > 100000000) {
		//	*iexit = 9;
		//	return;
		}

		if ( (*iexit == -2) && (it <= m) ){
			//use the ith point provided in set s
			for(j=0;j<n;j++){
				wt[i].x[j] = w[it].x[j];
				w1[j]      = w[it].x[j];
				if ( outlev >= 2 ){
					printf(" s(%d,%d)= %f\n",j,i, wt[i].x[j]);
				}
			}
			it = it + 1;
			viol = 0;
			/*---------------------------------*
			 * check whether the user point    *
			 * satisfies the box constraints   *
			 *---------------------------------*/
			for(j=0;j<n;j++){
				if ( ( w1[j] > ub[j] ) || ( w1[j] < lb[j] ) ){
					if ( outlev >= 2 ) printf(" the %dth initial point does not satisfy the box constraints\n",i);
					viol = 1;
				}
			}
			if ( viol ){
				for(j=0;j<n;j++){
					r = rand();
					r = r/(double)RAND_MAX;
					//random_num_(&r);
					w [i].x[j] = ( ub[j] - lb[j] ) * r + lb[j];
					w1[j]  = w [i].x[j];
					if ( outlev >= 2 ) printf(" s(%d,%d)= %f\n",j,i,w [i].x[j]);
				}
			}
		}else{
			for(j=0;j<n;j++){
				r = rand();
				r = r/(double)RAND_MAX;
				//random_num_(&r);
				w [i].x[j] = ( ub[j] - lb[j] ) * r + lb[j];
				w1[j]  = w [i].x[j];
				if ( outlev >= 2 ) printf(" -- s(%d,%d)= %f\n",j,i,w [i].x[j]);
			}
        }

		if( nc > 0 ) fconstr(n,nc,w1,c);

		for(j=0;j<nc;j++){
			if ( ( c[j] > cu[j] ) || ( c[j] < cl[j] ) ){
			   if ( outlev >= 2 )printf(" the initial point does not satisfy general constraints\n");
			   viol_count = viol_count + 1;
			   viol_tot   = viol_tot   + 1;
			   goto _90;
			}
		}

		f = funct(n,w1);

		//printf("done computing f value !\n");

		*nftot = *nftot + 1;

		w[i].f = f;

		if ( *iexit == -2 ){
			insert (wt,i,m,n);
		}else{
			insert (w, i,m,n);
		}

        if ( outlev >= 0 ) printf("  %d | %d | %11.4e | %11.4e | %9.2e\n",
								  0,*nftot,w[0].f,w[i].f,0.0);

	}

	if ( *iexit == -1 ){
		viol = 0;
		/*---------------------------------*
		 * check whether the user point    *
		 * satisfies the box constraints   *
		 *---------------------------------*/
		for(i=0;i<n;i++){
			if ( (xin[i] > ub[i]) || (xin[i] < lb[i]) ){
				if ( outlev >= 2 ) printf(" the initial point does not satisfy box constraints\n");
				viol = 1;
			}
		}
		if(!viol){
			/*---------------------------------*
			 * check whether the user point    *
			 * satisfies the gen. constraints  *
			 *---------------------------------*/
			if ( nc > 0 ) fconstr(n,nc,xin,c);
			for(j=0;j<nc;j++){
				if ( (c[j] > cu[j]) || (c[j] < cl[j]) ){
				   if ( outlev >= 2 )printf(" the initial point does not satisfy general constraints\n");
				   viol_count = viol_count + 1;
				   viol_tot   = viol_tot   + 1;
				   viol       = 1;
				}
			}

			/*---------------------------------*
			 * then add it to s                *
			 *---------------------------------*/
			f = funct (n,xin);
			*nftot = *nftot + 1;

			i_max = 0;
			for(i=1;i<m;i++) i_max = (w[i_max].f > w[i].f ? i_max : i);

			for(j=0;j<n;j++) w[i_max].x[j] = xin[j];
			w[i_max].f = f;
			
			insert (w,i_max,m,n);
            if ( outlev >= 0 ) printf("  %d | %d | %11.4e | %11.4e | %9.2e\n",
								  0,*nftot,w[0].f,w[m-1].f,0.0);
		}
	}

	f_max = w[0].f;
	f_min = w[0].f;
	for(i=1;i<m;i++){
		f_max = (f_max > w[i].f ? f_max : w[i].f);
		f_min = (f_min < w[i].f ? f_min : w[i].f);
	}

	while ( f_max-f_min <= 1.0e+0*tol){
		viol_count = 0;
_100:	1;
		//call cpu_time(timegen)
		//timegen = timegen - time
		//if( (timegen > maxtgen) || (viol_tot > 100000000) ) return;
		if( (viol_tot > 100000000) ) return;
		for(j=0;j<n;j++){
			r = rand();
			r = r/(double)RAND_MAX;
			//random_num_(&r);
			w1[j] = ( ub[j] - lb[j] ) * r + lb[j];
		}

		if( nc > 0 ) fconstr(n,nc,w1,c);

		for(j=0;j<nc;j++){
			if ( (c[j] > cu[j]) || (c[j] < cl[j]) ){
			   if ( outlev >= 2 )printf(" the initial point does not satisfy general constraints\n");
			   viol_count = viol_count + 1;
			   viol_tot   = viol_tot   + 1;
			   goto _100;
			}
		}

		f = funct (n,w1);
		*nftot = *nftot + 1;

		if( f<f_max ){
			i_max = 0;
			for(i=1;i<m;i++) i_max = (w[i_max].f > w[i].f ? i_max : i);

			for(j=0;j<n;j++) w[i_max].x[j] = w1[j];
			w[i_max].f = f;
			
			if ( *iexit == -2 ){
				insert (wt,i_max,m,n);
			}else{
				insert (w, i_max,m,n);
			}
			f_max = w[0].f;
			f_min = w[0].f;
			for(i=1;i<m;i++){
				f_max = (f_max > w[i].f ? f_max : w[i].f);
				f_min = (f_min < w[i].f ? f_min : w[i].f);
			}
            if ( outlev >= 0 ) printf("  %d | %d | %11.4e | %11.4e | %9.2e\n",
								  0,*nftot,w[0].f,w[m-1].f,0.0);
		}
	}

	if ( *iexit == -2 ){
		for(i=0;i<m;i++){
			for(j=0;j<n;j++) w[i].x[j] = wt[i].x[j];
			w[i].f = wt[i].f;
		}
	}

	free(w1);
	for(i=0;i<m;i++) free(wt[i].x);
	free(wt);

	return;
}

void insert(workingset *w, int i, int m, int n){
	double ftemp;
	double *xtemp;
	int j,l,h;

	xtemp = (double *)malloc(n*sizeof(double));

	for(j=0;j<n;j++) xtemp[j] = w[i].x[j];
	ftemp = w[i].f;

	for(j=0;j<=i-1;j++){
		if ( w[j].f > w[i].f){
			for(l=i;l>=j+1;l--){
				for(h=0;h<n;h++) w[l].x[h] = w[l-1].x[h];
				w[l].f = w[l-1].f;
			}
		   w[j].f = ftemp;
		   for(h=0;h<n;h++) w[j].x[h] = xtemp[h];
		   goto _TERMINATE;
		}
	}

_TERMINATE:

	free(xtemp);
	return;
}

