#============================================================================================
#    ACRS - Derivative-Free Adaptive Controlled Random Search algorithm for 
#           bound constrained global optimization problems 
#    Copyright (C) 2011  G.Liuzzi, S.Lucidi
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    G. Liuzzi, S. Lucidi, F. Parasiliti, M. Villani. Multiobjective optimization techniques 
#    for the design of induction motors, IEEE Transactions on Magnetics, 39(3): 1261-1264 (2003)
#    DOI: 10.1109/TMAG.2003.810193
#
#    L. Cirio, S. Lucidi, F. Parasiliti, M. Villani. A global optimization approach for 
#    the synchronous motors design by finite element analysis, International Journal of 
#    Applied Electromagnetics and Mechanics, 16(1): 13-27 (2002)
#============================================================================================
import math
sqrt = math.sqrt
import time
import random


#------------------------------------------------------------------------------------
class ACRS:
	def __init__(self,myfunc,lb,ub):
		self.funct = myfunc
		self.lb    = lb
		self.ub    = ub
		self.n	   = len(lb)

	#------------------------------------------------------------------------------------
	def INSERT(self,FVAL,S,i,m,n):

		XTEMP = [None] * n

		for k in range(0,n):
			XTEMP[k] = S[k][i]
		ftemp=FVAL[i]

		for j in range(0,i):
			if ( FVAL[j] > FVAL[i]):
				for l in range(i,j,-1):
					for h in range(0,n): 
						S[h][l] = S[h][l-1]
					FVAL[l] = FVAL[l-1]
				FVAL[j] = ftemp;
				for h in range(0,n):
					S[h][j] = XTEMP[h]
				break 

	def GENERATE_SET_S(self,maxtgen,n,m,nftot,outlev,tol,S,FVAL,errflag,xin):

		W1 = [None]*n
		tbeg = time.clock()
		ST   = [[None for j in range(0,m)] for i in range(0,n)]
		for i in range(0,m):
			FVAL[i] = 1.e+30
		if((errflag==0) or (errflag==-1)): 
			for i in range(0,n):
				for j in range(0,m):
					S[i][j] = None

		it = 0
		for i in range(0,m):
			timegen = time.clock()
			timegen = timegen - tbeg
			if(timegen > maxtgen):
				errflag = 9
				return nftot, errflag

			if((errflag == -2) and (it <= M-1)):
				#use the ith point provided in set S
				for j in range(0,n):
					ST[j][i] = S[j][it]
					W1[j] = S[j][it]
					if ( outlev >= 2 ):
						print(' S(',j,i,')= ' , ST[j][i])
				it += 1
				viol = False
				#---------------------------------
				# check whether the user point
				# satisfies the box constraints
				#---------------------------------
				for j in range(0,n):
					if ( ( W1[j] > self.ub[j] ) or ( W1[j] < self.lb[j] ) ):

						if ( outlev >= 2 ):
							print(' the',i,'th initial point does not satisfy box constraints')

						viol = True

				if(viol):
					for j in range(0,n):
						r = random.random()
						S[j][it]  = ( self.ub[j] - self.lb[j] ) * r + self.lb[j]
						ST[j][i] = S[j][it]
					
						W1[j]  = S[j][it]
						if ( outlev >= 2 ):
							print(' S(',j,i,')= ' , ST[j][i])
			else:
				for j in range(0,n):
					r = random.random()
					S[j][i] = ( self.ub[j] - self.lb[j] ) * r + self.lb[j]
					W1[j]  = S[j][i]
					if ( outlev >= 2 ):
						print(' S(',j,i,')= ' , S[j][i])

	#--------------------------------------------------------------------------
			f = self.funct( W1 )
			nftot += 1

			FVAL[i] = f

			format3020 = '  %6d | %6d | %11.4e | %11.4e | %9.2e'

			if(errflag == -2):
				self.INSERT (FVAL,ST,i,m,n)
			else:
				self.INSERT (FVAL,S,i,m,n)
		
			diff = FVAL[i] - FVAL[0]

			if(outlev >= 0):
				print(format3020 % (0,nftot,FVAL[0],FVAL[i],diff))

	#---------------------------------------------------------------------------
	# If errflag < 0 then we have to add the user point to the working set
	#---------------------------------------------------------------------------
		if(errflag == -1):
			viol = False
			#---------------------------------
			# check whether the user point
			# satisfies the box constraints
			#---------------------------------
			for i in range(0,n):
				if ( ( xin[i] > self.ub[i] ) or ( xin[i] < self.lb[i] ) ):

					if ( outlev >= 2 ):
						print(' the initial point does not satisfy box constraints')

					viol = True

			if(not viol):
				#---------------------------------
				# if the user point is feasible
				# then add it to S
				#---------------------------------
				f = self.funct ( xin )
				nftot += 1

				I_MAX = 0
				for i in range(1,n):
					if FVAL[I_MAX] < FVAL[i]:
						I_MAX = i
				for j in range(0,n):
					S[j][I_MAX] = xin[j]
				FVAL[I_MAX] = f
				self.INSERT (FVAL,S,I_MAX,M,N)
				if(outlev >= 0):
					print(format3020 % (0,nftot,FVAL[0],FVAL[m],0.0))

		F_MAX = max(FVAL)
		F_MIN = min(FVAL)

		while (F_MAX-F_MIN <= tol):
			timegen = time.clock()
			timegen = timegen - tbeg
			if(timegen > maxtgen):
				return nftot, errflag
			for j in range(0,n):
				r = random.random()
				W1[j] = ( self.ub[j] - self.lb[j] ) * r + self.lb[j]

	#--------------------------------------------------------------------------
			f = self.funct ( W1 )
			nftot += 1

			if(f<F_MAX):
				I_MAX = 0
				for i in range(1,n):
					if FVAL[I_MAX] < FVAL[i]:
						I_MAX = i

				for j in range(0,n):
					S[j][I_MAX] = W1[j]
				FVAL[I_MAX] = f
				if(errflag == -2):
					self.INSERT (FVAL,ST,I_MAX,M,N)
				else:
					self.INSERT (FVAL,S,I_MAX,M,N)
				F_MAX = max(FVAL)
				F_MIN = min(FVAL)
				if(outlev >= 0):
					print(format3020 % (0,nftot,FVAL[0],FVAL[m],0.0))

		if(errflag == -2):
			for i in range(0,m):
				for j in range(0,n):
					S[j][i] = ST[j][i]

		return nftot, errflag

	#------------------------------------------------------------------------------------
	# This subroutine implements the base price iteration that:
	# - random selection of N+1 points of S
	# - computing of weigthed centroid
	# - weigthed reflection
	# - updating of S, if necessary
	#------------------------------------------------------------------------------------
	def PRICE_BASE_ITERATION(self,N,M,S,FVAL,nftot,outlev,iexit,omega,diff,diff_initial,proj,collis_count,rej_count,viol_count):
		p         = 2.0
		COLLISION = False
		VIOL      = False
		F_MIN     = FVAL[0]
		F_MAX     = FVAL[M-1]
		IND       = [None] * (N+1)
		FI        = [None] * (N+1)
		cent      = [None] * N
		X_trial   = [None] * N
	#==============================================================
	# Step 2 : Choose at random  n+1  points  over  S and determine
	#          the centroid
	#==============================================================
	# Random choice of n+1 points

		if ( outlev >= 2 ):
			print(' Random choice of %d points among %d' % (n+1, m))

	#--------------------------------------------------------------------------
		rr = random.random()
		r  = (p**rr-1)/(p-1)

		IND[0]       = int(math.floor(M*r))

		for i in range(1,N+1):

			halt = False
			while not halt:
				rr = random.random()
				r = (p**rr-1)/(p-1)

				IND[i]       = int(math.floor(M*r))

				halt = True
				for j in range(0,i):

					if ( IND[i] == IND[j] ):

						if ( not COLLISION ): 
							collis_count = 0
						COLLISION     = True
						collis_count += 1

						if ( collis_count > 1000 * N ):
							iexit = 4
							return nftot, iexit, collis_count, rej_count, viol_count, omega, diff, diff_initial

						halt = False

	#--------------------------------------------------------------------------

		COLLISION = False

	#---------------------------------------------------------------------------
	# Determine the trial point according to the IMPROVED+MODIFIED Price scheme:
	#---------------------------------------------------------------------------
	#
	# - Determine the weighted centroid

		mass_i = 0
		for i in range(1,N+1):
			if FVAL[mass_i] < FVAL[IND[i]]:
				mass_i = IND[i] 

		if ( rej_count >= 10 * M ): 
			omega = omega / 2.0

		phi = omega * ( ( diff * diff ) /  diff_initial )

		if  (( diff_initial >= 1.e+4 ) and ( phi < 1.e-6 ) ):
			diff_initial = 10.0 * omega * diff
			if ( outlev >= 1 ):
				print(' Modified diff_initial determining the centroid ')

		sumc = 0.0
		for j in range(0,N+1):
			sumc += ( 1.0 / ( FVAL[IND[j]] - F_MIN + phi ) )

		for i in range(0,N):
			cent[i] = 0.0
			for j in range(1,N+1):
				cent[i] += S[i][IND[j]] / ( FVAL[IND[j]] - F_MIN + phi )
			cent[i] = cent[i] / sumc

	# - Determine the trial point by a weighted reflection

		f_w = 0.0
		for j in range(1,N+1):
			f_w += FVAL[IND[j]] / ( FVAL[IND[j]] - F_MIN + phi )
		f_w = f_w / sumc

		alpha = 1.0 - ( ( FVAL[mass_i] - f_w ) / ( F_MAX - F_MIN + phi ) )

		for i in range(0,N):
			X_trial[i] = ( 1.0 + alpha ) * cent[i] - alpha * S[i][mass_i]

		if ( outlev >= 2 ):
			for i in range(0,N):
				print(' X_trial(%d) = %f' % (i, X_trial[i]))

	# Check the consistency of the new trial point with the box constraints
		if ( proj ):

			for i in range(0,N):

				if ( ( X_trial[i] > self.ub[i] ) ):

					X_trial[i] = self.ub[i]

					if ( outlev >= 2 ):
						print(' trial point does not satisfy UB constraints: projected ')

				if ( ( X_trial[i] < self.lb[i] ) ):

					X_trial[i] = self.lb[i]

					if ( outlev >= 2 ):
						print(' trial point does not satisfy LB constraints: projected ')

		else:
			for i in range(0,N):
				if ( ( X_trial[i] > self.ub[i] + 1.e-12) or ( X_trial[i] < self.lb[i] - 1.e-12) ):

					if ( outlev >= 2 ):
						print(' the trial point does not satisfy box constraints')

					if ( not VIOL ): 
						viol_count = 0
					VIOL        = True
					viol_count += 1

					if ( viol_count > 1000 * N ):
						iexit = 5
						return nftot, iexit, collis_count, rej_count, viol_count, omega, diff, diff_initial

					return nftot, iexit, collis_count, rej_count, viol_count, omega, diff, diff_initial

		VIOL = False

		f = self.funct(X_trial)
		nftot += 1

		if ( outlev >= 2 ):
			print(' f_trial = %f' % f)

	#======================================================================
	# Step 3-4 : Comparison of f(X_trial) with f_max
	#======================================================================

		if ( f >= F_MAX ):
			if ( outlev >= 1 ):
				print(' trial point rejected ')

			rej_count += 1

		else:

			rej_count = 0

			if ( outlev >= 1 ):
				print(' trial point accepted ')

			for i in range(0,N):
				S[i][M-1] = X_trial[i]
			FVAL[M-1] = f

			self.INSERT (FVAL,S,M-1,M,N)

			diff = FVAL[M-1] - FVAL[0]

		return nftot, iexit, collis_count, rej_count, viol_count, omega, diff, diff_initial

	def optimize(self,S,FVAL,cpu_limit,maxtgen,maxiter,prnlev,tol,iflag,xsol,diff_initial,omega):
	#---------------------------------------------------------------------------------
	#	S(N,M)		: (ON ENTRY)
	#				  if IEXIT > -2 Unspecified.
	#				  if IEXIT = -2, -3 initial working set.
	#				  (ON EXIT) The last working set used by ACRS
	#
	#	FVAL(M)		: (ON ENTRY)
	#				  if IEXIT > -3 Unspecified.
	#				  if IEXIT = -3 then FVAL(i) = f(S(1:n,i)).
	#				  (ON EXIT) The last working set values used by ACRS
	#
	#	CPU_LIMIT	: (ON ENTRY) The maximum allowed cpu time in seconds. Unchanged on exit
	#
	#	MAXTGEN		: (ON ENTRY) The maximum allowed cpu time for initial working
	#	              set generation. MAXTGEN should be <= CPU_LIMIT. If MAXTGEN > CPU_LIMIT then
	#	              CPU_LIMIT is used instead. (ON EXIT) Unchanged
	#
	#	MAXITER		: (ON ENTRY) The maximum allowed number of iterations
	#
	#	PRNLEV		: (ON ENTRY) The printing level. (ON EXIT) Unchanged
	#				  If PRNLEV < -1 then OUTLEV is automatically reset to -1
	#				  If PRNLEV >  2 then OUTLEV is automatically reset to  2
	#				  PRNLEV = -1 -->  NO OUTPUT
	#                 PRNLEV =  0 -->  CONCISE OUTPUT
	#                 PRNLEV = +1 -->  VERBOSE
	#                 PRNLEV = +2 -->  VERY VERBOSE (for debug purpose only)
	#
	#	TOL			: (ON ENTRY) Tolerance in the stopping criterion. (ON EXIT) Unchanged
	#
	#	iflag		: (ON ENTRY)
	#				  0 --> Generate set S from scratch
	#				 -1 --> Generate set S then add point specified by XSOL to it
	#				 -2 --> Do not generate set S. Use the one provided
	#				 -3 --> Do not generate set S. Continue from previous set
	#	iexit		  (ON EXIT) An integer specifying the type of output.
	#				  0 --> Normal termination
	#				  Abnormal or early termination signals
	#				  1 --> N <= 0
	#				  2 --> maximum number of iterations reached
	#				  3 --> exceeded CPU time limit
	#				  4 --> unable to generate n+1 different points
	#				  5 --> exceeded number of constraints violations
	#				  6 --> exceeded max. number of f. evaluations
	#				  7 --> some var has an improper lower bound
	#				  8 --> some var has an improper upper bound
	#				  9 --> constraints are too tight to be easily satisfied
	#				 10 --> provided working set dimension too small
	#
	#	XSOL		: (ON ENTRY)
	#				  if IEXIT ==  0, unspecified
	#				  if IEXIT == -1, specifies the point to be added to the working set
	#				  (ON EXIT) The computed solution.
	#
	#	FOUT		: (ON ENTRY) Unspecified. (ON EXIT) The minimum computed function value.
	#
	#	ITER		: (ON ENTRY) Unspecified. (ON EXIT) The number of iterations to convergence
	#
	#	NFTOT		: (ON ENTRY) Unspecified. (ON EXIT) The number of objective function evaluations to convergence.
	#
	#---------------------------------------------------------------------------------

		N      = self.n
		M      = len(S[0])
		iexit  = iflag
		tbegin = time.clock()

		if N <= 0:
			iexit = 1
		

		if (iexit <= 0) and (M <= N+1):
			iexit = 10

		if iexit <= 0:
			if( prnlev < -1 ): 
				outlev = -1
			elif( prnlev >  2 ): 
				outlev =  2
			else:
				outlev = prnlev

			for I in range(0,N):
				if(self.lb[I] <= -1.e+6):
					iexit = 7
					break
				if(self.ub[I] >= 1.e+6):
					iexit = 8
					break

		if(iexit <= 0):

	#========================================================
	# Initializations of state variables
	#========================================================

			proj       = True
			if(iexit > -3): 
				omega      = float(N)

			ITER       = 0
			nftot      = 0

			viol       = False
			collision  = False
			rej_count  = 0
			viol_count = 0
			collis_count = 0

	#========================================================
	# Step 0 : Determine the initial set S = {x_1, ... ,x_m}
	#          and evaluate f at each point
	#========================================================
			if ( outlev >= 1 ):
				print('++++++++++++ Determine the initial set S +++++++++ ')

			if(iexit > -3):
				nftot, errflag = self.GENERATE_SET_S(maxtgen,N,M,nftot,outlev,tol,S,FVAL,iexit,xsol)
			
			timeset = time.clock()
			timeset = timeset - tbegin

			if(iexit <= 0):
				diff		 = FVAL[M-1] - FVAL[0]
				if(iexit > -3): 
					diff_initial = diff

				iexit = 0

				if ( diff_initial > 1.e+10 ):
					diff_initial = 1.e+10
					if ( outlev >= 1 ):
						print(' Modified diff_initial ')

	#========================================================
	# Main loop of the algorithm
	#========================================================
	 			format3000 = '   Iter      NF        FMIN          FMAX        DIFF'
				format3010 = '-----------------------------------------------------------'
				format3020 = '  %6d | %6d | %11.4e | %11.4e | %9.2e'
				format2000 = '%10d'
				format2003 = '  Total time = %10.2f seconds'
				format2004 = '  Final function value = %22.14e'
				format2005 = '  Minimum point : '
				format2006 = '  tolerance used in the (global) stopping criterion = %10.2e'
				format2007 = '   X _%2d   = %22.14e'
				format2008 = ' ********************* FINAL RESULTS ********************'
				format2009 = '  number of initial points used (m) = %6d'
				format2010 = '  dimension of the problem (n) = %3d'
				format2014 = '  Setup time = %10.2f seconds\n  Solve time = %10.2f seconds\n  Total time = %10.2f seconds'

				while True:

					#---------------------
					# do some printing
					#---------------------
					if(outlev >= 0):
						if (ITER % 30) ==0:
							print('')
							print(format3000)
							print(format3010)
						print(format3020 % (ITER,nftot,FVAL[0],FVAL[M-1],diff))

						if( outlev >= 1 ):
							print(' ')
							print('****** Iteration %d ****** ' % ITER)

					#------------------------------------------
					# Stopping criterion: stops either because
					#		tolerance    reached -- iexit = 0
					#		maxiter      reached -- iexit = 2
					#		max cpu time reached -- iexit = 3
					#------------------------------------------
					timesol = time.clock()
					timesol = timesol - tbegin
					if ( diff < tol ):
						iexit = 0
						break
					if ( ITER >= maxiter ):
						iexit = 2
						break
					if ( timesol >= cpu_limit ):
						iexit = 3
						break

					#-----------------------------------
					# perform the basic Price iteration
					#-----------------------------------
					nftot, iexit, collis_count, rej_count, viol_count, omega, diff, diff_initial = self.PRICE_BASE_ITERATION(N,M,S,FVAL,nftot,outlev,iexit,omega,diff,diff_initial,proj,collis_count,rej_count,viol_count)
					#-----------------------------------
					# if there were errors then STOP
					#-----------------------------------
					if(iexit > 0):
						break

					ITER += 1

	#----------------------------------------------------------------
	# Normal exit
	#----------------------------------------------------------------

			timetot = time.clock()
			timetot = timetot - tbegin
			timesol = timetot - timeset

			if(outlev >= -1):

				if iexit == 0:

	# Normal exit. Printing the final results

					print(format2008)
					print(' Modified Improved Price Algorithm ')
					print(format2010 % N)
					print(format2009 % M)
					print(format2006 % tol)

					print(format2004 % FVAL[0])

					print(' # iterations : %d' % ITER)
					print(' # function evaluations : %d' % nftot)
					print(format2014 % (timeset,timesol,timetot))
					print(format2005)
					for i in range(0,N):
						print(format2007 % (i , S[i][0]))
						xsol[i] = S[i][0]
					FOUT      = FVAL[0]
					print('')

					print(format2008)

					return xsol,FOUT,ITER,nftot,iexit,diff_initial,omega

				else:
	# Error exit
					if iexit == 1:
						print(' *** ERROR:  n  must be > 0        *** ')
						print(' *** Please, correct and resubmit. *** ')

					if iexit == 2:
						print(' *** WARNING: maximum number of iterations ***   (%d)' % ITER)
						print(format2003 % timetot)

					if iexit == 3:
						print(' *** WARNING: exceeded CPU time limit ***')
						print(format2003 % timetot)

					if iexit == 4:
						print(' *** ERROR: unable to generate n+1 different points ***')

					if iexit == 5:
						print(' *** ERROR: exceeded number of constraints violations ***')

					if iexit == 6:
						print(' *** WARNING: exceeded max number f.evaluations ***')

					if iexit == 7:
						print(' *** ERROR: some var has an improper lower bound ***')
						print(' *** Please, correct and resubmit.               ***')

					if iexit == 8:
						print(' *** ERROR: some var has an improper upper bound ***')
						print(' *** Please, correct and resubmit.               ***')

					if iexit == 9:
						print(' *** ERROR: constraints are too tight to be easily satisfied ***')
						print(' *** Please, correct and resubmit.                           ***')

					if iexit == 10:
						print(' *** ERROR: working set dimension too small ***')

				return xsol,FOUT,ITER,nftot,iexit,diff_initial,omega
