import problem
lb, ub = problem.setbounds()

from ACRS import *
acrs = ACRS(problem.funct,lb,ub)
n = len(lb)
m = 25*n

S    = [[None for j in range(0,m)] for i in range(0,n)]
FVAL = [None] * m
xsol = [None] * n

tol       = 1.e-6
maxiter   = 10000
cpu_limit = 3600.0
maxtgen   = 3600.0
diff_init = 0.0
omega     = 0.0
prnlev    = 0
iflag     = 0

#ACRS is a probabilistic code, hence we randomly set the random seed
#random.seed()


acrs.optimize(S,FVAL,cpu_limit,maxtgen,maxiter,prnlev,tol,iflag,xsol,diff_init,omega)
