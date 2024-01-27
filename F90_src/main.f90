program acrs_main
  implicit none

  INTEGER :: n, m, maxiter, prnlev, iout, iter, nftot, iexit
  DOUBLE PRECISION :: tol, fstoptol, fout, diff_initial, omega
  REAL :: cpu_limit, maxtgen
  DOUBLE PRECISION, ALLOCATABLE :: lb(:), ub(:), xout(:)
  DOUBLE PRECISION, ALLOCATABLE :: S(:,:), FVAL(:), xsol(:)

  EXTERNAL :: funct, funct_acrs
  
  call setdim(n)

  m = 25*n
  allocate(lb(n),ub(n),xout(n))
  allocate(S(n,m), FVAL(m), xsol(n))

  call setbounds(n,lb,ub)

  tol=  1.d-6
  fstoptol = 0.d0
  maxiter = 10000

  cpu_limit = 3600.0
  maxtgen   = 3600.0
  prnlev    = 0
  iout      = 6
  iexit     = 0

  call ACRS(N,LB,UB,M,S,FVAL,FUNCT_ACRS,cpu_limit,maxtgen,maxiter,prnlev,tol,IOUT, 	&
  	    XSOL,FOUT,iter,nftot,iexit,diff_initial,omega)

end program acrs_main

subroutine funct_acrs(x,n,f)
  implicit none
  integer :: n
  double precision :: x(n), f

  call funct(n,x,f)

  return
end subroutine funct_acrs
