!*************************************
!	POWELL FUNCTION
!*************************************
subroutine setdim(n)
	!set problem dimension n
	implicit none
	integer		:: n

	n = 4

	return
end subroutine setdim

!*************************************
subroutine setbounds(n,lb,ub)
	!set problem dimension n
	implicit none
	integer		:: n
	real*8		:: lb(n), ub(n)

	lb = - 10.d0
	ub =   10.d0

	return
end subroutine setbounds

!*************************************
SUBROUTINE FUNCT(N,X,F)
	implicit none
	integer		:: n
	real*8		:: x(n), f
	real*8		:: A, B, C, D

	A=X(1)+10.D0*X(2)
	B=X(3)-X(4)
	C=X(2)-2.D0*X(3)
	D=X(1)-X(4)
	F=A**2+5.D0*B**2+C**4+10.D0*D**4

	RETURN
END SUBROUTINE FUNCT


