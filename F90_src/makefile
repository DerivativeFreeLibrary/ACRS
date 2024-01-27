
FC = gfortran
CC = gcc
RM = rm -f

#FFLAGS = -O3
FFLAGS = -g


OBJS = main.o ACRS.o problem.o

all: acrs_f

acrs_f:  $(OBJS)
	$(FC) -o acrs_f $(OBJS)
	mv acrs_f ..

.SUFFIXES : .f90 .o

.f90.o: $* ; $(FC) $(FFLAGS) -c $*.f90

clean: 
	$(RM) *.o
	$(RM) *.mod
	$(RM) acrs_f

