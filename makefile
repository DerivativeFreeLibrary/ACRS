
FC = gfortran
CC = gcc
RM = rm -f

#FFLAGS = -O3
FFLAGS = -g

OBJS = main.o ACRS.o problem.o

all: fver cver

fver:
	(cd ./F90_src ; make ; cd ..)

cver:
	(cd ./C_src ; make ; cd ..)

clean: 
	(cd ./C_src ; make clean ; cd ..)
	(cd ./F90_src ; make clean ; cd ..)
	$(RM) acrs_f
	$(RM) acrs_c
