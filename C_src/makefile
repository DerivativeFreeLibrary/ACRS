
FC = gfortran
CC = gcc
RM = rm -f

#FFLAGS = -O3
FFLAGS = -g


OBJS = main.o ACRS.o problem.o

all: acrs_c

acrs_c:  $(OBJS)
	$(CC) -lm -o acrs_c $(OBJS)
	mv acrs_c ..

.SUFFIXES : .c   .o
	    
.c.o:   $* ; $(CC) $(FFLAGS) -c $*.c

clean: 
	$(RM) *.o
	$(RM) *.mod
	$(RM) acrs_c

