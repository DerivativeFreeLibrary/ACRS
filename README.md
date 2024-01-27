ACRS - A Controlled Random Search algorithm
-----------------------------------------------------------
 How to use the derivative-free optimizer ACRS for 
 bound constrained global optimization
-----------------------------------------------------------
 The package provides a FORTRAN90, C, PYTHON and MATLAB version of the code.

0. Gunzip and untar the archive in a folder on your computer by
   issuing in a directory of your choice (ex. curdir) the command

   $> tar -xvf ACRS.tar.gz

   FORTRAN VERSION OF THE CODE:

   Edit file curdir/F90_src/problem.f90 to define your own objective function.
   In particular, modify the subroutines 
   setdim    : which sets problem dimension
   setbounds : which sets upper and lower bounds on the variables
   funct     : which defines the objective function

   C VERSION OF THE CODE:

   Edit file curdir/C_src/problem.c to define your own objective function.
   In particular, modify the subroutines 
   setdim    : which sets problem dimension
   setbounds : which sets upper and lower bounds on the variables
   funct     : which defines the objective function

   PYTHON VERSION OF THE CODE:
   
   Edit file curdir/PYTHON/problem.py to define your own objective function.
   In particular, modify functions 
   setbounds : which sets upper and lower bounds on the variables
   funct     : which defines the objective function

1. At command prompt in curdir execute 

     $> make
 
   which will create the executables 'acrs_c' (using C sources) and
   'acrs_f' (using FORTRAN90 sources)

2. execute

     $> ./acrs_f
     $> ./acrs_c
	 $> python PYTHON/main.py

   MATLAB VERSION OF THE CODE

   In the folder curdir/Matlab there is a matlab version (acrs.m) of the
   code. File main.m contains an example usage of ACRS for Matlab to 
   minimize the funcion defined in powell.m
