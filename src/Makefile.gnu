FC = gfortran
FCFLAGS = -O2 -Wall -fopenmp -g -fbounds-check -fdefault-integer-8
BLAS = -lblas
LAPACK = -llapack
include ./Makefile.in
