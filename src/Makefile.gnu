FC = gfortran
FCFLAGS = -O2 -Wall -fopenmp -g -fbounds-check
BLAS = -lblas
LAPACK = -llapack
include ./Makefile.in
