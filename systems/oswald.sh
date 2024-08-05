
#!/bin/bash

# System-specific settings for belize.math.uwaterloo.ca

CC=icc
CXX=icpc
LD=mpiicpc
MPICXX=mpiicpc

# System-specific compiler flags
SYSTEM_CFLAGS="-Ofast"
SYSTEM_LDFLAGS=

# Compiler flags for debugging
DEBUG_CFLAGS="-g -DBZ_DEBUG"
DEBUG_LDFLAGS=

# Compiler flags for optimization
OPTIM_CFLAGS="-O3 -fp-model fast=2 -ftz"
OPTIM_LDFLAGS=$OPTIM_CFLAGS

# Compiler flags for extra optimization, such as -ip -ipo on icc
EXTRA_OPTIM_CFLAGS="-ip -ipo"
EXTRA_OPTIM_LDFLAGS=$EXTRA_OPTIM_CFLAGS

# Library names/locations/flags for MPI-compilation.  This will
# probably not be necessary on systems with a working mpicc
# alias
#MPICXX=mpic++
#MPI_CFLAGS=
#MPI_LIB=
#MPI_LIBDIR=
#MPI_INCDIR=

# Library names/locations for LAPACK
LAPACK_LIB="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
LAPACK_LIBDIR=
LAPACK_INCDIR=

# Library locations for blitz; leave blank to use system-installed
# or compiled-with-this-package version
BLITZ_LIBDIR=
BLITZ_INCDIR=

# Library locations for fftw
FFTW_LIBDIR=
FFTW_INCDIR=

# Library locations for UMFPACK
UMF_INCDIR=
UMF_LIBDIR=

# Location/library for BLAS
BLAS_LIB=$LAPACK_LIB
BLAS_LIBDIR=
BLAS_INCDIR=
