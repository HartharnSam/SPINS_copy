#!/bin/bash

# Helper script to generate an appropriate system.mk file, based on a 
# configuration script in this directory

# Read in the appropriate system script.  If none is specified on the
# command line, guess based on the hostname

if [ $# -gt 0 ]; then
   if [ -f $1 ]; then
      echo Reading system-specific variables from $1
      source $1
   elif [ -f $1.sh ]; then
      echo Reading system-specific variables from $1.sh
      source $1.sh
   else 
      echo Neither $1 nor $1.sh are found!
      exit 1
   fi
else
   echo Guessing that system-specific variables are in `hostname -s`.sh
   source `hostname -s`.sh || (echo "But they're not"; exit 1)
fi

# Now, generate ../src/system.mk 

echo "Generating ../src/system.mk"
cat > ../src/system.mk << EOF
# System-specific settings for `hostname -s`

# Compiler settings
CC:=$CC
CXX:=$CXX
LD:=$LD
MPICXX:=$MPICXX

# Compiler flags
SYSTEM_CFLAGS:=$SYSTEM_CFLAGS
SYSTEM_LDFLAGS:=$SYSTEM_LDFLAGS
SYSTEM_CXXFLAGS:=$SYSTEM_CXXFLAGS

# Debugging flags
DEBUG_CFLAGS:=$DEBUG_CFLAGS
DEBUG_LDFLAGSL=$DEBUG_LDFLAGS

# Optimization flags
OPTIM_CFLAGS:=$OPTIM_CFLAGS
OPTIM_LDFLAGS:=$OPTIM_LDFLAGS

# Extra-slow optimization flags
EXTRA_OPTIM_CFLAGS:=$EXTRA_OPTIM_CFLAGS
EXTRA_OPTIM_LDFLAGS:=$EXTRA_OPTIM_LDFLAGS

# Library names
# MPI
MPI_CFLAGS:=$MPI_CFLAGS
MPI_LIB:=$MPI_LIB
MPI_LIBDIR:=$MPI_LIBDIR
MPI_INCDIR:=$MPI_INCDIR

# LAPACK
LAPACK_LIB:=$LAPACK_LIB
LAPACK_LIBDIR:=$LAPACK_LIBDIR
LAPACK_INCDIR:=$LAPACK_INCDIR

# Blitz
BLITZ_LIBDIR:=$BLITZ_LIBDIR
BLITZ_INCDIR:=$BLITZ_INCDIR

# fftw
FFTW_LIBDIR:=$FFTW_LIBDIR
FFTW_INCDIR:=$FFTW_INCDIR

# UMFPack
UMF_LIBDIR:=$UMF_LIBDIR
UMF_INCDIR:=$UMF_INCDIR

# BLAS
BLAS_LIB:=$BLAS_LIB
BLAS_LIBDIR:=$BLAS_LIBDIR
BLAS_INCDIR:=$BLAS_INCDIR

EOF


echo "Sanity check - building trivial MPI executable"
rm -f mpitest.o mpitest.x mpitest.cpp
cat > mpitest.cpp << EOF
#include <mpi.h>
int main(int argc, char ** argv) {
   MPI_Init(&argc, &argv);
   MPI_Finalize();
   return 0;
}
EOF
echo "Executing: $MPICXX $SYSTEM_CFLAGS $SYSTEM_CXXFLAGS $MPI_CFLAGS $MPI_INCDIR -c -o mpitest.o ./mpitest.cpp"
$MPICXX $SYSTEM_CFLAGS $SYSTEM_CXXFLAGS $MPI_CFLAGS $MPI_INCDIR -c -o mpitest.o ./mpitest.cpp || echo "WARNING: MPI compilation failed!" && echo "MPI compilation success"
echo "Executing: $LD $SYSTEM_LDFLAGS $MPI_LIBDIR -o mpitest.x mpitest.o $MPI_LIB -lstdc++" 
$LD $SYSTEM_LDFLAGS $MPI_LIBDIR -o mpitest.x mpitest.o $MPI_LIB -lstdc++ || echo "WARNING: MPI linking failed!" && echo "MPI linking success"
rm -f mpitest.o mpitest.x mpitest.cpp
