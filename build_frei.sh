#!/bin/bash

echo
echo "------------------------------------------------------------"

if [ ! -d "build" ]; then
    echo
    echo "Creating build directory"
    mkdir build
    echo
    echo "------------------------------------------------------------"
fi

##################################################
###   Debug Build                              ###
##################################################

echo
echo "Building debug version of Frei..."
echo

PATH_TO_CBLAS_DIR="/usr/include"
PATH_TO_BLAS_LIBS="/usr/lib64"
PATH_TO_LAPACKE_INCLUDE_DIR="/usr/include"
PATH_TO_LIBGFORTRAN="/usr/lib64"
PATH_TO_LAPACK_BINARIES="/usr/lib64"

if test -f "build/frei_dbg-build.log" ; then
    mv build/frei_dbg-build.log build/frei_dbg-build.log.old
fi

chpl -o build/frei_dbg                                     \
     --warnings                                            \
     --warn-unstable                                       \
     -I$PATH_TO_CBLAS_DIR                                  \
     -L$PATH_TO_BLAS_LIBS -lcblas                          \
     -I$PATH_TO_LAPACKE_INCLUDE_DIR                        \
     -L$PATH_TO_LIBGFORTRAN -lgfortran                     \
     -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas \
     --main-module "FREI" src/FREI.chpl                    \
                          src/fr.chpl                      \
                          src/temporal.chpl                \
                          src/output.chpl                  \
                          src/boundary.chpl                \
                          src/frmesh.chpl                  \
                          src/mapping.chpl                 \
                          src/mesh.chpl                    \
                          src/gmesh.chpl                   \
                          src/riemann.chpl                 \
                          src/flux.chpl                    \
                          src/correction.chpl              \
                          src/limiter.chpl                 \
                          src/projection.chpl              \
                          src/quadrature.chpl              \
                          src/interpolation.chpl           \
                          src/polynomials.chpl             \
                          src/sourceterm.chpl              \
                          src/init.chpl                    \
                          src/ringleb.chpl                 \
                          src/config.chpl                  \
                          src/input.chpl                   \
                          src/parameters.chpl              \
                          src/testing.chpl                 \
    2>&1 | tee build/frei_dbg-build.log
if [ $? -eq 0 ]; then
  echo -e "\nSuccess"
else
  echo -e "\nFailed"
fi
echo
echo "------------------------------------------------------------"

##################################################
###   Optimized / Production Build             ###
##################################################

echo
echo "Building Optimized version of Frei..."
echo

PATH_TO_CBLAS_DIR="/usr/include"
PATH_TO_BLAS_LIBS="/usr/lib64"
PATH_TO_LAPACKE_INCLUDE_DIR="/usr/include"
PATH_TO_LIBGFORTRAN="/usr/lib64"
PATH_TO_LAPACK_BINARIES="/usr/lib64"

if test -f "build/frei_opt-build.log" ; then
    mv build/frei_opt-build.log build/frei_opt-build.log.old
fi

chpl -o build/frei_opt                                     \
     --fast                                                \
     -I$PATH_TO_CBLAS_DIR                                  \
     -L$PATH_TO_BLAS_LIBS -lcblas                          \
     -I$PATH_TO_LAPACKE_INCLUDE_DIR                        \
     -L$PATH_TO_LIBGFORTRAN -lgfortran                     \
     -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas \
     --main-module "FREI" src/FREI.chpl                    \
                          src/fr.chpl                      \
                          src/temporal.chpl                \
                          src/output.chpl                  \
                          src/boundary.chpl                \
                          src/frmesh.chpl                  \
                          src/mapping.chpl                 \
                          src/mesh.chpl                    \
                          src/gmesh.chpl                   \
                          src/riemann.chpl                 \
                          src/flux.chpl                    \
                          src/correction.chpl              \
                          src/limiter.chpl                 \
                          src/projection.chpl              \
                          src/quadrature.chpl              \
                          src/interpolation.chpl           \
                          src/polynomials.chpl             \
                          src/sourceterm.chpl              \
                          src/init.chpl                    \
                          src/ringleb.chpl                 \
                          src/config.chpl                  \
                          src/input.chpl                   \
                          src/parameters.chpl              \
                          src/testing.chpl                 \
    2>&1 | tee build/frei_opt-build.log
if [ $? -eq 0 ]; then
  echo -e "\nSuccess"
else
  echo -e "\nFailed"
fi
echo
echo "------------------------------------------------------------"

##################################################
###   Optimized / Production Build             ###
##################################################

echo
echo "Building Optimized Intel MKL version of Frei..."
echo

PATH_TO_CBLAS_DIR="/opt/intel/oneapi/mkl/latest/include"
PATH_TO_BLAS_LIBS="/opt/intel/oneapi/mkl/latest/lib/intel64"
PATH_TO_LAPACKE_INCLUDE_DIR="/opt/intel/oneapi/mkl/latest/include"
PATH_TO_LIBGFORTRAN="/opt/intel/oneapi/mkl/latest/lib/intel64"
PATH_TO_LAPACK_BINARIES="/opt/intel/oneapi/mkl/latest/lib/intel64"

if test -f "frei_opt_mkl-build.log" ; then
    mv frei_opt_mkl-build.log frei_opt_mkl-build.log.old
fi

chpl -o build/frei_opt_mkl                                 \
     --fast                                                \
     --set blasImpl=mkl                                    \
     --set lapackImpl=mkl                                  \
     -I$PATH_TO_CBLAS_DIR                                  \
     -L$PATH_TO_BLAS_LIBS -lblas                           \
     -I$PATH_TO_LAPACKE_INCLUDE_DIR                        \
     -L$PATH_TO_LIBGFORTRAN -lgfortran                     \
     -L$PATH_TO_LAPACK_BINARIES -llapacke -llapack -lcblas \
     --main-module "FREI" src/FREI.chpl                    \
                          src/fr.chpl                      \
                          src/temporal.chpl                \
                          src/output.chpl                  \
                          src/boundary.chpl                \
                          src/frmesh.chpl                  \
                          src/mapping.chpl                 \
                          src/mesh.chpl                    \
                          src/gmesh.chpl                   \
                          src/riemann.chpl                 \
                          src/flux.chpl                    \
                          src/correction.chpl              \
                          src/limiter.chpl                 \
                          src/projection.chpl              \
                          src/quadrature.chpl              \
                          src/interpolation.chpl           \
                          src/polynomials.chpl             \
                          src/sourceterm.chpl              \
                          src/init.chpl                    \
                          src/ringleb.chpl                 \
                          src/config.chpl                  \
                          src/input.chpl                   \
                          src/parameters.chpl              \
                          src/testing.chpl                 \
    2>&1 | tee build/frei_opt-build.log
if [ $? -eq 0 ]; then
  echo -e "\nSuccess"
else
  echo -e "\nFailed"
fi
echo
echo "------------------------------------------------------------"
