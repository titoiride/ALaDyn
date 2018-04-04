#!/bin/bash

module purge

module load compilers/gcc-7.1.0
module load compilers/openmpi-2.1.1_gcc-7.1.0
module load boost_1_64_0_gcc7_1_0

export FC=/shared/software/compilers/gcc-7.1.0/bin/gfortran
export CC=/shared/software/compilers/gcc-7.1.0/bin/gcc
export CXX=/shared/software/compilers/gcc-7.1.0/bin/g++

mkdir -p build ; cd build
cmake .. -DFFTW_ROOT_DIR=/shared/software/project/aladyn/fftw -DBoost_NO_BOOST_CMAKE=ON
cmake --build . --target install
cd ..
