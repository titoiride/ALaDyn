#!/bin/bash 

module purge

module load compilers/gcc-4.9.0
module load compilers/openmpi-1.8.1_gcc-4.9.0_with_cuda6.5
module load boost_1_56_0_gcc4_9_0

mkdir -p build_gcc ; cd build_gcc ; cmake .. -DCMAKE_PREFIX_PATH=/shared/software/project/aladyn/fftw -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_C_COMPILER=/shared/software/compilers/gcc-4.9.0/bin/gcc -DCMAKE_CXX_COMPILER=/shared/software/compilers/gcc-4.9.0/bin/g++ ; cmake --build . --target install ; cd ..
#mkdir -p build_gcc ; cd build_gcc ; cmake .. -DCMAKE_PREFIX_PATH=/shared/software/project/aladyn/fftw -DOTHER_LINK_DIR=/shared/software/project/aladyn/fftw/lib -DOTHER_INCLUDE_DIR=/shared/software/project/aladyn/fftw/include/ -DCMAKE_LINKER=/shared/software/compilers/gcc-4.9.2/bin/gfortran  -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_C_COMPILER=/shared/software/compilers/gcc-4.9.2/bin/gcc -DCMAKE_CXX_COMPILER=/shared/software/compilers/gcc-4.9.2/bin/g++ -DCMAKE_Fortran_COMPILER=/shared/software/compilers/gcc-4.9.2/bin/gfortran ; cmake --build . --target install ; cd ..

