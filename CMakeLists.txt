cmake_minimum_required (VERSION 2.8.10)
project (ALaDyn)
enable_language (Fortran)
enable_language (CXX)

set (CMAKE_CXX_STANDARD 11)
set (MAJOR_VERSION 1)
set (MINOR_VERSION 0)
set (FIX_RELEASE 3)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the root directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set (CMAKE_C_FLAGS_RELEASE "-O3")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_C_FLAGS_DEBUG "-Og -g")
set (CMAKE_CXX_FLAGS_DEBUG "-Og -g")

if (Fortran_COMPILER_NAME MATCHES "gfortran*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-fdefault-real-8 -funroll-all-loops -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fdefault-real-8 -Wall -Wextra -fbacktrace -fbounds-check -Og -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-real-size 64 -axcore-avx-i,SSE4.2 -ipo -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-real-size 64 -O0 -g -real-size 64 -g -check all -fpe0 -warn -traceback -debug extended")
else ()
  message ("CMAKE_Fortran_COMPILER found but unsupported")
  message ("Let's try using it...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ()


# build executables
add_executable (ALaDyn "src/ALaDyn.F90")

# Require Boost, MPI and FFTW3 for this project
find_package(Boost COMPONENTS system filesystem REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(MPI REQUIRED)
#find_package(FFTW REQUIRED)
find_library(FFTW_LIBRARIES NAMES fftw3 libfftw3)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS})
link_directories(${MPI_Fortran_LIBRARIES} ${FFTW_LIBRARIES} ${Boost_LIBRARY_DIR})
include_directories(${MPI_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH} ${Boost_INCLUDE_DIR} ${ALaDyn_INCLUDE_DIRS} ${PROJECT_BINARY_DIR})

target_link_libraries(ALaDyn 
	precision_def 
	system_utilities 
	cpp_folder_tree 
	find_last_addr 
	mpi_var 
	phys_param 
	grid_and_particles 
	code_util 
	util 
	struct_def 
	grid_param 
	control_bunch_input 
	ionize 
	pstruct_data 
	fstruct_data 
	all_param 
	parallel 
	pic_rutil 
	grid_fields 
	der_lib 
	pdf_moments 
	pic_in 
	fft_lib 
	pic_out 
	pic_dump 
	pic_evolve_in_time 
	particles 
	pwfa_bunch_field_calculation 
	read_input 
	pwfa_output_addons
	${MPI_Fortran_LIBRARIES} 
	${MPI_LIBRARIES} 
	${FFTW_LIBRARIES} 
	${Boost_LIBRARIES}
) 

# install executables and scripts
install (TARGETS ALaDyn DESTINATION ${CMAKE_SOURCE_DIR}/bin/)

add_library(cpp_folder_tree src/cpp_folder_tree.cpp)
add_library(find_last_addr src/find_last_addr.cpp)
add_library(precision_def src/precision_def.F90)
add_library(mpi_var src/mpi_var.f90)
add_library(phys_param src/phys_param.f90)
add_library(grid_and_particles src/grid_and_particles.f90)
add_library(code_util src/code_util.f90)
add_library(system_utilities src/system_utilities.f90)
add_library(util src/util.f90)
add_library(fft_lib src/fft_lib.F90)
add_library(struct_def src/struct_def.f90)
add_library(grid_param src/grid_param.f90)
add_library(control_bunch_input src/control_bunch_input.f90)
add_library(ionize src/ionize.f90)
add_library(pstruct_data src/pstruct_data.F90)
add_library(fstruct_data src/fstruct_data.f90)
add_library(all_param src/all_param.f90)
add_library(particles src/particles.F90)
add_library(parallel src/parallel.F90)
add_library(pic_rutil src/pic_rutil.f90)
add_library(der_lib src/der_lib.f90)
add_library(grid_fields src/grid_fields.f90)
add_library(pdf_moments src/pdf_moments.f90)
add_library(pic_in src/pic_in.f90)
add_library(pic_out src/pic_out.f90)
add_library(pic_dump src/pic_dump.f90)
add_library(pic_evolve_in_time src/pic_evolve_in_time.f90)
add_library(pwfa_bunch_field_calculation src/pwfa_bunch_field_calculation.F90)
add_library(read_input src/read_input.f90)
add_library(pwfa_output_addons src/pwfa_output_addons.f90)

add_dependencies(mpi_var precision_def)
add_dependencies(phys_param precision_def)
add_dependencies(grid_and_particles precision_def)
add_dependencies(code_util precision_def)
add_dependencies(system_utilities precision_def mpi_var)
add_dependencies(util precision_def code_util grid_and_particles phys_param mpi_var)
add_dependencies(fft_lib precision_def)
add_dependencies(struct_def precision_def)
add_dependencies(grid_param precision_def struct_def)
add_dependencies(control_bunch_input precision_def)
add_dependencies(ionize precision_def)
add_dependencies(pstruct_data precision_def struct_def)
add_dependencies(fstruct_data precision_def)
add_dependencies(all_param precision_def control_bunch_input ionize mpi_var phys_param grid_and_particles code_util grid_param)
add_dependencies(particles precision_def all_param pstruct_data fstruct_data)
add_dependencies(parallel precision_def all_param pstruct_data fstruct_data fft_lib)
add_dependencies(pic_rutil precision_def pstruct_data fstruct_data util all_param parallel)
add_dependencies(der_lib precision_def util)
add_dependencies(grid_fields precision_def util der_lib all_param fft_lib parallel fstruct_data)
add_dependencies(pdf_moments precision_def all_param fstruct_data pstruct_data pic_rutil grid_fields)
add_dependencies(pic_in precision_def pic_rutil fft_lib grid_fields particles)
add_dependencies(pic_out precision_def all_param pstruct_data fstruct_data parallel)
add_dependencies(pic_dump precision_def all_param pstruct_data fstruct_data parallel)
add_dependencies(pic_evolve_in_time precision_def pic_rutil particles grid_fields)
add_dependencies(pwfa_bunch_field_calculation precision_def pwfa_bunch_field_calculation pic_in pic_evolve_in_time particles pic_rutil fft_lib grid_fields all_param pstruct_data fstruct_data)
add_dependencies(read_input precision_def control_bunch_input grid_and_particles phys_param code_util mpi_var)
add_dependencies(pwfa_output_addons precision_def pic_in pic_evolve_in_time read_input pdf_moments system_utilities pic_out pic_dump)
