To build on Windows, you must have Intel Visual Fortran 16 installed. The Microsoft Visual Studio 2015  solution then should automatically work.

To build on Mac, you should have brew installed. Then, by running  
`brew install fftw gcc boost`  
the dependencies will be satisfied. Build with `make brew`

To build on Linux, install `boost`, `fftw`, an MPI library (and a Fortran compiler if missing) and then adapt the makefile. A plain and simple `make` could also work, if you're lucky.

You can find `make` solutions already existing for many HPC systems available in Italy and Europe.

Debug and profiling `make` solutions, ready or useful to adapt to your system, are also available.
