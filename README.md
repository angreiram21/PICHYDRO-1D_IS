Compile and link:
mkdir build
cd build
gfortran -c ../src/hydro_3D_IS.f ../src/cornelius.f90 -mcmodel=large
gfortran *.o -o hydro


ene0 = 146 MeV/fm3
n0   = 0.15 fm-3

ecrit = 5*ene0
ecrit is defined at the beginning of subroutine picrel
