Compile and link:
mkdir build
cd build
gfortran -ffixed-line-length-none -c ../src/hydro_1D_IS.f ../src/cornelius.f90 -mcmodel=large
gfortran -mcmodel=large *.o -o hydro



ene0 = 146 MeV/fm3
n0   = 0.15 fm-3

ecrit = 5*ene0
ecrit is defined at the beginning of subroutine picrel

SMASH freezeout EoS controls:
- file: `input/smash_fo.cfg`
- keys:
  - `use_nq_from_za 0|1`
  - `Z <value>`
  - `A <value>`
- behavior:
  - `use_nq_from_za=0`: uses `nQ=0` in SMASH lookup
  - `use_nq_from_za=1`: uses `nQ=(Z/A)*nB` in SMASH lookup
