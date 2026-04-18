Compile and link:
mkdir build
cd build
gfortran -ffixed-line-length-none -c ../src/ic_smash.f90 ../src/fluidization_smash.f90 ../src/hydro_1D_IS.f ../src/cornelius.f90 -mcmodel=large
gfortran -mcmodel=large *.o -o hydro

Initialization source auto-detection in `readinput`:
1. `../input/SMASH_IC.dat` (preferred for vHLLE comparisons)
2. `../input/SMASH_IC.oscar`
3. `../input/smash_particles.dat`
4. `../input/smash_ic.dat`
5. `../input/fort.100` (legacy fallback)

SMASH OSCAR2013Extended format accepted:
`t x y z mass p0 px py pz pdg ID charge ncoll form_time xsecfac proc_id_origin proc_type_origin time_last_coll pdg_mother1 pdg_mother2 baryon_number strangeness`

For `SMASH_IC.oscar`, a spectator filter is enabled by default:
- remove particles with `abs(baryon_number)>=1`,
- nucleon PDG (`2212` or `2112`, including anti-particles),
- `ncoll==0`,
- and `|y| >= 1.5` where `y = 0.5*ln((p0+pz)/(p0-pz))`.

SMASH particle format accepted (header from SMASH 3.2.1):
`tau x y eta mt px py Rap pdg charge baryon_number strangeness`

Legacy particle format also accepted:
`tau eta_s x y px py pz E baryon`

Particles from constant-`tau` are converted to `(t,z)` and propagated to a
common Cartesian start time `t=t0` before fluidization. If `t0<=0`, the code
uses the maximum particle `t` found in the file.

Automatic mapping diagnostics are printed at initialization:
- particle counts pre/post mapping,
- `sum(E), sum(Px), sum(Py), sum(Pz), sum(B)` pre/post,
- energy-weighted centroids `<x>_E, <y>_E, <z>_E, <t>_E` pre/post.



ene0 = 146 MeV/fm3
n0   = 0.15 fm-3

ecrit = 5*ene0
ecrit is defined at the beginning of subroutine picrel
