Compile and link:
mkdir build
cd build
gfortran -ffixed-line-length-none -c ../src/ic_smash.f90 ../src/fluidization_smash.f90 ../src/hydro_1D_IS.f ../src/cornelius.f90 -mcmodel=large
gfortran -mcmodel=large *.o -o hydro

`input/inhydro` now accepts labeled entries:
- `name = values`
- `name: values`
- plain values (legacy-compatible)
- comments starting with `#` or `!` are ignored

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
common Cartesian start time before fluidization. In the current setup this
`target_t` is selected directly as `max(t_i)` from the particle list.

If `SMASH_IC.dat` (or OSCAR) contains multiple events in the same file
(`event ... start/end` blocks), the code now:
- fluidizes each event independently (including its own auto `target_t`);
- accumulates conserved fields (`T00`, momentum densities, `N0`);
- averages those fields over the number of events before hydrodynamic closure.

Optional `inhydro` flag for SMASH Milne input:
- add one extra line after `npart npartz` with a single value `eta_t0_cut`.
- if `eta_t0_cut > 0`, the per-event synchronization time is estimated as
  `target_t = tau0 * cosh(eta_t0_cut)` (using `tau0` from the first column of
  `SMASH_IC.dat`), instead of `max(t_i)`.
- if missing or `eta_t0_cut <= 0`, the previous `max(t_i)` logic is used.

Safety mechanism for auto-selected `target_t`:
- when `t0<=0`, raw `target_t=max(t_i)` is clamped to `[0.05, 3.0] fm/c`
  (configurable in `fluidization_options`);
- a yellow warning is printed if clamping is applied.
- particles with `ti>target_t` are dropped before fluidization, with a yellow
  warning including dropped fraction.
- optional hard stop is available through
  `stop_if_large_time_drop=.true.` and
  `max_time_drop_fraction` (default `0.05`).

Automatic mapping diagnostics are printed at initialization:
- particle counts pre/post mapping,
- `sum(E), sum(Px), sum(Py), sum(Pz), sum(B)` pre/post,
- energy-weighted centroids `<x>_E, <y>_E, <z>_E, <t>_E` pre/post.



ene0 = 146 MeV/fm3
n0   = 0.15 fm-3

ecrit = 5*ene0
ecrit is defined at the beginning of subroutine picrel
