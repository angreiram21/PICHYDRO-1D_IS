module ic_smash
  implicit none
  private

  public :: smash_ic_options
  public :: read_smash_ic
  public :: map_smash_ic_to_hydro
  public :: enforce_smash_ic_constraints

  type :: smash_ic_options
    character(len=256) :: filename = '../input/smash_ic.dat'
    logical :: has_header = .true.
    logical :: convert_gev_to_mev = .true.
    logical :: clamp_superluminal = .true.
    real(8) :: min_nb = 0.0d0
    real(8) :: min_e = 0.0d0
  end type smash_ic_options

contains

  subroutine read_smash_ic(opts, ibar, jbar, kbar, nq, ijbar, dxdydz, &
                           e, capn, pressure, u, v, w, xmom, ymom, zmom, &
                           n_loaded, ierr)
    type(smash_ic_options), intent(in) :: opts
    integer, intent(in) :: ibar, jbar, kbar, nq, ijbar
    real(8), intent(in) :: dxdydz

    real(8), intent(inout) :: e(*)
    real(8), intent(inout) :: capn(*)
    real(8), intent(inout) :: pressure(*)
    real(8), intent(inout) :: u(*)
    real(8), intent(inout) :: v(*)
    real(8), intent(inout) :: w(*)
    real(8), intent(inout) :: xmom(*)
    real(8), intent(inout) :: ymom(*)
    real(8), intent(inout) :: zmom(*)

    integer, intent(out) :: n_loaded
    integer, intent(out) :: ierr

    integer :: unit, ios
    integer :: i1, j1, k1, ijk
    real(8) :: eijk, nbijk, pressijk, wijk, uijk, vijk
    real(8) :: zmomijk, xmomijk, ymomijk
    real(8) :: unit_scale

    n_loaded = 0
    ierr = 0

    call zero_hydro_arrays(ibar, jbar, kbar, nq, ijbar, e, capn, pressure, &
                           u, v, w, xmom, ymom, zmom)

    unit = 199
    open(unit=unit, file=trim(opts%filename), status='old', action='read', &
         iostat=ios)
    if (ios /= 0) then
      ierr = 1001
      return
    end if

    if (opts%convert_gev_to_mev) then
      unit_scale = 1000.0d0
    else
      unit_scale = 1.0d0
    end if

    if (opts%has_header) then
      call skip_smash_header(unit, ios)
      if (ios /= 0) then
        ierr = 1002
        close(unit)
        return
      end if
    end if

    do
      read(unit, *, iostat=ios) i1, j1, k1, eijk, nbijk, pressijk, wijk, &
                                 uijk, vijk, zmomijk, xmomijk, ymomijk
      if (ios /= 0) exit

      if (i1 < 1 .or. i1 > ibar) cycle
      if (j1 < 1 .or. j1 > jbar) cycle
      if (k1 < 1 .or. k1 > kbar) cycle

      ijk = nq * (ijbar * (k1 - 1) + ibar * (j1 - 1) + i1 - 1) + 1

      e(ijk) = max(eijk * unit_scale, opts%min_e)
      capn(ijk) = max(nbijk, opts%min_nb)
      pressure(ijk) = max(pressijk * unit_scale, 0.0d0)

      u(ijk) = uijk
      v(ijk) = vijk
      w(ijk) = wijk

      xmom(ijk) = xmomijk * unit_scale
      ymom(ijk) = ymomijk * unit_scale
      zmom(ijk) = zmomijk * unit_scale

      n_loaded = n_loaded + 1
    end do

    close(unit)

    if (n_loaded <= 0) then
      ierr = 1003
      return
    end if

    call enforce_smash_ic_constraints(ibar, jbar, kbar, nq, ijbar, dxdydz, &
                                      opts%clamp_superluminal, e, capn, &
                                      pressure, u, v, w, xmom, ymom, zmom)
  end subroutine read_smash_ic


  subroutine map_smash_ic_to_hydro()
    ! TODO(you): Implement when SMASH input is not already on hydro mesh.
    !
    ! Suggested responsibilities:
    !  1) Coordinate transform (tau,eta,x,y) -> (t,x,y,z) at chosen tau0.
    !  2) Spatial interpolation / smearing onto (i,j,k).
    !  3) Conservative remap checks for E and baryon number.
  end subroutine map_smash_ic_to_hydro


  subroutine enforce_smash_ic_constraints(ibar, jbar, kbar, nq, ijbar, dxdydz, &
                                          clamp_superluminal, e, capn, pressure, &
                                          u, v, w, xmom, ymom, zmom)
    integer, intent(in) :: ibar, jbar, kbar, nq, ijbar
    real(8), intent(in) :: dxdydz
    logical, intent(in) :: clamp_superluminal

    real(8), intent(inout) :: e(*)
    real(8), intent(inout) :: capn(*)
    real(8), intent(inout) :: pressure(*)
    real(8), intent(inout) :: u(*)
    real(8), intent(inout) :: v(*)
    real(8), intent(inout) :: w(*)
    real(8), intent(inout) :: xmom(*)
    real(8), intent(inout) :: ymom(*)
    real(8), intent(inout) :: zmom(*)

    integer :: i, j, k, ijk
    real(8) :: v2, scale

    do k = 1, kbar
      do j = 1, jbar
        do i = 1, ibar
          ijk = nq * (ijbar * (k - 1) + ibar * (j - 1) + i - 1) + 1

          if (e(ijk) < 0.0d0) e(ijk) = 0.0d0
          if (capn(ijk) < 0.0d0) capn(ijk) = 0.0d0
          if (pressure(ijk) < 0.0d0) pressure(ijk) = 0.0d0

          if (clamp_superluminal) then
            v2 = u(ijk) * u(ijk) + v(ijk) * v(ijk) + w(ijk) * w(ijk)
            if (v2 >= 1.0d0) then
              scale = 0.999999d0 / sqrt(v2)
              u(ijk) = u(ijk) * scale
              v(ijk) = v(ijk) * scale
              w(ijk) = w(ijk) * scale
            end if
          end if
        end do
      end do
    end do

    call report_smash_ic_integrals(ibar, jbar, kbar, nq, ijbar, dxdydz, &
                                   e, capn, xmom, ymom, zmom)
  end subroutine enforce_smash_ic_constraints


  subroutine zero_hydro_arrays(ibar, jbar, kbar, nq, ijbar, e, capn, pressure, &
                               u, v, w, xmom, ymom, zmom)
    integer, intent(in) :: ibar, jbar, kbar, nq, ijbar
    real(8), intent(inout) :: e(*)
    real(8), intent(inout) :: capn(*)
    real(8), intent(inout) :: pressure(*)
    real(8), intent(inout) :: u(*)
    real(8), intent(inout) :: v(*)
    real(8), intent(inout) :: w(*)
    real(8), intent(inout) :: xmom(*)
    real(8), intent(inout) :: ymom(*)
    real(8), intent(inout) :: zmom(*)

    integer :: i, j, k, ijk

    do k = 1, kbar
      do j = 1, jbar
        do i = 1, ibar
          ijk = nq * (ijbar * (k - 1) + ibar * (j - 1) + i - 1) + 1
          e(ijk) = 0.0d0
          capn(ijk) = 0.0d0
          pressure(ijk) = 0.0d0
          u(ijk) = 0.0d0
          v(ijk) = 0.0d0
          w(ijk) = 0.0d0
          xmom(ijk) = 0.0d0
          ymom(ijk) = 0.0d0
          zmom(ijk) = 0.0d0
        end do
      end do
    end do
  end subroutine zero_hydro_arrays


  subroutine skip_smash_header(unit, ios)
    integer, intent(in) :: unit
    integer, intent(out) :: ios
    character(len=512) :: line

    ios = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) return
      if (len_trim(line) == 0) cycle
      if (line(1:1) /= '#') then
        backspace(unit)
        return
      end if
    end do
  end subroutine skip_smash_header


  subroutine report_smash_ic_integrals(ibar, jbar, kbar, nq, ijbar, dxdydz, &
                                       e, capn, xmom, ymom, zmom)
    integer, intent(in) :: ibar, jbar, kbar, nq, ijbar
    real(8), intent(in) :: dxdydz
    real(8), intent(in) :: e(*)
    real(8), intent(in) :: capn(*)
    real(8), intent(in) :: xmom(*)
    real(8), intent(in) :: ymom(*)
    real(8), intent(in) :: zmom(*)

    integer :: i, j, k, ijk
    real(8) :: sum_e, sum_nb, sum_px, sum_py, sum_pz

    sum_e = 0.0d0
    sum_nb = 0.0d0
    sum_px = 0.0d0
    sum_py = 0.0d0
    sum_pz = 0.0d0

    do k = 1, kbar
      do j = 1, jbar
        do i = 1, ibar
          ijk = nq * (ijbar * (k - 1) + ibar * (j - 1) + i - 1) + 1
          sum_e = sum_e + e(ijk)
          sum_nb = sum_nb + capn(ijk)
          sum_px = sum_px + xmom(ijk)
          sum_py = sum_py + ymom(ijk)
          sum_pz = sum_pz + zmom(ijk)
        end do
      end do
    end do

    write(*,'(A,1PE14.6)') 'SMASH IC check: sum(T00)*V = ', sum_e * dxdydz
    write(*,'(A,1PE14.6)') 'SMASH IC check: sum(N0)*V  = ', sum_nb * dxdydz
    write(*,'(A,1PE14.6)') 'SMASH IC check: sum(T0X)*V = ', sum_px * dxdydz
    write(*,'(A,1PE14.6)') 'SMASH IC check: sum(T0Y)*V = ', sum_py * dxdydz
    write(*,'(A,1PE14.6)') 'SMASH IC check: sum(T0Z)*V = ', sum_pz * dxdydz
  end subroutine report_smash_ic_integrals

end module ic_smash
