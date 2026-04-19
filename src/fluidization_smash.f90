module fluidization_smash
  implicit none
  private

  public :: fluidization_options
  public :: fluidize_smash_particles
  public :: finalize_hydro_from_tensors

  type :: fluidization_options
    character(len=256) :: particle_file = '../input/SMASH_IC.dat'
    logical :: has_header = .true.
    logical :: input_milne_tau_eta = .true.
    logical :: propagate_to_target_t = .true.
    logical :: drop_if_t_above_target = .true.
    logical :: filter_oscar_spectators = .true.
    logical :: clamp_auto_target_t = .true.
    logical :: stop_if_large_time_drop = .false.
    logical :: input_momentum_in_gev = .true.
    logical :: input_energy_in_gev = .true.
    logical :: use_gaussian_kernel = .true.
    real(8) :: sigma_x = 0.8d0
    real(8) :: sigma_y = 0.8d0
    real(8) :: sigma_z = 0.8d0
    real(8) :: target_t = -1.0d0
    real(8) :: auto_target_t_max = 3.0d0
    real(8) :: auto_target_t_min = 0.05d0
    real(8) :: max_time_drop_fraction = 0.05d0
    real(8) :: eta_t0_cut = -1.0d0
    real(8) :: spectator_abs_rapidity_cut = 1.5d0
    real(8) :: min_e = 0.0d0
    real(8) :: min_nb = 0.0d0
    integer :: kernel_half_width = 2
  end type fluidization_options

contains

  subroutine fluidize_smash_particles(opts, ibar, jbar, kbar, nq, ijbar, &
                                      dx, dy, dz, x0, y0, z0, dxdydz, &
                                      e, capn, pressure, u, v, w, xmom, ymom, zmom, &
                                      n_particles, ierr)
    type(fluidization_options), intent(in) :: opts
    integer, intent(in) :: ibar, jbar, kbar, nq, ijbar
    real(8), intent(in) :: dx, dy, dz, x0, y0, z0, dxdydz

    real(8), intent(inout) :: e(*)
    real(8), intent(inout) :: capn(*)
    real(8), intent(inout) :: pressure(*)
    real(8), intent(inout) :: u(*)
    real(8), intent(inout) :: v(*)
    real(8), intent(inout) :: w(*)
    real(8), intent(inout) :: xmom(*)
    real(8), intent(inout) :: ymom(*)
    real(8), intent(inout) :: zmom(*)

    integer, intent(out) :: n_particles
    integer, intent(out) :: ierr

    integer :: unit, ios
    integer :: n_dropped_time, n_dropped_zeroe, n_dropped_spectators
    integer :: n_prefilter
    integer :: n_events
    character(len=512) :: line
    real(8) :: c1, c2, c3, c4, px, py, pz, ep, baryon
    logical :: is_milne_line, is_oscar_line
    real(8) :: pdg_meta, ncoll_meta
    real(8) :: t_part, x_part, y_part, z_part
    real(8) :: vx, vy, vz, dt_prop
    real(8) :: target_t_eff, tiny, target_t_raw
    real(8) :: time_drop_fraction, escale_diag, pscale_diag
    real(8) :: sum_e_pre, sum_px_pre, sum_py_pre, sum_pz_pre, sum_b_pre
    real(8) :: sum_ex_pre, sum_ey_pre, sum_ez_pre, sum_et_pre
    real(8) :: sum_e_post, sum_px_post, sum_py_post, sum_pz_post, sum_b_post
    real(8) :: sum_ex_post, sum_ey_post, sum_ez_post, sum_et_post
    real(8) :: sum_s_pre, sum_s_post
    real(8) :: sum_e_fluid, sum_px_fluid, sum_py_fluid, sum_pz_fluid, sum_b_fluid, sum_s_fluid
    logical :: has_event_markers, in_event
    integer :: nbuf, ncap
    real(8), allocatable :: buf_t(:), buf_x(:), buf_y(:), buf_z(:)
    real(8), allocatable :: buf_px(:), buf_py(:), buf_pz(:), buf_e(:), buf_b(:)
    integer :: n_prefilter_evt, n_dropped_zeroe_evt, n_dropped_time_evt
    integer :: n_particles_evt
    real(8) :: t_max_evt, tau0_evt
    real(8) :: sum_e_pre_evt, sum_px_pre_evt, sum_py_pre_evt, sum_pz_pre_evt, sum_b_pre_evt
    real(8) :: sum_ex_pre_evt, sum_ey_pre_evt, sum_ez_pre_evt, sum_et_pre_evt
    real(8) :: sum_e_post_evt, sum_px_post_evt, sum_py_post_evt, sum_pz_post_evt, sum_b_post_evt
    real(8) :: sum_ex_post_evt, sum_ey_post_evt, sum_ez_post_evt, sum_et_post_evt
    real(8) :: sum_s_pre_evt, sum_s_post_evt

    n_particles = 0
    ierr = 0
    n_events = 0
    n_dropped_time = 0
    n_dropped_zeroe = 0
    n_dropped_spectators = 0
    n_prefilter = 0
    tiny = 1.0d-20
    has_event_markers = .false.
    in_event = .false.
    nbuf = 0
    ncap = 0
    sum_e_pre = 0.0d0
    sum_px_pre = 0.0d0
    sum_py_pre = 0.0d0
    sum_pz_pre = 0.0d0
    sum_b_pre = 0.0d0
    sum_ex_pre = 0.0d0
    sum_ey_pre = 0.0d0
    sum_ez_pre = 0.0d0
    sum_et_pre = 0.0d0
    sum_e_post = 0.0d0
    sum_px_post = 0.0d0
    sum_py_post = 0.0d0
    sum_pz_post = 0.0d0
    sum_b_post = 0.0d0
    sum_ex_post = 0.0d0
    sum_ey_post = 0.0d0
    sum_ez_post = 0.0d0
    sum_et_post = 0.0d0
    sum_s_pre = 0.0d0
    sum_s_post = 0.0d0

    if (opts%input_energy_in_gev) then
      escale_diag = 1000.0d0
    else
      escale_diag = 1.0d0
    end if
    if (opts%input_momentum_in_gev) then
      pscale_diag = 1000.0d0
    else
      pscale_diag = 1.0d0
    end if

    call zero_hydro_arrays(ibar, jbar, kbar, nq, ijbar, e, capn, pressure, &
                           u, v, w, xmom, ymom, zmom)

    unit = 209
    open(unit=unit, file=trim(opts%particle_file), status='old', action='read', &
         iostat=ios)
    if (ios /= 0) then
      ierr = 2001
      return
    end if

    call reset_event_state()

    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') then
        if (is_event_start_marker(line)) then
          has_event_markers = .true.
          if (in_event) then
            call flush_event(.true.)
            if (ierr /= 0) then
              close(unit)
              call cleanup_buffers()
              return
            end if
          end if
          in_event = .true.
          cycle
        end if
        if (is_event_end_marker(line)) then
          has_event_markers = .true.
          if (in_event) then
            call flush_event(.true.)
            if (ierr /= 0) then
              close(unit)
              call cleanup_buffers()
              return
            end if
          end if
          in_event = .false.
          cycle
        end if
        cycle
      end if

      if (.not.in_event) in_event = .true.

      call parse_smash_particle_line(line, ios, c1, c2, c3, c4, px, py, pz, ep, baryon, &
                                     is_milne_line, is_oscar_line, pdg_meta, ncoll_meta)
      if (ios /= 0) cycle

      if (opts%filter_oscar_spectators .and. is_oscar_line) then
        if (is_oscar_spectator(ep, pz, pdg_meta, baryon, ncoll_meta, &
                               opts%spectator_abs_rapidity_cut)) then
          n_dropped_spectators = n_dropped_spectators + 1
          cycle
        end if
      end if

      if (is_milne_line) then
        if (tau0_evt <= 0.0d0) tau0_evt = c1
        call milne_to_cartesian(c1, c2, c3, c4, t_part, x_part, y_part, z_part)
      else
        t_part = c1
        x_part = c2
        y_part = c3
        z_part = c4
      end if

      n_prefilter_evt = n_prefilter_evt + 1
      sum_e_pre_evt = sum_e_pre_evt + ep * escale_diag
      sum_px_pre_evt = sum_px_pre_evt + px * pscale_diag
      sum_py_pre_evt = sum_py_pre_evt + py * pscale_diag
      sum_pz_pre_evt = sum_pz_pre_evt + pz * pscale_diag
      sum_b_pre_evt = sum_b_pre_evt + baryon
      sum_ex_pre_evt = sum_ex_pre_evt + ep * escale_diag * x_part
      sum_ey_pre_evt = sum_ey_pre_evt + ep * escale_diag * y_part
      sum_ez_pre_evt = sum_ez_pre_evt + ep * escale_diag * z_part
      sum_et_pre_evt = sum_et_pre_evt + ep * escale_diag * t_part
      sum_s_pre_evt = sum_s_pre_evt + entropy_proxy_scalar(ep * escale_diag)

      if (ep <= tiny) then
        n_dropped_zeroe_evt = n_dropped_zeroe_evt + 1
        cycle
      end if

      call append_particle(t_part, x_part, y_part, z_part, px, py, pz, ep, baryon, ierr)
      if (ierr /= 0) then
        close(unit)
        call cleanup_buffers()
        return
      end if
      if (t_part > t_max_evt) t_max_evt = t_part
    end do

    if (in_event .or. (.not.has_event_markers .and. n_prefilter_evt > 0)) then
      call flush_event(.true.)
      if (ierr /= 0) then
        close(unit)
        call cleanup_buffers()
        return
      end if
    end if

    close(unit)
    call cleanup_buffers()

    if (n_events <= 0 .or. n_particles <= 0) then
      ierr = 2003
      return
    end if

    call scale_conserved_arrays(ibar, jbar, kbar, nq, ijbar, &
                                1.0d0 / dble(n_events), e, capn, xmom, ymom, zmom)

    call finalize_hydro_from_tensors(ibar, jbar, kbar, nq, ijbar, dxdydz, &
                                     opts%min_e, opts%min_nb, e, capn, pressure, &
                                     u, v, w, xmom, ymom, zmom)
    call compute_fluid_stage_totals(ibar, jbar, kbar, nq, ijbar, dxdydz, e, capn, pressure, &
                                    xmom, ymom, zmom, sum_e_fluid, sum_px_fluid, sum_py_fluid, &
                                    sum_pz_fluid, sum_b_fluid, sum_s_fluid)
    write(*,'(A,I10)') 'Fluidization: events averaged = ', n_events
    write(*,'(A,I10)') 'Fluidization: particles kept (sum over events) = ', n_particles
    if (n_dropped_time > 0) then
      write(*,'(A,I10)') 'Fluidization: particles dropped (t>target_t) = ', n_dropped_time
      if (n_prefilter > 0) then
        time_drop_fraction = dble(n_dropped_time) / dble(n_prefilter)
      else
        time_drop_fraction = 0.0d0
      end if
      call print_time_drop_warning(n_dropped_time, n_prefilter, time_drop_fraction, &
                                   opts%max_time_drop_fraction, opts%stop_if_large_time_drop, ierr)
      if (ierr /= 0) return
    end if
    if (n_dropped_zeroe > 0) then
      write(*,'(A,I10)') 'Fluidization: particles dropped (E<=0) = ', n_dropped_zeroe
    end if
    if (n_dropped_spectators > 0) then
      write(*,'(A,I10)') 'Fluidization: particles dropped (OSCAR spectators) = ', n_dropped_spectators
    end if
    call print_mapping_diagnostics(n_prefilter, n_particles, &
                                   sum_e_pre, sum_px_pre, sum_py_pre, sum_pz_pre, sum_b_pre, &
                                   sum_ex_pre, sum_ey_pre, sum_ez_pre, sum_et_pre, &
                                   sum_e_post, sum_px_post, sum_py_post, sum_pz_post, sum_b_post, &
                                   sum_ex_post, sum_ey_post, sum_ez_post, sum_et_post)
    call print_conservation_table(n_events, sum_e_pre, sum_px_pre, sum_py_pre, sum_pz_pre, sum_b_pre, &
                                  sum_s_pre, sum_e_post, sum_px_post, sum_py_post, sum_pz_post, &
                                  sum_b_post, sum_s_post, sum_e_fluid, sum_px_fluid, sum_py_fluid, &
                                  sum_pz_fluid, sum_b_fluid, sum_s_fluid)
    if (opts%propagate_to_target_t) then
      call print_yellow_mapping_warning(-1.0d0)
    end if
  contains
    subroutine reset_event_state()
      nbuf = 0
      n_prefilter_evt = 0
      n_dropped_zeroe_evt = 0
      n_dropped_time_evt = 0
      n_particles_evt = 0
      t_max_evt = -1.0d30
      tau0_evt = -1.0d0
      sum_e_pre_evt = 0.0d0
      sum_px_pre_evt = 0.0d0
      sum_py_pre_evt = 0.0d0
      sum_pz_pre_evt = 0.0d0
      sum_b_pre_evt = 0.0d0
      sum_ex_pre_evt = 0.0d0
      sum_ey_pre_evt = 0.0d0
      sum_ez_pre_evt = 0.0d0
      sum_et_pre_evt = 0.0d0
      sum_e_post_evt = 0.0d0
      sum_px_post_evt = 0.0d0
      sum_py_post_evt = 0.0d0
      sum_pz_post_evt = 0.0d0
      sum_b_post_evt = 0.0d0
      sum_ex_post_evt = 0.0d0
      sum_ey_post_evt = 0.0d0
      sum_ez_post_evt = 0.0d0
      sum_et_post_evt = 0.0d0
      sum_s_pre_evt = 0.0d0
      sum_s_post_evt = 0.0d0
    end subroutine reset_event_state

    subroutine flush_event(count_event)
      logical, intent(in) :: count_event
      integer :: ip

      if (count_event) n_events = n_events + 1
      if (n_prefilter_evt <= 0 .and. nbuf <= 0) then
        call reset_event_state()
        return
      end if

      target_t_eff = opts%target_t
      target_t_raw = target_t_eff
      if (opts%propagate_to_target_t .and. target_t_eff <= 0.0d0) then
        if (opts%eta_t0_cut > 0.0d0 .and. tau0_evt > tiny) then
          target_t_raw = tau0_evt * cosh(opts%eta_t0_cut)
          target_t_eff = target_t_raw
        else
          if (t_max_evt < -1.0d20) then
            ierr = 2004
            return
          end if
          target_t_raw = t_max_evt
          target_t_eff = t_max_evt
        end if
        if (target_t_eff < -1.0d20) then
          ierr = 2004
          return
        end if
        if (opts%clamp_auto_target_t) then
          target_t_eff = min(target_t_eff, opts%auto_target_t_max)
          target_t_eff = max(target_t_eff, opts%auto_target_t_min)
        end if
      end if

      do ip = 1, nbuf
        t_part = buf_t(ip)
        x_part = buf_x(ip)
        y_part = buf_y(ip)
        z_part = buf_z(ip)
        px = buf_px(ip)
        py = buf_py(ip)
        pz = buf_pz(ip)
        ep = buf_e(ip)
        baryon = buf_b(ip)

        if (opts%propagate_to_target_t) then
          vx = px / ep
          vy = py / ep
          vz = pz / ep
          if (vx*vx + vy*vy + vz*vz > 1.0d0) then
            call clamp_velocity(vx, vy, vz)
          end if
          dt_prop = target_t_eff - t_part
          if (dt_prop < 0.0d0 .and. opts%drop_if_t_above_target) then
            n_dropped_time_evt = n_dropped_time_evt + 1
            cycle
          end if
          x_part = x_part + vx * dt_prop
          y_part = y_part + vy * dt_prop
          z_part = z_part + vz * dt_prop
          t_part = target_t_eff
        end if

        sum_e_post_evt = sum_e_post_evt + ep * escale_diag
        sum_px_post_evt = sum_px_post_evt + px * pscale_diag
        sum_py_post_evt = sum_py_post_evt + py * pscale_diag
        sum_pz_post_evt = sum_pz_post_evt + pz * pscale_diag
        sum_b_post_evt = sum_b_post_evt + baryon
        sum_ex_post_evt = sum_ex_post_evt + ep * escale_diag * x_part
        sum_ey_post_evt = sum_ey_post_evt + ep * escale_diag * y_part
        sum_ez_post_evt = sum_ez_post_evt + ep * escale_diag * z_part
        sum_et_post_evt = sum_et_post_evt + ep * escale_diag * t_part
        sum_s_post_evt = sum_s_post_evt + entropy_proxy_scalar(ep * escale_diag)

        call deposit_particle(opts, ibar, jbar, kbar, nq, ijbar, dx, dy, dz, &
                              x0, y0, z0, x_part, y_part, z_part, px, py, pz, ep, baryon, &
                              e, capn, xmom, ymom, zmom)
        n_particles_evt = n_particles_evt + 1
      end do

      n_prefilter = n_prefilter + n_prefilter_evt
      n_dropped_zeroe = n_dropped_zeroe + n_dropped_zeroe_evt
      n_dropped_time = n_dropped_time + n_dropped_time_evt
      n_particles = n_particles + n_particles_evt
      sum_e_pre = sum_e_pre + sum_e_pre_evt
      sum_px_pre = sum_px_pre + sum_px_pre_evt
      sum_py_pre = sum_py_pre + sum_py_pre_evt
      sum_pz_pre = sum_pz_pre + sum_pz_pre_evt
      sum_b_pre = sum_b_pre + sum_b_pre_evt
      sum_ex_pre = sum_ex_pre + sum_ex_pre_evt
      sum_ey_pre = sum_ey_pre + sum_ey_pre_evt
      sum_ez_pre = sum_ez_pre + sum_ez_pre_evt
      sum_et_pre = sum_et_pre + sum_et_pre_evt
      sum_e_post = sum_e_post + sum_e_post_evt
      sum_px_post = sum_px_post + sum_px_post_evt
      sum_py_post = sum_py_post + sum_py_post_evt
      sum_pz_post = sum_pz_post + sum_pz_post_evt
      sum_b_post = sum_b_post + sum_b_post_evt
      sum_ex_post = sum_ex_post + sum_ex_post_evt
      sum_ey_post = sum_ey_post + sum_ey_post_evt
      sum_ez_post = sum_ez_post + sum_ez_post_evt
      sum_et_post = sum_et_post + sum_et_post_evt
      sum_s_pre = sum_s_pre + sum_s_pre_evt
      sum_s_post = sum_s_post + sum_s_post_evt

      if (count_event) then
        write(*,'(A,I6,A,I8,A,I8,A,I8,A,1PE12.4)') &
          'Fluidization event ', n_events, ': pre=', n_prefilter_evt, &
          ' kept=', n_particles_evt, ' drop_t=', n_dropped_time_evt, &
          ' target_t=', target_t_eff
        if (opts%eta_t0_cut > 0.0d0 .and. tau0_evt > tiny) then
          write(*,'(A,1PE12.4,A,1PE12.4)') '  t0 mode: tau0*cosh(eta_cut), tau0=', &
                                            tau0_evt, ' eta_cut=', opts%eta_t0_cut
        end if
      end if

      if (opts%propagate_to_target_t .and. opts%clamp_auto_target_t .and. opts%target_t <= 0.0d0) then
        if (abs(target_t_eff-target_t_raw) > 1.0d-12) then
          call print_target_t_clamp_warning(target_t_raw, target_t_eff, &
                                            opts%auto_target_t_min, opts%auto_target_t_max)
        end if
      end if
      if (n_dropped_time_evt > 0) then
        if (n_prefilter_evt > 0) then
          time_drop_fraction = dble(n_dropped_time_evt) / dble(n_prefilter_evt)
        else
          time_drop_fraction = 0.0d0
        end if
        call print_time_drop_warning(n_dropped_time_evt, n_prefilter_evt, time_drop_fraction, &
                                     opts%max_time_drop_fraction, opts%stop_if_large_time_drop, ierr)
        if (ierr /= 0) return
      end if

      call reset_event_state()
    end subroutine flush_event

    subroutine append_particle(ti, xi, yi, zi, pxi, pyi, pzi, ei, bi, ierr_loc)
      real(8), intent(in) :: ti, xi, yi, zi, pxi, pyi, pzi, ei, bi
      integer, intent(out) :: ierr_loc
      ierr_loc = 0
      if (nbuf >= ncap) then
        call grow_particle_buffer(ierr_loc)
        if (ierr_loc /= 0) return
      end if
      nbuf = nbuf + 1
      buf_t(nbuf) = ti
      buf_x(nbuf) = xi
      buf_y(nbuf) = yi
      buf_z(nbuf) = zi
      buf_px(nbuf) = pxi
      buf_py(nbuf) = pyi
      buf_pz(nbuf) = pzi
      buf_e(nbuf) = ei
      buf_b(nbuf) = bi
    end subroutine append_particle

    subroutine grow_particle_buffer(ierr_loc)
      integer, intent(out) :: ierr_loc
      integer :: new_cap

      ierr_loc = 0
      if (ncap <= 0) then
        new_cap = 2048
      else
        new_cap = 2 * ncap
      end if

      call resize_real8_array(buf_t, new_cap, ierr_loc); if (ierr_loc /= 0) return
      call resize_real8_array(buf_x, new_cap, ierr_loc); if (ierr_loc /= 0) return
      call resize_real8_array(buf_y, new_cap, ierr_loc); if (ierr_loc /= 0) return
      call resize_real8_array(buf_z, new_cap, ierr_loc); if (ierr_loc /= 0) return
      call resize_real8_array(buf_px, new_cap, ierr_loc); if (ierr_loc /= 0) return
      call resize_real8_array(buf_py, new_cap, ierr_loc); if (ierr_loc /= 0) return
      call resize_real8_array(buf_pz, new_cap, ierr_loc); if (ierr_loc /= 0) return
      call resize_real8_array(buf_e, new_cap, ierr_loc); if (ierr_loc /= 0) return
      call resize_real8_array(buf_b, new_cap, ierr_loc); if (ierr_loc /= 0) return
      ncap = new_cap
    end subroutine grow_particle_buffer

    subroutine cleanup_buffers()
      if (allocated(buf_t)) deallocate(buf_t)
      if (allocated(buf_x)) deallocate(buf_x)
      if (allocated(buf_y)) deallocate(buf_y)
      if (allocated(buf_z)) deallocate(buf_z)
      if (allocated(buf_px)) deallocate(buf_px)
      if (allocated(buf_py)) deallocate(buf_py)
      if (allocated(buf_pz)) deallocate(buf_pz)
      if (allocated(buf_e)) deallocate(buf_e)
      if (allocated(buf_b)) deallocate(buf_b)
      nbuf = 0
      ncap = 0
    end subroutine cleanup_buffers
  end subroutine fluidize_smash_particles


  subroutine finalize_hydro_from_tensors(ibar, jbar, kbar, nq, ijbar, dxdydz, &
                                         min_e, min_nb, e, capn, pressure, &
                                         u, v, w, xmom, ymom, zmom)
    integer, intent(in) :: ibar, jbar, kbar, nq, ijbar
    real(8), intent(in) :: dxdydz, min_e, min_nb

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
    real(8) :: v2, scale, tiny
    real(8) :: sum_e, sum_nb

    tiny = 1.0d-20
    sum_e = 0.0d0
    sum_nb = 0.0d0

    do k = 1, kbar
      do j = 1, jbar
        do i = 1, ibar
          ijk = nq * (ijbar * (k - 1) + ibar * (j - 1) + i - 1) + 1

          e(ijk) = max(e(ijk), min_e)
          capn(ijk) = max(capn(ijk), min_nb)

          if (e(ijk) > tiny) then
            u(ijk) = xmom(ijk) / e(ijk)
            v(ijk) = ymom(ijk) / e(ijk)
            w(ijk) = zmom(ijk) / e(ijk)
          else
            u(ijk) = 0.0d0
            v(ijk) = 0.0d0
            w(ijk) = 0.0d0
          end if

          v2 = u(ijk)*u(ijk) + v(ijk)*v(ijk) + w(ijk)*w(ijk)
          if (v2 >= 1.0d0) then
            scale = 0.999999d0 / sqrt(v2)
            u(ijk) = u(ijk) * scale
            v(ijk) = v(ijk) * scale
            w(ijk) = w(ijk) * scale
          end if

          ! Placeholder closure: replace with EoS call if needed.
          pressure(ijk) = e(ijk) / 3.0d0

          sum_e = sum_e + e(ijk)
          sum_nb = sum_nb + capn(ijk)
        end do
      end do
    end do

    write(*,'(A,1PE14.6)') 'Fluidization check: sum(T00)*V = ', sum_e * dxdydz
    write(*,'(A,1PE14.6)') 'Fluidization check: sum(N0)*V  = ', sum_nb * dxdydz
  end subroutine finalize_hydro_from_tensors


  subroutine deposit_particle(opts, ibar, jbar, kbar, nq, ijbar, dx, dy, dz, &
                              x0, y0, z0, xp, yp, zp, px, py, pz, ep, baryon, &
                              e, capn, xmom, ymom, zmom)
    type(fluidization_options), intent(in) :: opts
    integer, intent(in) :: ibar, jbar, kbar, nq, ijbar
    real(8), intent(in) :: dx, dy, dz, x0, y0, z0
    real(8), intent(in) :: xp, yp, zp, px, py, pz, ep, baryon

    real(8), intent(inout) :: e(*)
    real(8), intent(inout) :: capn(*)
    real(8), intent(inout) :: xmom(*)
    real(8), intent(inout) :: ymom(*)
    real(8), intent(inout) :: zmom(*)

    integer :: ic, jc, kc
    integer :: i, j, k, ijk
    integer :: di, dj, dk
    integer :: r
    real(8) :: wx, wy, wz, wtot
    real(8) :: dxp, dyp, dzp
    real(8) :: escale, pscale, inv_cell_volume, cell_volume

    if (opts%input_energy_in_gev) then
      escale = 1000.0d0
    else
      escale = 1.0d0
    end if

    if (opts%input_momentum_in_gev) then
      pscale = 1000.0d0
    else
      pscale = 1.0d0
    end if

    cell_volume = dx * dy * dz
    if (cell_volume > 0.0d0) then
      inv_cell_volume = 1.0d0 / cell_volume
    else
      inv_cell_volume = 0.0d0
    end if

    ic = int((xp - x0) / dx) + 1
    jc = int((yp - y0) / dy) + 1
    kc = int((zp - z0) / dz) + 1

    if (.not.opts%use_gaussian_kernel) then
      if (ic < 1 .or. ic > ibar) return
      if (jc < 1 .or. jc > jbar) return
      if (kc < 1 .or. kc > kbar) return
      ijk = nq * (ijbar * (kc - 1) + ibar * (jc - 1) + ic - 1) + 1
      e(ijk) = e(ijk) + ep * escale * inv_cell_volume
      capn(ijk) = capn(ijk) + baryon * inv_cell_volume
      xmom(ijk) = xmom(ijk) + px * pscale * inv_cell_volume
      ymom(ijk) = ymom(ijk) + py * pscale * inv_cell_volume
      zmom(ijk) = zmom(ijk) + pz * pscale * inv_cell_volume
      return
    end if

    r = max(opts%kernel_half_width, 0)

    wtot = 0.0d0
    do dk = -r, r
      k = kc + dk
      if (k < 1 .or. k > kbar) cycle
      dzp = (z0 + (dble(k) - 0.5d0) * dz) - zp
      wz = gaussian_1d(dzp, opts%sigma_z)

      do dj = -r, r
        j = jc + dj
        if (j < 1 .or. j > jbar) cycle
        dyp = (y0 + (dble(j) - 0.5d0) * dy) - yp
        wy = gaussian_1d(dyp, opts%sigma_y)

        do di = -r, r
          i = ic + di
          if (i < 1 .or. i > ibar) cycle
          dxp = (x0 + (dble(i) - 0.5d0) * dx) - xp
          wx = gaussian_1d(dxp, opts%sigma_x)
          wtot = wtot + wx * wy * wz
        end do
      end do
    end do

    if (wtot <= 0.0d0) return

    do dk = -r, r
      k = kc + dk
      if (k < 1 .or. k > kbar) cycle
      dzp = (z0 + (dble(k) - 0.5d0) * dz) - zp
      wz = gaussian_1d(dzp, opts%sigma_z)

      do dj = -r, r
        j = jc + dj
        if (j < 1 .or. j > jbar) cycle
        dyp = (y0 + (dble(j) - 0.5d0) * dy) - yp
        wy = gaussian_1d(dyp, opts%sigma_y)

        do di = -r, r
          i = ic + di
          if (i < 1 .or. i > ibar) cycle
          dxp = (x0 + (dble(i) - 0.5d0) * dx) - xp
          wx = gaussian_1d(dxp, opts%sigma_x)
          ijk = nq * (ijbar * (k - 1) + ibar * (j - 1) + i - 1) + 1

          wx = (wx * wy * wz) / wtot
          e(ijk) = e(ijk) + wx * ep * escale * inv_cell_volume
          capn(ijk) = capn(ijk) + wx * baryon * inv_cell_volume
          xmom(ijk) = xmom(ijk) + wx * px * pscale * inv_cell_volume
          ymom(ijk) = ymom(ijk) + wx * py * pscale * inv_cell_volume
          zmom(ijk) = zmom(ijk) + wx * pz * pscale * inv_cell_volume
        end do
      end do
    end do
  end subroutine deposit_particle


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


  subroutine scale_conserved_arrays(ibar, jbar, kbar, nq, ijbar, factor, e, capn, xmom, ymom, zmom)
    integer, intent(in) :: ibar, jbar, kbar, nq, ijbar
    real(8), intent(in) :: factor
    real(8), intent(inout) :: e(*), capn(*), xmom(*), ymom(*), zmom(*)
    integer :: i, j, k, ijk

    do k = 1, kbar
      do j = 1, jbar
        do i = 1, ibar
          ijk = nq * (ijbar * (k - 1) + ibar * (j - 1) + i - 1) + 1
          e(ijk) = factor * e(ijk)
          capn(ijk) = factor * capn(ijk)
          xmom(ijk) = factor * xmom(ijk)
          ymom(ijk) = factor * ymom(ijk)
          zmom(ijk) = factor * zmom(ijk)
        end do
      end do
    end do
  end subroutine scale_conserved_arrays


  subroutine compute_fluid_stage_totals(ibar, jbar, kbar, nq, ijbar, dxdydz, e, capn, pressure, &
                                        xmom, ymom, zmom, se, spx, spy, spz, sb, ss)
    integer, intent(in) :: ibar, jbar, kbar, nq, ijbar
    real(8), intent(in) :: dxdydz
    real(8), intent(in) :: e(*), capn(*), pressure(*), xmom(*), ymom(*), zmom(*)
    real(8), intent(out) :: se, spx, spy, spz, sb, ss
    integer :: i, j, k, ijk

    se = 0.0d0
    spx = 0.0d0
    spy = 0.0d0
    spz = 0.0d0
    sb = 0.0d0
    ss = 0.0d0

    do k = 1, kbar
      do j = 1, jbar
        do i = 1, ibar
          ijk = nq * (ijbar * (k - 1) + ibar * (j - 1) + i - 1) + 1
          se = se + e(ijk) * dxdydz
          spx = spx + xmom(ijk) * dxdydz
          spy = spy + ymom(ijk) * dxdydz
          spz = spz + zmom(ijk) * dxdydz
          sb = sb + capn(ijk) * dxdydz
          ss = ss + entropy_proxy_scalar(max(e(ijk) + pressure(ijk), 0.0d0)) * dxdydz
        end do
      end do
    end do
  end subroutine compute_fluid_stage_totals


  real(8) function entropy_proxy_scalar(x)
    real(8), intent(in) :: x
    real(8), parameter :: tiny = 1.0d-20
    if (x <= tiny) then
      entropy_proxy_scalar = 0.0d0
    else
      entropy_proxy_scalar = x ** 0.75d0
    end if
  end function entropy_proxy_scalar


  subroutine print_conservation_table(n_events, e_pre, px_pre, py_pre, pz_pre, b_pre, s_pre, &
                                      e_post, px_post, py_post, pz_post, b_post, s_post, &
                                      e_fl, px_fl, py_fl, pz_fl, b_fl, s_fl)
    integer, intent(in) :: n_events
    real(8), intent(in) :: e_pre, px_pre, py_pre, pz_pre, b_pre, s_pre
    real(8), intent(in) :: e_post, px_post, py_post, pz_post, b_post, s_post
    real(8), intent(in) :: e_fl, px_fl, py_fl, pz_fl, b_fl, s_fl
    real(8) :: norm, a_e_pre, a_px_pre, a_py_pre, a_pz_pre, a_b_pre, a_s_pre
    real(8) :: a_e_post, a_px_post, a_py_post, a_pz_post, a_b_post, a_s_post
    real(8) :: a_e_fl, a_px_fl, a_py_fl, a_pz_fl, a_b_fl, a_s_fl
    character(len=*), parameter :: c_blue = achar(27)//'[34m'
    character(len=*), parameter :: c_reset = achar(27)//'[0m'

    norm = 1.0d0
    if (n_events > 1) norm = dble(n_events)

    a_e_pre = e_pre / norm
    a_px_pre = px_pre / norm
    a_py_pre = py_pre / norm
    a_pz_pre = pz_pre / norm
    a_b_pre = b_pre / norm
    a_s_pre = s_pre / norm

    a_e_post = e_post / norm
    a_px_post = px_post / norm
    a_py_post = py_post / norm
    a_pz_post = pz_post / norm
    a_b_post = b_post / norm
    a_s_post = s_post / norm

    ! Fluid stage is already per-event averaged when n_events>1.
    a_e_fl = e_fl
    a_px_fl = px_fl
    a_py_fl = py_fl
    a_pz_fl = pz_fl
    a_b_fl = b_fl
    a_s_fl = s_fl

    write(*,'(A)') c_blue//'================ Conservation Check (Blue) ================'
    if (n_events > 1) then
      write(*,'(A,I8)') 'Averaged over events: ', n_events
    else
      write(*,'(A)') 'Single-event initialization'
    end if
    write(*,'(A)') 'Stage                  E_tot        Px_tot       Py_tot       Pz_tot       B_tot        S_proxy'
    write(*,'(A,6(1X,1PE12.4))') 'Before mapping     ', a_e_pre, a_px_pre, a_py_pre, a_pz_pre, a_b_pre, a_s_pre
    write(*,'(A,6(1X,1PE12.4))') 'After mapping      ', a_e_post, a_px_post, a_py_post, a_pz_post, a_b_post, a_s_post
    write(*,'(A,6(1X,1PE12.4))') 'After fluidization ', a_e_fl, a_px_fl, a_py_fl, a_pz_fl, a_b_fl, a_s_fl
    write(*,'(A)') '============================================================'//c_reset
  end subroutine print_conservation_table


  subroutine resize_real8_array(arr, new_size, ierr)
    real(8), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: new_size
    integer, intent(out) :: ierr
    integer :: old_size
    real(8), allocatable :: tmp(:)

    ierr = 0
    old_size = 0
    if (allocated(arr)) old_size = size(arr)

    allocate(tmp(new_size), stat=ierr)
    if (ierr /= 0) then
      ierr = 2006
      return
    end if
    tmp = 0.0d0
    if (old_size > 0) tmp(1:old_size) = arr(1:old_size)
    if (allocated(arr)) deallocate(arr)
    call move_alloc(tmp, arr)
  end subroutine resize_real8_array


  logical function is_event_start_marker(line)
    character(len=*), intent(in) :: line
    is_event_start_marker = index(line, '# event') > 0 .and. &
                            (index(line, 'start') > 0 .or. index(line, ' out ') > 0)
  end function is_event_start_marker


  logical function is_event_end_marker(line)
    character(len=*), intent(in) :: line
    is_event_end_marker = index(line, '# event') > 0 .and. index(line, 'end') > 0
  end function is_event_end_marker


  subroutine skip_header(unit, ios)
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
  end subroutine skip_header


  subroutine compute_max_particle_time(unit, has_header, t_max_input, ios)
    integer, intent(in) :: unit
    logical, intent(in) :: has_header
    real(8), intent(out) :: t_max_input
    integer, intent(out) :: ios

    character(len=512) :: line
    real(8) :: c1, c2, c3, c4, px, py, pz, ep, baryon
    logical :: is_milne_line, is_oscar_line
    real(8) :: pdg_meta, ncoll_meta
    real(8) :: t_part, x_part, y_part, z_part

    t_max_input = -1.0d30
    rewind(unit)
    if (has_header) then
      call skip_header(unit, ios)
      if (ios /= 0) return
    end if

    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      call parse_smash_particle_line(line, ios, c1, c2, c3, c4, px, py, pz, ep, baryon, &
                                     is_milne_line, is_oscar_line, pdg_meta, ncoll_meta)
      if (ios /= 0) cycle
      if (is_milne_line) then
        call milne_to_cartesian(c1, c2, c3, c4, t_part, x_part, y_part, z_part)
      else
        t_part = c1
      end if
      if (t_part > t_max_input) t_max_input = t_part
    end do

    if (t_max_input < -1.0d20) then
      ios = 1
    else
      ios = 0
    end if
  end subroutine compute_max_particle_time


  subroutine milne_to_cartesian(tau, eta_s, x_in, y_in, t_out, x_out, y_out, z_out)
    real(8), intent(in) :: tau, eta_s, x_in, y_in
    real(8), intent(out) :: t_out, x_out, y_out, z_out

    t_out = tau * cosh(eta_s)
    z_out = tau * sinh(eta_s)
    x_out = x_in
    y_out = y_in
  end subroutine milne_to_cartesian


  subroutine clamp_velocity(vx, vy, vz)
    real(8), intent(inout) :: vx, vy, vz
    real(8) :: v2, scale

    v2 = vx*vx + vy*vy + vz*vz
    if (v2 >= 1.0d0) then
      scale = 0.999999d0 / sqrt(v2)
      vx = vx * scale
      vy = vy * scale
      vz = vz * scale
    end if
  end subroutine clamp_velocity


  subroutine parse_smash_particle_line(line, ios, c1, c2, c3, c4, px, py, pz, ep, baryon, &
                                       is_milne_line, is_oscar_line, pdg_meta, ncoll_meta)
    character(len=*), intent(in) :: line
    integer, intent(out) :: ios
    real(8), intent(out) :: c1, c2, c3, c4, px, py, pz, ep, baryon
    logical, intent(out) :: is_milne_line, is_oscar_line
    real(8), intent(out) :: pdg_meta, ncoll_meta

    real(8) :: tau_i, x_i, y_i, eta_i, mt_i, px_i, py_i, rap_i
    real(8) :: pdg_i, charge_i, baryon_i, strange_i
    real(8) :: t_i, z_i, mass_i, e_i, pz_i, id_i
    real(8) :: ncoll_i, form_i, xsecfac_i, procid_i, proctype_i
    real(8) :: tlast_i, mother1_i, mother2_i

    ios = 0
    is_milne_line = .false.
    is_oscar_line = .false.
    pdg_meta = 0.0d0
    ncoll_meta = -1.0d0

    ! OSCAR2013Extended SMASH_IC format:
    ! t x y z mass p0 px py pz pdg ID charge ncoll form_time xsecfac
    ! proc_id_origin proc_type_origin time_last_coll pdg_mother1 pdg_mother2
    ! baryon_number strangeness
    read(line, *, iostat=ios) t_i, x_i, y_i, z_i, mass_i, e_i, px_i, py_i, pz_i, &
                              pdg_i, id_i, charge_i, ncoll_i, form_i, xsecfac_i, &
                              procid_i, proctype_i, tlast_i, mother1_i, mother2_i, &
                              baryon_i, strange_i
    if (ios == 0) then
      c1 = t_i
      c2 = x_i
      c3 = y_i
      c4 = z_i
      px = px_i
      py = py_i
      pz = pz_i
      ep = e_i
      baryon = baryon_i
      pdg_meta = pdg_i
      ncoll_meta = ncoll_i
      is_oscar_line = .true.
      return
    end if

    ! SMASH header-described format:
    ! tau x y eta mt px py Rap pdg charge baryon_number strangeness
    read(line, *, iostat=ios) tau_i, x_i, y_i, eta_i, mt_i, px_i, py_i, rap_i, &
                              pdg_i, charge_i, baryon_i, strange_i
    if (ios == 0) then
      c1 = tau_i
      c2 = eta_i
      c3 = x_i
      c4 = y_i
      px = px_i
      py = py_i
      pz = mt_i * sinh(rap_i)
      ep = mt_i * cosh(rap_i)
      baryon = baryon_i
      pdg_meta = pdg_i
      is_milne_line = .true.
      return
    end if

    ! Backward-compatible fallback:
    ! tau eta x y px py pz E baryon
    ios = 0
    read(line, *, iostat=ios) c1, c2, c3, c4, px, py, pz, ep, baryon
    if (ios == 0) is_milne_line = .true.
  end subroutine parse_smash_particle_line


  logical function is_oscar_spectator(ep, pz, pdg, baryon, ncoll, ycut)
    real(8), intent(in) :: ep, pz, pdg, baryon, ncoll, ycut
    real(8) :: y

    is_oscar_spectator = .false.
    if (abs(baryon) < 0.5d0) return
    if (ncoll > 0.5d0) return

    if (abs(pdg) /= 2212.0d0 .and. abs(pdg) /= 2112.0d0) return

    y = momentum_rapidity(ep, pz)
    if (abs(y) >= ycut) is_oscar_spectator = .true.
  end function is_oscar_spectator


  real(8) function momentum_rapidity(ep, pz)
    real(8), intent(in) :: ep, pz
    real(8), parameter :: tiny = 1.0d-12
    real(8) :: num, den

    num = ep + pz
    den = ep - pz
    if (num <= tiny .or. den <= tiny) then
      if (pz >= 0.0d0) then
        momentum_rapidity = 1.0d3
      else
        momentum_rapidity = -1.0d3
      end if
      return
    end if
    momentum_rapidity = 0.5d0 * log(num / den)
  end function momentum_rapidity


  subroutine print_mapping_diagnostics(n_pre, n_post, &
                                       e_pre, px_pre, py_pre, pz_pre, b_pre, &
                                       ex_pre, ey_pre, ez_pre, et_pre, &
                                       e_post, px_post, py_post, pz_post, b_post, &
                                       ex_post, ey_post, ez_post, et_post)
    integer, intent(in) :: n_pre, n_post
    real(8), intent(in) :: e_pre, px_pre, py_pre, pz_pre, b_pre
    real(8), intent(in) :: ex_pre, ey_pre, ez_pre, et_pre
    real(8), intent(in) :: e_post, px_post, py_post, pz_post, b_post
    real(8), intent(in) :: ex_post, ey_post, ez_post, et_post
    real(8), parameter :: tiny = 1.0d-20
    real(8) :: xpre, ypre, zpre, tpre, xpost, ypost, zpost, tpost

    write(*,'(A)') '--- Fluidization mapping diagnostics (pre -> post) ---'
    write(*,'(A,I10,A,I10)') 'Particles: ', n_pre, ' -> ', n_post
    write(*,'(A,1PE14.6,A,1PE14.6)') 'sum(E): ', e_pre, ' -> ', e_post
    write(*,'(A,1PE14.6,A,1PE14.6)') 'sum(Px): ', px_pre, ' -> ', px_post
    write(*,'(A,1PE14.6,A,1PE14.6)') 'sum(Py): ', py_pre, ' -> ', py_post
    write(*,'(A,1PE14.6,A,1PE14.6)') 'sum(Pz): ', pz_pre, ' -> ', pz_post
    write(*,'(A,1PE14.6,A,1PE14.6)') 'sum(B): ', b_pre, ' -> ', b_post

    xpre = 0.0d0
    ypre = 0.0d0
    zpre = 0.0d0
    tpre = 0.0d0
    xpost = 0.0d0
    ypost = 0.0d0
    zpost = 0.0d0
    tpost = 0.0d0
    if (abs(e_pre) > tiny) then
      xpre = ex_pre / e_pre
      ypre = ey_pre / e_pre
      zpre = ez_pre / e_pre
      tpre = et_pre / e_pre
    end if
    if (abs(e_post) > tiny) then
      xpost = ex_post / e_post
      ypost = ey_post / e_post
      zpost = ez_post / e_post
      tpost = et_post / e_post
    end if

    write(*,'(A,1PE14.6,A,1PE14.6)') '<x>_E: ', xpre, ' -> ', xpost
    write(*,'(A,1PE14.6,A,1PE14.6)') '<y>_E: ', ypre, ' -> ', ypost
    write(*,'(A,1PE14.6,A,1PE14.6)') '<z>_E: ', zpre, ' -> ', zpost
    write(*,'(A,1PE14.6,A,1PE14.6)') '<t>_E: ', tpre, ' -> ', tpost
    write(*,'(A)') '------------------------------------------------------'
  end subroutine print_mapping_diagnostics


  subroutine print_yellow_mapping_warning(target_t_eff)
    real(8), intent(in) :: target_t_eff
    character(len=*), parameter :: c_yellow = achar(27)//'[33m'
    character(len=*), parameter :: c_reset = achar(27)//'[0m'

    write(*,'(A)') c_yellow// &
      'WARNING: SMASH->hydro mapping to constant t can alter the spatial structure of the initial state.'
    if (target_t_eff > 0.0d0) then
      write(*,'(A,1PE14.6)') 'WARNING: target t used for synchronization [fm/c] = ', target_t_eff
    else
      write(*,'(A)') 'WARNING: in multi-event mode, each event may use a different target t.'
    end if
    write(*,'(A)') 'WARNING: review pre/post mapping diagnostics above before physics comparisons.'//c_reset
  end subroutine print_yellow_mapping_warning


  subroutine print_target_t_clamp_warning(t_raw, t_used, tmin, tmax)
    real(8), intent(in) :: t_raw, t_used, tmin, tmax
    character(len=*), parameter :: c_yellow = achar(27)//'[33m'
    character(len=*), parameter :: c_reset = achar(27)//'[0m'

    write(*,'(A)') c_yellow// &
      'WARNING: auto target_t from particle max(t) was clamped for hydro safety.'
    write(*,'(A,1PE14.6,A,1PE14.6)') 'WARNING: raw/used target_t [fm/c] = ', t_raw, ' -> ', t_used
    write(*,'(A,1PE14.6,A,1PE14.6)') 'WARNING: clamp window [tmin,tmax] [fm/c] = ', tmin, ' / ', tmax
    write(*,'(A)') 'WARNING: set target_t explicitly if you need a different initialization time.'//c_reset
  end subroutine print_target_t_clamp_warning


  subroutine print_time_drop_warning(n_drop, n_total, frac, frac_lim, stop_on_limit, ierr)
    integer, intent(in) :: n_drop, n_total
    real(8), intent(in) :: frac, frac_lim
    logical, intent(in) :: stop_on_limit
    integer, intent(inout) :: ierr
    character(len=*), parameter :: c_yellow = achar(27)//'[33m'
    character(len=*), parameter :: c_red = achar(27)//'[31m'
    character(len=*), parameter :: c_reset = achar(27)//'[0m'

    write(*,'(A)') c_yellow// &
      'WARNING: particles with ti>target_t were removed before fluidization.'
    write(*,'(A,I10,A,I10,A,F8.3,A)') 'WARNING: dropped/considered = ', n_drop, ' / ', n_total, &
      ' (', 100.0d0*frac, '%)'
    write(*,'(A,F8.3,A)') 'WARNING: configured drop-fraction limit = ', 100.0d0*frac_lim, '%'//c_reset

    if (stop_on_limit .and. frac > frac_lim) then
      write(*,'(A)') c_red// &
        'ERROR: time-drop fraction exceeds configured limit; aborting initialization.'//c_reset
      ierr = 2005
    end if
  end subroutine print_time_drop_warning


  real(8) function gaussian_1d(dx, sigma)
    real(8), intent(in) :: dx, sigma
    real(8), parameter :: tiny = 1.0d-12

    if (sigma <= tiny) then
      if (abs(dx) <= tiny) then
        gaussian_1d = 1.0d0
      else
        gaussian_1d = 0.0d0
      end if
      return
    end if

    gaussian_1d = exp(-0.5d0 * (dx / sigma) * (dx / sigma))
  end function gaussian_1d

end module fluidization_smash
