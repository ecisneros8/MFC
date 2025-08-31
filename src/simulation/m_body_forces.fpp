#:include 'macros.fpp'

module m_body_forces

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion

    use m_nvtx

#ifdef MFC_OpenACC
    use openacc
#endif

    implicit none

    private; 
    public :: s_compute_body_forces_rhs, &
              s_initialize_body_forces_module, &
              s_finalize_body_forces_module

    integer, parameter :: spbf_num_freq = 8
    real(wp), allocatable, dimension(:) :: freq, phase
    real(wp), allocatable, dimension(:, :, :) :: rhoM

contains

    !> This subroutine inializes the module global array of mixture
    !! densities in each grid cell
    impure subroutine s_initialize_body_forces_module

        integer :: f !< frequency iterator
	real(wp) :: rf, rp
    	real(wp), parameter :: pi = 4._wp * atan(1._wp)
    
        ! Simulation is at least 2D
        if (n > 0) then
            ! Simulation is 3D
            if (p > 0) then
                @:ALLOCATE (rhoM(-buff_size:buff_size + m, &
                    -buff_size:buff_size + n, &
                    -buff_size:buff_size + p))
                ! Simulation is 2D
            else
                @:ALLOCATE (rhoM(-buff_size:buff_size + m, &
                    -buff_size:buff_size + n, &
                    0:0))
            end if
            ! Simulation is 1D
        else
            @:ALLOCATE (rhoM(-buff_size:buff_size + m, &
                0:0, &
                0:0))
        end if

	! Initialize the parameters for the body force with spatial dependence
	if (spatial_bf) then
		@:ALLOCATE(freq(spbf_num_freq), phase(spbf_num_freq))
		write(*, *) '---', spbf_amp, spbf_freq, spbf_xc, spbf_yc, spbf_conv_vel, spbf_sigma
		do f = 1, spbf_num_freq
	   	   call random_number(rf)
	   	   call random_number(rp)
	  	   rf = rf - 0.5 * 1._wp
	   	   rp = 2._wp * pi * rp 
	   	   freq(f) = 0.25 * spbf_freq * (f * 1._wp + rf) * pi
		   phase(f) = rp
		enddo
		write(*, *) freq(:)
		write(*, *) phase(:)
	endif
    end subroutine s_initialize_body_forces_module

    !> This subroutine computes the acceleration at time t
    subroutine s_compute_acceleration(t)

        real(wp), intent(in) :: t

        if (m > 0) then
            accel_bf(1) = g_x + k_x*sin(w_x*t - p_x)
            if (n > 0) then
                accel_bf(2) = g_y + k_y*sin(w_y*t - p_y)
                if (p > 0) then
                    accel_bf(3) = g_z + k_z*sin(w_z*t - p_z)
                end if
            end if
        end if

        $:GPU_UPDATE(device='[accel_bf]')

    end subroutine s_compute_acceleration

    !> This routine applies the body force of Wei & Freund (JFM, 2005)
    subroutine s_compute_body_force_with_spatial_support(t)


        real(wp), intent(in) :: t
	real(wp) :: support !< spatial support
	real(wp) :: theta_x, theta_y, pre_fac !< auxiliary variables
	integer :: f !< frequency iterator
	integer :: i, j, k, l !< standard iterators

	$:GPU_PARALLEL_LOOP(collapse=4)
	do l = 0, p
	   do k = 0, n
	      do j = 0, m
	         support = exp(-spbf_sigma * &
		 	   ((x_cc(j) - spbf_xc)**2 + (y_cc(k) - spbf_yc)**2))
		 spatial_bf_x(j, k, l) = 0._wp
		 spatial_bf_y(j, k, l) = 0._wp
		 do f = 1, spbf_num_freq
		    pre_fac = (freq(f) / spbf_conv_vel)
		    theta_x = pre_fac * &
		    	      (x_cc(j) - spbf_xc - spbf_conv_vel * t) + &
			      phase(f)
		    theta_y = (freq(f) / spbf_conv_vel) * y_cc(k) + phase(f)
		    spatial_bf_x(j, k, l) = spatial_bf_x(j, k, l) + &
		     		      	    spbf_amp * support * ( &
					    pre_fac * cos(theta_x) * sin(theta_y) - &
					    2 * spbf_sigma * (x_cc(j) - spbf_xc) * &
					    sin(theta_x) * sin(theta_y))
		     spatial_bf_y(j, k, l) = spatial_bf_y(j, k, l) - &
		     		     	     spbf_amp * support * ( &
					     pre_fac * sin(theta_x) * cos(theta_y) - &
					     2 * spbf_sigma * (y_cc(k) - spbf_yc) * &
					     sin(theta_x) * sin(theta_y))
		 end do
	      end do
	   end do
	end do
    end subroutine s_compute_body_force_with_spatial_support

    !> This subroutine calculates the mixture density at each cell
    !! center
    !! param q_cons_vf Conservative variable
    subroutine s_compute_mixture_density(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer :: i, j, k, l !< standard iterators

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rhoM(j, k, l) = 0._wp
                    do i = 1, num_fluids
                        rhoM(j, k, l) = rhoM(j, k, l) + &
                                        q_cons_vf(contxb + i - 1)%sf(j, k, l)
                    end do
                end do
            end do
        end do

    end subroutine s_compute_mixture_density

    !> This subroutine calculates the source term due to body forces
    !! so the system can be advanced in time
    !! @param q_cons_vf Conservative variables
    !! @param q_prim_vf Primitive variables
    subroutine s_compute_body_forces_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        integer :: i, j, k, l !< Loop variables

	if (bf_x .or. bf_y .or. bf_z) then
	        call s_compute_acceleration(mytime)
	endif

	if (spatial_bf) then
		call s_compute_body_force_with_spatial_support(mytime)
	endif

        call s_compute_mixture_density(q_cons_vf)

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = momxb, E_idx
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do

	if (spatial_bf) then

	   $:GPU_PARALLEL_LOOP(collapse=3)
	   do l = 0, p
	      do k = 0, n
	      	 do j = 0, m
		    rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
		    			      	rhoM(j, k, l) * spatial_bf_x(j, k, l)
		    rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
		    		   	       	    rhoM(j, k, l) * spatial_bf_y(j, k, l)
		    ! write(*, *) j, k, l, spatial_bf_x(j, k, l)
		 enddo
	      enddo
	   enddo
	endif

        if (bf_x) then ! x-direction body forces

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                    rhoM(j, k, l)*accel_bf(1)
                        rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                    q_cons_vf(momxb)%sf(j, k, l)*accel_bf(1)
                    end do
                end do
            end do
        end if

        if (bf_y) then ! y-direction body forces

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                        rhoM(j, k, l)*accel_bf(2)
                        rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                    q_cons_vf(momxb + 1)%sf(j, k, l)*accel_bf(2)
                    end do
                end do
            end do
        end if

        if (bf_z) then ! z-direction body forces

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(momxe)%sf(j, k, l) = rhs_vf(momxe)%sf(j, k, l) + &
                                                    (rhoM(j, k, l))*accel_bf(3)
                        rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                    q_cons_vf(momxe)%sf(j, k, l)*accel_bf(3)
                    end do
                end do
            end do

        end if

    end subroutine s_compute_body_forces_rhs

    impure subroutine s_finalize_body_forces_module

        @:DEALLOCATE(rhoM)

    end subroutine s_finalize_body_forces_module

end module m_body_forces
