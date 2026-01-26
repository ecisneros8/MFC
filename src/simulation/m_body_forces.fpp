#:include 'macros.fpp'

module m_body_forces

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion

    use m_nvtx

! $:USE_GPU_MODULE()

    implicit none

    private; 
    public :: s_compute_body_forces_rhs, &
              s_initialize_body_forces_module, &
              s_finalize_body_forces_module

    integer, parameter :: spbf_num_freq = 8
    real(wp) :: spbf_amp
    real(wp) :: spbf_xc
    real(wp) :: spbf_yc
    real(wp) :: spbf_conv_vel
    real(wp) :: spbf_sigma
    real(wp), allocatable, dimension(:) :: freq, phase
    real(wp), allocatable, dimension(:, :, :) :: rhoM
    $:GPU_DECLARE(create='[spbf_amp,spbf_xc,spbf_yc,spbf_conv_vel,spbf_sigma,freq,phase,rhoM]')

contains

    !> This subroutine initializes the module global array of mixture
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

	if (bf_spatial_support) then
	   call s_initialize_body_force_with_spatial_support
	endif
    end subroutine s_initialize_body_forces_module

    !> This subroutine initializes a body force with spatial
    !! support presented in Wei & Freund (JFM, 2005)
    subroutine s_initialize_body_force_with_spatial_support

    	   integer :: f !< frequency iterator

	   spbf_amp = spatial_bf%amp
	   spbf_xc = spatial_bf%x_centroid
	   spbf_yc = spatial_bf%y_centroid
	   spbf_conv_vel = spatial_bf%conv_vel
	   spbf_sigma = spatial_bf%sigma

    	   @:ALLOCATE(freq(spbf_num_freq), phase(spbf_num_freq))
    	   @:PREFER_GPU(freq)
    	   @:PREFER_GPU(phase)
    	   do f = 1, spbf_num_freq
	      freq(f) = spatial_bf%freq(f)
	      phase(f) = spatial_bf%phase(f)
	   enddo
	   $:GPU_UPDATE(device='[freq,phase]')

    end subroutine s_initialize_body_force_with_spatial_support

    !> This subroutine computes the acceleration at time t
    subroutine s_compute_acceleration(t)

        real(wp), intent(in) :: t

        #:for DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (bf_${XYZ}$) then
                accel_bf(${DIR}$) = g_${XYZ}$+k_${XYZ}$*sin(w_${XYZ}$*t - p_${XYZ}$)
            end if
        #:endfor

        $:GPU_UPDATE(device='[accel_bf]')

    end subroutine s_compute_acceleration

    !> This routine applies the body force of Wei & Freund (JFM, 2005)
    subroutine s_compute_body_force_with_spatial_support(t)

        real(wp), intent(in) :: t
	real(wp) :: support !< spatial support
	real(wp) :: theta_x, theta_y, pre_fac !< auxiliary variables
	integer :: f !< frequency iterator
	integer :: i, j, k, l !< standard iterators

	$:GPU_PARALLEL_LOOP(private='[support,theta_x,theta_y,pre_fac,f,j,k,l]', collapse=3)
	do l = 0, p
	   do k = 0, n
	      do j = 0, m
	         support = exp(-spbf_sigma * &
		 	   ((x_cc(j) - spbf_xc)**2 + (y_cc(k) - spbf_yc)**2))
		 spbf_source_x(j, k, l) = 0._wp
		 spbf_source_y(j, k, l) = 0._wp
		 do f = 1, spbf_num_freq
		    pre_fac = (freq(f) / spbf_conv_vel)
		    theta_x = pre_fac * &
		    	      (x_cc(j) - spbf_xc - spbf_conv_vel * t) + &
			      phase(f)
		    theta_y = pre_fac * (y_cc(k) - spbf_yc) + phase(f)
		    spbf_source_x(j, k, l) = spbf_source_x(j, k, l) + &
		     		      	    spbf_amp * support * ( &
					    pre_fac * sin(theta_x) * cos(theta_y) - &
					    2 * spbf_sigma * (y_cc(k) - spbf_yc) * &
					    sin(theta_x) * sin(theta_y))
		     spbf_source_y(j, k, l) = spbf_source_y(j, k, l) - &
		     		     	     spbf_amp * support * ( &
					     pre_fac * cos(theta_x) * sin(theta_y) - &
					     2 * spbf_sigma * (x_cc(j) - spbf_xc) * &
					     sin(theta_x) * sin(theta_y))
		 end do
	      end do
	   end do
	end do
	$:END_GPU_PARALLEL_LOOP()
    end subroutine s_compute_body_force_with_spatial_support

    !> This subroutine calculates the mixture density at each cell
    !! center
    !! param q_cons_vf Conservative variable
    subroutine s_compute_mixture_density(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer :: i, j, k, l !< standard iterators

        $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
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
        $:END_GPU_PARALLEL_LOOP()

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

	if (bf_spatial_support) then
		call s_compute_body_force_with_spatial_support(mytime)
	endif

        call s_compute_mixture_density(q_cons_vf)

        $:GPU_PARALLEL_LOOP(private='[i,j,k,l]', collapse=4)
        do i = momxb, E_idx
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

	if (bf_spatial_support) then

	   $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
	   do l = 0, p
	      do k = 0, n
	      	 do j = 0, m
		    rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
		    			      	spbf_source_x(j, k, l)
		    rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
		    		   	       	    spbf_source_y(j, k, l)
		    ! write(*, *) j, k, l, spatial_bf_x(j, k, l)
		 enddo
	      enddo
	   enddo
	   $:END_GPU_PARALLEL_LOOP()
	endif

        if (bf_x) then ! x-direction body forces

            $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
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
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (bf_y) then ! y-direction body forces

            $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
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
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (bf_z) then ! z-direction body forces

            $:GPU_PARALLEL_LOOP(private='[j,k,l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(momxe)%sf(j, k, l) = rhs_vf(momxe)%sf(j, k, l) + &
                                                    rhoM(j, k, l)*accel_bf(3)
                        rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                    q_cons_vf(momxe)%sf(j, k, l)*accel_bf(3)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end if

    end subroutine s_compute_body_forces_rhs

    impure subroutine s_finalize_body_forces_module

        @:DEALLOCATE(rhoM)

        if (bf_spatial_support) then
            @:DEALLOCATE(freq, phase)
        end if

    end subroutine s_finalize_body_forces_module

end module m_body_forces
