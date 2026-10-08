module adjoint_module
    implicit none

    integer :: source_model
    integer :: initial_condition
    integer :: mdata !!number of data points
    integer :: obs_index = 1 !! current observation used to build representer

    logical :: pseudo_1d_experiment !! .true.: pseudo-1D experiment; .false.: full-2D

    double precision :: eps_1d !! width of the heat kernel for pseudo-1D experiment
    double precision :: eps_2d !! width of the heat kernel for full-2D experiment
    double precision, dimension(:), allocatable :: xm, ym, tm !! spatial and temp. positions

    double precision :: beta, x0, y0

    double precision, dimension(:), allocatable :: dm !! data observation
    double precision, dimension(:), allocatable :: W_eps !! data weights
    double precision :: tfinal !! Final time for simulation

    
end module adjoint_module


double precision function dirac_delta(r)
    use adjoint_module, only : eps_1d, eps_2d, pseudo_1d_experiment
    implicit none

    double precision r, pi

    pi = 4.d0*atan(1.d0)

    if (pseudo_1d_experiment) then
        dirac_delta = exp(-r**2/(4*eps_1d))/sqrt(4*pi*eps_1d) ! 1D delta function approx. using heat kernel

    else
        dirac_delta = exp(-r**2/(4*eps_2d))/(4*pi*eps_2d) ! 2D delta function approx. using heat kernel
    
    endif
end function dirac_delta


double precision function dirac_delta_time(s)
    !! Temporal delta approximation. Time is one-dimensional regardless of
    !! the spatial dimension of the experiment, so the 1D heat-kernel
    !! normalization 1/sqrt(4*pi*eps) is always used (the 2D prefactor
    !! 1/(4*pi*eps) would inject a total impulse of 1/sqrt(4*pi*eps)
    !! instead of unity when integrated over time).
    use adjoint_module, only : eps_1d, eps_2d, pseudo_1d_experiment
    implicit none

    double precision s, pi, eps

    pi = 4.d0*atan(1.d0)

    if (pseudo_1d_experiment) then
        eps = eps_1d
    else
        eps = eps_2d
    endif

    dirac_delta_time = exp(-s**2/(4*eps))/sqrt(4*pi*eps)

end function dirac_delta_time



subroutine data(m)
    use adjoint_module
    implicit none

    integer m
    mdata = m
    
    !if (.not. allocated(xm)) allocate(xm(m), ym(m), tm(m))


    if (allocated(xm)) deallocate(xm, ym, tm, dm, W_eps)
    allocate(xm(m), ym(m), tm(m), dm(m), W_eps(m))

end subroutine data

subroutine set_observation_index(idx)
    use adjoint_module, only : obs_index, mdata
    implicit none

    integer, intent(in) :: idx

    if (idx < 1 .or. idx > mdata) then
        write(6,*) 'set_observation_index: invalid idx = ', idx
        stop
    endif

    obs_index = idx

end subroutine set_observation_index
