module smoke2d_module
    implicit none
      

    double precision :: pi, pi2

    !! Domain parameters
    double precision, dimension(2) :: latitude, longitude

    !! Smoke source term data
    double precision :: radius_sphere

    !! Observation data
    double precision, dimension(:), allocatable :: xm, ym, tm

    !! Delta function parameter
    double precision :: eps
    integer, parameter :: n_points = 1
    !double precision, allocatable :: uv(:,:)

end module smoke2d_module

subroutine generate_data(m)
    use smoke2d_module
    implicit none

    integer m

    allocate(xm(m),ym(m),tm(m))

    if (m > 1) then
        write(6,*) "Only one observations allowed!"
        stop
    endif

    xm(1) = -110.1
    ym(1) = 43.5
    tm(1) = 1.0

end subroutine generate_data    

double precision function delta(r)
    use smoke2d_module, only : pi, eps
    implicit none

    double precision r

    delta = exp(-r**2/(4*eps))/(4*pi*eps)    
    

end function  delta


