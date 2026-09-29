double precision function adjoint_initial(xp,yp)
    use adjoint_module, only : beta, x0, y0, initial_condition
    implicit none

    double precision :: xp, yp
    double precision :: r2, q0

    q0 = 0

    if (initial_condition .eq. 0) then 
        !! Zero initial condition
        q0 = 0
    else if (initial_condition .eq. 1) then
        !! Gaussian initial condition
        r2 = (xp - x0)**2 + (yp - y0)**2
        q0 = exp(-beta*r2)
    endif

    adjoint_initial = q0

end function adjoint_initial