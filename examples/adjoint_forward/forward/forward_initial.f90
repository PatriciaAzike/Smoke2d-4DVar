double precision function forward_initial(xp,yp)
!!forward_initial(xp,yp,mx,my,mbc,maux,aux,t)
    use forward_module, only : beta, initial_condition, x0, y0
    implicit none

    double precision :: xp, yp
    double precision :: r2, q0,t
    integer mbc, mx, my, maux

    !!double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

    integer i, j

    q0 = 0

    if (initial_condition .eq. 0) then 
        !! Zero initial condition
        q0 = 0
    else if (initial_condition .eq. 1) then
        !! Gaussian initial condition
        r2 = (xp - x0)**2 + (yp - y0)**2
        q0 = exp(-beta*r2)

    !!else if (initial_condition .eq. 2) then
    !!    do i = 1-mbc,mx+mbc
    !!        do j = 1-mbc,my+mbc
    !!            q0 = aux(i,j,3)
    !!        enddo
    !!    enddo
    endif

    forward_initial = q0

end function forward_initial