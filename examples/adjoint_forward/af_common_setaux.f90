subroutine af_common_setaux(maxmx,maxmy,mbc,mx,my, & 
    xlower,ylower,dx,dy,maux,aux,psi)
    implicit none

    integer mbc, mx, my, maux, maxmx, maxmy
    double precision xlower, ylower, dx, dy
    double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

    integer i, j
    double precision xll, yll, psi
    logical fwave

    fwave = .true.

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc

            if (fwave) then

                !! coordinates of lower left corner of grid cell:
                xll = xlower + (i-0.5d0)*dx
                yll = ylower + (j-1)*dy
                aux(i,j,1) = (psi(xll, yll+dy) - psi(xll,yll)) / dy

                xll = xlower + (i-1)*dx
                yll = ylower + (j-0.5d0)*dy
                aux(i,j,2) = -(psi(xll+dx, yll) - psi(xll,yll)) / dx

            else
                xll = xlower + (i-1)*dx
                yll = ylower + (j-1)*dy

                !! difference stream function psi to get normal velocities:
                aux(i,j,1) = (psi(xll, yll+dy) - psi(xll,yll)) / dy
                aux(i,j,2) = -(psi(xll+dx, yll) - psi(xll,yll)) / dx

            endif
        enddo
    enddo

    return
end  subroutine af_common_setaux



subroutine forward_setaux(maxmx,maxmy,mbc,mx,my, & 
    xlower,ylower,dx,dy,maux,aux)   
    implicit none

    external forward_psi
    integer mbc, mx, my, maux, maxmx, maxmy
    double precision xlower, ylower, dx, dy
    double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

    double precision forward_psi

    call af_common_setaux(maxmx,maxmy,mbc,mx,my, & 
                          xlower,ylower,dx,dy,maux,aux, & 
                          forward_psi)

    return
end subroutine forward_setaux


subroutine adjoint_setaux(maxmx,maxmy,mbc,mx,my, & 
    xlower,ylower,dx,dy,maux,aux)   
    implicit none

    external adjoint_psi
    integer mbc, mx, my, maux, maxmx, maxmy
    double precision xlower, ylower, dx, dy
    double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

    double precision adjoint_psi


    call af_common_setaux(maxmx,maxmy,mbc,mx,my, & 
                          xlower,ylower,dx,dy,maux,aux, &
                          adjoint_psi)

    return
end subroutine adjoint_setaux


