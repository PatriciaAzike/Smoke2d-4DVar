SUBROUTINE smoke3d_update_gauge (blockno, mx,my,mbc,meqn, & 
    xlower,ylower, dx,dy,q,xc,yc,qvar)

    implicit none

    integer :: blockno, mx, my, mbc, meqn
    double precision :: xlower, ylower, dx, dy
    double precision :: q(1-mbc:mx+mbc, 1-mbc:my+mbc,meqn)
    double precision :: qvar

    !! local variables:
    double precision :: xc,yc, a00, a01, a10, a11
    double precision :: quad(0:1,0:1)
    integer icell, jcell, ii, jj

    !! get dual cell containing (xc,yc)
    icell =  int((xc-(xlower + dx/2))/dx) + 1
    jcell =  int((yc-(ylower + dy/2))/dy) + 1

    do ii = 0,1
        do jj = 0,1
            quad(ii,jj) = q(icell+ii,jcell+jj,1)
        end do
    end do

    a00 = quad(0,0)
    a01 = quad(1,0) - quad(0,0)
    a10 = quad(0,1) - quad(0,0)
    a11 = quad(1,1) - quad(1,0) - quad(0,1) + quad(0,0)

    qvar = a00 + a01*xc + a10*yc + a11*xc*yc

!!    icell =  int((xc-xlow_long)/dx) + 1
!!    jcell =  int((yc-ylow_lat)/dy) + 1
!!
!!    xcenter  = xlow_long + (icell-0.5)*dx
!!    ycenter  = ylow_lat + (jcell-0.5)*dy
!!    xoff   = (x_long-xcenter)/dx
!!    yoff   = (y_lat-ycenter)/dy

!!    !! # Dry cell ? 
!!    h(1) = q(1,iindex,jindex)
!!    h(2) = q(1,iindex+1,jindex)
!!    h(3) = q(1,iindex,jindex+1)
!!    h(4) = q(1,iindex+1,jindex+1)

!!    !! Linear interpolation between four cells
!!    mq = 1
!!    qvar = (1.d0 - xoff) * (1.d0 - yoff) &
!!            * q(icell,jcell,mq)  &
!!            + xoff*(1.d0 - yoff) * q(icell+1,jcell,mq)  &
!!            + (1.d0 - xoff) * yoff * q(icell,jcell+1,mq)  &
!!            + xoff * yoff * q(icell+1,jcell+1,mq)
!!

    if (abs(qvar) < 1d-99) then
        qvar = 0.d0
    endif

END SUBROUTINE smoke3d_update_gauge
