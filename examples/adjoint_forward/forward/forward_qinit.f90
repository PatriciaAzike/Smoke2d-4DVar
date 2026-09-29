subroutine forward_qinit(maxmx,maxmy,meqn,mbc, & 
   mx,my,xlower,ylower,dx,dy,q,maux,aux)

    !!use forward_module
    use forward_module, only: W_i
    implicit none

    integer maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision xlower, ylower, dx, dy
    double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
    double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)
    integer i, j, mq

    double precision xc, yc
    double precision q0, forward_initial

    do mq = 1,meqn
        do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc
                xc = xlower + (i-0.5d0)*dx
                yc = ylower + (j-0.5d0)*dy

                q0 = forward_initial(xc,yc)
                !!q0 = (1/W_i) * forward_initial(xc,yc,mx,my,mbc,maux,aux,0)
                q(i,j,mq) = q0
            enddo
        enddo
    enddo

    return
end subroutine forward_qinit