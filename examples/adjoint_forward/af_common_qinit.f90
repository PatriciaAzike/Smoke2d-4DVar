subroutine af_common_qinit(maxmx,maxmy,meqn,mbc, & 
   mx,my,xlower,ylower,dx,dy,q,maux,aux, & 
   initial_condition,q0_initial)
    use adjoint_module, only : pseudo_1d
    implicit none

    external q0_initial
    integer maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision xlower, ylower, dx, dy
    double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
    double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)
    integer initial_condition

    double precision xc, yc, xlow,ylow,w
    double precision q0, q0_initial
    integer i, j, mq
    integer blockno, fc2d_clawpack46_get_block

    blockno = fc2d_clawpack46_get_block()

    do mq = 1,meqn
        do i = 1-mbc,mx+mbc
            xc = xlower + (i-0.5d0)*dx
            xlow = xlower + (i-1)*dx
            do j = 1-mbc,my+mbc
                yc = ylower + (j-0.5d0)*dy
                ylow = ylower + (j-1)*dy

                if (initial_condition .eq. 2) then
                    if (pseudo_1d .eq. 1) then
                        if (abs(xc - 0.5d0) .lt. 0.25d0) then 
                            w = 1
                        else 
                            w = 0
                        endif
                    elseif (pseudo_1d .eq. 2) then
                        call cellave2(blockno,xlow,ylow,dx,dy,w)
                    endif

                    q0 = w
                else
                    q0 = q0_initial(xc,yc)
                endif
                q(i,j,mq) = q0
         enddo
      enddo
    enddo

    return
end subroutine af_common_qinit


subroutine forward_qinit(maxmx,maxmy,meqn,mbc, & 
   mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use forward_module, only : initial_condition
    implicit none

    external forward_initial
    integer maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision xlower, ylower, dx, dy
    double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
    double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

    double precision forward_initial

    call af_common_qinit(maxmx,maxmy,meqn,mbc, & 
                         mx,my,xlower,ylower,dx,dy,q, & 
                         maux,aux,initial_condition, &
                         forward_initial)

end subroutine forward_qinit


subroutine  adjoint_qinit(maxmx,maxmy,meqn,mbc, & 
   mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use  adjoint_module, only : initial_condition
    implicit none

    external adjoint_initial
    integer maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision xlower, ylower, dx, dy
    double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
    double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

    double precision adjoint_initial

    call af_common_qinit(maxmx,maxmy,meqn,mbc, & 
                         mx,my,xlower,ylower,dx,dy,q, & 
                         maux,aux,initial_condition, & 
                         adjoint_initial)

end subroutine  adjoint_qinit
