subroutine adjoint_clawpack46_qinit(maxmx,maxmy,meqn,mbc, & 
   mx,my,xlower,ylower,dx,dy,q,maux,aux)

   use adjoint_module, only : initial_condition
   implicit none

   integer maxmx, maxmy, meqn, mbc, mx, my, maux
   double precision xlower, ylower, dx, dy
   double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
   double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

   integer i, j, mq
   double precision xc, yc, xlow,ylow,w

   double precision q0, adjoint_initial
   integer blockno, fc2d_clawpack46_get_block

   blockno = fc2d_clawpack46_get_block()

   do mq = 1,meqn
        do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            xc = xlower + (i-0.5d0)*dx
            do j = 1-mbc,my+mbc 
                yc = ylower + (j-0.5d0)*dy
                ylow = ylower + (j-1)*dy

                if (initial_condition .eq. 2)
                    call cellave2(blockno,xlow,ylow,dx,dy,w)
                    q0 = w
                else
                    q0 = adjoint_initial(xc,yc)
                endif
                q(i,j,mq) = q0
            end do
        end do
    end do

    return
end subroutine adjoint_clawpack46_qinit