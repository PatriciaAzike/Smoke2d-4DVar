subroutine adjoint_src2(maxmx,maxmy,meqn,mbc,mx,my, & 
    xlower,ylower,dx,dy,q,maux,aux,t,dt)
!!     =====================================================

!!       # Set initial conditions for q.
!!       # Sample scalar equation with data that is piecewise constant with
!!       # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
!!       #     0.1  otherwise

   use adjoint_module
   implicit none

   integer maxmx, maxmy, meqn, mbc, mx, my, maux
   double precision xlower, ylower, dx, dy,t, dt
   double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
   double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)
   double precision S, adjoint_source

   integer i, j, mq
   double precision xc,yc

   do mq = 1,meqn
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            
            !! # Cell center, in brick units     
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            S = adjoint_source(mq,xc,yc,t)
            q(i,j,mq) = q(i,j,mq) + dt*S
            
         enddo
      enddo
   enddo

   return
end subroutine adjoint_src2