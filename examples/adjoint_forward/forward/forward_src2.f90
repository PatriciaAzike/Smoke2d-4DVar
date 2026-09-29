subroutine forward_src2(maxmx,maxmy,meqn,mbc,mx,my, & 
    xlower,ylower,dx,dy,r,maux,aux,t,dt)
!!     =====================================================

!!       # Set initial conditions for q.
!!       # Sample scalar equation with data that is piecewise constant with
!!       # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
!!       #     0.1  otherwise
!!
!!       # After the overlap exchange, the aux routine contains the adjoint data

   use forward_module,only : W_f
   implicit none

   integer maxmx, maxmy, meqn, mbc, mx, my, maux
   double precision xlower, ylower, dx, dy,t, dt
   double precision r(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
   double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)
   double precision S, alpha

   integer i, j, mq

   do mq = 1,meqn
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            !! Let alpha be the adjoint solution stored in the aux array.
            alpha = aux(i,j,2+mq)   !! Adjoint solution
            S = (1.d0/W_f) *alpha

            r(i,j,mq) = r(i,j,mq) + dt*S

         enddo
      enddo
   enddo

   return
end subroutine forward_src2


