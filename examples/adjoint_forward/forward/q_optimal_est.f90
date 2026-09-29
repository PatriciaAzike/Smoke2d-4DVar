subroutine updateq_from_representer(mx,my,mbc,meqn,q_rep,q_model,beta_j)

!!     =====================================================

!!       # q_hat = q_F + beta*R_m
!!       # We update q_F here because we want to access glob from C
   
   implicit none

   integer meqn, mbc, mx, my,k
   double precision :: q_rep(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
   double precision :: q_model(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
   
   integer i, j, mq

   double precision :: beta_j


   do mq = 1, meqn
      do j = 1-mbc, my+mbc
         do i = 1-mbc, mx+mbc
            q_model(i,j,mq) = q_model(i,j,mq) + beta_j*q_rep(i,j,mq)
         end do
      end do
   end do

   return
end subroutine updateq_from_representer


subroutine updateq_from_representer_point(i,j,mx,my,mbc,meqn,qvar,q_model,beta_j)

!!     =====================================================

!!       # q_hat(i,j) = q_F(x,y,t) + beta*R_m(ij)
!!       # Pointwise AMR update using an interpolated representer value.

   implicit none

   integer i, j, meqn, mbc, mx, my
   double precision :: qvar(meqn)
   double precision :: q_model(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)

   integer mq

   double precision :: beta_j


   do mq = 1, meqn
      q_model(i,j,mq) = q_model(i,j,mq) + beta_j*qvar(mq)
   end do


   return
end subroutine updateq_from_representer_point
