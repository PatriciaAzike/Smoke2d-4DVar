subroutine rpn2_forward(ixy,maxm, meqn,mwaves,mbc, & 
   mx,ql,qr,auxl,auxr,wave,s, amdq,apdq,maux)
   implicit none

   integer maxm, mbc,mwaves,meqn,mx, maux
   integer ixy

   double precision wave(1-mbc:maxm+mbc, meqn, mwaves)   
   double precision    s(1-mbc:maxm+mbc, mwaves)
   double precision   ql(1-mbc:maxm+mbc, meqn)
   double precision   qr(1-mbc:maxm+mbc, meqn)
   double precision amdq(1-mbc:maxm+mbc, meqn)
   double precision apdq(1-mbc:maxm+mbc, meqn)
   double precision auxl(1-mbc:maxm+mbc,maux)
   double precision auxr(1-mbc:maxm+mbc,maux)


   integer i, idir
   double precision qll,qrr
   double precision urrot, ulrot, g, uhat

#if 0
   integer m,iface

   iface = ixy
   do i = 2-mbc, mx+mbc

      do m = 1,meqn
         wave(i,m,1) = ql(i,m) - qr(i-1,m)
      enddo
 
      s(i,1) = auxl(i,iface)

      do m = 1,meqn
         amdq(i,m) = min(s(i,1), 0.d0) * wave(i,m,1)
         !write(*,*) 'min(s(i,1), 0.d0) = ', min(s(i,1), 0.d0)
         !write(*,*) 'max(s(i,1), 0.d0) = ', max(s(i,1), 0.d0)
         apdq(i,m) = max(s(i,1), 0.d0) * wave(i,m,1)
      enddo

   enddo
#endif

   idir = ixy-1
   do i = 2-mbc, mx+mbc
      !!g = auxl(i,6+idir)  !! Edge length. 
      !!We are on a unit square so this is redundant. Hence g=1
      g = 1 !! 

       
      urrot = g * auxl(i,1+idir)
      ulrot = g * auxl(i,1+idir)

      qrr = ql(i,1)
      qll = qr(i-1,1)

      !! # Use Roe-average values         
      uhat = (ulrot + urrot)/2.d0

      if (uhat .ge. 0) then
         amdq(i,1) = 0.d0
         apdq(i,1) = urrot*qrr - ulrot*qll
      else
         amdq(i,1) = urrot*qrr - ulrot*qll
         apdq(i,1) = 0.d0
      endif
      wave(i,1,1) = urrot*qrr - ulrot*qll
      s(i,1) = uhat

   enddo


   return
end subroutine rpn2_forward
