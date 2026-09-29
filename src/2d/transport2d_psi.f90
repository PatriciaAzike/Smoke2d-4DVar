double precision function psi(xd,yd,zd,t)
   use smoke3d_module, only : pi
   implicit none

   double precision xd, yd, zd, t

   double precision l, th, lp, scale, kappa, tfinal


   kappa = 2.0
   tfinal = 5.0

   call map2polar(xd,yd,zd,l,th)

   lp = l - 2*pi*t/tfinal

   psi = kappa*sin(lp)**2*cos(th)**2*cos(pi*t/Tfinal) & 
      - 2*pi*sin(th)/Tfinal

   !! # Sign difference from Benchmark problem
   scale = 1.1e-3
   psi = scale*psi
   


end function psi

subroutine map2polar(x,y,z,lambda,th)
   implicit none

   double precision x,y,z,th,lambda
   double precision r

   double precision pi, pi2
   common /compi/ pi, pi2

   r = sqrt(x**2 + y**2 + z**2)

   th = asin(z/r)

   lambda = atan2(y,x)
   if (lambda .lt. 0) then
      lambda = lambda + 2*pi
   endif

end subroutine map2polar


