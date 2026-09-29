subroutine rpt2_forward(ixy,maxm,meqn, mwaves,mbc,mx, & 
  ql,qr,aux1,aux2,aux3,imp,asdq, bmasdq,bpasdq, maux)
      
  implicit none

  integer ixy, maxm, meqn,mwaves,mbc,mx,imp, maux

  double precision     ql(1-mbc:maxm+mbc, meqn)
  double precision     qr(1-mbc:maxm+mbc, meqn)
  double precision   asdq(1-mbc:maxm+mbc, meqn)
  double precision bmasdq(1-mbc:maxm+mbc, meqn)
  double precision bpasdq(1-mbc:maxm+mbc, meqn)
  double precision   aux1(1-mbc:maxm+mbc, maux)
  double precision   aux2(1-mbc:maxm+mbc, maux)
  double precision   aux3(1-mbc:maxm+mbc, maux)


  integer i, i1, idir
  double precision vrrot, vlrot, g, vhat

  !! # Direction of the normal solve
  !! # ixy : 1  Normal solve is in x
  !! # ixy : 2 Normal solve is in y
  
  idir = 2 - ixy

  do i = 2-mbc, mx+mbc
      i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq

      !! # -----------------------------------------
      !! # Lower faces - cell centered velocities
      !! # -----------------------------------------
           
      !! # 6-7    Edge lengths (x-face, y-face)
      g = aux2(i1,6+idir)

      !! # left-right : 2,3
      !! # top-bottom : 4,5         
      vrrot = g*aux2(i,  2 + 2*idir)   !! Left edge of right cell
      vlrot = g*aux2(i-1,3 + 2*idir)   !! Right edge of left cell

      !! # Get an edge value
      vhat = (vrrot + vlrot)/2.0

      bmasdq(i,1) = min(vhat,0.d0)*asdq(i,1)

      !! # -----------------------------------------
      !! # Upper faces - cell centered velocities
      !! # -----------------------------------------

      g = aux3(i1,6+idir)

      !! # left-right : 2,3
      !! # top-bottom : 4,5         
      vrrot = g*aux3(i,  2 + 2*idir)   !! Left edge of right cell
      vlrot = g*aux3(i-1,3 + 2*idir)   !! Right edge of left cell

      vhat = (vrrot + vlrot)/2.0

      bpasdq(i,1) = max(vhat,0.d0)*asdq(i,1)

    enddo


    return
end subroutine rpt2_forward
