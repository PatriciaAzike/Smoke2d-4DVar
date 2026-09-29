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
  !! # ixy : 2  Normal solve is in y
  
  idir = 2 - ixy

  do i = 2-mbc, mx+mbc
      i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq

      !! # Get an edge value
      bmasdq(i,1) = min(aux2(i1,idir+1),0.d0)*asdq(i,1)

      bpasdq(i,1) = max(aux3(i1,idir+1),0.d0)*asdq(i,1)
  enddo


    return
end subroutine rpt2_forward
