double precision function adjoint_psi(x,y,z)
    implicit none

    double precision x,y,z,r

    double precision pi, pi2
    common /compi/ pi, pi2

    double precision forward_psi


    adjoint_psi = -forward_psi(x,y,z)

    !!r = sqrt((x-1.d0)**2 + (y-1.d0)**2)

    !! # Rigid body rotation
    !!adjoint_psi = -(4.d0/3.d0)*r**3!r**2 !-y -- 1D example


    !! # Eventually, we want the adjoint problem to solve the forward problem, but in 
    !! # reverse. So we should negate the velocity below.            
    !! # Filament formation (negative for clockwise rotation)
    !!adjoint_psi = y - x!!(4.d0/3.d0)*r**3


    return
end function adjoint_psi

subroutine adjoint_get_psi_vel(xd1,xd2,ds,vn,t)
    implicit none

    double precision xd1(3),xd2(3), ds, vn, adjoint_psi,t

    vn = (adjoint_psi(xd1(1),xd1(2),xd1(3)) - & 
          adjoint_psi(xd2(1),xd2(2),xd2(3)))/ds

end subroutine adjoint_get_psi_vel
