double precision function forward_psi(x,y,z)
    use adjoint_module, only : pseudo_1d
    implicit none

    double precision x,y,z,r

    double precision pi, pi2
    common /compi/ pi, pi2

    r = sqrt((x-1.d0)**2 + (y-1.d0)**2)

    !! # Rigid body rotation
    if (pseudo_1d .eq. 1) then
        forward_psi = y                !! pseudo-1d example
    elseif (pseudo_1d .eq. 2) then
        forward_psi = (4.d0/3.d0)*r**3 !! full 2D example
    endif

    !! # Filament formation (negative for clockwise rotation)
    !!forward_psi = -(y - x)!!(4.d0/3.d0)*r**3

    return
end function forward_psi

subroutine forward_get_psi_vel(xd1,xd2,ds,vn,t)
    implicit none

    double precision xd1(3),xd2(3), ds, vn, forward_psi,t

    vn = (forward_psi(xd1(1),xd1(2),xd1(3)) - & 
          forward_psi(xd2(1),xd2(2),xd2(3)))/ds

end subroutine forward_get_psi_vel
