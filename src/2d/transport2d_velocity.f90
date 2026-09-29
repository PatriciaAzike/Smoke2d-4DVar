subroutine transport2d_velocity_components_cart(x,y,t,vcart)
    implicit none

    double precision x,y,t,vcart(3)
    double precision t1(3), t2(3)
    integer k
    logical :: is_cart

    double precision :: u(2)

    !! Initialize u to avoid compiler warnings
    u(1) = 0
    u(2) = 0

    call latlong_velocity(x,y,t, u,vcart,is_cart)
    


    !!if (time_dependent_velocity) then
    if (.not. is_cart) then

            !! # Velocity components are given in spherical components
            !! # and must be converted to Cartesian
            call map_covariant_basis(x, y, t1,t2)

            !! # Normalize basis vectors so that the speed is correct
            call transport_normalize_vector(t1)
            call transport_normalize_vector(t2)

            do k = 1,3
                vcart(k) = u(1)*t1(k) + u(2)*t2(k)
            enddo
            !!write(6,*) 'Stopping in transport2d velocity.f90'
            !!stop

    endif
    !!endif 

end subroutine transport2d_velocity_components_cart

subroutine transport_normalize_vector(v)
    implicit none

    double precision v(3), vn
    integer k

    vn = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
    vn = sqrt(vn)

    do k = 1,3
        v(k) = v(k)/vn
    enddo

end subroutine transport_normalize_vector

!! # ------------------------------------------------------------
!! # Called from set_aux to set velocity components in 
!! # Cartesian coordinates.
!! # 
!! # Inputs
!! # 
!! # (x,y)  : Coordinates in [0,1]x[0,1] (mapped to unit sphere)
!! #     t  : Time (s)
!! #  vcart : Velocity (u,v,w) relative to Cartesian basis.
!! # ------------------------------------------------------------

subroutine transport_velocity(x,y,t,vcart)
    implicit none

    double precision x,y,t, vcart(3)

    call transport2d_velocity_components_cart(x,y,t,vcart)
end subroutine transport_velocity






