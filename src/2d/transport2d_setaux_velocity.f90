subroutine transport2d_setaux_velocity(blockno, mx,my,mbc, & 
   dx,dy,xlower,ylower, t, xp,yp,zp, & 
   xnormals,ynormals,surfnormals, aux, maux)
    implicit none

    integer mx,my,mbc,maux,blockno
    double precision dx,dy, xlower,ylower, t
    double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
    double precision xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
    double precision ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)

    double precision surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)

    double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
    double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
    double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    double precision xc,yc,xc1, yc1
    double precision nv(3), vel(3), vdotn, transport_map_dot

    double precision nl(3), nr(3), nb(3), nt(3)
    double precision urrot, ulrot, ubrot, utrot      

    integer i,j, k



    !! # Cell-centered velocities : entries (4,5,6)
    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc

            !! # Brick units
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy


            !! # Map unit square in region specified by lat/long coordinates
            !! # to fraction unit square defined over sphere.
            call map_brick2unitsq(blockno,xc,yc,xc1,yc1)
            
            !! # Get velocity in unit square spherical coordinates.
            call transport_velocity(xc1,yc1,t, vel)


            !! # subtract out normal components
            do k = 1,3
                nv(k) = surfnormals(i,j,k)
            enddo

            vdotn = transport_map_dot(vel,nv)

            !! # Subtract out component in the normal direction
            do k = 1,3
                vel(k) = vel(k) - vdotn*nv(k)
            end do

            do k = 1,3
                nl(k)  = xnormals(i,  j,  k)
                nr(k)  = xnormals(i+1,j,  k)
                nb(k)  = ynormals(i,  j,  k)
                nt(k)  = ynormals(i,  j+1,k)
            enddo

            ulrot = transport_map_dot(nl,vel)
            urrot = transport_map_dot(nr,vel)
            ubrot = transport_map_dot(nb,vel)
            utrot = transport_map_dot(nt,vel)

            

            aux(i,j,2) = ulrot
            aux(i,j,3) = urrot
            aux(i,j,4) = ubrot
            aux(i,j,5) = utrot

            !print *, 'aux(i,j,2) = ', aux(i,j,2)
            !print *, 'stopping in setaux_velocity'
            !stop

         enddo
    enddo

end subroutine transport2d_setaux_velocity




double precision function transport_map_dot(u,v)
    implicit none

    double precision u(3),v(3)

    transport_map_dot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

end
