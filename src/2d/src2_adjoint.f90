subroutine src2_adjoint(maxmx,maxmy,meqn,mbc,mx,my, & 
    xlower,ylower,dx,dy,q,maux,aux,t,dt)
    use smoke3d_module
    implicit none

    integer maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision xlower, ylower, dx, dy, t, dt
    double precision   q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
    double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

    integer fc2d_clawpack46_get_block
    integer blockno      

    double precision xc,yc, xp,yp,zp
    double precision x_long, y_lat
    double precision S, delta, r, rmax, rmin
    integer i, j

    blockno = fc2d_clawpack46_get_block()

    rmax = 0
    rmin = 200
    do j = 1-mbc,my+mbc
        do i = 1-mbc,mx+mbc

            !! # Cell center, in brick units     
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            !! # Return physical location (in meters)
            call map_brick2latlong(blockno,xc,yc,x_long,y_lat)
            call transport_mapc2m_latlong(blockno,xc,yc,xp,yp,zp)

            r  = sqrt((x_long - xm(1))**2 + (y_lat - ym(1))**2)
            rmax = max(r,rmax)
            rmin = min(r,rmin)            
            S = delta(r)*delta(t - tm(1))  !! Should be a delta function source term

            !! Multiply by fraction w of cell in fire disk
            q(i,j,1) = q(i,j,1) + dt*S
        enddo
    enddo
    !!write(6,*) 'rmax = ', rmax
    !!write(6,*) 'rmin = ', rmin
    !!write(6,*) xm(1), ym(1)
    !!stop

    return
end subroutine src2_adjoint
