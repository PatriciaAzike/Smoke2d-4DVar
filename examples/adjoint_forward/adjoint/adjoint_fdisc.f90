double precision function  fdisc(blockno,xc,yc)
    use adjoint_module,only : x0, y0
    implicit none

    integer blockno
    double precision xc,yc

    double precision r

    r = sqrt((xc-x0)**2 + (yc-y0)**2)

    if (blockno > 0) stop

    fdisc = r-0.25d0
end function fdisc
