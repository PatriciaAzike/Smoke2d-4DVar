double precision function  forward_source(mq,xc,yc,t)
    use forward_module
    implicit none

    integer mq
    double precision xc, yc,t
    double precision S

    S = 5

    forward_source = S

    return

end function forward_source