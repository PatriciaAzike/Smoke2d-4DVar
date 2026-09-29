double precision function  adjoint_source(mq,xc,yc,t)
    use adjoint_module, only : xm, ym, tm, mdata,tfinal,pseudo_1d, obs_index
    implicit none

    integer mq
    double precision xc, yc,t
    double precision dirac_delta, dirac_delta_time, r, S

    double precision theta

    double precision pi, pi2
    common /compi/ pi, pi2

    theta = pi/4.d0

#if 1
    if (mq > mdata) then
        write(6,*) "adjoint_source: mq > mdata"
        stop
    endif
#endif

    if (.not. allocated(xm)) then
        write(6,*) "adjoint_source: data not allocated"
        stop
    endif

    if (pseudo_1d .eq. 1) then
        r = abs(xc - xm(obs_index))                      !! Euclidean distance for pseudo-1d experiment
    elseif (pseudo_1d .eq. 2) then
        r = sqrt((xc-xm(obs_index))**2 + (yc-ym(obs_index))**2) !! Euclidean distance for full 2D experiment
    endif
    !write(6,*) 'obs_index = ', obs_index, ' xm = ', xm(obs_index), ' ym = ', ym(obs_index)
    
    !r = (xc - xm(mq))*cos(theta) + (yc - ym(mq))*sin(theta)
    
    S = dirac_delta(r) * dirac_delta_time(tfinal-t-tm(obs_index))
    
    adjoint_source = S

    return

end function adjoint_source