subroutine adjoint_src2(maxmx,maxmy,meqn,mbc,mx,my, &
    xlower,ylower,dx,dy,q,maux,aux,tau,dt)
!! Apply an observation as an instantaneous event and integrate backwards.
!!
!! The adjoint advances in tau = T - t, where t is physical time.
!! Observation m therefore occurs at tau_m = T - t_m.  When a step crosses
!! tau_m, add the complete spatial heat kernel once:
!!
!!   alpha(t_m^-) = alpha(t_m^+) + G_eps(x - x_m).
!!
!! No factor of dt and no temporal Gaussian appear here.  With output times
!! chosen to include tau_m, the event is deposited at the end of the step that
!! lands on tau_m, before that state is written.

    use adjoint_module, only : xm, ym, tm, tfinal, pseudo_1d_experiment, obs_index
    implicit none

    integer maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision xlower,ylower,dx,dy,tau,dt
    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

    integer i,j,mq
    double precision xc,yc,r,tau_m,tol,dirac_delta
    logical crosses_event

    tau_m = tfinal - tm(obs_index)

    !! tol is far above round-off and far below any dt.
    tol = 1.d-12

    !! Each step owns (tau, tau+dt].  Comparing both ends with the same point,
    !! tau_m - tol, catches a step that lands a round-off short of tau_m and
    !! applies the event exactly once.
    crosses_event = (tau .lt. tau_m - tol) .and. &
                    (tau + dt .ge. tau_m - tol)

    if (.not. crosses_event) return

    do mq = 1,meqn
        do i = 1-mbc,mx+mbc
            xc = xlower + (i-0.5d0)*dx
            do j = 1-mbc,my+mbc
                yc = ylower + (j-0.5d0)*dy

                if (pseudo_1d_experiment) then
                    r = abs(xc-xm(obs_index))
                else
                    r = sqrt((xc-xm(obs_index))**2 + &
                             (yc-ym(obs_index))**2)
                endif

                q(i,j,mq) = q(i,j,mq) + dirac_delta(r)
            enddo
        enddo
    enddo

end subroutine adjoint_src2
