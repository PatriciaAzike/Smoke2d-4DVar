subroutine adjoint_src2(maxmx,maxmy,meqn,mbc,mx,my, &
    xlower,ylower,dx,dy,q,maux,aux,t,dt)
!! Apply an observation as an instantaneous event and integrate backwards.
!!
!! The adjoint advances in tau = T - t_phys.  Observation m therefore occurs
!! at tau_m = T - t_m.  When a step crosses tau_m, add the complete spatial
!! heat kernel once:
!!
!!   alpha(t_m^-) = alpha(t_m^+) + G_eps(x - x_m).
!!
!! No factor of dt and no temporal Gaussian appear here.  With output times
!! chosen to include tau_m, the event is deposited at the end of the step that
!! lands on tau_m, before that state is written.

    use adjoint_module, only : xm, ym, tm, tfinal, pseudo_1d, obs_index
    implicit none

    integer maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision xlower,ylower,dx,dy,t,dt
    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

    integer i,j,mq
    double precision xc,yc,r,event_time,tol,dirac_delta
    logical crosses_event

    event_time = tfinal - tm(obs_index)
    tol = 100.d0*epsilon(1.d0)*max(1.d0,abs(event_time),abs(t),abs(dt))
    crosses_event = (t .lt. event_time) .and. &
                    (t + dt .ge. event_time - tol)

    if (.not. crosses_event) return

    do mq = 1,meqn
        do i = 1-mbc,mx+mbc
            xc = xlower + (i-0.5d0)*dx
            do j = 1-mbc,my+mbc
                yc = ylower + (j-0.5d0)*dy

                if (pseudo_1d .eq. 1) then
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
