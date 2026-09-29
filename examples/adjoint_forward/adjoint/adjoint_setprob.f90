subroutine adjoint_setprob
    use adjoint_module, only : beta, x0, y0, initial_condition, mdata, &
                               pseudo_1d,eps_1d,eps_2d,xm,ym,tm,dm, &
                               W_eps,tfinal,obs_index
    implicit none

    double precision pi, pi2
    common /compi/ pi, pi2

    integer example
    common /com_ex/ example

    integer i

    pi = 4.d0*atan(1.d0)
    pi2 = 2*pi

    open(10,file='setprob.data')
    read(10,*) example
    read(10,*) mdata
    read(10,*) initial_condition
    read(10,*) eps_1d
    read(10,*) eps_2d
    read(10,*) beta
    read(10,*) x0
    read(10,*) y0
    read(10,*) tfinal

    call data(mdata)

    do i = 1, mdata
        read(10,*) xm(i)
        read(10,*) ym(i)
        read(10,*) tm(i)
        read(10,*) dm(i)
        read(10,*) W_eps(i)
    end do

    read(10,*) pseudo_1d

    obs_index = 1

    close(10)
    

end subroutine adjoint_setprob
