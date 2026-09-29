subroutine forward_setprob
    use forward_module
    implicit none

    double precision pi, pi2
    common /compi/ pi, pi2

    integer example
    common /com_ex/ example

    pi = 4.d0*atan(1.d0)
    pi2 = 2*pi

    open(10,file='setprob.data')
    read(10,*) example
    read(10,*) initial_condition
    read(10,*) W_f
    read(10,*) W_i
    read(10,*) beta
    read(10,*) x0
    read(10,*) y0
    close(10)

    
end subroutine forward_setprob
