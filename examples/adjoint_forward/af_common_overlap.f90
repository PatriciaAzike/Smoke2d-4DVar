subroutine setaux_from_adjoint(mx,my,mbc,meqn,maux,adjoint,aux_representer)
    implicit none
    integer meqn, mbc, mx, my, maux
    double precision       adjoint(1-mbc:mx+mbc,1-mbc:my+mbc, meqn)
    double precision   aux_representer(1-mbc:mx+mbc,1-mbc:my+mbc, maux)
    
    integer i,j

    !! Copy solution into forward aux array so we can access it from the
    !! source term routine
    !! First two entries in aux array are the velocity
    do i = 1-mbc,mx+mbc
        do j = 1-mbc, my+mbc
            aux_representer(i,j,3) = adjoint(i,j,1)
        end do
    end do

    return
end subroutine setaux_from_adjoint


subroutine setq_from_adjoint(mx,my,mbc,meqn,adjoint,q_representer)
    use forward_module, only : W_i
    implicit none
    integer meqn, mbc, mx, my
    double precision       adjoint(1-mbc:mx+mbc,1-mbc:my+mbc, meqn)
    double precision     q_representer(1-mbc:mx+mbc,1-mbc:my+mbc, meqn)
    
    integer i,j

    !! Copy adjoint solution into forward array 
    do i = 1-mbc,mx+mbc
        do j = 1-mbc, my+mbc
            q_representer(i,j,1) = (1.d0/W_i) * adjoint(i,j,1)
        end do
    end do

    return
end subroutine setq_from_adjoint


subroutine setaux_from_adjoint_point(i,j,mx,my,mbc,meqn,maux,qvar,aux_representer)
    implicit none
    integer i, j, meqn, mbc, mx, my, maux
    double precision       qvar(meqn)
    double precision   aux_representer(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

    integer mq

    !! Store one interpolated adjoint value at one target cell center.
    do mq = 1, meqn
        aux_representer(i,j,2+mq) = qvar(mq)
    end do

    return
end subroutine setaux_from_adjoint_point


subroutine setq_from_adjoint_point(i,j,mx,my,mbc,meqn,qvar,q_representer)
    use forward_module, only : W_i
    implicit none
    integer i, j, meqn, mbc, mx, my
    double precision       qvar(meqn)
    double precision     q_representer(1-mbc:mx+mbc,1-mbc:my+mbc, meqn)

    integer mq

    !! Initialize the representer at one target cell center using one
    !! interpolated adjoint value.
    do mq = 1, meqn
        q_representer(i,j,mq) = (1.d0/W_i) * qvar(mq)
    end do

    return
end subroutine setq_from_adjoint_point
