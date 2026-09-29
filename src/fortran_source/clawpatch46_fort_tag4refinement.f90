subroutine clawpatch46_fort_tag4refinement(mx,my,mbc,&
    meqn,xlower,ylower,dx,dy,blockno,q, &
    refine_threshold, init_flag,tag_patch)
    implicit none

    integer mx,my,mz,mbc,meqn,tag_patch,init_flag
    integer blockno
    double precision xlower,ylower,dx,dy,dz
    double precision refine_threshold
    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    integer i,j,mq
    double precision qmin, qmax

    tag_patch = 0

    !! Refine based only on first variable in system.
    mq = 1
    qmin = q(1,1,mq)
    qmax = q(1,1,mq)

    do j = 1-mbc,my+mbc
        do i = 1-mbc,mx+mbc
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            if (qmax - qmin .gt. refine_threshold) then
                tag_patch = 1
                return
            endif
        enddo
    enddo
    

end subroutine clawpatch46_fort_tag4refinement
