!!We tag for coarsening if this coarsened patch isn't tagged for refinement
subroutine clawpatch46_fort_tag4coarsening(mx,my,mbc,meqn,&
    xlower,ylower,dx,dy,blockno, q0, q1, q2, q3,&
    coarsen_threshold, init_flag, tag_patch)
    implicit none

    integer mx,my,mbc,meqn,tag_patch
    integer blockno, init_flag
    double precision xlower(0:3),ylower(0:3),dx,dy
    double precision coarsen_threshold
    double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    integer mq
    double precision qmin, qmax

    !! Assume that we will coarsen a family unless we find a grid
    !!that doesn't pass the coarsening test.
    tag_patch = 1
    mq = 1
    qmin = q0(1,1,mq)
    qmax = q0(1,1,mq)


    !! If we find that (qmax-qmin > coarsen_threshold) on any
    !! grid, we return immediately, since the family will then
    !! not be coarsened.

    call get_minmax(mx,my,mbc,meqn,mq,q0,qmin,qmax,&
        coarsen_threshold,tag_patch)
    if (tag_patch == 0) return

    call get_minmax(mx,my,mbc,meqn,mq,q1,qmin,qmax,&
        coarsen_threshold,tag_patch)
    if (tag_patch == 0) return

    call get_minmax(mx,my,mbc,meqn,mq,q2,qmin,qmax,&
        coarsen_threshold,tag_patch)
    if (tag_patch == 0) return

    call get_minmax(mx,my,mbc,meqn,mq,q3,qmin,qmax,&
        coarsen_threshold,tag_patch)

end subroutine clawpatch46_fort_tag4coarsening

subroutine get_minmax(mx,my,mbc,meqn,mq,q,&
    qmin,qmax,coarsen_threshold,tag_patch)

    implicit none
    integer mx,my,mbc,meqn,mq,tag_patch
    double precision coarsen_threshold
    double precision qmin,qmax
    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    integer i,j

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            if (qmax - qmin .gt. coarsen_threshold) then
            !! We won't coarsen this family because at least one
            !! grid fails the coarsening test.
                tag_patch = 0
                return
            endif
        enddo
    enddo
    
end subroutine get_minmax
