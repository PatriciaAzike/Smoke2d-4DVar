module forward_module
    implicit none

    integer :: source_model
    integer :: initial_condition
    double precision :: W_f !! weight associated with the err in the forcing
    double precision :: W_i !! weight associated with the err in the initial cond.
    
    double precision :: beta, x0, y0
    double precision, dimension(:), allocatable :: beta_rep_value


end module forward_module




subroutine set_beta_values(beta,m)
    !! here we set beta which is obtained from the gauges
    use forward_module, only : beta_rep_value
    implicit none

    integer, intent(in) :: m
    double precision, intent(in) :: beta(m)

    integer :: i

    if (.not. allocated(beta_rep_value)) then
        allocate(beta_rep_value(m))
    end if

    do i = 1, m
        beta_rep_value(i) = beta(i)
    end do

end subroutine set_beta_values
