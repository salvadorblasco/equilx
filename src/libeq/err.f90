! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>

subroutine var_conc(C, P, varB, varC, NCP, NEQ)
    ! Compute the variance of the concentrations.
    implicit none
    double precision, intent(in),  dimension(NCP+NEQ)  :: C
    double precision, intent(out), dimension(NCP+NEQ)  :: varC
    double precision, intent(out), dimension(NEQ)      :: varB
    integer,          intent(in),  dimension(NEQ, NCP) :: P
    integer,          intent(in)                       :: NCP     ! number of components
    integer,          intent(in)                       :: NEQ     ! number of equilibria

    double precision, dimension(NCP+NEQ, NEQ) :: dlogcdlogb
    integer i

    call calc_dlogcdlogb(C, P, NCP, NEQ, dlogcdlogb(1:NCP,:))
    call calc_extdlogcdlogb(dlogcdlogb, P, NCP, NEQ)

    ! calculate the errors
    do i=1,NCP+NEQ
        varC(i) = C(i)**2*sum(varB * dlogcdlogb(i,:)**2)
    end do
end subroutine var_conc


subroutine calc_extdlogcdlogb(dlogcdlogb, P, NCP, NEQ)
    double precision, intent(inout), dimension(NCP+NEQ, NEQ) :: dlogcdlogb
    integer,          intent(in),    dimension(NEQ, NCP)     :: P
    integer,          intent(in)                             :: NCP     ! number of components
    integer,          intent(in)                             :: NEQ     ! number of equilibria

    integer i
    do i=1,NEQ
        dlogcdlogb(i,:) = sum(P(i,:) * dlogcdlogb(:,i))
        dlogcdlogb(i,i) = dlogcdlogb(i,i) + 1.0
    end do
end subroutine calc_extdlogcdlogb


subroutine calc_dlogcdlogb(C, P, NCP, NEQ, dlcdlb)
    double precision, intent(in),  dimension(NCP+NEQ)  :: C
    integer,          intent(in),  dimension(NEQ, NCP) :: P
    integer,          intent(in)                       :: NCP     ! number of components
    integer,          intent(in)                       :: NEQ     ! number of equilibria
    double precision, intent(out), dimension(NCP, NEQ) :: dlcdlb  ! the result

    double precision, dimension(NCP, NCP) :: Am
    double precision, dimension(NCP, NEQ) :: bm
    integer i, j

    ! For DGELSY 
    integer, dimension(NCP) :: ipiv   ! pivot indices
    integer :: lwork, info
    double precision, allocatable, dimension(:) :: work

    ! compute A matrix
    call amatrix(Am, C, P, NCP, NEQ)

    ! compute b matrix
    do i=1,NCP
        do j=1,NEQ
            bm(i,j) = -P(j,i) * C(NCP+j)
        end do
    end do

    ! solve the system
    lwork = -1
    allocate(work(1))
    call dsysv("U", NCP, NEQ, Am, NCP, ipiv, bm, NCP, work, lwork, info)
    lwork = nint(work(1))
    deallocate(work)

    allocate(work(lwork))
    call dsysv("U", NCP, NEQ, Am, NCP, ipiv, bm, NCP, work, lwork, info)
    deallocate(work)

    dlcdlb = bm
end subroutine calc_dlogcdlogb
