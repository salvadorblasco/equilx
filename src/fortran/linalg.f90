module linalg
    implicit none (type, external)
    public

contains
subroutine las_cond(A, B, x, cond, info)
    double precision, intent(in) :: A(:,:)
    double precision, intent(in) :: B(:)
    double precision, intent(out) :: x(:)
    double precision, intent(out) :: cond
    integer, intent(inout) :: info

    integer :: n
    double precision :: anorm, work(3*size(x))
    integer :: iwork(size(x))
    integer :: linfo
    double precision, allocatable :: tmpx(:,:)
    double precision :: dlansy

    external dpotrf
    external dpocon
    external dpotrs
    external dlansy

    n = size(A, dim=1)
    allocate(tmpx(n,1))

    ! 1-norm of the original matrix (needed by dpocon)
    anorm = dlansy('1', 'U', n, A, n, work)
    ! Cholesky factorisation

    call dpotrf('U', n, A, n, linfo)
    if (linfo /= 0) then
       cond = 0.0
       return
    end if

    ! Reciprocal condition-number estimate
    call dpocon('U', n, A, n, anorm, cond, work, iwork, info)
    ! info from dpocon is almost always 0; we ignore it

    tmpx(:,1) = b
    call dpotrs('U', n, 1, A, n, tmpx, n, info)
    x = tmpx(:,1)

    ! return condition, not reciprocal
    cond = 1.0/cond
end subroutine las_cond

end module linalg
