module linalg
    implicit none (type, external)
    public

    interface nnls
        module procedure nnls_vector, nnls_matrix
    end interface nnls

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


! Minimize ||A - C*E||_F subject to E >= 0, one wavelength at a time.
! info: 0 success, -1 invalid dimensions, -2 non-finite input,
!       1 iteration limit, 2 numerical rank failure, 3 LAPACK or numerical failure.
! On failure E is zero; failed_column identifies a solver failure, or 0 otherwise.
subroutine nnls_matrix(c, a, e, info, failed_column)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    double precision, intent(in) :: c(:,:), a(:,:)
    double precision, intent(out) :: e(:,:)
    integer, intent(out) :: info
    integer, intent(out), optional :: failed_column

    integer :: j, npt, nsp, nwl

    npt = size(c, 1)
    nsp = size(c, 2)
    nwl = size(a, 2)
    e = 0.0d0
    info = 0
    if (present(failed_column)) failed_column = 0
    if (npt < 1 .or. nsp < 1 .or. nwl < 1 .or. size(a, 1) /= npt .or. &
        size(e, 1) /= nsp .or. size(e, 2) /= nwl) then
        info = -1
        return
    end if
    if (.not. all(ieee_is_finite(c)) .or. .not. all(ieee_is_finite(a))) then
        info = -2
        return
    end if

    do j = 1, nwl
        call nnls_vector(c, a(:,j), e(:,j), info)
        if (info /= 0) then
            e = 0.0d0
            if (present(failed_column)) failed_column = j
            return
        end if
    end do
end subroutine nnls_matrix

! Minimize ||A*x - b||_2 subject to x >= 0. Status codes are as above.
subroutine nnls_vector(a, b, x, info)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    double precision, intent(in) :: a(:,:), b(:)
    double precision, intent(out) :: x(:)
    integer, intent(out) :: info

    integer :: m, n, i, j, k, nfree, rank, lapack_info, lwork, iterations
    integer, allocatable :: free_idx(:), jpvt(:)
    double precision :: bscale, alpha, candidate_alpha, tol, rcond, largest
    double precision :: query(1)
    double precision, allocatable :: c(:,:), column_scale(:), column_norm(:), rhs(:)
    double precision, allocatable :: y(:), z(:), dual(:), residual(:), work(:), ls_a(:,:), ls_b(:,:)
    logical, allocatable :: passive(:), rejected(:)
    external dgelsy

    m = size(a, 1)
    n = size(a, 2)
    x = 0.0d0
    info = 0
    if (m < 1 .or. n < 1 .or. size(b) /= m .or. size(x) /= n) then
        info = -1
        return
    end if
    if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) then
        info = -2
        return
    end if

    bscale = maxval(abs(b))
    if (bscale == 0.0d0) return

    allocate(c(m,n), column_scale(n), column_norm(n), rhs(m))
    allocate(y(n), z(n), dual(n), residual(m), free_idx(n), jpvt(n))
    allocate(passive(n), rejected(n), ls_a(m,n), ls_b(max(m,n),1))

    ! Scaling columns makes the dual stopping test independent of their units.
    do j = 1, n
        column_scale(j) = maxval(abs(a(:,j)))
        column_norm(j) = 1.0d0
        c(:,j) = 0.0d0
        if (column_scale(j) > 0.0d0) then
            c(:,j) = a(:,j) / column_scale(j)
            column_norm(j) = norm2(c(:,j))
            c(:,j) = c(:,j) / column_norm(j)
        end if
    end do
    rhs = b / bscale
    rcond = epsilon(1.0d0) * max(m,n)
    y = 0.0d0
    passive = .false.
    rejected = .false.
    iterations = 0

    ls_a = 0.0d0
    ls_b = 0.0d0
    jpvt = 0
    call dgelsy(m, n, 1, ls_a, m, ls_b, max(m,n), jpvt, rcond, rank, query, -1, lapack_info)
    if (lapack_info /= 0 .or. .not. ieee_is_finite(query(1))) then
        info = 3
        return
    end if
    lwork = max(1, ceiling(query(1)))
    allocate(work(lwork))

    do
        residual = rhs - matmul(c, y)
        dual = matmul(transpose(c), residual)
        ! Account for cancellation when large coefficients produce a small residual.
        tol = 10.0d0 * epsilon(1.0d0) * max(m,n) * (norm2(rhs) + sum(abs(y)))
        k = 0
        largest = tol
        do i = 1, n
            if (.not. passive(i) .and. .not. rejected(i) .and. dual(i) > largest) then
                k = i
                largest = dual(i)
            end if
        end do
        if (k == 0) then
            if (any(rejected .and. dual > tol)) then
                info = 2
                return
            end if
            exit
        end if

        passive(k) = .true.
        do
            iterations = iterations + 1
            if (iterations > 30 * n) then
                info = 1
                return
            end if

            nfree = 0
            do i = 1, n
                if (passive(i)) then
                    nfree = nfree + 1
                    free_idx(nfree) = i
                    ls_a(:,nfree) = c(:,i)
                end if
            end do
            ls_b = 0.0d0
            ls_b(1:m,1) = rhs
            jpvt = 0
            call dgelsy(m, nfree, 1, ls_a, m, ls_b, max(m,n), jpvt, rcond, rank, work, lwork, lapack_info)
            if (lapack_info /= 0) then
                info = 3
                return
            end if
            if (rank < nfree) then
                passive(k) = .false.
                rejected(k) = .true.
                exit
            end if

            z = 0.0d0
            do j = 1, nfree
                z(free_idx(j)) = ls_b(j,1)
            end do
            if (.not. all(ieee_is_finite(z))) then
                info = 3
                return
            end if
            if (all(z(free_idx(1:nfree)) > 0.0d0)) then
                y = z
                rejected = .false.
                exit
            end if

            alpha = 1.0d0
            j = 0
            do i = 1, n
                if (passive(i) .and. z(i) <= 0.0d0) then
                    if (y(i) == 0.0d0) then
                        candidate_alpha = 0.0d0
                    else
                        candidate_alpha = y(i) / (y(i) - z(i))
                    end if
                    if (candidate_alpha <= alpha) then
                        alpha = candidate_alpha
                        j = i
                    end if
                end if
            end do
            y = (1.0d0 - alpha) * y + alpha * z
            y(j) = 0.0d0
            where (passive .and. y <= 0.0d0)
                passive = .false.
                y = 0.0d0
            end where
            rejected = .false.
        end do
    end do

    do i = 1, n
        if (column_scale(i) > 0.0d0) x(i) = (y(i) * bscale / column_scale(i)) / column_norm(i)
    end do
    if (.not. all(ieee_is_finite(x))) then
        x = 0.0d0
        info = 3
    end if
end subroutine nnls_vector

end module linalg
