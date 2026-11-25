program test_linalg
    use linalg, only: nnls
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none (type, external)

    call test_positive()
    call test_boundary()
    call test_scaled_columns()
    call test_active_set_step()
    call test_underdetermined()
    call test_dependent_columns()
    call test_ill_conditioned_fit()
    call test_invalid_input()
    call test_matrix_fit()
    call test_matrix_invalid_input()
    call test_matrix_solver_failure()
    print *, "NNLS tests passed"

contains

subroutine assert_close(actual, expected, tolerance)
    double precision, intent(in) :: actual(:), expected(:), tolerance

    if (any(abs(actual - expected) > tolerance * max(1.0d0, maxval(abs(expected))))) then
        print *, "Expected:", expected
        print *, "Actual:  ", actual
        error stop "NNLS result mismatch"
    end if
end subroutine assert_close

subroutine test_positive()
    double precision :: a(3,2), b(3), x(2)
    integer :: info

    a(:,1) = [1.0d0, 0.0d0, 1.0d0]
    a(:,2) = [0.0d0, 1.0d0, 1.0d0]
    b = [2.0d0, 3.0d0, 5.0d0]
    call nnls(a, b, x, info)
    if (info /= 0) error stop "positive solution: unexpected status"
    call assert_close(x, [2.0d0, 3.0d0], 1.0d-12)
end subroutine test_positive

subroutine test_boundary()
    double precision :: a(2,2), b(2), x(2)
    integer :: info

    a(:,1) = [1.0d0, 0.0d0]
    a(:,2) = [1.0d0, 1.0d0]
    b = [0.0d0, 1.0d0]
    call nnls(a, b, x, info)
    if (info /= 0) error stop "boundary solution: unexpected status"
    call assert_close(x, [0.0d0, 0.5d0], 1.0d-12)

    b = [-1.0d0, -1.0d0]
    call nnls(a, b, x, info)
    if (info /= 0) error stop "zero solution: unexpected status"
    call assert_close(x, [0.0d0, 0.0d0], 1.0d-12)
end subroutine test_boundary

subroutine test_scaled_columns()
    double precision :: a(4,3), b(4), x(3)
    integer :: info, i

    a(:,1) = [2.0d0, 1.0d0, 0.0d0, 2.0d0]
    a(:,2) = [2.0d0, 4.0d0, 2.0d0, 1.0d0]
    a(:,3) = [3.0d0, 4.0d0, 2.0d0, 2.0d0]
    b = [2.0d0, 4.0d0, 5.0d0, 4.0d0]
    do i = 1, 2
        call nnls(a, b, x, info)
        if (info /= 0) error stop "scaled columns: unexpected status"
        call assert_close(x, [0.0d0, 2.0d0/41.0d0, 48.0d0/41.0d0] * 1.0d14**(i-1), 1.0d-11)
        a = a * 1.0d-14
    end do
end subroutine test_scaled_columns

subroutine test_active_set_step()
    double precision :: a(5,4), b(5), x(4)
    integer :: info

    a(:,1) = [6.0d0, 0.0d0, 4.0d0, 1.0d0, 6.0d0]
    a(:,2) = [7.0d0, 1.0d0, 5.0d0, 2.0d0, 5.0d0]
    a(:,3) = [4.0d0, 8.0d0, 4.0d0, 8.0d0, 3.0d0]
    a(:,4) = [5.0d0, 0.0d0, 5.0d0, 1.0d0, 5.0d0]
    b = [9.0d0, 0.0d0, 9.0d0, 9.0d0, 9.0d0]
    call nnls(a, b, x, info)
    if (info /= 0) error stop "active-set step: unexpected status"
    call assert_close(x, [144.0d0/11803.0d0, 0.0d0, 180.0d0/407.0d0, 17883.0d0/11803.0d0], 1.0d-11)
end subroutine test_active_set_step

subroutine test_underdetermined()
    double precision :: a(2,3), b(2), x(3)
    integer :: info

    a(:,1) = [1.0d0, 0.0d0]
    a(:,2) = [0.0d0, 1.0d0]
    a(:,3) = [1.0d0, 1.0d0]
    b = [1.0d0, 1.0d0]
    call nnls(a, b, x, info)
    if (info /= 0 .or. any(x < 0.0d0)) error stop "underdetermined: unexpected status or negative solution"
    call assert_close(matmul(a, x), b, 1.0d-12)
end subroutine test_underdetermined

subroutine test_dependent_columns()
    double precision :: a(2,3), b(2), x(3)
    integer :: info

    a(:,1) = [1.0d0, 2.0d0]
    a(:,2) = a(:,1)
    a(:,3) = 0.0d0
    b = [2.0d0, 4.0d0]
    call nnls(a, b, x, info)
    if (info /= 0 .or. any(x < 0.0d0)) error stop "dependent columns: unexpected status or negative solution"
    call assert_close(matmul(a, x), b, 1.0d-12)
    call assert_close(x(3:3), [0.0d0], 1.0d-12)
end subroutine test_dependent_columns

subroutine test_ill_conditioned_fit()
    double precision :: a(3,7), b(3), x(7)
    integer :: info

    a(:,1) = [6.0570872396751782d-4, -4.0901535006915072d-5, 3.8944420819454664d-4]
    a(:,2) = a(:,1)
    a(:,3) = [-0.057704510553853353d0, -0.064821639990299038d0, -0.074196332864847434d0]
    a(:,4) = [0.73076825532548728d0, -0.86182008115797859d0, 0.84918640004789059d0]
    a(:,5) = [-7.2138981422861592d0, 1.1031206782439984d0, -3.1817356361509397d0]
    a(:,6) = [-82.883649438018537d0, 46.231645604681091d0, -31.518763253993921d0]
    a(:,7) = [-642.30028632221399d0, -118.33868638270050d0, 72.803496253993921d0]
    b = [0.43647379342733128d0, 0.036924463373496419d0, -0.37448121890729036d0]
    call nnls(a, b, x, info)
    if (info /= 0 .or. any(x < 0.0d0)) error stop "ill-conditioned fit: unexpected status or negative solution"
    call assert_close(matmul(a, x), b, 1.0d-10)
end subroutine test_ill_conditioned_fit

subroutine test_invalid_input()
    double precision :: a(2,1), b(2), x(1)
    integer :: info

    a = 1.0d0
    b = 1.0d0
    call nnls(a, b(1:1), x, info)
    if (info /= -1) error stop "invalid dimensions: unexpected status"
    b(1) = ieee_value(0.0d0, ieee_quiet_nan)
    call nnls(a, b, x, info)
    if (info /= -2) error stop "non-finite input: unexpected status"
end subroutine test_invalid_input

subroutine test_matrix_fit()
    double precision :: c(3,2), a(3,4), e(2,4), residual(3), dual(2)
    integer :: info, failed_column, j

    c(:,1) = [1.0d0, 0.0d0, 1.0d0]
    c(:,2) = [1.0d0, 1.0d0, 1.0d0]
    a(:,1) = [5.0d0, 3.0d0, 5.0d0]
    a(:,2) = [0.0d0, 1.0d0, 0.0d0]
    a(:,3) = 0.0d0
    a(:,4) = -1.0d0
    call nnls(c, a, e, info, failed_column)
    if (info /= 0 .or. failed_column /= 0) error stop "matrix fit: unexpected status"
    call assert_close(e(:,1), [2.0d0, 3.0d0], 1.0d-12)
    call assert_close(e(:,2), [0.0d0, 1.0d0/3.0d0], 1.0d-12)
    call assert_close(reshape(e(:,3:4), [4]), [0.0d0, 0.0d0, 0.0d0, 0.0d0], 1.0d-12)
    do j = 1, size(a, 2)
        residual = a(:,j) - matmul(c, e(:,j))
        dual = matmul(transpose(c), residual)
        if (any(e(:,j) < 0.0d0) .or. any(dual > 1.0d-12) .or. &
            any(abs(e(:,j) * dual) > 1.0d-12)) error stop "matrix fit: KKT conditions not satisfied"
    end do
end subroutine test_matrix_fit

subroutine test_matrix_invalid_input()
    double precision :: c(2,1), a(2,2), e(1,2), wrong_e(2,2)
    integer :: info, failed_column

    c = 1.0d0
    a = 1.0d0
    call nnls(c, a(1:1,:), e, info, failed_column)
    if (info /= -1 .or. failed_column /= 0 .or. any(e /= 0.0d0)) error stop "matrix dimensions"
    call nnls(c, a, wrong_e, info)
    if (info /= -1 .or. any(wrong_e /= 0.0d0)) error stop "matrix output dimensions"
    call nnls(c, a(:,1:0), e(:,1:0), info)
    if (info /= -1) error stop "zero wavelengths"
    a(1,2) = ieee_value(0.0d0, ieee_quiet_nan)
    call nnls(c, a, e, info, failed_column)
    if (info /= -2 .or. failed_column /= 0 .or. any(e /= 0.0d0)) error stop "matrix non-finite input"
end subroutine test_matrix_invalid_input

subroutine test_matrix_solver_failure()
    double precision :: c(1,1), a(1,2), e(1,2)
    integer :: info, failed_column

    c = 1.0d-200
    a(1,:) = [1.0d-200, 1.0d200]
    call nnls(c, a, e, info, failed_column)
    if (info /= 3 .or. failed_column /= 2 .or. any(e /= 0.0d0)) error stop "matrix solver failure"
end subroutine test_matrix_solver_failure

end program test_linalg
