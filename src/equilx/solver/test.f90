module simple1
    implicit none

    integer, PARAMETER :: NPT = 5
    integer, PARAMETER :: NCP = 2
    integer, PARAMETER :: NEQ = 1
    double precision, dimension(NPT, NCP) :: T = transpose(reshape( &
        (/1.0, 0.0, 0.7, 0.3, 0.5, 0.5,0.3,0.7,0.0,1.0/), &
        (/NCP, NPT/)))
    integer, dimension(NEQ, NCP) :: P = reshape((/1,1/), (/NEQ, NCP/))
    DOUBLE PRECISION, DIMENSION(NEQ) :: B = 1.0
    double precision, dimension(NPT,NCP+NEQ) :: tC = transpose(reshape( &
        (/1.000000000, 0.000000000, 0.000000000, &
          0.588819442, 0.188819442, 0.111180558, &
          0.366025404, 0.366025404, 0.133974596, &
          0.188819442, 0.588819442, 0.111180558, &
          0.000000000, 1.000000000, 0.000000000/), (/NCP+NEQ,NPT/)))
    double precision, dimension(NPT, NCP, NCP) :: tJ = reshape( &
         (/ 1.        ,0.        , &
            0.        ,0.        , &
            0.7       ,0.11118056, &
            0.11118056,0.3       , &
            0.5       ,0.1339746 , &
            0.1339746 ,0.5       , &
            0.3       ,0.11118056, &
            0.11118056,0.7       , &
            0.        ,0.        , &
            0.        ,1. /), (/ NPT, NCP, NCP /), order=[3,2,1])
end module simple1

module simple2
    implicit none

    integer, PARAMETER :: NPT = 5
    integer, PARAMETER :: NCP = 2
    integer, PARAMETER :: NEQ = 1
    integer, PARAMETER :: NSL = 1
    double precision, dimension(NPT, NCP) :: T = transpose(reshape( &
        (/1.0, 0.0, 0.7, 0.3, 0.5, 0.5,0.3,0.7,0.0,1.0/), &
        (/NCP, NPT/)))
    integer, dimension(NEQ, NCP) :: P = reshape((/1,1/), (/NEQ, NCP/))
    integer, dimension(NSL, NCP) :: Q = reshape((/1,1/), (/NSL, NCP/))
    DOUBLE PRECISION, DIMENSION(NEQ) :: B = 1.0
    DOUBLE PRECISION, DIMENSION(NSL) :: KS = 0.12
    double precision, dimension(NPT,NCP+NEQ) :: tC = transpose(reshape( &
     (/1.000000000, 0.000000000, 0.000000000, & 
       0.588819442, 0.188819442, 0.111180558, &
       0.346410162, 0.346410162, 0.120000000, &
       0.188819442, 0.588819442, 0.111180558, &
       0.000000000, 1.000000000, 0.000000000/), (/NCP+NEQ,NPT/)))
    double precision, dimension(NPT,NSL) :: tCS = transpose(reshape( &
      (/0.000000000, 0.000000000, 0.033589838, 0.000000000, 0.000000000/), &
      (/NSL, NPT/)))
end module simple2


module hexaprotic
    implicit none

    integer, PARAMETER :: NPT = 89
    integer, PARAMETER :: NCP = 2
    integer, PARAMETER :: NEQ = 7
    double precision, dimension(NPT) :: V = (/  0.00 , 0.12 , 0.24 , 0.36 , 0.48 , 0.60 , 0.72 , 0.84 , 0.96 , 1.08 , 1.20 , 1.32, &
    1.44 , 1.56 , 1.68 , 1.80 , 1.92 , 2.04 , 2.16 , 2.28 , 2.40 , 2.52 , 2.64 , 2.76 , 2.88 , 3.00 , 3.12 , 3.24 , 3.36 , 3.48, &
    3.60 , 3.72 , 3.84 , 3.96 , 4.08 , 4.20 , 4.32 , 4.44 , 4.56 , 4.68 , 4.80 , 4.92 , 5.04 , 5.16 , 5.28 , 5.40 , 5.52 , 5.64, &
    5.76 , 5.88 , 6.00 , 6.12 , 6.24 , 6.36 , 6.48 , 6.60 , 6.72 , 6.84 , 6.96 , 7.08 , 7.20 , 7.32 , 7.44 , 7.56 , 7.68 , 7.80, &
    7.92 , 8.04 , 8.16 , 8.28 , 8.40 , 8.52 , 8.64 , 8.76 , 8.88 , 9.00 , 9.12 , 9.24 , 9.36 , 9.48 , 9.72 , 10.08, 10.32, 10.56,&
    10.80, 10.92, 11.28, 11.64, 12.00 /)
    double precision :: V0 = 25
    double precision, dimension(NPT, NCP) :: T
    DOUBLE PRECISION, DIMENSION(NEQ) :: B 
    DOUBLE PRECISION, DIMENSION(NEQ) :: logB = (/ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 /) 
    integer, dimension(NEQ, NCP) :: P = reshape((/ 1,1,1,1,1,1,0, 1,2,3,4,5,6,-1 /), &
                                                shape=(/ NEQ, NCP /))
contains

subroutine initialize()
    integer i
    do i=1,NPT
        T(i,1) = 5.0E-3*V0 / (V(i) + V0)
        T(i,2) = (30.0E-3*V0 - 0.1* V(i))/ (V(i) + V0)
    end do
    B = 10.0**logB
end subroutine initialize
end module hexaprotic


module print_utils
    implicit none
    private
    public :: print_matrix

    interface print_matrix
        module procedure print_real_matrix
        module procedure print_integer_matrix
    end interface

contains
    subroutine print_real_matrix(A, name)
        real(kind=8), intent(in) :: A(:,:)
        character(len=*), intent(in), optional :: name
        
        integer :: i, j, rows, cols
        character(len=40) :: data_fmt
        
        rows = size(A, 1)
        cols = size(A, 2)
        
        if (present(name)) then
            print '(/,A)', trim(name)
        end if
        
        ! Column headers
        write(*, '(A6)', advance='no') '     '
        do j = 1, cols
            write(*, '(A9)', advance='no') '     Col'//trim(adjustl(itoa(j)))
        end do
        print *
        
        ! Separator
        write(*, '(A)') repeat('-', 6 + cols*9)
        
        ! Data - one write per row with repeated format
        write(data_fmt, '(A,I0,A)') '(', cols, '(D14.6))'
        
        do i = 1, rows
            write(*, '(I5,A1)', advance='no') i, '|'
            write(*, data_fmt) A(i, :)
        end do
        print *
    end subroutine print_real_matrix


    subroutine print_integer_matrix(A, name)
        integer, intent(in) :: A(:,:)
        character(len=*), intent(in), optional :: name
        
        integer :: i, j, rows, cols
        character(len=40) :: data_fmt
        
        rows = size(A, 1)
        cols = size(A, 2)
        
        if (present(name)) then
            print '(/,A)', trim(name)
        end if
        
        ! Column headers
        write(*, '(A6)', advance='no') '     '
        do j = 1, cols
            write(*, '(A9)', advance='no') '     Col'//trim(adjustl(itoa(j)))
        end do
        print *
        
        ! Separator
        write(*, '(A)') repeat('-', 6 + cols*9)
        
        ! Data
        write(data_fmt, '(A,I0,A)') '(', cols, '(I10))'
        
        do i = 1, rows
            write(*, '(I5,A1)', advance='no') i, '|'
            write(*, data_fmt) A(i, :)
        end do
        print *
    end subroutine print_integer_matrix


    pure function itoa(n) result(str)
        integer, intent(in) :: n
        character(len=12) :: str
        write(str, '(I0)') n
        str = adjustl(str)
    end function itoa

end module print_utils


program test
    implicit none
    call test_cpluss1()
    call test_cpluss3()
    call test_mb1()
    call test_pcf1()
    call test_mb2()
    call test_si2()
    call test_jac1()
    call test_lsolve()
    call test_scale()
    call test_nr1()
    call test_solve_point1()
    call test_solve1()
    call test_solve3()
contains


subroutine test_nr1()
    use simple1
    implicit none
    double precision, dimension(NCP+NEQ) :: C
    double precision, PARAMETER :: TOLERANCE = 1D-7
    integer i, CTRL
    logical passed
    integer, PARAMETER :: POINT=5
    passed = .true.
    C = 0.99*tC(POINT,:)
    C(1) = 1e-12
    call nr(C, T(POINT,:), B, P, CTRL, NCP, NEQ)
    if (CTRL .ne. 0) then
        if (passed) call print_result("scale", "failed")
        passed = .false.
        print *, "  control value returned ", CTRL
        return
    end if
    do i = 1,NCP+NEQ
        if (abs(tC(POINT, i) - C(i))  > TOLERANCE) then
            passed = .false.
            print "(A I2 2(A,F10.4))", "  SP=", i, " calc=", C(i), "true=", tC(POINT, i) 
        end if
    end do
    if (passed) then
        call print_result("nr1", "passed")
    end if
end subroutine test_nr1


subroutine test_pcf1()
    use simple1
    double precision, dimension(NCP+NEQ) :: C
    logical passed
    integer i, j
    double precision :: test
    double precision, PARAMETER :: TOLERANCE = 0.95
    passed = .true.
    do i=1,NPT
        do j=1,NCP
            if (T(i,j) > 0.0) then
                C(j) = 0.01*T(i,j)
            else
                C(j) = 1D-12
            end if
        end do
        call pcf(C, T(i,:), B, P, 1D-9, NCP, NEQ)
        do j=1,NCP
            test = abs(C(j) - tC(i,j))/(C(j) + tC(i,j))
            if (test == 1.0) then
                test = abs(C(j) - tC(i,j))
            end if
            if (test > TOLERANCE) then
                passed = .false.
                print *, "test pcf1 .......... failed at ", i, j, tC(i, j), C(j), test
            end if
        end do
    end do
    if (passed) then
        call print_result("pcf1", "passed")
    end if
end subroutine test_pcf1


subroutine test_solve_point1()
    use simple1
    implicit none
    double precision, dimension(NCP+NEQ) :: C
    double precision, PARAMETER :: TOLERANCE = 1D-8
    integer i, CTRL
    logical passed
    integer, PARAMETER :: POINT=5
    passed = .true.
    do i = 1,NCP                 ! initialize first point
        if (T(POINT,i) > 0.0) then
            C(i) = 0.01*T(POINT,i)
        else
            C(i) = 1D-12
        end if
    end do
    call solve_point(C, T(POINT,:), B, P, CTRL, .true., NCP, NEQ)
    if (CTRL .ne. 0) then
        !         ....|....|....|....v|....|
        print *, "test solve_point1 .. failed -> control value returned ", CTRL
        return
    end if
    do i = 1,NCP+NEQ
        if (abs(C(i) - tC(POINT, i)) > TOLERANCE) then
            passed = .false.
            print *, "test solve_point1 .. failed at ", i, tC(POINT, i), C(i)
        end if
    end do
    if (passed) then
        call print_result("solve_point1", "passed")
    end if
end subroutine test_solve_point1


subroutine test_solve1()
    use simple1
    implicit none
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, PARAMETER :: TOLERANCE = 1D-7
    integer i, j, CTRL
    logical passed
    passed = .true.
    C(:,:) = 0.0
    do i = 1,NCP                 ! initialize first point
        if (T(1,i) > 0.0) then
            C(1,i) = 0.01*T(1,i)
        else
            C(1,i) = 1D-12
        end if
    end do
    call solve(C, T, B, P, CTRL, NPT, NCP, NEQ)
    do i=1,NPT
        do j=1,NCP+NEQ
            if (abs(C(i,j) - tC(i,j)) > TOLERANCE) then
                if (passed) call print_result("solve1", "failed")
                passed = .false.
                print "(A,I4,2(A,E10.4))","   CP", i, "  ", tC(i,j)," != ", C(i,j)
            end if
        end do
    end do
    if (passed) call print_result("solve1", "passed")
end subroutine test_solve1


subroutine test_solve3()
    use hexaprotic
    use print_utils
    implicit none
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-8
    integer i, j, CTRL
    logical passed
    passed = .true.
    call initialize()
    C(:,:) = 0.0
    do i = 1,NCP                 ! initialize first point
        if (T(1,i) > 0.0) then
            C(1,i) = 0.01*T(1,i)
        else
            C(1,i) = 1D-10
        end if
    end do
    !call solve(C, T, B, P, CTRL, NPT, NCP, NEQ)
    call solve(C(1:10,:), T(1:10,:), B, P, CTRL, 10, NCP, NEQ)
    if (CTRL .ne. 0) then
        if (passed) then
            passed = .false.
            call print_result("solve3", "failed")
        end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,10
        call mb_soluble(F, C(i,:), P, T(i,:), NCP, NEQ)
        do j=1,NCP
            if (abs(F(j)) > TOLERANCE) then
                if (passed) then
                    passed = .false.
                    print *, "test solve3 ........ failed"
                end if
                print "(A,I3,A,I3,A,E10.4)", "   F(", i, ",", j, ")=", F(j)
            end if
        end do
    end do
    if (passed) then
        call print_result("solve3", "passed")
    end if
end subroutine test_solve3


subroutine test_cpluss1()
    use simple1
    implicit none
    double precision, dimension(NCP+NEQ) :: C
    double precision, PARAMETER:: TOL = 1e-8
    integer i, j
    logical passed
    passed = .true.
    do i=1,NPT
        C(1:NCP) = tC(i,1:NCP)
        call cpluss(C, B, P, NCP, NEQ)
        do j=1,(NCP+NEQ)
            if (abs(tC(i,j)-C(j)) > TOL) then
                print *, "test cpluss1 ....... failed at ", i, tC(i,j), C(j)
                passed = .false.
            end if
        end do
    end do
    if (passed) call print_result("cpluss1", "passed")
end subroutine test_cpluss1


subroutine test_cpluss3()
    use hexaprotic
    implicit none
    double precision, dimension(NCP+NEQ) :: C
    double precision, dimension(NCP+NEQ) :: tC = (/4.4289E-22, 3.7200E-03, 1.6476E-14, 6.1291E-09, 2.2800E-05, & 
                                                   8.4818E-04, 3.1553E-03, 1.1738E-03, 4.5651E-12/)
    double precision, PARAMETER:: TOL = 1e-6
    integer j
    logical passed
    call initialize()
    passed = .true.
    C(1:2) = (/ 4.4289D-22, 3.7200D-03 /)

    call cpluss(C, B, P, NCP, NEQ)
    do j=1,(NCP+NEQ)
        if (abs(tC(j)-C(j)) > TOL) then
            print *, "test cpluss3 ....... failed at ", j, tC(j), C(j)
            passed = .false.
        end if
    end do
    if (passed) call print_result("cpluss3", "passed")
end subroutine test_cpluss3


subroutine test_mb1()
    use simple1
    use print_utils
    implicit none
    double precision, dimension(NPT,NCP) :: F
    double precision, PARAMETER:: TOL = 1e-7
    logical passed
    integer i,j
    passed = .true.
    do i=1,NPT
        call mb_soluble(F(i,:), tC(i,:), P, T(i,:), NCP, NEQ)
        do j=1,NCP
            if (abs(F(i,j)) > TOL) then
                print *, "test mb1 ........... failed at (", i, ",", j, ")"
                passed = .false.
            end if
        end do
    end do
    if (passed) call print_result("mb1", "passed")
end subroutine test_mb1


subroutine test_mb2()
    use simple2
    implicit none
    double precision, dimension(NPT,NCP) :: F
    double precision, PARAMETER:: TOL = 1e-7
    integer i,j
    logical passed
    passed = .true.
    do i=1,NPT
        call mb_solids(F(i,:), tC(i,:), P, T(i,:), Q, tCS(i,:), NCP, NEQ, NSL)
        do j=1,NCP
            if (abs(F(i,j)) > TOL) then
                print *, "test mb2 ........... failed at (", i, ",", j, ")", F(i,j)
                passed = .false.
            end if
        end do
    end do
    if (passed) call print_result("mb2", "passed")
end subroutine test_mb2


subroutine test_si2()
    use simple2
    implicit none
    double precision, dimension(NPT, NSL) :: SI
    double precision, dimension(NPT) :: tSI = (/0.0, 0.92650466, 1.0, 0.92650466, 0.0 /)
    double precision, PARAMETER:: TOL = 1e-7
    logical passed
    integer i
    passed = .true.
    do i=1,NPT
        call saturation_product(tC(i,:), Q, NCP, NSL, SI(i,:))
        SI(i,1) = SI(i,1) / KS(1)
        if (abs(SI(i,1)-tSI(i)) > TOL) then
            print *, "test si2 ........... failed at", i, SI(i,1)
            passed = .false.
        end if
    end do
    if (passed) call print_result("si2", "passed")
end subroutine test_si2


subroutine test_jac1()
    use simple1
    implicit none
    DOUBLE PRECISION, DIMENSION(NCP, NCP) :: Jac
    logical passed
    integer i, j, k
    double precision, PARAMETER:: TOL = 1D-7
    passed = .true.
    do i=1,NPT
        call jacobian_soluble(Jac, tC(i,:), P, NCP, NEQ)
        do j=1,NCP
            do k=1,NCP
                !print *, i,j,k, tJ(i,j,k)
                if (abs(Jac(j,k)-tJ(i,j,k)) > TOL) then
                    print *, "test jac1 .......... failed at ", i, j, k, Jac(j,k), tJ(i,j,k)
                    passed = .false.
                end if
            end do
        end do
    end do
    if (passed) call print_result("jacsl", "passed")
end subroutine test_jac1


subroutine test_lsolve()
    double PRECISION, dimension(5,5) :: A = reshape( &
     (/0.05877024, 0.28411459, 0.40437397, 0.06793575, 0.03309539, &
       0.5086304 , 0.16324818, 0.41410313, 0.92346009, 0.22741105, &
       0.19705674, 0.08176619, 0.95449347, 0.14038565, 0.51148029, &
       0.43147132, 0.07351833, 0.67603311, 0.53159246, 0.1574564 , &
       0.50511729, 0.63200965, 0.53917495, 0.35968986, 0.79263917/), (/5,5/), order=[2,1])
    double PRECISION, dimension(5) :: b = (/0.20798025, 0.47081858, 0.30973142, 0.46441424, 0.38044028/)
    double PRECISION, dimension(5) :: truex = (/0.45757758,  0.17643674,  0.31928777,  0.1409372 , -0.23345483/)
    double PRECISION, dimension(5) :: calcx
    integer CTRL, i
    logical passed
    double precision, PARAMETER:: TOL = 1D-7
    passed = .true.
    CTRL = 0

    call lsolve(A, b, calcx, 5, CTRL)
    if (CTRL .ne. 0) then
        print *, "test lsolve ........ failed -> CTRL must be 0 (is ", CTRL, ")"
        passed = .false.
    end if
    do i=1,5
        if (abs(calcx(i)-truex(i)) > TOL) then
            print *, "test lsolve ........ failed at ", i, calcx(i), truex(i)
            passed = .false.
        end if
    end do
    if (passed) call print_result("lsolve", "passed")
end subroutine test_lsolve


subroutine test_scale()
    double PRECISION, dimension(5,5) :: A = reshape( &
     (/0.05877024, 0.28411459, 0.40437397, 0.06793575, 0.03309539, &
       0.5086304 , 0.16324818, 0.41410313, 0.92346009, 0.22741105, &
       0.19705674, 0.08176619, 0.95449347, 0.14038565, 0.51148029, &
       0.43147132, 0.07351833, 0.67603311, 0.53159246, 0.1574564 , &
       0.50511729, 0.63200965, 0.53917495, 0.35968986, 0.79263917/), (/5,5/), order=[2,1])
    double PRECISION, dimension(5,5) :: sA = reshape( &
     (/ 1.        , 2.90061861, 1.70733202, 0.38435296, 0.15333841, &
        5.19277386, 1.        , 1.04905403, 3.13476182, 0.63219292, &
        0.83200529, 0.20713958, 1.        , 0.19708194, 0.58803701, &
        2.44109004, 0.24956406, 0.94905653, 1.        , 0.24256779, &
        2.34032269, 1.75695957, 0.61987691, 0.55411643, 1.        /), (/5,5/), order=[2,1])
    double PRECISION, dimension(5) :: b = (/0.20798025, 0.47081858, 0.30973142, 0.46441424, 0.38044028/)
    double PRECISION, dimension(5) :: sb = (/ 0.85791326, 1.16527764, 0.31702884, 0.6369658 , 0.42731558 /)
    double PRECISION, dimension(5) :: sdiag = (/ 0.05877024, 0.16324818, 0.95449347, 0.53159246, 0.79263917 /)
    double PRECISION, dimension(5) :: diag
    integer i, j
    logical passed
    double precision, PARAMETER:: TOL = 1D-6

    passed = .true.
    do i=1,5
        diag(i) = A(i,i)
    end do
    call scale(A, b, diag, 5)
    do i=1,5
        do j=1,5
            if (abs(A(i,j)-sA(i,j)) > TOL) then
                print *, "test scale ......... failed for A(", i,",",j,")=", A(i,j)," should be ", sA(i,j)
                passed = .false.
            end if
        end do
    end do
    do i=1,5
        if (abs(b(i)-sb(i)) > TOL) then
            print *, "test scale ......... failed for b(", i,")=", b(i)," should be ", sb(i)
            passed = .false.
        end if
        if (abs(diag(i)-sdiag(i)) > TOL) then
            print *, "test scale ......... failed for diag(", i,")=", diag(i)," should be ", sdiag(i)
            passed = .false.
        end if
    end do
    if (passed) call print_result("scale", "passed")
end subroutine test_scale


subroutine print_result(test_name, result_tag)
    character(len=*), intent(in) :: result_tag, test_name
    integer, parameter :: total_width = 28
    
    character(len=1) :: filler
    integer :: left_len, right_len, fill_len
    
    filler = "."
    left_len  = len_trim(test_name)
    right_len = len_trim(result_tag)
    
    ! Space + filler + Space = 2 extra characters
    fill_len = total_width - left_len - right_len - 2
    
    if (fill_len < 0) fill_len = 0   ! prevent negative repeat
    
    print "(A, A, A, A, A, A)", &
          "test ", &
          trim(test_name),      " ", &                    ! left + space
          repeat(filler, fill_len), &                ! filling
          " ", trim(result_tag)                           ! space + right
end subroutine

! from collections import namedtuple
! 
! import numpy as np
! 
! 
! data = namedtuple('data', ['analyticalc', 'stoichiometry', 'log_beta',
!                            'solid_stoichiometry', 'log_ks',
!                            'real_freec', 'real_jac', 'gfun'])
! 
! def synthetic_simple1_wsld():
!     analyticalc = np.array([[1.0, 0.0],[0.7, 0.3], [0.5, 0.5],[0.3,0.7],[0.0,1.0]])
!     stoichiometry = np.array([[1,1]], dtype=int)
!     log_beta = np.log10(np.array([1.0]))
!     solid_stoichiometry = np.array([[1,1]], dtype=int)
!     log_ks =   np.log10(np.array([0.12]))
!     #                       c1           c2           c3           s1
!     real_freec = np.array([[1.000000000, 0.000000000, 0.000000000, 0.000000000],
!                            [0.588819442, 0.188819442, 0.111180558, 0.000000000],
!                            [0.346410162, 0.346410162, 0.120000000, 0.033589838],
!                            [0.188819442, 0.588819442, 0.111180558, 0.000000000],
!                            [0.000000000, 1.000000000, 0.000000000, 0.000000000]])
!     real_jac = np.array([[[1.        , 0.        , 0.],
!                          [0.        , 0.        , 0.],
!                          [0.        , 0.        , 0.]],
!                         [[0.7       , 0.11118056, 0.],
!                          [0.11118056, 0.3       , 0.],
!                          [0.        , 0.        , 0.]],
!                         [[0.46641016, 0.12      , 1.],
!                          [0.12      , 0.46641016, 1.],
!                          [-2.80722e-09,-2.80722e-09,0.]],
!                         [[0.3       , 0.11118056, 0.],
!                          [0.11118056, 0.7       , 0.],
!                          [0.        , 0.        , 0.]],
!                         [[0.        , 0.        , 0.],
!                          [0.        , 1.        , 0.],
!                          [0.        , 0.        , 0.]]])
!     gfun = np.array([[0.0, 0.0, 2.80722e-9, 0.0, 0.0]]).T
!     return data(analyticalc, stoichiometry, log_beta, solid_stoichiometry, log_ks,
!                 real_freec, real_jac, gfun)
! 
! 
! def synthetic_simple1_wosld():
!     analyticalc = np.array([[1.0, 0.0],[0.7, 0.3], [0.5, 0.5],[0.3,0.7],[0.0,1.0]])
!     stoichiometry = np.array([[1,1]], dtype=int)
!     log_beta = np.log10(np.array([1.0]))
!     solid_stoichiometry = np.array([[]], dtype=int)
!     log_ks = np.array([])
!     #                       c1           c2           c3         
!     real_freec = np.array([[1.000000000, 0.000000000, 0.000000000],
!                            [0.588819442, 0.188819442, 0.111180558],
!                            [0.366025404, 0.366025404, 0.133974596],
!                            [0.188819442, 0.588819442, 0.111180558],
!                            [0.000000000, 1.000000000, 0.000000000]])
!     real_jac = np.array([[[1.        , 0.        ],
!                          [0.        , 0.        ]],
!                         [[0.7       , 0.11118056],
!                          [0.11118056, 0.3       ]],
!                         [[0.5       , 0.1339746 ],
!                          [0.1339746 , 0.5       ]],
!                         [[0.3       , 0.11118056],
!                          [0.11118056, 0.7       ]],
!                         [[0.        , 0.        ],
!                          [0.        , 1.        ]]])
!     gfun = np.array([[]]).T
!     return data(analyticalc, stoichiometry, log_beta, solid_stoichiometry, log_ks,
!                 real_freec, real_jac, gfun)
! 
end program test

