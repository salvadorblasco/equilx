! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>

! testing for the library

program test

    implicit none
    real t_start, t_stop

    call CPU_TIME(t_start)

    call test_cpluss1()
    call test_cpluss2()
    call test_cpluss3()
    call test_mb1()
    call test_mb2()
    call test_pcf1()
    call test_pcf2()
    call test_pcf3()
    call test_si2()
    call test_jac1()
    call test_jac2()
    call test_lsolve()
    call test_sysolve()
    call test_scale()
    call test_nr1()
    call test_nr2()
    call test_solve_point1()
    call test_solve1()
    call test_case_hexaprotic()
    call test_case_phosphoric()
    call test_case_molybdate()
    call test_case_aluminium()
    call test_case_large()
    call test_case_gallic()
    call test_case_coedta()
    call test_case_caco3()

    call CPU_TIME(t_stop)
    print "(A,F10.6,A)", "done in ", (t_stop-t_start), "s"
contains


subroutine test_case_caco3()
    use calcium_carbonate
    implicit none
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NCP+NSL) :: F
    double precision, dimension(NSL) :: SI
    double precision, dimension(NPT,NSL)     :: CS
    double precision, PARAMETER :: TOLERANCE = 1D-12
    integer i, j, CTRL
    logical passed
    character(len=*), parameter :: test_name = "case Ca carbonate"

    passed = .true.
    call initialize()
    CS(:,:) = 0.0
    C = 0.0
    call initial_guess(C, T, NPT, NCP, NEQ)

    call eq_0s0(C, T, B, P, CS, KS, Q, CTRL, NPT, NCP, NEQ, NSL)

    if (CTRL .ne. 0) then
        if (passed) then
            passed = .false.
            call print_result(test_name, "failed")
       end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,NPT
        call mb_solids(F(1:NCP), C(i,:), P, T(i,:), Q, CS(i,:), NCP, NEQ, NSL)
        call saturation_product(C(i,:), Q, NCP, NSL, SI)
        SI = SI / KS
        do j=1,NSL
            if (SI(j)>1.0) then
                F(NCP+j) = 1.0-SI(j)
            else
                F(NCP+j) = 0.0
            end if
        end do
        if (any(F > TOLERANCE)) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,E10.4)", "   max F =", maxval(F)
        end if
    end do
    if (passed) call print_result(test_name, "passed")
end subroutine test_case_caco3


subroutine test_case_coedta()
    use cobalt_edta
    implicit none
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NCP+NSL) :: F
    double precision, dimension(NSL) :: SI
    double precision, dimension(NPT,NSL)     :: CS
    double precision, PARAMETER :: TOLERANCE = 1D-12
    integer i, j, CTRL
    logical passed
    character(len=*), parameter :: test_name = "case cobalt/edta"

    passed = .true.
    call initialize()
    CS(:,:) = 0.0
    C = 0.0
    call initial_guess(C, T, NPT, NCP, NEQ)

    call eq_0s0(C, T, B, P, CS, KS, Q, CTRL, NPT, NCP, NEQ, NSL)

    if (CTRL .ne. 0) then
        if (passed) then
            passed = .false.
            call print_result(test_name, "failed")
       end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,NPT
        call mb_solids(F(1:NCP), C(i,:), P, T(i,:), Q, CS(i,:), NCP, NEQ, NSL)
        call saturation_product(C(i,:), Q, NCP, NSL, SI)
        SI = SI / KS
        do j=1,NSL
            if (SI(j)>1.0) then
                F(NCP+j) = 1.0-SI(j)
            else
                F(NCP+j) = 0.0
            end if
        end do
        if (any(F > TOLERANCE)) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,E10.4)", "   max F =", maxval(F)
        end if
    end do
    if (passed) call print_result(test_name, "passed")
end subroutine test_case_coedta


subroutine test_case_gallic()
    use gallic
    double precision, dimension(NCP+NEQ) :: C
    double precision, dimension(NEQ)        :: fB  ! the frozen equilibrium constants
    integer,          dimension(NEQ, NCP-1) :: fP  ! The frozen stoichiometry array
    double precision, dimension(NEQ+NCP-1)  :: fC  ! frozen initial estimates (in) and refined values (out)
    double precision, dimension(NCP-1)      :: fT  ! frozen analytical concentration values
    double precision :: Cref
    integer, parameter :: REF=1
    integer j, CTRL
    logical passed
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-14
    character(len=*), parameter :: test_name = "case gallic"
    passed = .true.

    call initialize()

    Cref = tC(1)

    call freeze_stoich(P, fP, REF, NCP, NEQ)
    call freeze_component(Cref, B, fB, T, fT, P, REF, NCP, NEQ)

    fC(1) = 1e-5
    fC(2) = 1e-5

    call solve_point(fC, fT, fB, fP, CTRL, .true., NCP-1, NEQ)
    do j=1,NCP+NEQ
        if (j < REF) C(j) = fC(j)
        if (j > REF) C(j) = fC(j-1)
    end do
    C(REF) = Cref
    T(REF) = Cref + sum(P(:,REF)*C(NCP+1:NCP+NEQ))

    call mb_soluble(F, C, P, T, NCP, NEQ)
    do j=1,NCP
        if (abs(F(j)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I3,A,E10.4)", "   F(",  j, ")=", F(j)
        end if
    end do

    if (passed) call print_result(test_name, "passed")
end subroutine test_case_gallic


subroutine test_case_large()
    use large_fixture
    implicit none
    double precision, dimension(NPT, NCP+NEQ) :: C
    !!double precision, dimension(NCP+NEQ) :: C
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-10
    integer i, j, CTRL
    logical passed
    character(len=*), parameter :: test_name = "case p53/Cu/Zn/Im"
    passed = .true.
    call initialize()
    C = 0.0
    !C(1:4) = (/ 1.3372024240083939e-22, 4.8139229216387732e-07, 2.1914186611768873e-06, 3.1268176522999525e-06 /)
    !call nr(C, T(9,:), B, P, ctrl, NCP, NEQ)

    call initial_guess(C, T, NPT, NCP, NEQ)
    call eq_000(C, T, B, P, CTRL, NPT, NCP, NEQ)
    if (CTRL .ne. 0) then
        if (passed) then
            passed = .false.
            call print_result(test_name, "failed")
        end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,NPT
        call mb_soluble(F, C(i,:), P, T(i,:), NCP, NEQ)
        do j=1,NCP
            if (abs(F(j)) > TOLERANCE) then
                if (passed) call print_result(test_name, "failed")
                passed = .false.
                print "(A,I3,A,I3,A,E10.4)", "   F(", i, ",", j, ")=", F(j)
            end if
        end do
    end do
    if (passed) then
        call print_result(test_name, "passed")
    end if

end subroutine test_case_large

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
        if (passed) call print_result("nr1", "failed")
        passed = .false.
        print *, "  control value returned ", CTRL
        return
    end if
    do i = 1,NCP+NEQ
        if (abs(tC(POINT, i) - C(i))  > TOLERANCE) then
            if (passed) call print_result("nr1", "failed")
            passed = .false.
            print "(A I2 2(A,F10.4))", "  SP=", i, " calc=", C(i), "true=", tC(POINT, i) 
        end if
    end do
    if (passed) then
        call print_result("nr1", "passed")
    end if
end subroutine test_nr1


subroutine test_nr2()
    use simple2
    implicit none
    double precision, dimension(NCP+NEQ) :: C
    double precision, dimension(NSL)     :: CS
    double precision, PARAMETER :: TOLERANCE = 1D-6
    integer i, CTRL, NX
    logical passed
    logical, dimension(NSL) :: do_solids = (/ .true. /)
    integer, PARAMETER :: POINT=3
    passed = .true.
    C = 0.9*tC(POINT,:)
    CS = .9*tCS(POINT,:)
    CTRL=0
    NX = NCP + 1
    call nr_solids(C, CS, T(POINT,:), B, P, KS, Q, ctrl, do_solids, NCP, NEQ, NSL, NX)
    if (CTRL .ne. 0) then
        if (passed) call print_result("nr2", "failed")
        passed = .false.
        print *, "  control value returned ", CTRL
        return
    end if
    do i = 1,NCP+NEQ
        if (abs(tC(POINT, i) - C(i))  > TOLERANCE) then
            if (passed) call print_result("nr2", "failed")
            passed = .false.
            print "(A,I2,2(A,E12.6))", "  SP=", i, " calc=", C(i), "true=", tC(POINT, i) 
        end if
    end do
    do i = 1,NSL
        if (abs(tCS(POINT, i) - CS(i))  > TOLERANCE) then
            if (passed) call print_result("nr2", "failed")
            passed = .false.
            print "(A,I2,2(A,E12.6))", "  SL=", i, " calc=", CS(i), "true=", tCS(POINT, i) 
        end if
    end do
    if (passed) then
        call print_result("nr2", "passed")
    end if
end subroutine test_nr2


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
        call pcf(C, T(i,:), B, P, 1E-9, NCP, NEQ)
        do j=1,NCP
            test = abs(C(j) - tC(i,j))/(C(j) + tC(i,j))
            if (test == 1.0) then
                test = abs(C(j) - tC(i,j))
            end if
            if (test > TOLERANCE) then
                if (passed) call print_result("pcf1", "failed")
                passed = .false.
                print "(A,2I3, 3E12.4)", "  ", i, j, tC(i, j), C(j), test
            end if
        end do
    end do
    if (passed) then
        call print_result("pcf1", "passed")
    end if
end subroutine test_pcf1


subroutine test_pcf2()
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
        call pcf(C, T(i,:), B, P, 1E-9, NCP, NEQ)
        do j=1,NCP
            test = abs(C(j) - tC(i,j))/(C(j) + tC(i,j))
            if (test == 1.0) then
                test = abs(C(j) - tC(i,j))
            end if
            if (test > TOLERANCE) then
                if (passed) call print_result("pcf2", "failed")
                passed = .false.
                print "(A,2I3, 3E12.4)", "  ", i, j, tC(i, j), C(j), test
            end if
        end do
    end do
    if (passed) then
        call print_result("pcf2", "passed")
    end if
end subroutine test_pcf2


subroutine test_pcf3()
    use hexaprotic
    implicit none
    double precision, dimension(NCP+NEQ) :: C
    double precision, dimension(NCP) :: fprev, fpost
    logical passed
    integer i, j
    double precision, PARAMETER :: TOLERANCE = 0.95
    passed = .true.
    call initialize()
    do i=1,NPT
    !do i=2,2
        do j=1,NCP
            if (T(i,j) > 0.0) then
                C(j) = 0.01*T(i,j)
            else
                C(j) = 1D-12
            end if
        end do

        call cpluss(C, B, P, NCP, NEQ)
        call mb_soluble(fprev, C, P, T, NCP, NEQ)
        call pcf(C, T(i,:), B, P, 1E-9, NCP, NEQ)
        call mb_soluble(fpost, C, P, T, NCP, NEQ)
        do j=1,NCP
            if (abs(fpost(j)) > abs(fprev(j))) then
                if (passed) call print_result("pcf3", "failed")
                passed = .false.
                print "(A,2I3,E12.4,A,E12.4)", "  ", i, j, fprev(j), ">",  fpost(j)
            end if
        end do
    end do
    if (passed) then
        call print_result("pcf3", "passed")
    end if
end subroutine test_pcf3


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
        if (passed) call print_result("solve_point1", "failed")
        passed = .false.
        print "(A, I3)", "  control value returned ", CTRL
        return
    end if
    do i = 1,NCP+NEQ
        if (abs(C(i) - tC(POINT, i)) > TOLERANCE) then
            if (passed) call print_result("solve_point1", "failed")
            passed = .false.
            print "(A,I3, 2E12.6)", "  ", i, tC(POINT, i), C(i)
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
    double precision, PARAMETER :: TOLERANCE = 1D-6
    integer i, j, CTRL
    logical passed
    passed = .true.
    call initial_guess(C, T, NPT, NCP, NEQ)
    call eq_000(C, T, B, P, CTRL, NPT, NCP, NEQ)
    do i=1,NPT
        do j=1,NCP+NEQ
            if (abs(C(i,j) - tC(i,j)) > TOLERANCE) then
                if (passed) call print_result("solve1", "failed")
                passed = .false.
                print "(A,I4,2(A,E14.8))","   CP", i, "  ", tC(i,j)," != ", C(i,j)
            end if
        end do
    end do
    if (passed) call print_result("solve1", "passed")
end subroutine test_solve1


subroutine test_case_hexaprotic()
    use hexaprotic
    implicit none
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-10
    integer i, j, CTRL
    logical passed
    character(len=*), parameter :: test_name = "case hexaprotic"
    passed = .true.
    call initialize()
    call initial_guess(C, T, NPT, NCP, NEQ)
    call eq_000(C, T, B, P, CTRL, NPT, NCP, NEQ)
    if (CTRL .ne. 0) then
        if (passed) then
            passed = .false.
            call print_result(test_name, "failed")
        end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,NPT
        call mb_soluble(F, C(i,:), P, T(i,:), NCP, NEQ)
        do j=1,NCP
            if (abs(F(j)) > TOLERANCE) then
                if (passed) call print_result(test_name, "failed")
                passed = .false.
                print "(A,I3,A,I3,A,E10.4)", "   F(", i, ",", j, ")=", F(j)
            end if
        end do
    end do
    if (passed) then
        call print_result("solve3", "passed")
    end if
end subroutine test_case_hexaprotic


subroutine test_case_phosphoric()
    use phosphoric
    implicit none
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-10
    integer i, j, CTRL
    logical passed
    character(len=*), parameter :: test_name = "case phosphoric"

    passed = .true.
    call initialize()
    C(1,1) = 1e-10
    call eq_f00(C, T, B, P, 2, CTRL, NPT, NCP, NEQ)
    if (CTRL .ne. 0) then
        if (passed) then
            passed = .false.
            call print_result(test_name, "failed")
        end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,NPT
        call mb_soluble(F, C(i,:), P, T(i,:), NCP, NEQ)
        do j=1,NCP
            if (abs(F(j)) > TOLERANCE) then
                if (passed) call print_result(test_name, "failed")
                passed = .false.
                print "(A,I3,A,I3,A,E10.4)", "   F(", i, ",", j, ")=", F(j)
            end if
        end do
    end do
    if (passed) then
        call print_result(test_name, "passed")
    end if
end subroutine test_case_phosphoric


subroutine test_case_molybdate()
    use molybdate
    implicit none
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-15
    integer i, j, CTRL
    logical passed
    character(len=*), parameter :: test_name = "case molybdate"

    passed = .true.
    call initialize()
    call initial_guess(C, T, NPT, NCP, NEQ)

    call eq_00i(C, T, B, P, Iref, Iback, charge, DHA, DHB, DHP, CTRL, NPT, NCP, NEQ)
    if (CTRL .ne. 0) then
        if (passed) then
            passed = .false.
            call print_result(test_name, "failed")
        end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,NPT
        call mb_soluble(F, C(i,:), P, T(i,:), NCP, NEQ)
        do j=1,NCP
            if (abs(F(j)) > TOLERANCE) then
                if (passed) call print_result(test_name, "failed")
                passed = .false.
                print "(A,I3,A,I3,A,E10.4)", "   F(", i, ",", j, ")=", F(j)
            end if
        end do
    end do
    if (passed) then
        call print_result(test_name, "passed")
    end if
end subroutine test_case_molybdate


subroutine test_case_aluminium()
    use aluminium
    implicit none
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NPT,NSL)     :: CS
    double precision, PARAMETER :: TOLERANCE = 1D-9
    integer i, j, CTRL
    logical passed
    character(len=*), parameter :: test_name = "case aluminium"

    passed = .true.
    call initialize()
    CS(:,:) = 0.0
    C = 0.0
    C(1,1) = 1e-10
    C(:,2) = T(:,2)

    call eq_fs0(C, T, B, P, CS, KS, Q, 2, CTRL, NPT, NEQ, NCP, NSL)
    if (CTRL .ne. 0) then
        if (passed) then
            passed = .false.
            call print_result(test_name, "failed")
       end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,NPT
        !call mb_solids(F, C(i,:), P, T(i,:), Q, CS(i,:), NCP, NEQ, NSL)
        do j=1,NCP+NEQ
            if (abs(C(i,j)-tC(i,j)) > TOLERANCE) then
                if (passed) call print_result(test_name, "failed")
                passed = .false.
                print "(2(A,I3),2(A,E10.4))", "   C(", i, ",", j, ")=", C(i,j), " should be ", tC(i,j)
            end if
        end do
        do j=1,NSL
            if (abs(CS(i,j)-tCS(i,j)) > TOLERANCE) then
                if (passed) call print_result(test_name, "failed")
                passed = .false.
                print "(2(A,I3),2(A,E10.4))", "  CS(", i, ",", j, ")=", CS(i,j), " should be ", tCS(i,j)
            end if
        end do
    end do
    if (passed) then
        call print_result(test_name, "passed")
    end if
end subroutine test_case_aluminium


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
                if (passed) call print_result("cpluss1", "failed")
                passed = .false.
                print "(A,I3,2E12.6)", "  it", i, tC(i,j), C(j)
            end if
        end do
    end do
    if (passed) call print_result("cpluss1", "passed")
end subroutine test_cpluss1


subroutine test_cpluss2()
    use simple2
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
                print *, "test cpluss2 ....... failed at ", i, tC(i,j), C(j)
                passed = .false.
            end if
        end do
    end do
    if (passed) call print_result("cpluss2", "passed")
end subroutine test_cpluss2


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
                    if (passed) call print_result("jac1", "failed")
                    passed = .false.
                    print "(A,3I4, 2E12.4)", "  ", i, j, k, Jac(j,k), tJ(i,j,k)
                end if
            end do
        end do
    end do
    if (passed) call print_result("jac1", "passed")
end subroutine test_jac1


subroutine test_jac2()
    use simple2
    implicit none
    DOUBLE PRECISION, DIMENSION(NCP+NSL, NCP+NSL) :: Jac
    logical passed
    integer i, j, k, nx
    double precision, PARAMETER:: TOL = 1D-6
    logical, dimension(NSL) :: do_solids = (/.true./)
    passed = .true.
    nx = NCP + 1
    do i=1,NPT
        call jacobian_solid(Jac, tC(i,:), KS, P, Q, NCP, NEQ, NSL, nx, do_solids)
        do j=1,NCP+NSL
            do k=1,NCP+NSL
                if (abs(Jac(j,k)-tJ(i,j,k)) > TOL) then
                    if (passed) call print_result("jac2", "failed")
                    passed = .false.
                    print "(A,3I4,2(A,E12.4))", "  ", i, j, k, " Jcalc=",Jac(j,k), " should be ", tJ(i,j,k)
                end if
            end do
        end do
    end do
    if (passed) call print_result("jac2", "passed")
end subroutine test_jac2


subroutine test_lsolve()
    character(len=*), parameter :: test_name = "lesolve"
    integer, parameter :: N=5
    double PRECISION, dimension(N,N) :: A = reshape( &
     (/0.05877024, 0.28411459, 0.40437397, 0.06793575, 0.03309539, &
       0.5086304 , 0.16324818, 0.41410313, 0.92346009, 0.22741105, &
       0.19705674, 0.08176619, 0.95449347, 0.14038565, 0.51148029, &
       0.43147132, 0.07351833, 0.67603311, 0.53159246, 0.1574564 , &
       0.50511729, 0.63200965, 0.53917495, 0.35968986, 0.79263917/), (/N,N/), order=[2,1])
    double PRECISION, dimension(N) :: b = (/0.20798025, 0.47081858, 0.30973142, 0.46441424, 0.38044028/)
    double PRECISION, dimension(N) :: truex = (/0.45757758,  0.17643674,  0.31928777,  0.1409372 , -0.23345483/)
    double PRECISION, dimension(N) :: calcx
    integer CTRL, i
    logical passed
    double precision, PARAMETER:: TOL = 1D-7
    passed = .true.
    CTRL = 0

    call lesolve(A, b, calcx, N, CTRL)
    if (CTRL .ne. 0) then
        call print_result(test_name, "failed")
        print *, "  CTRL must be 0 (is ", CTRL, ")"
        passed = .false.
    end if
    do i=1,N
        if (abs(calcx(i)-truex(i)) > TOL) then
        call print_result(test_name, "failed")
            print "(A,I3,2F10.4)", "  ", i, calcx(i), truex(i)
            passed = .false.
        end if
    end do
    if (passed) call print_result(test_name, "passed")
end subroutine test_lsolve


subroutine test_sysolve()
    implicit none
    integer, parameter :: N=4
    double PRECISION, dimension(N,N) :: A = reshape( &
      (/ 1.08099127E-03, 1.08098930E-03, 1.07835948E-03, 1.08281530E-03, &
         1.08098930E-03, 1.08098953E-03, 1.07835948E-03, 1.08280751E-03, &
         1.07835948E-03, 1.07835948E-03, 1.07923939E-03, 1.07877450E-03, &
         1.08281530E-03, 1.08280751E-03, 1.07877450E-03, 1.08745856E-03 /), (/N,N/), order=[2,1])
    double PRECISION, dimension(N) :: b = (/ 8.09912257E-05, 8.09893566E-05, 7.92393462E-05, 2.93463972E-04 /)
    double PRECISION, dimension(N) :: calcx
    integer CTRL, i
    logical passed
    double precision, PARAMETER:: TOL = 1D-11
    character(len=*), parameter :: test_name = "sysolve"
    passed = .true.
    CTRL = 0

    call sysolve(A, b, calcx, 4, CTRL)

    if (CTRL .ne. 0) then
        if (passed) call print_result(test_name, "failed")
        print *, "  CTRL must be 0 (is ", CTRL, ")"
        passed = .false.
    end if
    do i=1,N
        if (any(abs(matmul(A,calcx) - b) > TOL)) then
            if (passed) call print_result(test_name, "failed")
            print *, "  ", (matmul(A,calcx) - b)
            passed = .false.
        end if
    end do
    if (passed) call print_result(test_name, "passed")
end subroutine test_sysolve


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
    integer, parameter :: total_width = 30
    
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


subroutine initial_guess(C, T, NPT, NCP, NEQ)
    integer, intent(in) :: NPT
    integer, intent(in) :: NCP
    integer, intent(in) :: NEQ
    double precision, intent(out), dimension(NPT, NCP+NEQ) :: C
    double precision, intent(in),  dimension(NPT, NCP) :: T
    integer i
    C(:,:) = 0.0
    do i = 1,NCP                 ! initialize first point
        if (T(1,i) > 0.0) then
            C(1,i) = 0.01*T(1,i)
        else
            C(1,i) = 1D-10
        end if
    end do
end subroutine initial_guess

end program test

