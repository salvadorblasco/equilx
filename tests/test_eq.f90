! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>

! testing for the library

program test
    use eq
    use solids, only: saturation_product
    use mb, only: mb_soluble, mb_solids
    use mcpluss, only: cpluss
    use mnr, only: jacobian_soluble, jacobian_solid, sysolve, lesolve

    implicit none (type, external)

    logical :: passed
    real :: t_start, t_stop

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
    call test_nr1()
    call test_nr2()
    call test_solve_point1()
    call test_case_hexaprotic()
    call test_case_phosphoric()
    call test_case_molybdate()
    call test_case_aluminium()
    call test_case_large()
    call test_case_gallic()
    call test_case_coedta()
    call test_case_caco3()
    call test_special_phosphoric()

    call CPU_TIME(t_stop)
    print "(A,F10.6,A)", "done in ", (t_stop-t_start), "s"

contains

subroutine test_case_caco3()
    use calcium_carbonate
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NCP+NSL) :: F
    double precision, dimension(NSL) :: SI
    double precision, dimension(NPT,NSL)     :: CS
    double precision, PARAMETER :: TOLERANCE = 1D-12
    integer :: i, j, CTRL
    character(len=*), parameter :: test_name = "case Ca carbonate"

    passed = .true.
    call initialize()
    CS(:,:) = 0.0
    C = 0.0
    call initial_guess(C, T, NPT, NCP, NEQ)

    call eq_0s0(C, T, B, P, CS, KS, Q, CTRL, NPT, NCP, NEQ, NSL)

    if (CTRL /= 0) then
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
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NCP+NSL) :: F
    double precision, dimension(NSL) :: SI
    double precision, dimension(NPT,NSL)     :: CS
    double precision, PARAMETER :: TOLERANCE = 1D-12
    integer :: i, j, CTRL
    character(len=*), parameter :: test_name = "case cobalt/edta"

    passed = .true.
    call initialize()
    CS(:,:) = 0.0
    C = 0.0
    call initial_guess(C, T, NPT, NCP, NEQ)

    call eq_0s0(C, T, B, P, CS, KS, Q, CTRL, NPT, NCP, NEQ, NSL)

    if (CTRL /= 0) then
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
    integer :: j, CTRL
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
    double precision, dimension(NPT, NCP+NEQ) :: C
    !!double precision, dimension(NCP+NEQ) :: C
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-10
    integer :: i, j, CTRL
    character(len=*), parameter :: test_name = "case p53/Cu/Zn/Im"
    passed = .true.
    call initialize()
    C = 0.0
    !C(1:4) = (/ 1.3372024240083939e-22, 4.8139229216387732e-07, 2.1914186611768873e-06, 3.1268176522999525e-06 /)
    !call nr(C, T(9,:), B, P, ctrl, NCP, NEQ)

    call initial_guess(C, T, NPT, NCP, NEQ)
    call eq_000(C, T, B, P, CTRL, NPT, NCP, NEQ)
    if (CTRL /= 0) then
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
    double precision, dimension(NCP+NEQ) :: C
    double precision, PARAMETER :: TOLERANCE = 1D-7
    integer :: i, CTRL
    integer, PARAMETER :: POINT=5
    passed = .true.
    C = 0.99*tC(POINT,:)
    C(1) = 1e-12
    call nr(C, T(POINT,:), B, P, CTRL, NCP, NEQ)
    if (CTRL /= 0) then
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
    double precision, dimension(NCP+NEQ) :: C
    double precision, dimension(NSL)     :: CS
    double precision, PARAMETER :: TOLERANCE = 1D-6
    integer :: i, CTRL, NX
    logical, dimension(NSL) :: do_solids
    integer, PARAMETER :: POINT=3
    passed = .true.
    do_solids = [ .true. ]
    C = 0.9*tC(POINT,:)
    CS = .9*tCS(POINT,:)
    CTRL=0
    NX = NCP + 1
    call nr_solids(C, CS, T(POINT,:), B, P, KS, Q, ctrl, do_solids, NCP, NEQ, NSL, NX)
    if (CTRL /= 0) then
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
    integer :: i, j
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
        call pcf(C, DBLE(T(i,:)), B, P, 1D-9, NCP, NEQ)
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
    integer :: i, j
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
        call pcf(C, DBLE(T(i,:)), B, P, 1D-9, NCP, NEQ)
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
    double precision, dimension(NCP+NEQ) :: C
    double precision, dimension(NCP) :: fprev, fpost
    integer :: i, j
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
        call pcf(C, T(i,:), B, P, 1D-9, NCP, NEQ)
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
    double precision, dimension(NCP+NEQ) :: C
    double precision, PARAMETER :: TOLERANCE = 1D-8
    integer :: i, CTRL
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
    if (CTRL /= 0) then
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


subroutine test_case_hexaprotic()
    use hexaprotic
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-15
    integer :: i, j, CTRL
    character(len=*), parameter :: test_name = "case hexaprotic"
    passed = .true.
    call initialize()
    call initial_guess(C, T, NPT, NCP, NEQ)
    call eq_000(C, T, B, P, CTRL, NPT, NCP, NEQ)
    if (CTRL /= 0) then
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
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-10
    integer :: i, j, CTRL
    character(len=*), parameter :: test_name = "case phosphoric"

    passed = .true.
    call initialize()
    C(1,1) = 1e-10
    call eq_f00(C, T, B, P, 2, CTRL, NPT, NCP, NEQ)
    if (CTRL /= 0) then
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
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-15
    integer :: i, j, CTRL
    character(len=*), parameter :: test_name = "case molybdate"

    passed = .true.
    call initialize()
    call initial_guess(C, T, NPT, NCP, NEQ)

    call eq_00i(C, T, B, P, Iref, Iback, charge, DHA, DHB, DHP, CTRL, NPT, NCP, NEQ)
    if (CTRL /= 0) then
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
    double precision, dimension(NPT, NCP+NEQ) :: C
    double precision, dimension(NPT,NSL)     :: CS
    double precision, dimension(NCP+NSL)     :: F
    double precision, PARAMETER :: TOLERANCE = 1D-9
    integer :: i, j, CTRL
    character(len=*), parameter :: test_name = "case aluminium"

    passed = .true.
    call initialize()
    CS(:,:) = 0.0
    C = 0.0
    C(1,1) = 1e-10
    C(:,2) = T(:,2)

    call eq_fs0(C, T, B, P, CS, KS, Q, 2, CTRL, NPT, NEQ, NCP, NSL)
    if (CTRL /= 0) then
        if (passed) then
            passed = .false.
            call print_result(test_name, "failed")
       end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,NPT
        call mb_solids(F, C(i,:), P, T(i,:), Q, CS(i,:), NCP, NEQ, NSL)
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
    double precision, dimension(NCP+NEQ) :: C
    double precision, PARAMETER:: TOL = 1e-8
    integer :: i, j
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
    double precision, dimension(NCP+NEQ) :: C
    double precision, PARAMETER:: TOL = 1e-8
    integer :: i, j
    passed = .true.
    do i=1,NPT
        C(1:NCP) = tC(i,1:NCP)
        call cpluss(C, B, P, NCP, NEQ)
        do j=1,(NCP+NEQ)
            if (abs(tC(i,j)-C(j)) > TOL) then
                if (passed) call print_result("cpluss2", "failed")
                print *, "  at ", i, tC(i,j), C(j)
                passed = .false.
            end if
        end do
    end do
    if (passed) call print_result("cpluss2", "passed")
end subroutine test_cpluss2


subroutine test_cpluss3()
    use hexaprotic
    double precision, dimension(NCP+NEQ) :: C
    double precision, dimension(NCP+NEQ) :: tC
    double precision, PARAMETER:: TOL = 1e-6
    integer :: j
    call initialize()
    passed = .true.
    C(1:2) = [ 4.4289D-22, 3.7200D-03 ]
    tC = [4.4289E-22, 3.7200E-03, 1.6476E-14, 6.1291E-09, 2.2800E-05, 8.4818E-04, 3.1553E-03, 1.1738E-03, 4.5651E-12]

    call cpluss(C, B, P, NCP, NEQ)
    do j=1,(NCP+NEQ)
        if (abs(tC(j)-C(j)) > TOL) then
            if (passed) call print_result("cpluss3", "failed")
            print *, "  at ", j, tC(j), C(j)
            passed = .false.
        end if
    end do
    if (passed) call print_result("cpluss3", "passed")
end subroutine test_cpluss3


subroutine test_mb1()
    use simple1
    double precision, dimension(NPT,NCP) :: F
    double precision, PARAMETER:: TOL = 1e-7
    integer :: i,j
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
    double precision, dimension(NPT,NCP) :: F
    double precision, PARAMETER:: TOL = 1e-7
    integer :: i,j
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
    double precision, dimension(NPT, NSL) :: SI
    double precision, dimension(NPT) :: tSI
    double precision, PARAMETER:: TOL = 1e-7
    integer :: i
    passed = .true.
    tSI = [0.0, 0.92650466, 1.0, 0.92650466, 0.0 ]
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
    DOUBLE PRECISION, DIMENSION(NCP, NCP) :: Jac
    integer :: i, j, k
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
    DOUBLE PRECISION, DIMENSION(NCP+NSL, NCP+NSL) :: Jac
    integer :: i, j, k, nx
    double precision, PARAMETER:: TOL = 1D-6
    logical, dimension(NSL) :: do_solids
    do_solids = [.true.]
    passed = .true.
    nx = NCP + 1
    do i=1,NPT
        call jacobian_solid(Jac, tC(i,:), tCS(i,:), KS, P, Q, NCP, NEQ, NSL, nx, do_solids)
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
    double PRECISION, dimension(N,N) :: A
    double PRECISION, dimension(N) :: b
    double PRECISION, dimension(N) :: truex
    double PRECISION, dimension(N) :: calcx
    integer :: CTRL, i
    double precision, PARAMETER:: TOL = 1D-7
    passed = .true.
    CTRL = 0
    A = reshape([0.05877024, 0.28411459, 0.40437397, 0.06793575, 0.03309539, &
                 0.5086304 , 0.16324818, 0.41410313, 0.92346009, 0.22741105, &
                 0.19705674, 0.08176619, 0.95449347, 0.14038565, 0.51148029, &
                 0.43147132, 0.07351833, 0.67603311, 0.53159246, 0.1574564 , &
                 0.50511729, 0.63200965, 0.53917495, 0.35968986, 0.79263917], [N,N], order=[2,1])
    b = [0.20798025, 0.47081858, 0.30973142, 0.46441424, 0.38044028]
    truex = [0.45757758,  0.17643674,  0.31928777,  0.1409372 , -0.23345483]

    call lesolve(A, b, calcx, N, CTRL)
    if (CTRL /= 0) then
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
    integer, parameter :: N=4
    double PRECISION, dimension(N,N) :: A
    double PRECISION, dimension(N) :: b
    double PRECISION, dimension(N) :: calcx
    integer :: CTRL, i
    double precision, PARAMETER:: TOL = 1D-11
    character(len=*), parameter :: test_name = "sysolve"
    passed = .true.
    CTRL = 0

    A = reshape([1.08099127E-03, 1.08098930E-03, 1.07835948E-03, 1.08281530E-03, &
                 1.08098930E-03, 1.08098953E-03, 1.07835948E-03, 1.08280751E-03, &
                 1.07835948E-03, 1.07835948E-03, 1.07923939E-03, 1.07877450E-03, &
                 1.08281530E-03, 1.08280751E-03, 1.07877450E-03, 1.08745856E-03], [N,N], order=[2,1])
    b = [8.09912257E-05, 8.09893566E-05, 7.92393462E-05, 2.93463972E-04]

    call sysolve(A, b, calcx, 4, CTRL)

    if (CTRL /= 0) then
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
end subroutine print_result


subroutine test_special_phosphoric()
    integer, parameter :: NCP = 2
    integer, parameter :: NEQ = 4
    integer, parameter :: NPT=50
    double precision, dimension(NCP) :: F
    double precision, PARAMETER :: TOLERANCE = 1D-10
    integer :: i, j, CTRL
    character(len=*), parameter :: test_name = "case special phosph."

    double precision :: freec(50,6), analc(50,2), beta(4)
    integer :: stoich(4,2)

    freec = 0.0
    freec(:,1) = [ 1.00000000e-14, 9.80392157e-15, 9.61538462e-15, 9.43396226e-15, 9.25925926e-15, 9.09090909e-15, &
                   8.92857143e-15, 8.77192982e-15, 8.62068966e-15, 8.47457627e-15, 8.33333333e-15, 8.19672131e-15, &
                   8.06451613e-15, 7.93650794e-15, 7.81250000e-15, 7.69230769e-15, 7.57575758e-15, 7.46268657e-15, &
                   7.35294118e-15, 7.24637681e-15, 7.14285714e-15, 7.04225352e-15, 6.94444444e-15, 6.84931507e-15, &
                   6.75675676e-15, 6.66666667e-15, 6.57894737e-15, 6.49350649e-15, 6.41025641e-15, 6.32911392e-15, &
                   6.25000000e-15, 6.17283951e-15, 6.09756098e-15, 6.02409639e-15, 5.95238095e-15, 5.88235294e-15, &
                   5.81395349e-15, 5.74712644e-15, 5.68181818e-15, 5.61797753e-15, 5.55555556e-15, 5.49450549e-15, &
                   5.43478261e-15, 5.37634409e-15, 5.31914894e-15, 5.26315789e-15, 5.20833333e-15, 5.15463918e-15, &
                   5.10204082e-15, 5.05050505e-15 ]
    freec(:,2) = [ 5.00000000e-11, 4.70588235e-11, 4.42307692e-11, 4.15094340e-11, 3.88888889e-11, 3.63636364e-11, &
                   3.39285714e-11, 3.15789474e-11, 2.93103448e-11, 2.71186441e-11, 2.50000000e-11, 2.29508197e-11, &
                   2.09677419e-11, 1.90476190e-11, 1.71875000e-11, 1.53846154e-11, 1.36363636e-11, 1.19402985e-11, &
                   1.02941176e-11, 8.69565217e-12, 7.14285714e-12, 5.63380282e-12, 4.16666667e-12, 2.73972603e-12, &
                   1.35135135e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, &
                   1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, &
                   1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, &
                   1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, 1.00000000e-12, &
                   1.00000000e-12, 1.00000000e-12 ]
    analc(:,1) = [ 1.00000000e-08, 9.80392157e-09, 9.61538462e-09, 9.43396226e-09, 9.25925926e-09, 9.09090909e-09, &
                   8.92857143e-09, 8.77192982e-09, 8.62068966e-09, 8.47457627e-09, 8.33333333e-09, 8.19672131e-09, &
                   8.06451613e-09, 7.93650794e-09, 7.81250000e-09, 7.69230769e-09, 7.57575758e-09, 7.46268657e-09, &
                   7.35294118e-09, 7.24637681e-09, 7.14285714e-09, 7.04225352e-09, 6.94444444e-09, 6.84931507e-09, &
                   6.75675676e-09, 6.66666667e-09, 6.57894737e-09, 6.49350649e-09, 6.41025641e-09, 6.32911392e-09, &
                   6.25000000e-09, 6.17283951e-09, 6.09756098e-09, 6.02409639e-09, 5.95238095e-09, 5.88235294e-09, &
                   5.81395349e-09, 5.74712644e-09, 5.68181818e-09, 5.61797753e-09, 5.55555556e-09, 5.49450549e-09, &
                   5.43478261e-09, 5.37634409e-09, 5.31914894e-09, 5.26315789e-09, 5.20833333e-09, 5.15463918e-09, &
                   5.10204082e-09, 5.05050505e-09 ]
    analc(:,2)= [ 5.00000000e-05, 4.70588235e-05, 4.42307692e-05, 4.15094340e-05, 3.88888889e-05, 3.63636364e-05, &
                  3.39285714e-05, 3.15789474e-05, 2.93103448e-05, 2.71186441e-05, 2.50000000e-05, 2.29508197e-05, &
                  2.09677419e-05, 1.90476190e-05, 1.71875000e-05, 1.53846154e-05, 1.36363636e-05, 1.19402985e-05, &
                  1.02941176e-05, 8.69565217e-06, 7.14285714e-06, 5.63380282e-06, 4.16666667e-06, 2.73972603e-06, &
                  1.35135135e-06, 0.00000000e+00, -1.31578947e-06, -2.59740260e-06, -3.84615385e-06, -5.06329114e-06, &
                  -6.25000000e-06, -7.40740741e-06, -8.53658537e-06, -9.63855422e-06, -1.07142857e-05, -1.17647059e-05, &
                  -1.27906977e-05, -1.37931034e-05, -1.47727273e-05, -1.57303371e-05, -1.66666667e-05, -1.75824176e-05, &
                  -1.84782609e-05, -1.93548387e-05, -2.02127660e-05, -2.10526316e-05, -2.18750000e-05, -2.26804124e-05, &
                  -2.34693878e-05, -2.42424242e-05 ]
    beta = [2.39883292e+12, 3.78442585e+19, 5.32108259e+21, 1.86208714e-14]
    stoich(:,1) = [1,1,1,0]
    stoich(:,2) = [1,2,3,-1]

    passed = .true.
    call eq_f00(freec, analc, beta, stoich, 2, CTRL, 50, 2, 4)
    if (CTRL /= 0) then
        if (passed) then
            passed = .false.
            call print_result(test_name, "failed")
        end if
       print *, "  control value returned ", CTRL
    end if
    do i=1,NPT
        call mb_soluble(F, freec(i,:), stoich, analc(i,:), NCP, NEQ)
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
end subroutine test_special_phosphoric


subroutine initial_guess(C, T, NPT, NCP, NEQ)
    integer, intent(in) :: NPT
    integer, intent(in) :: NCP
    integer, intent(in) :: NEQ
    double precision, intent(out), dimension(NPT, NCP+NEQ) :: C
    double precision, intent(in),  dimension(NPT, NCP) :: T
    integer :: i
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

