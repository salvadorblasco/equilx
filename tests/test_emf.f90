program test_libemf
    implicit none (type, external)
    real :: t_start, t_stop

    call CPU_TIME(t_start)

    print "(A)", "+-------------------------+"
    print "(A)", "|  Testing libemf module  |"
    print "(A)", "+-------------------------+"
    print *, ""

    !call test_gradient()
    !call test_hexaprotic1()
    call test_hexaprotic2()
    !call test_hexaprotic3()
    !call test_hexaprotic4()
    !call test_hexaprotic5()
    !call test_hexaprotic6()
    !call test_hexaprotic7()
    !call test_hexaprotic8()
    !call test_hexaprotic9()
    !call test_hexaprotic10()
    !call test_hexaprotic11()
    !call test_hexaprotic12()
    !call test_hexaprotic13()
    !call test_hexaprotic14()
    !call test_hexaprotic15()
    !call test_znedta1()

    call CPU_TIME(t_stop)
    print "(A,F10.6,A)", "done in ", (t_stop-t_start), "s"
contains

subroutine test_hexaprotic1()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    character(len=*), parameter :: test_name = "hexaprotic simu1"
    logical :: passed
    integer :: i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.02

    passed = .true.

    titr_flags = 0
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic1


subroutine test_hexaprotic2()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    character(len=*), parameter :: test_name = "hexaprotic free buret"
    logical :: passed
    integer :: i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .false.
    real, parameter :: TOLERANCE = 0.05

    passed = .true.

    titr_flags = 0
    titr_flags(2,2) = 1
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, -0.098 ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic2


subroutine test_hexaprotic3()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    character(len=*), parameter :: test_name = "hexaprotic constr buret"
    logical :: passed
    integer :: i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    real, parameter :: TOLERANCE = 0.05
    logical, parameter :: SILENT = .true.

    passed = .true.

    titr_flags = 0
    titr_flags(2,2) = 2
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, -0.098 ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic3


subroutine test_hexaprotic4()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic E0=+1mV fixed"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    real, parameter :: TOLERANCE = 0.20
    logical, parameter :: SILENT = .true.
    real, parameter :: new_emf0 = 406.0

    passed = .true.

    titr_flags = 0
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]

    call load_titration(titre1, remf1, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic4


subroutine test_hexaprotic5()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic E0=+1mV refined"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    real, parameter :: TOLERANCE = 0.05
    logical, parameter :: SILENT = .true.
    real, parameter :: new_emf0 = 406.0

    passed = .true.

    titr_flags = 0
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]

    call load_titration(titre1, remf1, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 1, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 1, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 1, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 1, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 1, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 1, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic5


subroutine test_hexaprotic6()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic E0=+1mV constrained"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    real, parameter :: TOLERANCE = 0.05
    logical, parameter :: SILENT = .true.
    real, parameter :: new_emf0 = 406.0

    passed = .true.

    titr_flags = 0
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]

    call load_titration(titre1, remf1, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 2, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 2, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 2, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 2, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 2, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, new_emf0, melectr_slope, &
                        titr_flags, 2, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic6


subroutine test_hexaprotic7()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic no junction"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(0) :: junc_flag
    integer, dimension(0) :: junc_comp
    real, dimension(0) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.

    passed = .true.

    titr_flags = 0
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 0, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 0, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > 0.02) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic7


subroutine test_hexaprotic8()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic junction constrained"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.05

    passed = .true.

    titr_flags = 0
    junc_flag = 2
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(30, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic8


subroutine test_hexaprotic9()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic H refined"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.05

    passed = .true.

    titr_flags = 0
    titr_flags(2,1) = 1
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic9


subroutine test_hexaprotic10()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic H constrained"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.06

    passed = .true.

    titr_flags = 0
    titr_flags(2,1) = 2
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic10


subroutine test_hexaprotic11()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic fix CA-2% ref H"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.13

    passed = .true.

    titr_flags = 0
    titr_flags(2,1) = 1
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)*0.98
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic11


subroutine test_hexaprotic12()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic free A, H"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.02

    passed = .true.

    titr_flags = 0
    titr_flags(:,1) = 1
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic12


subroutine test_hexaprotic13()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic A-2% free A, H"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.02

    passed = .true.

    titr_flags = 0
    titr_flags(:,1) = 1
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(1,1) = titr_init0(1,1)*0.98
    titr_params(2,1) = titr_init0(2,1)
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()
    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic13


subroutine test_hexaprotic14()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic fixed ionic"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.02
    real, dimension(2+3*7) :: dhparms
    integer, dimension(2+3*7) :: dhp_flags
    real, dimension(7) :: init_ions, buret_ions

    passed = .true.

    dhparms = 0.0
    dhparms(1) = 0.5115
    dhparms(2) = 1.5
    dhp_flags = 0

    titr_flags = 0
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()

    init_ions = 0.1
    buret_ions = 0.1
    call set_ionic_constant(charge, 0.1, init_ions, buret_ions, dhparms)

    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic14


subroutine test_hexaprotic15()
    use emf
    use case_hexaprotic
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic variable ionic"
    logical passed
    integer i
    integer, DIMENSION(7) :: rlbflags = [1,1,1,1,1,1,0]
    real, dimension(7) :: ini_logbeta = [ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 ]
    real, dimension(2,2) :: titr_params
    integer, dimension(2,2) :: titr_flags
    integer, dimension(1) :: junc_flag
    integer, dimension(1) :: junc_comp
    real, dimension(1) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.05
    real, dimension(2+3*7) :: dhparms
    integer, dimension(2+3*7) :: dhp_flags
    real, dimension(7) :: init_ions, buret_ions

    passed = .true.

    dhparms = 0.0
    dhparms(1) = 0.5115
    dhparms(2) = 1.5
    dhp_flags = 0

    titr_flags = 0
    junc_flag = 0
    junc_comp = 2
    junc_pot = electr_junctionA
    do i=1,6
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = ini_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, 20)
    else
        call init_module(mNPT, 6, mP, ini_logbeta, rlbflags, 1, output_unit)
    end if

    titr_params(:,1) = titr_init0(1,:)
    titr_params(:,2) = [ 0.0, buret ]
    call load_titration(titre1, remf1, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if

    call load_titration(titre1b, remf1b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(2,:)
    call load_titration(titre2, remf2, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre2b, remf2b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    titr_params(:,1) = titr_init0(3,:)
    call load_titration(titre3, remf3, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    call load_titration(titre3b, remf3b, titr_params, 2, initvol, melectr_emf0, melectr_slope, &
                        titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    if (SILENT) then
        close(20)
    end if

    call set_weight_gradient()

    init_ions = 0.1
    buret_ions = 0.1
    call set_ionic_variable(charge, 0.1, init_ions, buret_ions, dhparms)

    call fit_emf(20, temp)

    do i=1,6
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I1,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic15


subroutine test_gradient()
    use misc, only: derivative_1d
    implicit none
    character(len=*), parameter :: test_name = "gradient"
    logical passed
    real, dimension(89) :: titre = [ .00, .12, .24, .36, .48, .60, .72, .84, .96, 1.08, 1.20, 1.32, 1.44, 1.56, 1.68,&
        1.80, 1.92, 2.04, 2.16, 2.28, 2.40, 2.52, 2.64, 2.76, 2.88, 3.00, 3.12, 3.24, 3.36, 3.48, 3.60, 3.72, 3.84, &
        3.96, 4.08, 4.20, 4.32, 4.44, 4.56, 4.68, 4.80, 4.92, 5.04, 5.16, 5.28, 5.40, 5.52, 5.64, 5.76, 5.88, 6.00, &
        6.12, 6.24, 6.36, 6.48, 6.60, 6.72, 6.84, 6.96, 7.08, 7.20, 7.32, 7.44, 7.56, 7.68, 7.80, 7.92, 8.04, 8.16, &
        8.28, 8.40, 8.52, 8.64, 8.76, 8.88, 9.00, 9.12, 9.24, 9.36, 9.48, 9.72, 10.08, 10.32, 10.56, 10.80, 10.92, &
        11.28, 11.64, 12.00 ]
    real, dimension(89) :: emf = [ 265.6, 263.5, 261.4, 259.1, 256.7, 254.2, 251.6, 248.8, 245.9, 242.9, 239.7, 236.3, &
        232.8, 229.1, 225.2, 221.1, 216.9, 212.4, 207.8, 203.0, 197.9, 192.7, 187.2, 181.5, 175.5, 169.2, 162.4, 155.0, &
        146.7, 137.1, 125.6, 112.2, 98.0, 84.9, 73.6, 63.6, 54.2, 45.1, 35.7, 25.5, 14.0, .6, -14.0, -27.9, -40.1, -50.6, &
        -60.2, -69.4, -78.6, -88.4, -99.2, -111.6, -125.6, -139.7, -152.3, -163.0, -172.5, -181.1, -189.3, -197.2, -205.1, &
        -212.9, -220.4, -227.5, -234.0, -239.7, -244.8, -249.2, -253.0, -256.5, -259.5, -262.3, -264.8, -267.1, -269.1, &
        -271.0, -272.9, -274.5, -276.0, -277.5, -280.2, -283.7, -285.8, -287.6, -289.4, -290.2, -292.5, -294.6, -296.5 ]
    real, dimension(89) :: w = [ 10.78, 10.78, 10.75, 10.70, 10.67, 10.63, 10.58, 10.51, 10.48, 10.42, 10.33, 10.26, 10.20, 10.10,&
        10.00, 9.93, 9.82, 9.71, 9.64, 9.49, 9.38, 9.27, 9.12, 8.98, 8.81, 8.58, 8.25, 7.82, 7.19, 6.29, 5.31, 4.69, 4.76, 5.46, &
        6.27, 6.77, 7.01, 7.02, 6.71, 6.16, 5.35, 4.62, 4.52, 5.06, 5.90, 6.59, 6.93, 7.04, 6.88, 6.44, 5.77, 4.98, 4.59, 4.90, &
        5.73, 6.55, 7.12, 7.49, 7.68, 7.76, 7.77, 7.89, 8.09, 8.39, 8.83, 9.24, 9.61, 9.97, 10.17, 10.36, 10.50, 10.59, 10.68, &
        10.77, 10.83, 10.83, 10.88, 10.93, 10.94, 10.96, 10.99, 11.02, 11.04, 11.05, 11.05, 11.06, 11.07, 11.08, 11.08 ]
    real, parameter :: errV=0.003
    real, parameter :: errE=0.3
    real, parameter :: TOL=0.99
    real, dimension(89) :: dy
    real, dimension(89) :: calcw
    integer i

    passed = .true.

    call derivative_1d(89, titre, emf, dy)
    calcw = 1.0/(dy**2*errV**2 + errE**2)

    do i=1,89
    if (abs(calcw(i)-w(i)) > TOL) then
        if (passed) call print_result(test_name, 'failed')
        print '(A,2F6.2)', "  ", calcw(i), w(i)
        passed = .false.
    end if
    end do

    if (passed) call print_result(test_name, 'passed')

end subroutine test_gradient


subroutine print_result(test_name, result_tag)
    character(len=*), intent(in) :: result_tag, test_name
    integer, parameter :: total_width = 45
    
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


subroutine test_znedta1()
    use emf
    use case_znedta
    use iso_fortran_env, only: output_unit

    character(len=*), parameter :: test_name = "znedta refine beta, no jH"
    logical :: passed
    integer :: i
    integer, DIMENSION(15) :: rlbflags = [0,0,0,0,0,0,0,0,0,0,0,1,1,1,0]
    real, dimension(15) :: ini_logbeta
    real, dimension(3,2) :: titr_params
    integer, dimension(3,2) :: titr_flags
    integer, dimension(0) :: junc_flag
    integer, dimension(0) :: junc_comp
    real, dimension(0) :: junc_pot
    real :: rnd
    logical, parameter :: SILENT = .true.
    real, parameter :: TOLERANCE = 0.05

    passed = .true.

    titr_flags = 0
    !junc_flag = 0
    !junc_comp = 2
    !junc_pot = electr_junctionA
    ini_logbeta = real_logbeta
    do i=12,14
        call RANDOM_NUMBER(rnd)
        ini_logbeta(i) = real_logbeta(i) + (rnd-0.5)
    end do

    if (SILENT) then
        open(unit=20, file='/dev/null', status='old')
        call init_module(mNPT, 4, mP, ini_logbeta, rlbflags, 0, 20)
    else
        call init_module(mNPT, 4, mP, ini_logbeta, rlbflags, 0, output_unit)
    end if

    do i=1,4
        titr_params(:,1) = titr_init0(i,:)
        titr_params(:,2) = [ 0.0, 0.0, buret ]
        call load_titration(titre1, remf_nja(:,i), titr_params, 3, initvol, melectr_emf0, melectr_slope, &
                            titr_flags, 0, 0, error_emf, error_vol, junc_pot, junc_comp, junc_flag)
        if (trim(emf_diagnostic) /= 'no error') then
            if (passed) call print_result(test_name, "failed")
            print "(A,I2,2A)", "  loading titration", i, "-> ", emf_diagnostic
            passed=.false.
        end if
    end do

    if (SILENT) close(20)

    call set_weight_gradient()

    call set_ionic_ignored()

    call fit_emf(20, temp)

    do i=12,14
        if (abs(lbeta(i)-real_logbeta(i)) > TOLERANCE) then
            if (passed) call print_result(test_name, "failed")
            passed = .false.
            print "(A,I3,A, 2F10.3)", "  log beta(",i,") too large: ", lbeta(i), real_logbeta(i)
        end if
    end do

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_znedta1

end program test_libemf

