program test_libemf
    implicit none

    print "(A)", "+-------------------------+"
    print "(A)", "|  Testing libemf module  |"
    print "(A)", "+-------------------------+"
    print *, ""

    !call test_init()
    call test_gradient()
    call test_hexaprotic()
contains

!subroutine test_init()
!    ! Test initiation routines
!    use emf
!    implicit none
!    double precision, dimension(6) :: logbeta
!    integer, dimension(6) :: flgbeta
!    integer, dimension(6,3) :: P
!    real, dimension(5,2) :: data
!    real, dimension(3,2) :: tparms
!    integer, dimension(3,2) :: tflags
!
!    character(len=*), parameter :: test_name = "initiation"
!    logical passed
!    passed = .true.
!
!    call init_module(10, 2, 0, P, logbeta, flgbeta)
!    if (emf_diagnostic .ne. 0) then
!        if (passed) call print_result(test_name, "failed")
!        passed = .false.
!        print *, "  init_module returned ", emf_diagnostic
!    end if
!
!    call load_titration(data(:,1), data(:,2), tparms, 1.0, 300.0, 1.00, tflags, 0, 0)
!    if (emf_diagnostic .ne. 0) then
!        if (passed) call print_result(test_name, "failed")
!        passed = .false.
!        print *, "  load_titration returned ", emf_diagnostic
!    end if
!
!    call load_titration(data(:,1), data(:,2), tparms, 1.0, 300.0, 1.00, tflags, 0, 0)
!    if (emf_diagnostic .ne. 0) then
!        if (passed) call print_result(test_name, "failed")
!        passed = .false.
!        print *, "  load_titration returned ", emf_diagnostic
!    end if
!
!    call finalize()
!    if (passed) call print_result(test_name, "passed")
!end subroutine test_init


subroutine test_hexaprotic()
    use emf
    use iso_fortran_env, only: output_unit

    implicit none
    character(len=*), parameter :: test_name = "hexaprotic"
    logical passed
    integer :: unit, ios
    integer :: maxits, rNEQ, rNCP, aux
    integer i
    character(len=80) :: line
    real, allocatable, dimension(:) :: rlogb
    integer, allocatable, DIMENSION(:,:) :: rP
    integer, allocatable, DIMENSION(:) :: rlbflags
    real :: r_temp

    passed = .true.

    ! ---- open file ----
    open(newunit=unit, file='../../tests/simu1.sup', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file'
      stop
    end if

    ! header
    read(unit, '(A)', iostat=ios) line
    ! control numbers
    read(unit, '(4I5)', iostat=ios) maxits, rNEQ, rNCP, aux

    ! components names
    do i=1,rNCP
        read(unit, '(A)', iostat=ios) line
    end do
    ! temperature
    read(unit, '(F8.2)', iostat=ios) r_temp

    rNEQ = 7
    allocate(rlogb(rNEQ))
    allocate(rlbflags(rNEQ))
    allocate(rP(rNEQ, rNCP))

    do i=1,rNEQ
        read(unit, '(F9.4,3I5)', iostat=ios) rlogb(i), rP(i,:), rlbflags(i)
    end do

    read(unit, '(A)', iostat=ios) line    ! blank

    ! Open unit 99 to /dev/null (Unix/Linux) or 'NUL' (Windows)
    ! open(unit=99, file='/dev/null', status='unknown')
    ! call init_module(513, 6, rP, rlogb, rlbflags, output_unit)

    rlogb(1) = rlogb(1) + 0.5
    rlogb(2) = rlogb(2) - 0.5
    rlogb(3) = rlogb(3) + 0.5
    rlogb(4) = rlogb(4) - 0.5
    rlogb(5) = rlogb(5) + 0.5
    rlogb(6) = rlogb(6) - 0.5

    call init_module(513, 6, rP, rlogb, rlbflags, output_unit)

    call aux_read_suptitr(89, unit)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 1 -> " // emf_diagnostic
        passed=.false.
    end if
    call aux_read_suptitr(88, unit)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 2 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if
    call aux_read_suptitr(74, unit)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 3 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if
    call aux_read_suptitr(78, unit)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 4 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if
    call aux_read_suptitr(92, unit)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 5 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if
    call aux_read_suptitr(92, unit)
    if (trim(emf_diagnostic) /= 'no error') then
        if (passed) call print_result(test_name, "failed")
        print *, "  loading titration 6 -> " // emf_diagnostic
        passed=.false.
        print *, emf_diagnostic
    end if

    close(unit, iostat=ios)
    !close(99)

    call set_weight_gradient()
    call fit_emf(r_temp)

    call finalize()

    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic


subroutine aux_read_suptitr(n, unit)
    use emf
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: unit

    character(len=80) :: line
    integer ios
    integer dummy
    real, dimension(2,2) :: tparams
    integer, dimension(2,2) :: tflags
    real, dimension(n) :: titre
    real, dimension(n) :: emfr
    real :: initv, errv
    real :: emf0, erremf
    integer :: eact, e0flag
    integer i

    do i=1,NCP
        read(unit, '(I3,I2,2F9.5,I6,I5)', iostat=ios) dummy, dummy, tparams(i,:), tflags(i,:)
    end do

    read(unit, '(A)', iostat=ios) line    ! blank

    read(unit, '(F13.5,F10.5)', iostat=ios) initv, errv
    read(unit, '(I3,I2,2F9.5,2I6)', iostat=ios) dummy, eact, emf0, erremf, e0flag

    read(unit, '(A)', iostat=ios) line    ! blank
    do i=1,n
        read(unit, *, iostat=ios) titre(i), emfr(i)
    end do

    call load_titration(titre, emfr, tparams, eact, initv, emf0, 1.00, tflags, e0flag, 0, &
        erremf, errv)

    read(unit, '(A)', iostat=ios) line    ! blank
end subroutine aux_read_suptitr


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

end program test_libemf

