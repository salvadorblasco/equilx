program test_libemf
    implicit none

    print "(A)", "+-------------------------+"
    print "(A)", "|  Testing libemf module  |"
    print "(A)", "+-------------------------+"
    print *, ""

    !call test_init()
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
    implicit none
    character(len=*), parameter :: test_name = "load hexaprotic"
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
    read(unit, '(A)', iostat=ios) r_temp

    rNEQ = 7
    allocate(rlogb(rNEQ))
    allocate(rlbflags(rNEQ))
    allocate(rP(rNEQ, rNCP))

    do i=1,rNEQ
        read(unit, '(F9.4,3I5)', iostat=ios) rlogb(i), rP(i,:), rlbflags(i)
    end do

    read(unit, '(A)', iostat=ios) line    ! blank

    call init_module(513, 6, rP, rlogb, rlbflags)

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

    !call finalize()

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
        read(unit, '(F6.2,F12.2)', iostat=ios) titre(i), emfr(i)
    end do

    call load_titration(titre, emfr, tparams, eact, initv, emf0, 1.00, tflags, e0flag, 0, &
        erremf, errv)

    read(unit, '(A)', iostat=ios) line    ! blank
end subroutine aux_read_suptitr


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
