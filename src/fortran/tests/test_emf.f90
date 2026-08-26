program test_libemf
    implicit none

    print "(A)", "+-------------------------+"
    print "(A)", "|  Testing libemf module  |"
    print "(A)", "+-------------------------+"
    print *, ""

    call test_init()
    call test_hexaprotic()
contains

subroutine test_init()
    ! Test initiation routines
    use emf
    implicit none
    double precision, dimension(6) :: logbeta
    integer, dimension(6) :: flgbeta
    integer, dimension(6,3) :: P
    real, dimension(5,2) :: data
    real, dimension(3,2) :: tparms
    integer, dimension(3,2) :: tflags

    character(len=*), parameter :: test_name = "initiation"
    logical passed
    passed = .true.

    call init_module(10, 2, 0, P, logbeta, flgbeta)
    if (diagnostic .ne. 0) then
        if (passed) call print_result(test_name, "failed")
        passed = .false.
        print *, "  init_module returned ", diagnostic
    end if

    call load_titration(data(:,1), data(:,2), tparms, 1.0, 300.0, 1.00, tflags, 0, 0)
    if (diagnostic .ne. 0) then
        if (passed) call print_result(test_name, "failed")
        passed = .false.
        print *, "  load_titration returned ", diagnostic
    end if

    call load_titration(data(:,1), data(:,2), tparms, 1.0, 300.0, 1.00, tflags, 0, 0)
    if (diagnostic .ne. 0) then
        if (passed) call print_result(test_name, "failed")
        passed = .false.
        print *, "  load_titration returned ", diagnostic
    end if

    call finalize()
    if (passed) call print_result(test_name, "passed")
end subroutine test_init


subroutine test_hexaprotic()
    use emf
    implicit none
    integer, dimension(7, 2) :: P = reshape((/ 1,1,1,1,1,1,0, 1,2,3,4,5,6,-1 /), shape=(/ 7, 2 /))
    DOUBLE PRECISION, DIMENSION(7) :: logB = (/ 10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77 /) 
    real, dimension(89) :: titre = (/ 0.00, 0.12, 0.24, 0.36, 0.48, 0.60, 0.72, 0.84, 0.96, 1.08, 1.20, 1.32, &
    1.44 , 1.56 , 1.68 , 1.80 , 1.92 , 2.04 , 2.16 , 2.28 , 2.40 , 2.52 , 2.64, 2.76, 2.88, 3.00, 3.12, 3.24, 3.36, 3.48, &
    3.60 , 3.72 , 3.84 , 3.96 , 4.08 , 4.20 , 4.32 , 4.44 , 4.56 , 4.68 , 4.80, 4.92, 5.04, 5.16, 5.28, 5.40, 5.52, 5.64, &
    5.76 , 5.88 , 6.00 , 6.12 , 6.24 , 6.36 , 6.48 , 6.60 , 6.72 , 6.84 , 6.96, 7.08, 7.20, 7.32, 7.44, 7.56, 7.68, 7.80, &
    7.92 , 8.04 , 8.16 , 8.28 , 8.40 , 8.52 , 8.64 , 8.76 , 8.88 , 9.00 , 9.12, 9.24, 9.36, 9.48, 9.72,10.08,10.32, 10.56,&
    10.80, 10.92, 11.28, 11.64, 12.00 /)
    real, dimension(89) :: emf = (/ 265.56, 263.49, 261.36, 259.05, 256.69, 254.20, 251.60, 248.82, 245.92, &
        242.90, 239.71, 236.33, 232.79, 229.06, 225.21, 221.13, 216.87, 212.43, 207.82, 202.97, 197.94, 192.67, 187.23, &
        181.49, 175.52, 169.19, 162.44, 155.05, 146.71, 137.12, 125.65, 112.22, 98.02, 84.94, 73.64, 63.59, 54.18, 45.07, &
        35.66, 25.49, 13.95, 0.58, -14.03, -27.93, -40.06, -50.59, -60.18, -69.40, -78.63, -88.39, -99.16, -111.64, -125.61, &
        -139.75, -152.29, -163.00, -172.46, -181.10, -189.32, -197.25, -205.12, -212.87, -220.44, -227.54, -233.99, -239.73, &
        -244.75, -249.19, -253.04, -256.47, -259.54, -262.27, -264.75, -267.06, -269.13, -271.02, -272.86, -274.51, -276.05, &
        -277.53, -280.19, -283.68, -285.75, -287.65, -289.42, -290.19, -292.50, -294.57, -296.46 /)
    integer, dimension(2,2) :: tflags = reshape([2,0,0,0], [2,2])
    real, dimension(2,2) :: tparms = reshape([ 1.125, 0.75, 0.0, -0.1 ], [2,2])
    integer, dimension(7) :: flgbeta
    character(len=*), parameter :: test_name = "load hexaprotic"
    logical passed
    passed = .true.

    flgbeta = 0
    call init_module(89, 1, 1, P, logb, flgbeta)
    if (diagnostic .ne. 0) then
        if (passed) call print_result(test_name, "failed")
        passed = .false.
        print *, "  init_module returned ", diagnostic
    end if

    call load_titration(titre, emf, tparms, 25.0, 300.0, 1.00, tflags, 0, 0)
    if (diagnostic .ne. 0) then
        if (passed) call print_result(test_name, "failed")
        passed = .false.
        print *, "  load_titration returned ", diagnostic
    end if

    call finalize()
    if (passed) call print_result(test_name, "passed")
end subroutine test_hexaprotic


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
