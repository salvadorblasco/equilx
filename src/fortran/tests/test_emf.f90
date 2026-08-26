program test_libemf
    use emf

    double precision, dimension(6) :: logbeta
    integer, dimension(6) :: flag_beta
    integer, dimension(6,3) :: P

    print *, "Testing libemf module"

    call init_module(10, 2, 0, P, logbeta, flgbeta)

    print *, "> ", NPT, NTT

end program test_libemf
