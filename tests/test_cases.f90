! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>
!
! This file contains data for the unit tests

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
    double precision, dimension(NPT, NCP+NSL, NCP+NSL) :: tJ = reshape( &
      (/ 1.        , 0.        , 0., &
         0.        , 0.        , 0., &
         0.        , 0.        , 0., &
         0.7       , 0.11118056, 0., &
         0.11118056, 0.3       , 0., &
        -0.92650465,-0.92650465, 0., &
         0.46641016, 0.12      , 0.03358983851, &
         0.12      , 0.46641016, 0.03358983851, &
        -1.0       ,-1.0      ,  0., &
         0.3       , 0.11118056, 0., &
         0.11118056, 0.7       , 0., &
        -0.92650465,-0.92650465, 0., &
         0.        , 0.        , 0., &
         0.        , 1.        , 0., &
         0.        , 0.        , 0. /), (/ NPT, NCP+NSL, NCP+NSL /), order=[3,2,1])
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
    DOUBLE PRECISION, DIMENSION(NEQ) :: logB = [10.0, 18.0, 24.0, 28.0, 31.0, 33.0, -13.77]
    integer, dimension(NEQ, NCP) :: P = reshape([1,1,1,1,1,1,0, 1,2,3,4,5,6,-1], &
                                                shape=[NEQ, NCP])
contains

subroutine initialize()
    integer :: i
    do i=1,NPT
        T(i,1) = 5.0E-3*V0 / (V(i) + V0)
        T(i,2) = (30.0E-3*V0 - 0.1* V(i))/ (V(i) + V0)
    end do
    B = 10.0**logB
end subroutine initialize
end module hexaprotic


module molybdate
    implicit none
    integer, PARAMETER :: NPT = 50
    integer, PARAMETER :: NCP = 2
    integer, PARAMETER :: NEQ = 7
    double precision :: V0 = 25.0
    double precision :: Iref = 0.1
    double precision, dimension(NPT) :: V
    double precision, dimension(NPT, NCP) :: T
    DOUBLE PRECISION, DIMENSION(NEQ) :: B 
    DOUBLE PRECISION, DIMENSION(NEQ) :: logB = (/ -14.0, 4.47, 8.12, 51.93, 58.90, 64.45, 68.67 /) 
    integer, dimension(NEQ, NCP) :: P = reshape((/ 0,1,1,7,7,7,7,-1,1,2,8,9,10,11 /), &
                                                shape=(/ NEQ, NCP /))
    double PRECISION :: DHA = 0.51
    double PRECISION :: DHB = 1.5
    double precision, dimension(NEQ,3) :: DHP = reshape((/ -0.14, 0.57, 0.61, 1.38, 1.82, 2.74, 3.44, &
                                                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                                  -0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /), &
                                              shape=(/ NEQ,3 /))
    integer, DIMENSION(NCP+NEQ) :: charge
    double precision, dimension(NPT) :: Iback

contains

subroutine initialize()
    use mouter, only: full_charge
    integer :: i
    double precision, parameter :: initial_molybd = 0.001009
    double precision, parameter :: initial_proton = 0.006036
    double precision, parameter :: buret_proton = -0.09985
    double precision, parameter :: background_Ia = 0.20890
    double precision, parameter :: background_Ib = 0.09985
    do i=1,NPT
        V(i) = 2.0 * (i-1)/(NPT-1)
        T(i,1) =  initial_molybd*V0 / (V(i) + V0)
        T(i,2) = (initial_proton*V0 - buret_proton* V(i))/ (V(i) + V0)
        Iback(i) = background_Ia + background_Ib*V(i)
    end do
    B = 10.0**logB
    charge(1) = -2
    charge(2) =  1
    call full_charge(charge, P, NCP, NEQ)
end subroutine initialize

end module molybdate


module phosphoric
    implicit none
    integer, PARAMETER :: NPT = 50
    integer, PARAMETER :: NCP = 2
    integer, PARAMETER :: NEQ = 4
    double precision, dimension(NPT, NCP) :: T
    DOUBLE PRECISION, DIMENSION(NEQ) :: B 
    DOUBLE PRECISION, DIMENSION(NEQ) :: logB = (/ 12.38, 19.578, 21.726, -13.73 /) 
    integer, dimension(NEQ, NCP) :: P = reshape((/ 1,1,1,0, 1,2,3,-1 /), &
                                                shape=(/ NEQ, NCP /))
    double precision, dimension(NPT, NCP+NEQ) :: C

contains

subroutine initialize()
    integer i
    double precision, parameter :: initial_phos = 1e-8
    double precision, parameter :: initial_pH = 1.0
    double precision, parameter :: final_pH= 14.0
    C = 0.0
    do i=1,NPT
        T(i,1) =  initial_phos
        T(i,2) = 0.0
        C(i,2) = 10.0**(-(i-1)*(final_pH-initial_pH)/(NPT-1) - initial_pH)
    end do
    B = 10.0**logB
end subroutine initialize

end module phosphoric

module aluminium
    implicit none
    integer, PARAMETER :: NPT = 51
    integer, PARAMETER :: NCP = 2
    integer, PARAMETER :: NSL = 1
    integer, PARAMETER :: NEQ = 5
    double precision :: V0 = 1.0
    double precision, dimension(NPT) :: V
    double precision, dimension(NPT, NCP) :: T
    DOUBLE PRECISION, DIMENSION(NEQ) :: B 
    DOUBLE PRECISION, DIMENSION(NEQ) :: logB = (/ -5.0, -9.3, -15.0, -23.0, -14.0 /) 
    integer, dimension(NEQ, NCP) :: P = reshape((/ 1,1,1,1,0, -1,-2,-3,-4,-1 /), &
                                                shape=(/ NEQ, NCP /))
    integer, dimension(NSL, NCP) :: Q = reshape((/1,-3/), (/NSL, NCP/))
    DOUBLE PRECISION, DIMENSION(NSL) :: KS = (/ 10.0**5.5 /)
    double precision, dimension(NPT,NCP+NEQ):: tC = reshape( &
        (/ 9.989960e-09, 1.000000e-02, 9.989960e-12, 5.006840e-14, 9.989960e-18, 9.989960e-24, 1.000000e-12,&
           9.984051e-09, 6.309573e-03, 1.582365e-11, 1.256918e-13, 3.974722e-17, 6.299510e-23, 1.584893e-12,&
           9.974629e-09, 3.981072e-03, 2.505514e-11, 3.154255e-13, 1.580872e-16, 3.970971e-22, 2.511886e-12,&
           9.959559e-09, 2.511886e-03, 3.964972e-11, 7.911159e-13, 6.284057e-16, 2.501728e-21, 3.981072e-12,&
           9.935327e-09, 1.584893e-03, 6.268768e-11, 1.982358e-12, 2.495641e-15, 1.574643e-20, 6.309573e-12,&
           9.896070e-09, 1.000000e-03, 9.896070e-11, 4.959784e-12, 9.896070e-15, 9.896070e-20, 1.000000e-11,&
           9.831761e-09, 6.309573e-04, 1.558229e-10, 1.237745e-11, 3.914094e-14, 6.203421e-19, 1.584893e-11,&
           9.724817e-09, 3.981072e-04, 2.442764e-10, 3.075257e-11, 1.541280e-13, 3.871519e-18, 2.511886e-11,&
           9.543650e-09, 2.511886e-04, 3.799396e-10, 7.580791e-11, 6.021636e-13, 2.397257e-17, 3.981072e-11,&
           9.231057e-09, 1.584893e-04, 5.824403e-10, 1.841838e-10, 2.318737e-12, 1.463024e-16, 6.309573e-11,&
           8.687200e-09, 1.000000e-04, 8.687200e-10, 4.353914e-10, 8.687200e-12, 8.687200e-16, 1.000000e-10,&
           7.761784e-09, 6.309573e-05, 1.230160e-09, 9.771507e-10, 3.090022e-11, 4.897355e-15, 1.584893e-10,&
           6.316045e-09, 3.981072e-05, 1.586519e-09, 1.997309e-09, 1.001026e-10, 2.514463e-14, 2.511886e-10,&
           4.433496e-09, 2.511886e-05, 1.765007e-09, 3.521651e-09, 2.797347e-10, 1.113644e-13, 3.981072e-10,&
           1.258925e-09, 1.584893e-05, 7.943282e-10, 2.511886e-09, 3.162278e-10, 1.995262e-13, 6.309573e-10,&
           3.162278e-10, 1.000000e-05, 3.162278e-10, 1.584893e-09, 3.162278e-10, 3.162278e-13, 1.000000e-09,&
           7.943282e-11, 6.309573e-06, 1.258925e-10, 1.000000e-09, 3.162278e-10, 5.011872e-13, 1.584893e-09,&
           1.995262e-11, 3.981072e-06, 5.011872e-11, 6.309573e-10, 3.162278e-10, 7.943282e-13, 2.511886e-09,&
           5.011872e-12, 2.511886e-06, 1.995262e-11, 3.981072e-10, 3.162278e-10, 1.258925e-12, 3.981072e-09,&
           1.258925e-12, 1.584893e-06, 7.943282e-12, 2.511886e-10, 3.162278e-10, 1.995262e-12, 6.309573e-09,&
           3.162278e-13, 1.000000e-06, 3.162278e-12, 1.584893e-10, 3.162278e-10, 3.162278e-12, 1.000000e-08,&
           7.943282e-14, 6.309573e-07, 1.258925e-12, 1.000000e-10, 3.162278e-10, 5.011872e-12, 1.584893e-08,&
           1.995262e-14, 3.981072e-07, 5.011872e-13, 6.309573e-11, 3.162278e-10, 7.943282e-12, 2.511886e-08,&
           5.011872e-15, 2.511886e-07, 1.995262e-13, 3.981072e-11, 3.162278e-10, 1.258925e-11, 3.981072e-08,&
           1.258925e-15, 1.584893e-07, 7.943282e-14, 2.511886e-11, 3.162278e-10, 1.995262e-11, 6.309573e-08,&
           3.162278e-16, 1.000000e-07, 3.162278e-14, 1.584893e-11, 3.162278e-10, 3.162278e-11, 1.000000e-07,&
           7.943282e-17, 6.309573e-08, 1.258925e-14, 1.000000e-11, 3.162278e-10, 5.011872e-11, 1.584893e-07,&
           1.995262e-17, 3.981072e-08, 5.011872e-15, 6.309573e-12, 3.162278e-10, 7.943282e-11, 2.511886e-07,&
           5.011872e-18, 2.511886e-08, 1.995262e-15, 3.981072e-12, 3.162278e-10, 1.258925e-10, 3.981072e-07,&
           1.258925e-18, 1.584893e-08, 7.943282e-16, 2.511886e-12, 3.162278e-10, 1.995262e-10, 6.309573e-07,&
           3.162278e-19, 1.000000e-08, 3.162278e-16, 1.584893e-12, 3.162278e-10, 3.162278e-10, 1.000000e-06,&
           7.943282e-20, 6.309573e-09, 1.258925e-16, 1.000000e-12, 3.162278e-10, 5.011872e-10, 1.584893e-06,&
           1.995262e-20, 3.981072e-09, 5.011872e-17, 6.309573e-13, 3.162278e-10, 7.943282e-10, 2.511886e-06,&
           5.011872e-21, 2.511886e-09, 1.995262e-17, 3.981072e-13, 3.162278e-10, 1.258925e-09, 3.981072e-06,&
           1.258925e-21, 1.584893e-09, 7.943282e-18, 2.511886e-13, 3.162278e-10, 1.995262e-09, 6.309573e-06,&
           3.162278e-22, 1.000000e-09, 3.162278e-18, 1.584893e-13, 3.162278e-10, 3.162278e-09, 1.000000e-05,&
           7.943282e-23, 6.309573e-10, 1.258925e-18, 1.000000e-13, 3.162278e-10, 5.011872e-09, 1.584893e-05,&
           1.995262e-23, 3.981072e-10, 5.011872e-19, 6.309573e-14, 3.162278e-10, 7.943282e-09, 2.511886e-05,&
           3.883510e-24, 2.511886e-10, 1.546053e-19, 3.084782e-14, 2.450329e-10, 9.754936e-09, 3.981072e-05,&
           6.211126e-25, 1.584893e-10, 3.918956e-20, 1.239283e-14, 1.560164e-10, 9.843971e-09, 6.309573e-05,&
           9.900985e-26, 1.000000e-10, 9.900985e-21, 4.962247e-15, 9.900985e-11, 9.900985e-09, 1.000000e-04,&
           1.574956e-26, 6.309573e-11, 2.496136e-21, 1.982752e-15, 6.270011e-11, 9.937298e-09, 1.584893e-04,&
           2.501926e-27, 3.981072e-11, 6.284554e-22, 7.911784e-16, 3.965285e-11, 9.960346e-09, 2.511886e-04,&
           3.971097e-28, 2.511886e-11, 1.580922e-22, 3.154354e-16, 2.505593e-11, 9.974944e-09, 3.981072e-04,&
           6.299589e-29, 1.584893e-11, 3.974772e-23, 1.256933e-16, 1.582385e-11, 9.984176e-09, 6.309573e-04,&
           9.990010e-30, 1.000000e-11, 9.990010e-24, 5.006865e-17, 9.990010e-12, 9.990010e-09, 1.000000e-03,&
           1.583894e-30, 6.309573e-12, 2.510303e-24, 1.994004e-17, 6.305595e-12, 9.993694e-09, 1.584893e-03,&
           2.510887e-31, 3.981072e-12, 6.307063e-25, 7.940121e-18, 3.979487e-12, 9.996021e-09, 2.511886e-03,&
           3.980072e-32, 2.511886e-12, 1.584495e-25, 3.161484e-18, 2.511256e-12, 9.997489e-09, 3.981072e-03,&
           6.308574e-33, 1.584893e-12, 3.980441e-26, 1.258726e-18, 1.584642e-12, 9.998415e-09, 6.309573e-03,&
           9.999000e-34, 1.000000e-12, 9.999000e-27, 5.011371e-19, 9.999000e-13, 9.999000e-09, 1.000000e-02 /), &
           shape=(/ NPT, NCP+NEQ /), order=[2, 1])
       double precision, dimension(NPT, NSL) :: tSI = reshape( (/ 3.159103e-08, 1.256918e-07, 4.999157e-07, &
           1.987193e-06, 7.891911e-06, 3.129412e-05, 1.237745e-04, 4.873954e-04, 1.904209e-03, 7.332489e-03, &
           2.747134e-02, 9.771507e-02, 3.165521e-01, 8.845988e-01, 1.000000e+00, 1.000000e+00, 1.000000e+00, &
           1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, &
           1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, &
           1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, &
           7.748621e-01, 4.933673e-01, 3.130966e-01, 1.982752e-01, 1.253933e-01, 7.923379e-02, 5.003942e-02, &
           3.159119e-02, 1.994004e-02, 1.258424e-02, 7.941288e-03, 5.011078e-03, 3.161961e-03 /), shape=(/NPT,NSL/))
        double precision, dimension(NPT, NSL) :: tCS = reshape( (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
           0.0, 0.0, 0.0, 0.0, 0.0, 5.118433e-09, 7.466107e-09, 8.477946e-09, 8.981949e-09, 9.259442e-09, 9.421386e-09, &
           9.518642e-09, 9.577422e-09, 9.612212e-09, 9.631168e-09, 9.638620e-09, 9.636269e-09, 9.623641e-09, 9.598025e-09, &
           9.553897e-09, 9.481733e-09, 9.365959e-09, 9.181585e-09, 8.888813e-09, 8.424449e-09, 7.688259e-09, 6.521336e-09, &
           4.671800e-09, 1.740427e-09, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /), shape=(/NPT,NSL/))

contains

subroutine initialize()
    integer i
    double precision, parameter :: initial_alum = 1e-8
    double precision, parameter :: initial_pH = 2.0
    double precision, parameter :: final_pH= 12.0
    do i=1,NPT
        T(i,1) =  initial_alum
        T(i,2) = 10.0**(-(i-1)*(final_pH-initial_pH)/(NPT-1) - initial_pH)
    end do
    B = 10.0**logB
end subroutine initialize
end module aluminium


module large_fixture
    ! labels = "L Cu Im H"
    implicit none
    integer, PARAMETER :: NPT = 20
    integer, PARAMETER :: NCP = 4
    integer, PARAMETER :: NEQ = 25
    double precision, dimension(NPT, NCP) :: T
    DOUBLE PRECISION, DIMENSION(NEQ) :: B 
    DOUBLE PRECISION, DIMENSION(NEQ) :: logB = (/ 10.66, 20.07, 28.31, 35.67, 42.57, 48.44,  7.09,  4.20,  7.73, 10.60, &
                                                  13.00, 44.36, 40.54, 34.41, 26.97, 17.21, 25.18, 17.97, 50.43, 45.82, &
                                                  38.97, 36.30, 29.79, 22.28, -13.77 /)
    integer, dimension(NEQ, NCP) :: P = reshape([ 1, 0, 0,  1, &
                                                  1, 0, 0,  2, &
                                                  1, 0, 0,  3, &
                                                  1, 0, 0,  4, &
                                                  1, 0, 0,  5, &
                                                  1, 0, 0,  6, &
                                                  0, 0, 1,  1, &
                                                  0, 1, 1,  0, &
                                                  0, 1, 2,  0, &
                                                  0, 1, 3,  0, &
                                                  0, 1, 4,  0, &
                                                  1, 1, 0,  4, &
                                                  1, 1, 0,  3, &
                                                  1, 1, 0,  2, &
                                                  1, 1, 0,  1, &
                                                  1, 1, 0,  0, &
                                                  1, 2, 0,  0, &
                                                  1, 2, 0, -1, &
                                                  1, 1, 1,  4, &
                                                  1, 1, 1,  3, &
                                                  1, 1, 1,  2, &
                                                  1, 1, 1,  1, &
                                                  1, 2, 1,  0, &
                                                  1, 2, 1, -1, &
                                                  0, 0, 0, -1], shape=(/ NEQ, NCP /), order=[2, 1])

contains
subroutine initialize()
    integer i
    double precision, parameter :: initial_M = 0.001
    double precision, parameter :: initial_H = 0.005
    double precision, parameter :: final_H= -0.005
    T(:,1:3) = initial_M
    do i=1,NPT
        T(i,4) = (i-1)*(final_H-initial_H)/(NPT-1) + initial_H
    end do
    B = 10.0**logB
end subroutine initialize
end module large_fixture


module gallic
    ! case for aluminium with gallic acid. Values taken from Carrayrou (2004)
    implicit none
    integer, PARAMETER :: NCP = 3
    integer, PARAMETER :: NEQ = 14
    double precision, dimension(NCP) :: T = [1.58e-6, 1e-3, 1e-3]
    DOUBLE PRECISION, DIMENSION(NEQ) :: B 
    DOUBLE PRECISION, DIMENSION(NEQ) :: logB = (/ -14.0, -4.15, -12.59, -23.67, &
                                                  -4.93, -9.43, -21.98, -37.69, &
                                                  -22.65, -27.81, -32.87, -39.56, &
                                                  -20.25, -12.52 /) 
    integer, dimension(NEQ, NCP) :: P = reshape((/ -1, 0, 0, &
                                                   -1, 0, 1, &
                                                   -2, 0, 1, &
                                                   -3, 0, 1, &
                                                   -2, 1, 1, &
                                                   -3, 1, 1, &
                                                   -6, 1, 2, &
                                                   -9, 1, 3, &
                                                   -8, 2, 3, &
                                                   -9, 2, 3, &
                                                  -10, 2, 3, &
                                                  -11, 2, 3, &
                                                   -9, 4, 3, &
                                                   -5, 3, 1  /), &
                                                shape=(/ NEQ, NCP /), order=[2,1])
    double precision, dimension(NCP+NEQ) :: tC

contains

subroutine initialize()
    B = 10.0**logB
    tC(1) = 1.58e-6
    tC(2) = 2.03e-5
    tC(3) = 2.59e-7
end subroutine initialize

end module gallic


module cobalt_edta
    implicit none
    integer, PARAMETER :: NPT = 25
    integer, PARAMETER :: NCP = 3
    integer, PARAMETER :: NEQ = 10
    integer, PARAMETER :: NSL = 1
    double precision, parameter :: V0 = 10.0
    double precision, parameter :: final_titre = 1.25
    double precision, dimension(NPT) :: V
    double precision, dimension(NPT, NCP) :: T
    DOUBLE PRECISION, DIMENSION(NEQ) :: B 
    DOUBLE PRECISION, DIMENSION(NEQ) :: logB = [ -14.0, -12.2, -19.5, -31.7, 15.9, 3.1, 11.0, 17.3, 19.8, 21.7]
    integer, dimension(NEQ, NCP) :: P = reshape([ 0, 0,-1, & 
                                                  1, 0,-1, &
                                                  1, 0,-2, &
                                                  1, 0,-3, &
                                                  1, 1, 0, &
                                                  1, 1, 1, &
                                                  0, 1, 1, &
                                                  0, 1, 2, &
                                                  0, 1, 3, &
                                                  0, 1, 4 ], shape=[NEQ, NCP], order=[2,1])
    integer, dimension(NSL, NCP) :: Q = reshape([1,0,-2], [NSL, NCP])
    DOUBLE PRECISION, DIMENSION(NSL) :: KS = (/ 10.0**12.65 /)
contains
subroutine initialize()
    integer i
    double precision, parameter :: initial_co = 1.3e-3
    double precision, parameter :: initial_edta = 1.0e-3
    double precision, parameter :: initial_H= 4.0e-3
    double precision, parameter :: buret_H= -0.1
    do i=1,NPT
        V(i) = final_titre*(i-1)/(NPT-1)
        T(i,1) = initial_co/(V(i)+V0)
        T(i,3) = (initial_H+buret_H*V(i))/(V(i)+V0)
        T(i,2) = initial_edta/(V(i)+V0)
    end do
    B = 10.0**logB
end subroutine initialize
end module cobalt_edta


module calcium_carbonate
    implicit none
    integer, PARAMETER :: NPT = 20
    integer, PARAMETER :: NCP = 3
    integer, PARAMETER :: NEQ = 4
    integer, PARAMETER :: NSL = 2
    double precision, dimension(NPT) :: V
    double precision, dimension(NPT, NCP) :: T
    DOUBLE PRECISION, DIMENSION(NEQ) :: B 
    DOUBLE PRECISION, DIMENSION(NEQ) :: logB = [16.33, 22.68, -12.06, -14.0]
    integer, dimension(NEQ, NCP) :: P = reshape([ 1, 0, 1, & 
                                                  1, 0, 2, &
                                                  1, 1, 0, &
                                                  0, 0,-1 ], shape=[NEQ, NCP], order=[2,1])
    integer, dimension(NSL, NCP) :: Q = reshape([0, 1,-2, &
                                                 1, 1, 0], [NSL, NCP], order=[2, 1])
    DOUBLE PRECISION, DIMENSION(NSL) :: KS = [ 10.0**23.0, 10**(-8.4) ]
    ! charges=np.array([-3, 2, 1], dtype=int),
    ! component_labels=['CO3', 'Ca', 'H'])
contains
subroutine initialize()
    integer i
    double precision, parameter :: V0 = 10.0
    double precision, parameter :: final_titre = 3.8
    double precision, parameter :: initial_carb = 0.05
    double precision, parameter :: initial_ca = 0.1
    double precision, parameter :: initial_H= 0.2
    double precision, parameter :: buret_H= -0.1
    do i=1,NPT
        V(i) = final_titre*(i-1)/(NPT-1)
        T(i,1) = initial_carb/(V(i)+V0)
        T(i,2) = initial_ca/(V(i)+V0)
        T(i,3) = (initial_H+buret_H*V(i))/(V(i)+V0)
    end do
    B = 10.0**logB
end subroutine initialize
end module calcium_carbonate
