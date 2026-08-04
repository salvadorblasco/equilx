! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>


subroutine lsolve(A, b, x, n, CTRL)
  ! Solve system of lineal equations using LAPACK. This routine requires
  ! the -llapack flag at compilation and, probably also -lblas in some systems.
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: A   ! The nxn array of coefficients
  double precision, dimension(n),   intent(in) :: b   ! the n array of independent terms
  double precision, dimension(n),   intent(out) :: x  ! the n array where the result is returned
  integer, intent(inout) :: CTRL                      ! control value. 0 on no error.

  integer, dimension(n) :: ipiv   ! pivot indices
  integer :: info, nrhs, lda, ldb

  ! External procedures defined in LAPACK
  external dgesv

  ! Store b in x to prevent it from being overwritten by LAPACK
  x = b

  ! Set up solve by LAPACK routine DGESV
  nrhs = 1 ! number of right hand sides in b
  lda = n  ! leading dimension of a
  ldb = n  ! leading dimension of b

  call dgesv(n, nrhs, A, lda, ipiv, x, ldb, info)

  if (info /= 0) then
      CTRL = 99    ! unable to invert matrix
      return
  end if
end subroutine lsolve


subroutine jacobian_soluble(J, C, P, NCP, NEQ)
  ! calculate the jacobian df/dlogc of the mass balance not considering solids
  implicit none
  double precision, INTENT(INOUT), dimension(NCP, NCP) ::  J  ! jacobian
  double precision, INTENT(IN),    dimension(NCP+NEQ)  :: C   ! free concentrations
  integer,          INTENT(IN),    dimension(NEQ, NCP) :: P   ! The soluble stoichiometry array
  integer, intent(in) :: NCP                                  ! number of components
  integer, intent(in) :: NEQ                                  ! number of equilibria
  J(:,:) = 0.0
  call amatrix(J, C, P, NCP, NEQ)
end subroutine jacobian_soluble


subroutine jacobian_solid(JC, C, KS, P, Q, NCP, NEQ, NSL, do_solids, NX)
  implicit none
  integer, intent(in) :: NCP                                       ! number of components
  integer, intent(in) :: NEQ                                       ! number of equilibria
  integer, intent(in) :: NSL                                       ! number of solids
  integer,          intent(in),  dimension(NEQ, NCP)         :: P  ! The soluble stoichiometry array
  integer,          intent(in),  dimension(NSL, NCP)         :: Q  ! The solids stoichiometry array
  double precision, intent(in),  dimension(NEQ+NCP)          :: C  ! free concentation initial estimates (in) and refined values (out)
  double precision, intent(in),  dimension(NSL)              :: KS ! the solubility product
  double precision, INTENT(OUT), dimension(NCP+NSL, NCP+NSL) :: JC  ! jacobian
  logical,          intent(in),  dimension(NSL)              :: do_solids
  integer,          intent(out)                              :: NX

  integer i, j, k
  double PRECISION aux
  JC(:,:) = 0.0
  call jacobian_soluble(JC(1:NCP,1:NCP), C, P, NCP, NEQ)
  NX = NCP
  do i=1,NSL
    if (.not. do_solids(i)) cycle
    NX = NX + 1
    do j=1,NCP
      JC(j,NCP+i) = DBLE((Q(i,j)))
      aux = -Q(i,j)
      do k=1,NCP
        if (j .eq. k) then
            aux = aux * C(k)**(Q(i,k)-1)
        else
            aux = aux * C(k)**Q(i,k)
        end if
      end do
      JC(NCP+i,j) = aux/KS(i)
    end do
  end do
end subroutine jacobian_solid


subroutine amatrix(A, C, P, NCP, NEQ)
  integer, intent(in) :: NCP                             ! number of components
  integer, intent(in) :: NEQ                             ! numnber of equilibria
  integer, INTENT(IN), dimension(NEQ, NCP) :: P          ! The soluble stoichiometry array
  double precision, dimension(NCP+NEQ), intent(in) :: C  ! free concentrations
  double precision, dimension(NCP,NCP), intent(out) :: A ! output values of the A matrix
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  integer :: i, j
  do i=1,NCP
    do j=1,NCP
      if (i .eq. j) then
        A(i,j) = C(i)
      else
        A(i,j) = 0.0
      end if
      do k=1,NEQ
        A(i,j) = A(i,j) + P(k,i)*P(k,j)*C(NCP+k)
      end do
    end do
  end do
end subroutine amatrix


subroutine nr(C, T, B, P, ctrl, NCP, NEQ)
  ! Perform the actual Newton-Raphson iterations. The initial value is supplied
  ! in the C variable and it is replaced with the refined value upon exit.
  ! The free concentrations are refined in the logarithmic space.
  implicit none
  integer, intent(in) :: NCP                               ! number of components
  integer, intent(in) :: NEQ                               ! number of equilibria
  integer, intent(inout) :: ctrl                           ! ctrl -> control value. 
  integer, dimension(NEQ, NCP), intent(in) :: P            ! stoichiometry array
  double precision, dimension(NCP+NEQ), intent(inout) :: C ! initial estimates (in) and refined values (out)
  double precision, dimension(NCP), intent(in) :: T        ! analytical concentration values
  double precision, dimension(NEQ), intent(in) ::  B       ! the equilibrium constants

  integer its, i
  double precision, dimension(NCP)      :: X               ! natural log of free concentrations
  double precision, dimension(NCP)      :: fvec
  double precision, dimension(NCP, NCP) :: jac
  double precision, dimension(NCP)      :: dx
  logical check, converged
  integer MAXITS
  parameter (MAXITS=100)
  
  X(:) = LOG(C(1:NCP))

  call cpluss(C, B, P, NCP, NEQ)
  call mb_soluble(fvec, C, P, T, NCP, NEQ)
  dx(:) = 0.0

  call is_root(fvec, NCP, check)  ! check whether initial guess
  if (check) then                 ! is already a root
    CTRL=1
    return
  end if

  ! Main loopC(i) - 
  do its=1,MAXITS
    call jacobian_soluble(jac, C, P, NCP, NEQ)

    ! Solve linear system.
    call lsolve(jac, fvec, dx, NCP, CTRL)
    if (CTRL == 99) return   ! matrix could not be inverted

    do i=1,NCP
      X(i) = X(i) - dx(i)    ! dx comes sign inverted from lsolve
      C(i) = EXP(X(i))
    end do

    call cpluss(C, B, P, NCP, NEQ)
    call mb_soluble(fvec, C, P, T, NCP, NEQ)
    !call debug(its, fvec, NCP)

    ! test for convergence on relative values
    call is_converged(fvec, T, NCP, converged)
    if (converged) then
       ctrl = 0  ! normal exit
       return
    end if
  end do
  ctrl = 2 ! Exit on MAXITS
end subroutine nr


subroutine debug(i, f, n)
    integer, INTENT(IN) :: i, n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: f
    integer j
    DOUBLE PRECISION s
    s = 0.0
    do j=1,n
        s = s + f(j)**2
    end do
    print "(A, I6, A, E14.6)", "  it=",i, " resid=", s
end subroutine debug


subroutine nr_solids(C, CS, T, B, P, KS, Q, ctrl, do_solids, NCP, NEQ, NSL)
  ! Perform the actual Newton-Raphson iterations. The initial value is supplied
  ! in the X variable and it is replaced with the refined value upon exit.
  implicit none
  integer, intent(in)    :: NCP                            ! number of components
  integer, intent(in)    :: NEQ                            ! number of equilibria
  integer, intent(in)    :: NSL                            ! number of solids
  integer, intent(inout) :: ctrl                           ! control value
  integer,          intent(in),    dimension(NEQ, NCP) :: P         ! stoichiometry array
  double precision, intent(inout), dimension(NCP+NEQ)  :: C         ! initial estimates (in) and refined values (out)
  double precision, intent(inout), dimension(NSL)      :: CS        ! solid concentration initial estimates (in) and refined values (out)
  double precision, intent(in),    dimension(NCP)      :: T         ! analytical concentration values
  double precision, intent(in),    dimension(NEQ)      :: B         ! the equilibrium constants
  integer,          intent(in),    dimension(NSL, NCP) :: Q         ! the equilibrium constants
  double precision, intent(in),    dimension(NSL)      :: KS        ! the equilibrium constants
  logical,          intent(inout), dimension(NSL)      :: do_solids ! which solids are in the model

  integer its, i
  integer nx                                 ! number of refined variables
  double precision, dimension(NCP+NSL) :: X  ! 1:NCP -> free concentrations in log space
                                             ! NCP+1:NCP+NSL -> solid concentration in linear space
  double precision, dimension(NCP+NSL)          :: fvec
  double precision, dimension(NCP+NSL, NCP+NSL) :: jac
  double precision, dimension(NCP+NSL)          :: dx, temp_dx
  double precision, dimension(NSL)              :: SI
  logical check, converged
  integer MAXITS
  parameter (MAXITS=100)
  
  X(1:NCP) = LOG(C(1:NCP))  ! initialise C in log space
  X(NCP+1:NCP+NSL) = CS     ! and CS in linear space

  call cpluss(C, B, P, NCP, NEQ)
  call mb_solids(fvec, C, P, T, Q, CS, NCP, NEQ, NSL)
  call saturation_product(C, Q, NCP, NSL, SI)
  SI = SI / KS    ! convert saturation product to solubility index
  call update_solids(do_solids, SI, NSL)
  fvec(NCP+1:NCP+NSL) = merge(1.0 - SI, 0.0*fvec(NCP+1:NCP+NSL), do_solids)

  dx(:) = 0.0

  call is_root(fvec, NCP, check)  ! check whether initial guess
  if (check) then                 ! is already a root
    CTRL=1
    return
  end if

  ! Main loop
  do its=1,MAXITS
    call jacobian_solid(jac, C, KS, P, Q, NCP, NEQ, NSL, do_solids, nx)

    ! Solve linear system.
    call lsolve(jac, fvec, temp_dx, nx, CTRL)
    if (CTRL == 99) return
    dx(1:NCP) = -temp_dx(1:NCP)
    dx(NCP+1:NCP+NSL) = merge(-temp_dx(NCP+1:NCP+NSL), 0.0*temp_dx(NCP+1:NCP+NSL), do_solids)

    do i=1,NCP
      x(i) = x(i) + dx(i)
      C(i) = EXP(X(i))
    end do
    do i=1,NSL
        if (do_solids(i)) then
            x(i+NCP) = x(i+NCP) + dx(i+NCP)
            CS(i) = X(i+NCP)
        end if
    end do

    call cpluss(C, B, P, NCP, NEQ)
    call mb_solids(fvec, C, P, T, Q, CS, NCP, NEQ, NSL)
    call saturation_product(C, Q, NCP, NSL, SI)
    SI = SI / KS    ! convert saturation product to solubility index
    call update_solids(do_solids, SI, NSL)
    fvec(NCP+1:NCP+NSL) = merge(1.0 - SI, 0.0*fvec(NCP+1:NCP+NSL), do_solids)

    !call debug(its, fvec, NCP+NSL)

    ! test for convergence on relative values
    call is_converged(fvec, T, NCP, converged)
    if (converged) then
       ctrl = 0  ! normal exit
       return
    end if
  end do
  ctrl = 2 ! Exit on MAXITS
end subroutine nr_solids


subroutine is_converged(fvec, T, fsize, retval)
    ! test for convergence on relative values
    implicit none
    integer,          intent(in)                   :: fsize
    double precision, intent(in), dimension(fsize) :: fvec, T
    logical,          intent(inout)                :: retval

    double precision test
    double precision, parameter :: TEST_TOLERANCE = 1e-16
    double precision, parameter :: F_TOLERANCE = 1e-13
    integer i

    retval = .false.
    test = 0.
    do i=1,fsize
       if (T(i) .ne. 0.0) then
         test = test + (fvec(i) / T(i))**2
       end if
    end do
    if ((test <= TEST_TOLERANCE) .and. (maxval(abs(fvec)) <= F_TOLERANCE)) then
        retval = .true.
    end if
end subroutine is_converged


subroutine scale(M, V, D, N)
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N,N) :: M
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N) :: V
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N) :: D
    INTEGER, INTENT(IN) :: N
  
    integer i, j
    do i=1,N
      D(i) = abs(M(i,i))
    end do
    if (any(D .eq. 0.0)) then   ! if there is any diagonal value equal to zero
        D(:) = 1.0              ! do not scale
        return
    end if
    do i=1,N
      do j=1,N
        M(i, j) = M(i, j) / sqrt(D(i)*D(j))
      end do
      v(i) = v(i) / sqrt(D(i))
    end do
end subroutine scale


subroutine is_root(fvec, NCP, check)
  ! Test for initial guess being a root. Use more stringent test than simply TOLF.
  implicit none
  double precision, intent(in), dimension(NCP) :: fvec
  integer, intent(in) :: NCP     ! number of components
  logical, intent(inout) :: check

  DOUBLE PRECISION test
  integer i

  check = .false.
  test=0.
  do i=1,NCP
    if (abs(fvec(i)) .gt. test) test=abs(fvec(i))
  end do
  if (test .lt. 1D-25) then
    check=.true.
  end if
end subroutine is_root


subroutine update_solids(do_solids, si, NSL)
    ! update the solids that are accounted for in the model according to the procedure outlined
    ! in GEOCHEMICAL AND BIOGEOCHEMICAL REACTION MODELING, 2nd ed. CRAIG M. BETHKE 4.4.3p69
    implicit none
    double precision, intent(in),    dimension(NSL) :: si         ! the solubility index
    logical,          intent(inout), dimension(NSL) :: do_solids  ! the solids to be included
    integer,          intent(in) :: NSL                           ! the number of solids
    double precision, parameter :: SOLUBILITY_TOLERANCE = 1E-6

    integer idx

    if (any(si < 1.0 + SOLUBILITY_TOLERANCE)) then    ! Any undersaturated minerals in the system?
        idx = minloc(si, dim=1)                       ! remove most undersaturated component
        !print *, "removed ", idx
        do_solids(idx) = .false.
        return
    end if

    if (any(si > 1.0 - SOLUBILITY_TOLERANCE)) then    ! Any supersaturated mineral?
        idx = maxloc(si, dim=1)                       ! add the most saturated one to the model
        !print *, "added ", idx
        do_solids(idx) = .true.
        return
    end if
end subroutine update_solids

