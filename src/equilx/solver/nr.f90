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


subroutine jacobian_solid(J, C, SI, P, Q, NCP, NEQ, NSL)
  implicit none
  integer, intent(in) :: NCP                                         ! number of components
  integer, intent(in) :: NEQ                                         ! number of equilibria
  integer, intent(in) :: NSL                                         ! number of solids
  integer,          intent(in),    dimension(NEQ, NCP)         :: P  ! The soluble stoichiometry array
  integer,          intent(in),    dimension(NSL, NCP)         :: Q  ! The solids stoichiometry array
  double precision, intent(in),    dimension(NEQ+NCP)          :: C  ! free concentation initial estimates (in) and refined values (out)
  double precision, intent(in),    dimension(NSL)              :: SI ! the solubility index
  double precision, INTENT(INOUT), dimension(NCP+NSL, NCP+NSL) :: J  ! jacobian
  !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  integer i, k
  call jacobian_soluble(J(1:NCP,1:NCP), C, P, NCP, NEQ)
  J(1:NCP,NCP+1:) = DBLE(Q(:,:))
  J(NCP+1:NCP+NSL,NCP+1:NCP+NSL) = 0.0
  do i=1,NCP
    do k=1, NSL
       J(i,k) = -Q(k,i)*SI(i)/C(i)
    end do
  end do
end subroutine jacobian_solid


subroutine amatrix(A, C, P, NCP, NEQ)
  integer, intent(in) :: NCP                               ! number of components
  integer, intent(in) :: NEQ                               ! numnber of equilibria
  integer, INTENT(IN), dimension(NEQ, NCP) :: P            ! The soluble stoichiometry array
  double precision, dimension(NCP+NEQ), intent(in) :: C    ! free concentrations
  double precision, dimension(NCP,NCP), intent(inout) :: A ! output values of the A matrix
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
  double precision, dimension(NCP)      :: fvec, jdiag
  double precision, dimension(NCP, NCP) :: jac
  double precision, dimension(NCP)      :: dx
  double precision test
  logical check, converged
  integer MAXITS
  parameter (MAXITS=100)
  
  X(:) = LOG(C(1:NCP))

  call cpluss(C, B, P, NCP, NEQ)
  call mb_soluble(fvec, C, P, T, NCP, NEQ)
  dx(:) = 0.0

  call is_root(fvec, NCP, check)  ! check whether initial guess
  if (check) then             ! is already a root
    CTRL=1
    return
  end if

  ! Main loopC(i) - 
  do its=1,MAXITS
    call jacobian_soluble(jac, C, P, NCP, NEQ)

    ! Scale
    call scale(jac, fvec, jdiag, NCP)

    ! Solve linear system. Includes negation of fvec
    call lsolve(jac, fvec, dx, NCP, CTRL)
    if (CTRL == 99) return   ! matrix could not be inverted

    ! De-scale
    do i=1,NCP
      dx(i) = - dx(i) / sqrt(jdiag(i)) 
    end do

    do i=1,NCP
      X(i) = X(i) + dx(i)
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
    !print *, "DEBUG 4 > it 0", its, fvec
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
    print "(I4, E12.6)", i, s
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
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
  !f2py integer, intent(hide), depend(q) :: nsl = shape(q,0)

  integer its, i
  double precision, dimension(NCP+NSL) :: X  ! 1:NCP -> free concentrations in log space
                                             ! NCP+1:NCP+NSL -> solid concentration in linear space
  double precision, dimension(NCP)              :: fvec, jdiag
  double precision, dimension(NCP+NSL, NCP+NSL) :: jac
  double precision, dimension(NCP+NSL)          :: dx
  double precision, dimension(NSL)              :: SI
  logical check, converged
  integer MAXITS
  parameter (MAXITS=100)
  
  X(1:NCP) = LOG(C(1:NCP))  ! initialise C in log space
  X(NCP+1:NCP+NSL) = CS     ! and CS in linear space

  call cpluss(C, B, P, NEQ, NCP)
  call mb_solids(fvec, C, P, T, Q, CS, NCP, NEQ, NSL)

  dx(:) = 0.0

  call is_root(fvec, NCP, check)  ! check whether initial guess
  if (check) then             ! is already a root
    CTRL=1
    return
  end if

  ! Main loop
  do its=1,MAXITS
    call saturation_product(C, Q, NCP, NSL, SI)
    SI = SI / KS    ! convert saturation product to solubility index
    do i=1,NSL
        if (SI(i) < 1.0) do_solids(i) = .False.
        if (SI(i) > 1.0) do_solids(i) = .True.
    end do
    call jacobian_solid(jac, C, SI, P, Q, NCP, NEQ, NSL)

    ! Scale
    call scale(jac, fvec, jdiag, NCP)

    ! Solve linear system. Includes negation of fvec
    call lsolve(jac, fvec, dx, NCP, CTRL)
    if (CTRL == 99) return

    ! De-scale
    do i=1,NCP
      dx(i) = - dx(i) / sqrt(jdiag(i)) 
    end do

    do i=1,NCP
      x(i) = x(i) + dx(i)
      C(i) = EXP(X(i))
    end do
    do i=1,NSL
        if (do_solids(i)) then
            x(i) = x(i) + dx(i)
            if (X(NCP+i) > 0.0) then
                CS(i) = X(i)
            else
                CS(i) = 0.0
            end if
        end if
    end do
    call cpluss(C, B, P, NEQ, NCP)
    call mb_solids(fvec, C, P, T, Q, CS, NCP, NEQ, NSL)

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
    double precision, parameter :: F_TOLERANCE = 1e-10
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

