! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>
!
! Compile with $ f2py3 -llapack -c -m solver cpluss.o pcf.o mb.o eq.f90
!
! control codes
! 0 -> normal exit, no errors
! 1 -> unable to invert matrix

subroutine solve1_tit(C, T, B, P, CTRL, NPT, NEQ, NCP)
  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 1D. The first initial value must be provided in the first row
  ! of C. The following points take as starting value the value of 
  ! the last iteration.

  implicit none
  integer, intent(in) :: NPT     ! number of titration points
  integer, intent(in) :: NCP     ! number of components
  integer, intent(in) :: NEQ     ! number of equilibria
  integer, intent(inout) :: CTRL ! control value. 0 on no error, 1 otherwise
  integer, dimension(NEQ, NCP), intent(in) :: P ! The stoichiometry array
  double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C ! initial estimates (in) and refined values (out)
  double precision, dimension(NPT, NCP), intent(in) :: T ! analytical concentration values
  double precision, dimension(NEQ), intent(in) :: B  ! the equilibrium constants
  !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  integer :: i
  CTRL = 0
  do i=1,NPT
    if (i == 1) then
      call solve_point(C(i,:), T(i,:), B, P, CTRL, .true., NEQ, NCP)
    else
      C(i,:) = C(i-1,:)
      call solve_point(C(i,:), T(i,:), B, P, CTRL, .false., NEQ, NCP)
    end if
  end do
end subroutine solve1_tit

subroutine solve1_solids_tit(C, CS, T, B, P, KS, Q, CTRL, NPT, NEQ, NCP, NSL)
  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 1D and solids are present. The first initial value must be provided
  ! in the first row of C. The following points take as starting value the value of 
  ! the last iteration.

  implicit none
  integer, intent(in) :: NPT     ! number of titration points
  integer, intent(in) :: NCP     ! number of components
  integer, intent(in) :: NEQ     ! number of equilibria
  integer, intent(in) :: NSL     ! number of solids
  integer, intent(inout) :: CTRL ! control value. 0 on no error, 1 otherwise
  integer, dimension(NEQ, NCP), intent(in) :: P ! The soluble stoichiometry array
  integer, dimension(NSL, NCP), intent(in) :: Q ! The solids stoichiometry array
  double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C ! free concentation initial estimates (in) and refined values (out)
  double precision, dimension(NPT, NSL, intent(inout) :: CS     ! solid concentration initial estimates (in) and refined values (out)
  double precision, dimension(NPT, NCP), intent(in) :: T        ! analytical concentration values
  double precision, dimension(NEQ), intent(in) :: B  ! the equilibrium constants
  !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  integer :: i
  CTRL = 0

  do i=1,NPT
    if (i > 1) then
      C(i,:) = C(i-1,:)
      CS(i,:) = CS(i-1,:)
    end if
    if (i == 1) then
      call pcf(C(i,:), T(i,:), B, P, 1E-9, NCP, NEQ)
    else
      call pcf(C(i,:), T(i,:), B, P, 0.25, NCP, NEQ)
    end if
    call nr(C(i,:), T(i,:), B, P, CTRL, NEQ, NCP)
  end do
end subroutine solve1_solids_tit

subroutine solve2_tit(C, T, B, P, CTRL, NPT, NEQ, NCP)
  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 2D. The first initial value must be provided in the first row
  ! of C. The following points take as starting value the value of 
  ! the last iteration.

  implicit none
  integer, intent(in) :: NPT     ! number of titration points
  integer, intent(in) :: NCP     ! number of components
  integer, intent(in) :: NEQ     ! number of equilibria
  integer, intent(inout) :: CTRL ! control value. 0 on no error, 1 otherwise
  integer, dimension(NEQ, NCP), intent(in) :: P ! The stoichiometry array
  double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C ! initial estimates (in) and refined values (out)
  double precision, dimension(NPT, NCP), intent(in) :: T ! analytical concentration values
  double precision, dimension(NPT, NEQ), intent(in) :: B  ! the equilibrium constants

  !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  integer :: i

  CTRL = 0

  do i=1,NPT
    if (i == 1) then
      call solve_point(C(i,:), T(i,:), B(i,:), P, CTRL, .true., NEQ, NCP)
    else
      C(i,:) = C(i-1,:)
      call solve_point(C(i,:), T(i,:), B(i,:), P, CTRL, .false., NEQ, NCP)
    end if
  end do
end subroutine solve2_tit

subroutine solve2_solids_tit(C, CS, T, B, P, KS, Q, CTRL, NPT, NEQ, NCP)
end subroutine solve2_solids_tit


subroutine solve_point(C, T, B, P, CTRL, is_first, NEQ, NCP)
  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 1D. The first initial value must be provided in the first row
  ! of C. The following points take as starting value the value of 
  ! the last iteration.

  implicit none
  integer, intent(in) :: NCP     ! number of components
  integer, intent(in) :: NEQ     ! number of equilibria
  integer, intent(inout) :: CTRL ! control value. 0 on no error, 1 otherwise
  logical, intent(in) :: is_first   ! whther it is first point or not
  integer, dimension(NEQ, NCP), intent(in) :: P ! The stoichiometry array
  double precision, dimension(NEQ+NCP), intent(inout) :: C ! initial estimates (in) and refined values (out)
  double precision, dimension(NCP), intent(in) :: T ! analytical concentration values
  double precision, dimension(NEQ), intent(in) :: B  ! the equilibrium constants
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  CTRL = 0

  if (is_first) then
    call pcf(C, T, B, P, 1E-9, NCP, NEQ)
  else
    call pcf(C, T, B, P, 0.25, NCP, NEQ)
  end if
  call nr(C, T, B, P, CTRL, NEQ, NCP)
end subroutine solve_point


subroutine solve_point_solid(C, CS, T, B, P, KS, Q, CTRL, is_first, NEQ, NCP, NSL)
  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 1D. The first initial value must be provided in the first row
  ! of C. The following points take as starting value the value of 
  ! the last iteration.

  implicit none
  integer, intent(in) :: NCP     ! number of components
  integer, intent(in) :: NEQ     ! number of equilibria
  integer, intent(in) :: NSL     ! number of solids
  integer, intent(inout) :: CTRL ! control value. 0 on no error, 1 otherwise
  logical, intent(in) :: is_first   ! whther it is first point or not
  integer, dimension(NEQ, NCP), intent(in) :: P ! The stoichiometry array
  double precision, dimension(NEQ+NCP), intent(inout) :: C ! initial estimates (in) and refined values (out)
  double precision, dimension(NSL), intent(inout) :: CS     ! solid concentration initial estimates (in) and refined values (out)
  double precision, dimension(NCP), intent(in) :: T ! analytical concentration values
  double precision, dimension(NEQ), intent(in) :: B  ! the equilibrium constants
  double precision, dimension(NSL, NCP), intent(in) :: Q  ! the equilibrium constants
  double precision, dimension(NSL), intent(in) :: KS  ! the equilibrium constants
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  real, PARAMETER :: si_tolerance = 1e-6
  double PRECISION, dimension(NSL) :: si
  logical, dimension(NSL) :: do_solid
  integer :: i
  do_solid(:) = .false.

  call solve_point(C, T, B, P, CTRL, is_first, NEQ, NCP)
  do
    call saturation_product(C, Q, si)
    do i=1,NSL
      if (si(i)/KS(i) < 1.0 - si_tolerance) then
        do_solid(i) = .false.  ! remove from model
        cycle
      end if
      if (si(i)/KS(i) > 1.0 + si_tolerance) then
        do_solid(i) = .true.
        cycle
      end if
    end do
    if (.not. any(do_solid)) exit    ! no solids to do
    call nr_solids()
  end do
end subroutine solve_point_solid


subroutine lsolve(A, b, x, n, CTRL)
  ! Solve system of lineal equations using LAPACK

  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: A
  ! The nxn array of coefficients
  double precision, dimension(n), intent(in) :: b
  ! the n array of independent terms
  double precision, dimension(n), intent(out) :: x
  ! the n array where the result is returned
  integer, intent(inout) :: CTRL ! control value. 0 on no error, 1 otherwise

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
      CTRL = 1
      return
  end if
end subroutine lsolve


subroutine jacobian_soluble(J, C, P, NCP, NEQ)
  call amatrix(J, C, P, NCP, NEQ)
end subroutine jacobian_soluble

subroutine jacobian_solid(J, C, P, E, S)
end subroutine jacobian_solid


subroutine amatrix(A, C, P, NCP, NEQ)
  integer, intent(in) :: NCP    ! number of components
  integer, intent(in) :: NEQ    ! numnber of equilibria
  double precision, dimension(NCP+NEQ), intent(out) :: C   ! free concentrations
  double precision, dimension(NCP,NCP), intent(inout) :: A ! output values
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


subroutine nr(X, T, B, P, ctrl, NEQ, NCP)
  ! Perform the actual Newton-Raphson iterations. The initial value is supplied
  ! in the X variable and it is replaced with the refined value upon exit.
  implicit none
  integer, intent(in) :: NEQ, NCP
  integer, intent(inout) :: ctrl
  ! ctrl -> control value. 
  !         0 = normal exit
  !         2 = normal exit upon initial guess being a root
  !         99 = abnormal exit on MAXITS
  integer, dimension(NEQ, NCP), intent(in) :: P
  double precision, dimension(NCP+NEQ), intent(inout) :: X
  double precision, dimension(NCP), intent(in) :: T
  double precision, dimension(NEQ), intent(in) ::  B

  integer its, i, j
  double precision, dimension(NCP) :: fvec, jdiag
  double precision, dimension(NCP, NCP) ::  jac
  double precision :: f
  double precision, dimension(NCP) :: xold
  double precision, dimension(NCP) :: dx
  double precision test
  logical check
  integer MAXITS
  double precision relax
  parameter (MAXITS=100)
  
  f = 0.0
  do i=1,NCP
    xold(i) = X(i)
  end do
  call cpluss(X, B, P, NEQ, NCP)
  call fobj(fvec, X, P, T, NEQ, NCP)

  ! Test for initial guess being a root. Use more stringent test than simply TOLF.
  test=0.
  do i=1,NCP
    if (abs(fvec(i)) .gt. test) test=abs(fvec(i))
    dx(i) = 0.0
  end do
  if (test .lt. 1D-25) then
    check=.false.
    ctrl=2
    return
  end if

  ! Main loop
  do its=1,MAXITS
    call jacobian(jac, X, P, NEQ, NCP)

    ! Scale
    do i=1,NCP
      jdiag(i) = jac(i,i)
    end do
    do i=1,NCP
      do j=1,NCP
        jac(i, j) = jac(i, j) / sqrt(jdiag(i)*jdiag(j))
      end do
      fvec(i) = fvec(i) / sqrt(jdiag(i))
    end do

    ! Solve linear system. Includes negation of fvec
    call lsolve(jac, fvec, dx, NCP)

    ! De-scale
    do i=1,NCP
      dx(i) = - dx(i) / sqrt(jdiag(i)) 
    end do

    do i=1,NCP
      relax = 1.0 / max(1.0, -2*dx(i)/x(i))
      x(i) = x(i) + relax*dx(i)
    end do
    call cpluss(X, B, P, NEQ, NCP)
    call fobj(fvec, X, P, T, NEQ, NCP)

    ! test for convergence on relative values
    test = 0.
    do i=1,NCP
       if (T(i) .ne. 0.0) then
         test = test + (fvec(i) / T(i))**2
       end if
    end do
    if ((test <= 1e-16) .and. (maxval(abs(fvec)) <= 1e-9)) then
       ctrl = 0
       !DEBUG! print *, "Normal exit", its, sum(fvec**2), x(1:2)
       return
    end if
  end do
  ctrl=99
  !! print *, "Exit on MAXITS"
end subroutine nr


subroutine nr_solids(C, CS, T, B, P, KS, Q, ctrl, do_solids, NEQ, NCP, NSL)
  ! Perform the actual Newton-Raphson iterations. The initial value is supplied
  ! in the X variable and it is replaced with the refined value upon exit.
  implicit none
  integer, intent(in) :: NCP     ! number of components
  integer, intent(in) :: NEQ     ! number of equilibria
  integer, intent(in) :: NSL     ! number of solids
  integer, intent(inout) :: ctrl
  integer, dimension(NEQ, NCP), intent(in) :: P
  double precision, dimension(NCP+NEQ), intent(inout) :: C ! initial estimates (in) and refined values (out)
  double precision, dimension(NSL), intent(inout) :: CS    ! solid concentration initial estimates (in) and refined values (out)
  double precision, dimension(NCP), intent(in) :: T ! analytical concentration values
  double precision, dimension(NEQ), intent(in) :: B  ! the equilibrium constants
  double precision, dimension(NSL, NCP), intent(in) :: Q  ! the equilibrium constants
  double precision, dimension(NSL), intent(in) :: KS  ! the equilibrium constants
  logical, dimension(NSL), intent(in) :: do_solid
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
  !f2py integer, intent(hide), depend(q) :: nsl = shape(q,0)

  integer its, i, j
  double precision, dimension(NCP) :: fvec, jdiag
  double precision, dimension(NCP, NCP) ::  jac
  double precision :: f
  double precision, dimension(NCP) :: xold
  double precision, dimension(NCP) :: dx
  double precision test
  logical check
  integer MAXITS
  double precision relax
  parameter (MAXITS=100)
  
  f = 0.0
  do i=1,NCP
    xold(i) = X(i)
  end do
  call cpluss(C, B, P, NEQ, NCP)
  call mb_solids(fvec, C, P, T, Q, CS, NEQ, NCP, NSL)

  ! Test for initial guess being a root. Use more stringent test than simply TOLF.
  test=0.
  do i=1,NCP
    if (abs(fvec(i)) .gt. test) test=abs(fvec(i))
    dx(i) = 0.0
  end do
  if (test .lt. 1D-25) then
    check=.false.
    ctrl=2
    return
  end if

  ! Main loop
  do its=1,MAXITS
    call jacobian(jac, X, P, NEQ, NCP)

    ! Scale
    call scale(jac, fvec, NCP)

    ! Solve linear system. Includes negation of fvec
    call lsolve(jac, fvec, dx, NCP)

    ! De-scale
    do i=1,NCP
      dx(i) = - dx(i) / sqrt(jdiag(i)) 
    end do

    do i=1,NCP
      relax = 1.0 / max(1.0, -2*dx(i)/x(i))
      x(i) = x(i) + relax*dx(i)
    end do
    call cpluss(C, B, P, NEQ, NCP)
    call mb_solids(fvec, C, P, T, Q, CS, NEQ, NCP, NSL)

    ! test for convergence on relative values
    test = 0.
    do i=1,NCP
       if (T(i) .ne. 0.0) then
         test = test + (fvec(i) / T(i))**2
       end if
    end do
    if ((test <= 1e-16) .and. (maxval(abs(fvec)) <= 1e-9)) then
       ctrl = 0
       !DEBUG! print *, "Normal exit", its, sum(fvec**2), x(1:2)
       return
    end if
  end do
  ctrl=99
  !! print *, "Exit on MAXITS"
end subroutine nr_solids


subroutine scale(M, V, N)
  DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N,N) :: M
  DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N) :: V
  INTEGER, INTENT(IN) :: N

  DOUBLE PRECISION, DIMENSION(N) :: diag
  integer i, j
    do i=1,N
      diag(i) = M(i,i)
    end do
    do i=1,N
      do j=1,N
        M(i, j) = M(i, j) / sqrt(diag(i)*diag(j))
      end do
      v(i) = v(i) / sqrt(diag(i))
    end do
end subroutine scale
