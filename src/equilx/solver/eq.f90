! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>
!
! Compile with $ f2py3 -llapack -c -m solver cpluss.o pcf.o eq.f90
!
! control codes
! 0 -> normal exit, no errors
! 1 -> unable to invert matrix

subroutine solve1(C, T, B, P, CTRL, NPT, NEQ, NCP)
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
    if (i > 1) then
      C(i,:) = C(i-1,:)
    end if
    if (i == 1) then
      call pcf(C(i,:), T(i,:), B, P, 1E-9, NCP, NEQ)
    else
      call pcf(C(i,:), T(i,:), B, P, 0.25, NCP, NEQ)
    end if
    call nr(C(i,:), T(i,:), B, P, CTRL, NEQ, NCP)
  end do
end subroutine solve1

subroutine solve1_solids(C, CS, T, B, P, KS, Q, CTRL, NPT, NEQ, NCP)
end subroutine solve1_solids

subroutine solve2(C, T, B, P, CTRL, NPT, NEQ, NCP)
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
    if (i > 1) then
      C(i,:) = C(i-1,:)
    end if
    if (i == 1) then
      call pcf(C(i,:), T(i,:), B(i,:), P, 1E-9, NCP, NEQ)
    else
      call pcf(C(i,:), T(i,:), B(i,:), P, 0.25, NCP, NEQ)
    end if
    call nr(C(i,:), T(i,:), B(i,:), P, CTRL, NEQ, NCP)
  end do
end subroutine solve2

subroutine solve2_solids(C, CS, T, B, P, KS, Q, CTRL, NPT, NEQ, NCP)
end subroutine solve2_solids

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


subroutine fobj(F, C, P, T, E, S)
  ! Compute the objective function
  ! f_i = c_i + \sum_{j=1}^E { p_{ji}c_{j+S} } - T_i
  integer, intent(in) :: E, S
  double precision, dimension(S), intent(inout) :: F
  double precision, dimension(S), intent(in) :: T
  integer, dimension(E, S), intent(in):: P
  double precision, dimension(S+E), intent(in) :: C
  double precision AUX
  integer i, j
  DO i=1,S
    AUX = C(i)-T(i)
    DO j=1,E
      AUX = AUX+P(j,i)*C(j+S)
    END DO
    F(i)=AUX
  END DO
END subroutine fobj


subroutine jacobian(Jac, C, P, E, S)
  ! Compute the free concentrations for all the components.
  ! J_{ij} = \delta_{ij} + c_j^{-1} \sum_{k=1}^E {
  !        p_{ki}p_{kj}c_{k+S} }
  integer, intent(in) :: E, S
  integer, dimension(E, S), intent(in) :: P
  double precision, dimension(S+ E), intent(in) :: C
  double precision, dimension(S, S), intent(inout) :: Jac
  integer i, j, k
  double precision aux
  do i=1,S
    do j=1,S
      aux = 0.0
      do k=1,E
        aux = aux + P(k,i)*P(k,j)*C(k+s)
      end do
      Jac(i,j)=aux/C(j)
      if (j .eq. i) Jac(i,j) = 1.0 + Jac(i,j)
    end do
  end do
end subroutine jacobian


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
