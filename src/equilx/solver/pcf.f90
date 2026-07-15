! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! ROUTINE FOR POSITIVE CONTINUOUS FRACTION
!
! by Salvador Blasco <salvador.blasco@uv.es>
!
! Compile with $ f2py3 -c -m mfeq pcf.f90

subroutine pcf_tit(C, T, B, P, THRESH, NPT, NCP, NEQ)
  integer, intent(in) :: NPT    ! number of titration points
  integer, intent(in) :: NCP    ! number of components
  integer, intent(in) :: NEQ    ! numnber of equilibria
  double precision, dimension(NPT, NCP), intent(in) :: T
  double precision, dimension(NEQ), intent(in) :: B
  integer, dimension(NEQ, NCP), intent(in) :: P
  double precision, dimension(NPT,NCP+NEQ), intent(inout) :: C
  ! ==================== f2py directives ====================
  !f2py intent(inout) :: c
  !f2py intent(in)    :: t, b, p, thresh
  !f2py intent(hide), depend(t) :: npt = shape(t,0)
  !f2py intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py intent(hide), depend(p) :: neq = shape(p,0)
  !f2py depend(npt, ncp, neq) :: c
  !f2py depend(npt, ncp, neq) :: c
  ! =======================================================
  do n = 1, NPT
    call pcf(C(n,:), T(n,:), B, P, THRESH, NCP, NEQ)
    ! TODO calc R and copy only if result is better
    ! if (n < NPT) C(n+1,:) = C(n,:)
  end do
end subroutine pcf_tit

subroutine pcf(C, T, B, P, THRESH, NCP, NEQ)
  ! Implementation of the Positive Continuous Fraction (PCF) as described
  ! by Carrayrou (2017)
  integer, intent(in) :: NCP, NEQ                             ! number of species, number of equilibria
  double precision, dimension(NCP+NEQ), intent(inout) :: C
  double precision, dimension(NCP), intent(in) :: T
  double precision, dimension(NEQ), intent(in) :: B
  integer, dimension(NEQ, NCP), intent(in):: P
  real, intent(in) :: THRESH
  !f2py intent(in)    :: ncp, neq, t, b, p, thresh
  !f2py intent(inout) :: c
  !f2py depend(ncp, neq) :: c
  !f2py depend(ncp)      :: t
  !f2py depend(neq)      :: b
  !f2py depend(neq, ncp) :: p

  double precision, dimension(NCP) :: expn, convaux  !, DEBUG_F
  integer, dimension(NEQ+NCP, NCP) :: morel
  integer :: aux, i, j, guard
  logical :: converged
  double precision :: theta, sump, sumr, backtrack

  ! find exponents
  aux = 1
  do i=1,NCP
    do j=1,NEQ
      if ((P(j,i) < aux) .and. (P(j,i) > 0)) then
        aux = P(j,i)
      end if
    end do
    expn(i) = 1/aux
  end do

  ! build Morel table. REMOVE IN THE FUTURE
  do i=1,NCP
    do j=1,NCP
      if (i .eq. j )then
        morel(i,j) = 1
      else
        morel(i,j) = 0
      end if
    end do
    do j=1,NEQ
      morel(j+NCP, i) = P(j, i)
    end do
  end do

  ! main loop
  converged = .false.
  guard = 0
  backtrack = 1.0
  do while (.not. converged)
    guard = guard + 1
    call cpluss(C, B, P, NEQ, NCP)
    ! calculate sump and sumr
    do i=1,NCP
      if (T(i) .ge. 0.0) then
        sump = T(i)
        sumr = 0.0
      else
        sumr = -T(i)
        sump = 0.0
      end if
      do j=1,NCP+NEQ
        if (morel(j,i) > 0) then
          sumr = sumr + morel(j,i) * C(j)
        else
          sump = sump - morel(j,i) * C(j)
        end if
      end do
      if (sump < sumr) then
        theta = backtrack * (0.9 - 0.8 * sump / sumr)
      else
        theta = backtrack * (0.9 - 0.8 * sumr / sump)
      end if
      ! update concentrations
      C(i) = theta * C(i) * (sump/sumr)**expn(i) + (1-theta)*C(i)
      ! check convergence
      convaux(i) = abs(sump-sumr)/(sump + sump)
    end do

    converged = maxval(convaux) .lt. THRESH
    if (guard == 20 .or. guard == 40 .or. guard == 60) then
        backtrack = backtrack * 0.5
    end if
    if (guard >100) then
        converged = .true.
    end if
  end do
end subroutine pcf
