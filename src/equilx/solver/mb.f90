! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco
!
! This file contains subroutines for mass balance resolution
!
! Compile with $ f2py3  -c -m mb mb.f90

subroutine mb_soluble(F, C, P, T, NCP, NEQ)
  ! Compute the objective function not considering solids
  ! f_i = c_i + \sum_{j=1}^E { p_{ji}c_{j+S} } - T_i
  integer, intent(in) :: NEQ  ! number of equlibria
  integer, intent(in) :: NCP  ! number of components
  double precision, dimension(NCP), intent(inout) :: F   ! result output
  double precision, dimension(NCP), intent(in) :: T      ! the analytical concentration
  integer, dimension(NEQ, NCP), intent(in):: P           ! the soluble stoichiometry
  double precision, dimension(NCP+NEQ), intent(in) :: C  ! the free concentrations
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  double precision AUX
  integer i, j
  DO i=1,NCP
    AUX = C(i)-T(i)
    DO j=1,NEQ
      AUX = AUX+P(j,i)*C(j+NCP)
    END DO
    F(i)=AUX
  END DO
END subroutine mb_soluble


subroutine mb_solids(F, C, P, T, Q, CS, NCP, NEQ, NSL)
  ! Compute the objective function not considering solids
  ! f_i = c_i + \sum_{j=1}^E { p_{ji}c_{j+S} } - T_i + \sum_j^{solids} {q_{ij}S_j}
  integer, intent(in) :: NEQ  ! number of equlibria
  integer, intent(in) :: NCP  ! number of components
  double precision, dimension(NCP), intent(inout) :: F   ! result output
  double precision, dimension(NCP), intent(in) :: T      ! the analytical concentration
  integer, dimension(NEQ, NCP), intent(in):: P           ! the soluble stoichiometry
  integer, dimension(NSL, NCP), intent(in):: Q           ! the solid stoichiometry
  double precision, dimension(NCP+NEQ), intent(in) :: C  ! the free concentrations
  double precision, dimension(NSL), intent(in) :: CS     ! the solids "concentrations"
  !f2py depend(ncp, neq) :: c
  !f2py double precision, depend(nsl)  :: cs
  !f2py depend(ncp)      :: t
  !f2py depend(neq)      :: b
  !f2py depend(neq, ncp) :: p
  !f2py integer, intent(hide), depend(q) :: nsl = shape(Q,1)

  call mb_soluble(F, C, P, T, NCP, NEQ)
  do i=1,NCP
    do j=1,NSL
      F(i) = F(i) + CS(j)*Q(j,i)
    end do
  end do
END subroutine mb_solids
