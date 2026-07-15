! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco
!
! This file contains subroutines for mass balance resolution
!
! Compile with $ f2py3  -c -m mb mb.f90

subroutine mb_soluble(F, C, P, T, NEQ, NCP)
  ! Compute the objective function not considering solids
  ! f_i = c_i + \sum_{j=1}^E { p_{ji}c_{j+S} } - T_i
  integer, intent(in) :: NEQ  ! number of equlibria
  integer, intent(in) :: NCP  ! number of components
  double precision, dimension(NCP), intent(inout) :: F   ! result output
  double precision, dimension(NCP), intent(in) :: T      ! the analytical concentration
  integer, dimension(NEQ, NCP), intent(in):: P           ! the soluble stoichiometry
  double precision, dimension(NCP+NEQ), intent(in) :: C  ! the free concentrations
  !f2py depend(ncp, neq) :: c
  !f2py depend(ncp)      :: t
  !f2py depend(neq, ncp) :: p
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


subroutine mb_soluble_tit(F, C, P, T, NPT, NEQ, NCP)
  ! Compute the objective function not considering solids
  ! f_i = c_i + \sum_{j=1}^E { p_{ji}c_{j+S} } - T_i
  integer, intent(in) :: NPT  ! number of titration points
  integer, intent(in) :: NEQ  ! number of equlibria
  integer, intent(in) :: NCP  ! number of components
  double precision, dimension(NPT, NCP), intent(inout) :: F   ! result output
  double precision, dimension(NPT, NCP), intent(in) :: T      ! the analytical concentration
  integer, dimension(NEQ, NCP), intent(in):: P           ! the soluble stoichiometry
  double precision, dimension(NPT, NCP+NEQ), intent(in) :: C  ! the free concentrations
  !f2py depend(npt, ncp, neq) :: c
  !f2py depend(npt, ncp)      :: t
  !f2py depend(neq, ncp) :: p
  !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  integer i

  DO i=1,NPT
    call mb_soluble(F(i,:), C(i,:), P, T(i,:), NEQ, NCP)
  END DO
END subroutine mb_soluble_tit


subroutine mb_solids(F, C, P, T, Q, CS, NEQ, NCP, NSL)
  ! Compute the objective function not considering solids
  ! f_i = c_i + \sum_{j=1}^E { p_{ji}c_{j+S} } - T_i + \sum_j^{solids} {q_{ij}S_j}
  integer, intent(in) :: NEQ  ! number of equlibria
  integer, intent(in) :: NCP  ! number of components
  double precision, dimension(NCP), intent(inout) :: F   ! result output
  double precision, dimension(NCP), intent(in) :: T      ! the analytical concentration
  integer, dimension(NEQ, NCP), intent(in):: P           ! the soluble stoichiometry
  integer, dimension(NEQ, NSL), intent(in):: Q           ! the solid stoichiometry
  double precision, dimension(NCP+NEQ), intent(in) :: C  ! the free concentrations
  double precision, dimension(NSL), intent(in) :: CS     ! the solids "concentrations"
  !f2py depend(ncp, neq) :: c
  !f2py double precision, depend(nsl)  :: cs
  !f2py depend(ncp)      :: t
  !f2py depend(neq)      :: b
  !f2py depend(neq, ncp) :: p
  !f2py integer, intent(hide), depend(q) :: nsl = shape(Q,1)

  double precision AUX

  call mb_soluble(F, C, P, T, NEQ, NCP)

  do i=1,NCP
    aux = 0.0
    do j=1,NSL
      aux = aux + CS(i)*Q(j,i)
    end do
    F(i) = F(i) + aux
  end do
END subroutine mb_solids


subroutine mb_solids_tit(F, C, P, T, Q, CS, NPT, NEQ, NCP, NSL)
  ! Compute the objective function not considering solids
  ! f_i = c_i + \sum_{j=1}^E { p_{ji}c_{j+S} } - T_i + \sum_j^{solids} {q_{ij}S_j}
  integer, intent(in) :: NPT  ! number of titration points
  integer, intent(in) :: NEQ  ! number of equlibria
  integer, intent(in) :: NCP  ! number of components
  double precision, dimension(NPT, NCP), intent(inout) :: F   ! result output
  double precision, dimension(NPT, NCP), intent(in) :: T      ! the analytical concentration
  integer, dimension(NEQ, NCP), intent(in):: P           ! the soluble stoichiometry
  integer, dimension(NEQ, NSL), intent(in):: Q           ! the solid stoichiometry
  double precision, dimension(NPT, NCP+NEQ), intent(in) :: C  ! the free concentrations
  double precision, dimension(NPT, NSL), intent(in) :: CS     ! the solids "concentrations"
  !f2py depend(ncp, neq) :: c
  !f2py double precision, depend(nsl)  :: cs
  !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
  !f2py depend(ncp)      :: t
  !f2py depend(neq)      :: b
  !f2py depend(neq, ncp) :: p

  integer i
  DO i=1,NPT
    call mb_solids(F(i,:), C(i,:), P, T(i,:), Q, CS(i,:), NEQ, NCP, NSL)
  END DO
END subroutine mb_solids_tit
