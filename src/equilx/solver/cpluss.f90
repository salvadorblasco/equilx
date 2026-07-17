! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco

subroutine cpluss(C, B, P, NCP, NEQ)
  ! Compute the free concentrations for all the components.
  integer,          intent(in) :: NEQ    ! number of equilibria
  integer,          intent(in) :: NCP    ! number of components
  double precision, intent(inout), dimension(NCP+NEQ)  :: C  ! free concentrations
  integer,          intent(in),    dimension(NEQ, NCP) :: P  ! stoichiometry
  double precision, intent(in),    dimension(NEQ)      :: B  ! equilibrium constants

  double precision AUX
  integer i, j

  do i=1,NEQ
    AUX = 1D0
    do j=1,NCP
      AUX = AUX*C(j)**P(i,j)
    end do
    C(i+NCP) = B(i)*AUX
  end do
end subroutine cpluss
