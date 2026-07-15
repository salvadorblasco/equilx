
! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco

subroutine cpluss_tit(C, B, P, N, E, S)
  ! Compute the free concentrations for all the components.
  integer N, E, S
  double precision C(N, S+E)
  integer P(E, S)
  double precision B(E)
  intent(in) N, E, S, B, P
  intent(inout) C
  integer i
  DO i=1,N
    call cpluss(C(i,:), B, P, E, S)
  END DO
END subroutine cpluss_tit


subroutine cpluss(C, B, P, E, S)
  ! Compute the free concentrations for all the components.
  integer, intent(in) :: E
  integer, intent(in) :: S
  double precision, dimension(S+E), intent(inout) :: C
  integer, dimension(E, S), intent(in):: P
  double precision, dimension(E) :: B
  double precision AUX
  integer i, j
  do i=1,E
    AUX = 1D0
    do j=1,S
      AUX = AUX*C(j)**P(i,j)
    end do
    C(i+S) = B(i)*AUX
  end do
end subroutine cpluss
