! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! Compile with
! $ uv run f2py -c -m si si.f90
!  OR
! $ make si

subroutine saturation_product(C, Q, NCP, NSL, SI)
  integer, intent(in) :: NCP, NSL                  ! number of components, number of solids
  double precision, dimension(NCP), intent(in) :: C
  integer, dimension(NSL, NCP), intent(in):: Q
  double precision, DIMENSION(NSL), INTENT(INOUT) :: SI
  !f2py intent(hide), depend(Q) :: NCP=shape(Q,2), NSL=shape(Q,1)

  double PRECISION :: aux

  SI(:) = 0.0
  do i=1,NSL
    aux = 1.0
    do j=1,NCP
      aux = aux*C(j)**Q(i,j)
    end do
    SI(i) = aux
  end do
end subroutine saturation_product


subroutine saturation_product_tit(C, Q, NPT, NCP, NSL, SI)
  integer, intent(in) :: NPT, NCP, NSL          ! nomber of titration points, number of species, number of solids
  double precision, dimension(NPT,NCP), intent(in) :: C
  integer, dimension(NSL, NCP), intent(in):: Q
  double precision, DIMENSION(NPT, NSL), INTENT(OUT) :: SI
  !f2py intent(hide), depend(Q) :: NCP=shape(Q,2), NSL=shape(Q,1)

  do n=1,NPT
    call saturation_product(C(n,:), Q, NCP, NSL, SI(n,:))
  end do
end
