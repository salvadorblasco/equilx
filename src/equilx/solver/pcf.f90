! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! ROUTINE FOR POSITIVE CONTINUOUS FRACTION
!
! by Salvador Blasco <salvador.blasco@uv.es>
!
! Compile with $ f2py3 -c -m mfeq pcf.f90


subroutine pcf(C, T, B, P, THRESH, NCP, NEQ)
  ! Implementation of the Positive Continuous Fraction (PCF) as described
  ! by Carrayrou (2017)
  integer, intent(in) :: NCP                                ! number of species
  integer, intent(in) :: NEQ                                ! number of equilibria
  double precision, dimension(NCP+NEQ), intent(inout) :: C  ! the free concentrations
  double precision, dimension(NCP), intent(in) :: T         ! the analytical concentration
  double precision, dimension(NEQ), intent(in) :: B         ! equilibrium constants
  integer, dimension(NEQ, NCP), intent(in):: P              ! stoichiometry matrix
  double precision, intent(in) :: THRESH                                ! acceptable threshold
  !f2py intent(in)       :: ncp, neq, t, b, p, thresh
  !f2py intent(inout)    :: c
  !f2py depend(ncp, neq) :: c
  !f2py depend(ncp)      :: t
  !f2py depend(neq)      :: b
  !f2py depend(neq, ncp) :: p

  double precision, dimension(NCP) :: expn, convaux, convbest  !, DEBUG_F
  integer :: aux, i, j, iter
  logical :: converged
  double precision, dimension(NCP) :: theta
  double precision :: sump, sumr, aux_theta
  double precision, parameter :: THETA_ALPHA = 0.9
  double precision, parameter :: THETA_BETA = 0.8
  integer, parameter :: MAXITS = 5000

  do i=1,NCP                                        ! find the smallest strictly positive
    aux = maxval(P(:,i))                            ! stoichiometric coefficient that 
    aux = max(aux, 1)                               ! the jth component presents in Morel's Table.
    do j=1,NEQ
      if ((P(j,i) < aux) .and. (P(j,i) > 0)) then
        aux = P(j,i)
      end if
    end do
    expn(i) = 1.0/aux
  end do

  ! main loop
  converged = .false.
  iter = 0
  convbest(:) = 1.1
  do while (.not. converged)
    iter = iter + 1
    call cpluss(C, B, P, NCP, NEQ)

    ! calculate sump and sumr
    do i=1,NCP
      if (T(i) .ge. 0.0) then
        sump = T(i)
        sumr = C(i)
      else
        sump = 0.0
        sumr = abs(T(i)) + C(i)
      end if
      do j=1,NEQ
        if (P(j,i) > 0) then
          sumr = sumr + P(j,i) * C(j+NCP)
        else
          sump = sump + abs(P(j,i)) * C(j+NCP)
        end if
      end do
      aux_theta = THETA_ALPHA - THETA_BETA * min(sump/sumr, sumr/sump)
      ! check convergence
      convaux(i) = abs(sump-sumr)
      if ((sumr > 0.0 .and. sump > 0.0) .and. ((sumr + sump) > 0.0)) then   ! if one of the components is zero the division
          convaux(i) = convaux(i)/(sumr+sump)                               ! will be 1 thus ruining the value of convergence testing
      end if                                                                ! if the sum of both is zero, also problem
      if (convaux(i) < convbest(i)) then
          convbest(i) = convaux(i)
          theta(i) = max(aux_theta, theta(i))
        else
          theta(i) = min(aux_theta, theta(i))
      end if
      ! update concentrations
      C(i) = theta(i) * C(i) * (sump/sumr)**expn(i) + (1-theta(i))*C(i)
    end do
      !print *, iter, convaux, theta
    converged = (maxval(convaux) .lt. THRESH)
    !if (mod(iter,10) .eq. 0) print "(I5, 4E12.4)", iter, convaux, theta

    if (iter > MAXITS) then     ! exit on max iterations
        converged = .true.
    end if
  end do
  ! update C before finishing
  call cpluss(C, B, P, NCP, NEQ)
end subroutine pcf
