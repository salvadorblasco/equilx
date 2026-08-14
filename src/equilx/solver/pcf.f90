! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! ROUTINE FOR POSITIVE CONTINUOUS FRACTION
!
! by Salvador Blasco <salvador.blasco@uv.es>


subroutine pcf(C, T, B, P, THRESH, NCP, NEQ)
  ! Implementation of the Positive Continuous Fraction (PCF) as described
  ! by Carrayrou (2004) and Marinoni (2017).
  integer,          intent(in) :: NCP                        ! number of species
  integer,          intent(in) :: NEQ                        ! number of equilibria
  double precision, intent(in) :: THRESH                     ! acceptable threshold
  double precision, intent(inout), dimension(NCP+NEQ)  :: C  ! the free concentrations
  double precision, intent(in),    dimension(NCP)      :: T  ! the analytical concentration
  double precision, intent(in),    dimension(NEQ)      :: B  ! equilibrium constants
  integer,          intent(in),    dimension(NEQ, NCP) :: P  ! stoichiometry matrix

  double precision, dimension(NCP) :: expn, convaux
  integer :: i, iter
  double precision :: theta, check_convergence
  double precision            :: sump, sumr, calc_theta
  integer,          parameter :: MAXITS = 100


  call find_exponent(P, NCP, NEQ, expn)

  do iter=1,MAXITS
    call cpluss(C, B, P, NCP, NEQ)

    ! calculate sump and sumr
    do i=1,NCP
        call calc_sum(i, T, C, P, sump, sumr, NCP, NEQ)
        convaux(i) = check_convergence(sumr, sump)
        theta = calc_theta(sumr, sump)
        C(i) = theta * C(i) * (sump/sumr)**expn(i) + (1-theta)*C(i)
    end do
    !print *, iter, convaux
    if (maxval(convaux) .lt. THRESH) exit ! converged
  end do
  ! update C before finishing
  call cpluss(C, B, P, NCP, NEQ)
end subroutine pcf


subroutine pcf_solids(C, CS, T, B, KS, P, Q, THRESH, NCP, NEQ, NSL)
    ! Implementation of the Positive Continuous Fraction (PCF) as described
    ! by Carrayrou (2017) when solids are present.
    implicit none
    integer,          intent(in) :: NCP                        ! number of species
    integer,          intent(in) :: NEQ                        ! number of equilibria
    integer,          intent(in) :: NSL                        ! number of solids
    double precision, intent(in) :: THRESH                     ! acceptable threshold
    double precision, intent(inout), dimension(NCP+NEQ)  :: C  ! the free concentrations
    double precision, intent(inout), dimension(NSL)      :: CS  ! the free concentrations
    double precision, intent(in),    dimension(NCP)      :: T  ! the analytical concentration
    double precision, intent(in),    dimension(NEQ)      :: B  ! equilibrium constants
    integer,          intent(in),    dimension(NEQ, NCP) :: P  ! stoichiometry matrix
    double precision, intent(in),    dimension(NSL)      :: KS  ! equilibrium constants
    integer,          intent(in),    dimension(NSL, NCP) :: Q  ! solid stoichiometry matrix
  
    double precision, dimension(NCP) :: expn, convaux
    integer :: i, j, iter
    double precision :: theta, auxd
    double precision            :: sump, sumr, calc_theta, check_convergence
    integer,          parameter :: MAXITS = 500
    double precision, DIMENSION(NSL) :: SI
    integer :: current_solid, comp_control
    double precision, parameter :: SI_THRESH = 5.0
  
    call find_exponent(P, NCP, NEQ, expn)
    current_solid = -1         ! no current solid
    comp_control = -1          ! no component controlling

    do iter=1,MAXITS
        ! I. check whether there is any solid to precipitate
        call saturation_product(C, Q, NCP, NSL, SI)
        SI = SI / KS
        if (any(SI > SI_THRESH)) then
            current_solid = maxloc(SI, dim=1)    ! tag the solid to be precipitated

            ! choose component that will control this solid
            comp_control = maxloc(Q(current_solid,:), dim=1)
        else
            ! no solids to precipitate
            if (current_solid > 0) then
                if (CS(current_solid) < 0.0) CS(current_solid) = 0.0
            end if
            current_solid = -1         ! no current solid
            comp_control = -1          ! no component controlling
        end if
        
        !print *, iter, C(1), C(2), CS, SI

        call cpluss(C, B, P, NCP, NEQ)
        ! calculate sump and sumr
        do i=1,NCP
            if (i .eq. comp_control) then
                sump = 1.0
                sumr = SI(current_solid)**(1.0/Q(current_solid, comp_control))
            else   
                !call calc_sum(i, T, C, P, sump, sumr, NCP, NEQ)
                call calc_sum_solids(i, T, C, CS, P, Q, sump, sumr, NCP, NEQ, NSL)
            end if
            convaux(i) = check_convergence(sumr, sump)

            ! update concentrations
            theta = calc_theta(sump, sumr)
            if (i .eq. comp_control) then
                auxd = C(comp_control)
                auxd = auxd + sum(P(:, comp_control)*C(NCP+1:NCP+NEQ))
                do j=1,NSL
                    if (j .eq. current_solid) cycle
                    auxd = auxd + Q(j, comp_control)*CS(j)
                end do
                CS(current_solid) = (T(comp_control) - auxd)/Q(current_solid, comp_control)
            end if
            C(i) = theta * C(i) * (sump/sumr)**expn(i) + (1-theta)*C(i)
        end do
        !print *, iter, SI, CS, current_solid, comp_control, convaux
        if (maxval(convaux) .lt. THRESH) exit  ! converged
    end do
    !print *, "END"
    !print "(A, 3E14.6)", "refined values: ", C(1:2), CS
    !print *, "should be      : 1.258925E-15	1.584893E-07  9.638620E-09"
end subroutine pcf_solids


subroutine find_exponent(P, NCP, NEQ, expn)
    ! find the smallest strictly positive
    ! stoichiometric coefficient that
    ! the jth component presents in Morel's Table. 

    implicit none
    integer,          intent(in) :: NCP                        ! number of species
    integer,          intent(in) :: NEQ                        ! number of equilibria
    integer,          intent(in),    dimension(NEQ, NCP) :: P  ! stoichiometry matrix
    double precision, intent(out),  dimension(NCP) :: expn

    integer, dimension(NCP) :: aux

    aux = maxval(abs(P), dim=1)
    expn = 1.0/aux
end subroutine find_exponent


subroutine calc_sum(comp, T, C, P, sump, sumr, NCP, NEQ)
    implicit none
    integer,          intent(in) :: NCP                        ! number of species
    integer,          intent(in) :: NEQ                        ! number of equilibria
    integer, intent(in) :: comp
    double precision, intent(out) :: sump, sumr
    DOUBLE PRECISION, INTENT(IN), DIMENSION(NCP+NEQ) :: C
    integer,          intent(in),    dimension(NEQ, NCP) :: P  ! stoichiometry matrix
    double precision, intent(in),    dimension(NCP)      :: T  ! the analytical concentration

    integer j
    if (T(comp) .ge. 0.0) then
      sump = T(comp)
      sumr = C(comp)
    else
      sump = 0.0
      sumr = abs(T(comp)) + C(comp)
    end if
    do j=1,NEQ
      if (P(j,comp) > 0) then
        sumr = sumr + P(j,comp) * C(j+NCP)
      else
        sump = sump + abs(P(j,comp)) * C(j+NCP)
      end if
    end do
end subroutine calc_sum


subroutine calc_sum_solids(comp, T, C, CS, P, Q, sump, sumr, NCP, NEQ, NSL)
    implicit none
    integer,          intent(in) :: NCP                        ! number of species
    integer,          intent(in) :: NEQ                        ! number of equilibria
    integer,          intent(in) :: NSL                        ! number of solids
    integer, intent(in) :: comp
    double precision, intent(out) :: sump, sumr
    DOUBLE PRECISION, INTENT(IN), DIMENSION(NCP+NEQ) :: C
    DOUBLE PRECISION, INTENT(IN), DIMENSION(NSL) :: CS
    integer,          intent(in),    dimension(NEQ, NCP) :: P  ! stoichiometry matrix
    integer,          intent(in),    dimension(NSL, NCP) :: Q  ! stoichiometry matrix
    double precision, intent(in),    dimension(NCP)      :: T  ! the analytical concentration

    integer j
    if (T(comp) .ge. 0.0) then
      sump = T(comp)
      sumr = C(comp)
    else
      sump = 0.0
      sumr = abs(T(comp)) + C(comp)
    end if
    do j=1,NEQ
      if (P(j,comp) > 0) then
        sumr = sumr + P(j,comp) * C(j+NCP)
      else
        sump = sump + abs(P(j,comp)) * C(j+NCP)
      end if
    end do
    do j=1,NSL
      if (CS(j) < 0.0) cycle        ! CS might go negative; do not include in balance or else everythin breaks
      if (Q(j,comp) > 0) then
        sumr = sumr + Q(j,comp) * CS(j)
      else
        sump = sump + abs(Q(j,comp)) * CS(j)
      end if
  end do
end subroutine calc_sum_solids

function calc_theta(sumr, sump) result(theta)
    double precision, intent(in)       :: sump, sumr
    double precision :: theta
    double precision, parameter :: THETA_ALPHA = 0.9
    double precision, parameter :: THETA_BETA = 0.8
    if (sumr > sump) then
        theta = THETA_ALPHA - THETA_BETA * sump/sumr
    else
        theta = THETA_ALPHA - THETA_BETA * sumr/sump
    end if
end function calc_theta


function check_convergence(sumr, sump) result(conv)
    DOUBLE PRECISION :: conv
    double precision, intent(in)       :: sump, sumr
    conv = abs(sump-sumr)
    if ((sumr + sump) > 0.0) then    ! if one of the components is zero the division
        conv = conv/(sumr+sump)      ! will be 1 thus ruining the value of convergence testing
    end if                           ! if the sum of both is zero, also problem
end function check_convergence
