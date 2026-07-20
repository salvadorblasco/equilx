! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>
!
! Compile with $ f2py -llapack -c -m solver cpluss.o pcf.o mb.o eq.f90
!
! control codes
! 0  -> normal exit, no errors
! 1  -> normal exit upon initial guess being a root
! 2  -> abnormal exit on MAXITS
! 3  -> too many iterations on solids
! 99 -> unable to invert matrix
! 98 -> solids could not be solved
! 97 -> invalid reference component for simulation was provided


! ====================================================================================================
! PUBLIC PYTHON MODULE INTERFACE
! ====================================================================================================


subroutine solve(C, T, B, P, CTRL, NPT, NEQ, NCP)
  ! Solvereal equilibrium. The first initial value must be provided in the first row of C. The following 
  ! points take as starting value the value of the last iteration.
  implicit none
  integer, intent(in) :: NPT     ! number of titration points
  integer, intent(in) :: NCP     ! number of components
  integer, intent(in) :: NEQ     ! number of equilibria
  integer, intent(inout) :: CTRL ! control value. 0 on no error, 1 otherwise
  integer, dimension(NEQ, NCP), intent(in) :: P                 ! The stoichiometry array
  double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C ! initial estimates (in) and refined values (out)
  double precision, dimension(NPT, NCP),     intent(in) :: T    ! analytical concentration values
  double precision, dimension(NEQ),          intent(in) :: B    ! the equilibrium constants
  !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  integer i
  CTRL = 0
  do i=1,NPT
    if (i == 1) then
      call solve_point(C(i,:), T(i,:), B, P, CTRL, .true., NEQ, NCP)
      return ! DEBUG
    else
      C(i,:) = C(i-1,:)
      call solve_point(C(i,:), T(i,:), B, P, CTRL, .false., NEQ, NCP)
    end if
    if (CTRL > 50) exit  ! error unable to continue
  end do
end subroutine solve


subroutine solve_solids(C, CS, T, B, P, KS, Q, CTRL, NPT, NEQ, NCP, NSL)
  ! Solve equilibrium when solids are present. The first initial value must be provided
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
  double precision, dimension(NPT, NSL), intent(inout) :: CS     ! solid concentration initial estimates (in) and refined values (out)
  double precision, dimension(NPT, NCP), intent(in) :: T        ! analytical concentration values
  double precision, dimension(NEQ), intent(in) :: B  ! the equilibrium constants
  double precision, dimension(NSL), intent(in) :: KS ! the solubility constant
  !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  integer :: i
  CTRL = 0
  do i=1,NPT
    if (i == 1) then
      call solve_point_solid(C(i,:), CS(i,:), T(i,:), B, P, KS, Q, CTRL, .true., NEQ, NCP, NSL)
    else
      C(i,:) = C(i-1,:)
      CS(i,:) = CS(i-1,:)
      call solve_point_solid(C(i,:), CS(i,:), T(i,:), B, P, KS, Q, CTRL, .false., NEQ, NCP, NSL)
    end if
    if (CTRL > 50) exit  ! error unable to continue
  end do
end subroutine solve_solids


subroutine solve_distribution(C, T, B, P, CTRL, REF, NPT, NEQ, NCP)
    implicit none
    integer, intent(in) :: NPT     ! number of titration points
    integer, intent(in) :: NCP     ! number of components
    integer, intent(in) :: NEQ     ! number of equilibria
    integer, intent(in) :: REF     ! component of reference
    integer, intent(inout) :: CTRL ! control value. 0 on no error, 1 otherwise
    integer, dimension(NEQ, NCP), intent(in) :: P ! The stoichiometry array
    double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C ! initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NCP), intent(in) :: T ! analytical concentration values
    double precision, dimension(NEQ), intent(in) :: B  ! the equilibrium constants
    !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
    !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
    !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
  
    double precision, dimension(NPT, NEQ)       :: frozenB  ! the frozen equilibrium constants
    integer,          dimension(NEQ, NCP-1)     :: frozenP ! The frozen stoichiometry array
    double precision, dimension(NPT, NEQ+NCP-1) :: frozenC ! frozen initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NCP-1)     :: frozenT ! frozen analytical concentration values
    integer i, j
    CTRL = 0
  
    call freeze_component(B, frozenB, T, frozenT, P, frozenP, REF, NPT, NEQ, NCP)
  
    do i=1,NPT
        if (i == 1) then
            do j=1,NCP     ! copy first point to frozenC
                if (j < REF) frozenC(i,j) = C(i,j)
                if (j > REF) frozenC(i,j-1) = C(i,j)
            end do
          call solve_point(frozenC(i,:), frozenT(i,:), frozenB(i,:), frozenP, CTRL, .true., NEQ, NCP)
      else
          frozenC(i,:) = frozenC(i-1,:)
          call solve_point(frozenC(i,:), frozenT(i,:), frozenB(i,:), frozenP, CTRL, .false., NEQ, NCP)
      end if
      if (CTRL > 50) exit  ! error unable to continue
      do j=1,NCP+NEQ
        if (j < REF) C(i,j) = frozenC(i,j)
        if (j > REF) C(i,j) = frozenC(i,j-1)
        if (j == REF) C(i,REF) = T(i,REF)
      end do
    end do
end subroutine solve_distribution


subroutine solve_distribution_solids(C, CS, T, B, P, KS, Q, CTRL, REF, NPT, NEQ, NCP, NSL)
    ! Perform Newton-Raphson iterations per point when the beta matrix
    ! is 1D and solids are present. The first initial value must be provided
    ! in the first row of C. The following points take as starting value the value of 
    ! the last iteration.

    implicit none
    integer, intent(in) :: NPT     ! number of titration points
    integer, intent(in) :: NCP     ! number of components
    integer, intent(in) :: NEQ     ! number of equilibria
    integer, intent(in) :: NSL     ! number of solids
    integer, intent(in) :: REF     ! component of reference
    integer, intent(inout) :: CTRL ! control value. 0 on no error, 1 otherwise
    integer, dimension(NEQ, NCP), intent(in) :: P ! The soluble stoichiometry array
    integer, dimension(NSL, NCP), intent(in) :: Q ! The solids stoichiometry array
    double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C ! free concentation initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NSL), intent(inout) :: CS     ! solid concentration initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NCP), intent(in) :: T        ! analytical concentration values
    double precision, dimension(NEQ), intent(in) :: B  ! the equilibrium constants
    double precision, dimension(NSL), intent(in) :: KS ! the solubility constant
    !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
    !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
    !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
  
    double precision, dimension(NPT, NEQ)       :: fB  ! the frozen equilibrium constants
    double precision, dimension(NPT, NSL)       :: fKS ! the frozen equilibrium constants
    integer,          dimension(NSL, NCP-1)     :: fQ  ! the frozen equilibrium constants
    integer,          dimension(NEQ, NCP-1)     :: fP  ! The frozen stoichiometry array
    double precision, dimension(NPT, NEQ+NCP-1) :: fC  ! frozen initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NCP-1)     :: fT  ! frozen analytical concentration values
    integer i, j
    CTRL = 0
  
    ! Freeze component
    call freeze_component(B, fB, T, fT, P, fP, REF, NPT, NEQ, NCP)
    call freeze_component_solid(T(:,REF), Q, fQ, KS, fKS, REF, NPT, NCP, NSL)
  
    ! solve point-wise
    do i=1,NPT
        if (i == 1) then
            do j=1,NCP     ! copy first point to frozenC
                if (j < REF) fC(i,j) = C(i,j)
                if (j > REF) fC(i,j-1) = C(i,j)
            end do
          call solve_point_solid(fC(i,:), CS(i,:), fT(i,:), fB(i,:), fP, fKS(i,:), fQ, CTRL, .true., NEQ, NCP, NSL)
      else
          fC(i,:) = fC(i-1,:)
          call solve_point_solid(fC(i,:), CS(i,:), fT(i,:), fB(i,:), fP, fKS(i,:), fQ, CTRL, .false., NEQ, NCP, NSL)
      end if
      if (CTRL > 50) exit  ! error unable to continue
      do j=1,NCP+NEQ
        if (j < REF) C(i,j) = fC(i,j)
        if (j > REF) C(i,j) = fC(i,j-1)
        if (j == REF) C(i,REF) = T(i,REF)
      end do
    end do
end subroutine solve_distribution_solids

! ====================================================================================================
! END PUBLIC PYTHON MODULE INTERFACE
! ====================================================================================================


subroutine freeze_component(B, fB, T, fT, P, fP, REF, NPT, NEQ, NCP)
    integer, intent(in) :: NPT     ! number of titration points
    integer, intent(in) :: NCP     ! number of components
    integer, intent(in) :: NEQ     ! number of equilibria
    integer, intent(in) :: REF     ! component of reference
    integer, dimension(NEQ, NCP), intent(in) :: P ! The stoichiometry array
    double precision, dimension(NPT, NCP), intent(in) :: T ! analytical concentration values
    double precision, dimension(NEQ), intent(in) :: B  ! the equilibrium constants
    double precision, intent(inout), dimension(NPT, NEQ) :: fB  ! the frozen equilibrium constants
    integer, INTENT(INOUT), dimension(NEQ, NCP-1) :: fP ! The stoichiometry array
    double precision,INTENT(INOUT),  dimension(NPT, NCP-1):: fT ! analytical concentration values

    INTEGER i, j
    do i=1,NCP
        if (i<REF) fP(:,i) = P(:,i)
        if (i>REF) fP(:,i-1) = P(:,i)
    end do
    do i=1,NPT
        do j=1,NEQ
            fB(i,j) = B(j)*T(i,REF)**P(j,REF)
        end do
        do j=1,NCP
            if (j<REF) fT(i,j) = T(i,j)
            if (j>REF) fT(i,j-1) = T(i,j)
        end do
    end do
end subroutine freeze_component


subroutine freeze_component_solid(CR, Q, fQ, KS, fKS, REF, NPT, NCP, NSL)
    integer, intent(in) :: NPT     ! number of titration points
    integer, intent(in) :: NSL     ! number of solids
    integer, intent(in) :: REF     ! component of reference
    integer, intent(in) :: NCP     ! number of components
    double precision, dimension(NPT), intent(in) :: CR ! concentration of the reference component
    double precision, INTENT(IN), dimension(NSL) :: KS ! solubility constants
    double precision, INTENT(INOUT), dimension(NPT,NSL) :: fKS ! frozen solubility constants
    integer, dimension(NSL, NCP), intent(in) :: Q ! The solids stoichiometry array
    integer, dimension(NSL, NCP-1), intent(inout) :: fQ ! The solids stoichiometry array

    INTEGER i, j

    do i=1,NCP
        if (i<REF) fQ(:,i) = Q(:,i)
        if (i>REF) fQ(:,i-1) = Q(:,i)
    end do
    do i=1,NPT
        do j=1,NSL
            fKS(i,j) = KS(j)*CR(i)**Q(j, REF)
        end do
    end do
end subroutine freeze_component_solid


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

  DOUBLE PRECISION, parameter :: PCF_TOLERANCE = 0.25
  DOUBLE PRECISION, parameter :: PCF_TOLERANCE_FIRST = 1D-9
  CTRL = 0

  print *, " DEBUG 1 > ", c
  if (is_first) then
    call pcf(C, T, B, P, PCF_TOLERANCE_FIRST, NCP, NEQ)
  else
    call pcf(C, T, B, P, PCF_TOLERANCE, NCP, NEQ)
  end if
  print *, " DEBUG 2 > ", c
  call nr(C, T, B, P, CTRL, NEQ, NCP)
  print *, " DEBUG 2 > ", c
end subroutine solve_point


subroutine solve_point_solid(C, CS, T, B, P, KS, Q, CTRL, is_first, NEQ, NCP, NSL)
  !f2py intent(in) :: NEQ, NCP, NSL
  !f2py intent(in) :: T, B, P, KS, Q
  !f2py double precision, dimension(NEQ+NCP), intent(inout) :: C
  !f2py intent(inout) :: CS, CTRL
  !f2py logical intent(in) :: is_first
  !f2py integer, intent(hide), depend(P) :: ncp = shape(P, 1)
  !f2py integer, intent(hide), depend(P) :: neq = shape(P, 0)
  !f2py integer, intent(hide), depend(KS) :: nsl = size(KS)

  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 1D. The first initial value must be provided in the first row
  ! of C. The following points take as starting value the value of 
  ! the last iteration.
  implicit none
  integer,         intent(in) :: NCP                               ! number of components
  integer,         intent(in) :: NEQ                               ! number of equilibria
  integer,         intent(in) :: NSL                               ! number of solids
  integer,         intent(inout) :: CTRL                           ! control value. 0 on no error, 1 otherwise
  logical,         intent(in) :: is_first                          ! whther it is first point or not
  integer,         dimension(NEQ, NCP), intent(in) :: P            ! The stoichiometry array
  double precision, dimension(NEQ+NCP), intent(inout) :: C ! initial estimates (in) and refined values (out)
  double precision, dimension(NSL), intent(inout) :: CS    ! solid concentration initial estimates (in) and refined values (out)
  double precision, dimension(NCP), intent(in) :: T        ! analytical concentration values
  double precision, dimension(NEQ), intent(in) :: B        ! the equilibrium constants
  integer, dimension(NSL, NCP), intent(in) :: Q   ! the equilibrium constants
  double precision, dimension(NSL), intent(in) :: KS       ! the precipitation constants

  DOUBLE PRECISION, PARAMETER :: si_tolerance = 1e-6
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
    call nr_solids(C, CS, T, B, P, KS, Q, ctrl, do_solid, NEQ, NCP, NSL)
  end do
end subroutine solve_point_solid
