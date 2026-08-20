! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>
!
! Compile with $ f2py -llapack -c -m solver cpluss.o pcf.o mb.o eq.f90
! 
! Simple titration solving:
!
!                 activity=1.0   activity!=1.0
! no solids          solve         solve_vis
! with solids     solve_solids   solve_solids_vis
!
! control codes
! 0  -> normal exit, no errors
! 1  -> normal exit upon initial guess being a root
! 2  -> abnormal exit on MAXITS
! 3  -> too many iterations on solids
! 4  -> outer point routine exited on MAXITS
! 5  -> solid solver exited on MAXITS
! 6  -> the solid solver found a concentration equal to zero (caller must remove solid and recall)
! 99 -> unable to invert matrix
! 98 -> solids could not be solved
! 97 -> invalid reference component for simulation was provided


! ====================================================================================================
! PUBLIC PYTHON MODULE INTERFACE
! ====================================================================================================


subroutine solve(C, T, B, P, CTRL, NPT, NCP, NEQ)
    ! Solve equilibrium titration when no solids are present and ionic strengh is constant.
    ! The first initial value must be provided in C(1,1:NCP). The following 
    ! points take as starting value the value of the last iteration.
    implicit none
    integer, intent(in) :: NPT                                    ! number of titration points
    integer, intent(in) :: NCP                                    ! number of components
    integer, intent(in) :: NEQ                                    ! number of equilibria
    integer, intent(inout) :: CTRL                                ! control value. See return values in header of this file.
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
        if (i .eq. 1) then
            call solve_point(C(1,:), T(1,:), B, P, CTRL, .true., NCP, NEQ)
        else
            C(i,:) = C(i-1,:)
            call solve_point(C(i,:), T(i,:), B, P, CTRL, .false., NCP, NEQ)
        end if
        if (CTRL > 50) exit  ! error unable to continue
    end do
end subroutine solve


subroutine solve_vis(C, T, Bref, Iref, Iback, P, charge, DHA, DHB, DHP, CTRL, NPT, NCP, NEQ)
    ! Solve equilibrium titration when no solids are present and ionic strengh is variable.
    ! The first initial value must be provided in C(1,1:NCP). The following 
    ! points take as starting value the value of the last iteration.
    implicit none
    integer,          intent(in)    :: NPT                               ! number of titration points
    integer,          intent(in)    :: NCP                               ! number of components
    integer,          intent(in)    :: NEQ                               ! number of equilibria
    integer,          intent(inout) :: CTRL                              ! control value. See return values in header of this file.
    double precision, intent(inout), dimension(NPT, NEQ+NCP) :: C        ! initial estimates (in) and refined values (out)
    double precision, intent(in),    dimension(NPT, NCP)     :: T        ! analytical concentration values
    double precision, intent(in),    dimension(NEQ)          :: Bref     ! the equilibrium constants
    DOUBLE PRECISION, INTENT(IN)                             :: Iref     ! reference ionic strength
    DOUBLE PRECISION, INTENT(IN),    dimension(NPT)          :: Iback    ! background ionic strength
    integer,          intent(in),    dimension(NEQ, NCP)     :: P        ! The stoichiometry array
    integer,          INTENT(IN),    DIMENSION(NCP+NEQ)      :: charge   ! charge for every species
    DOUBLE PRECISION, INTENT(IN)                             :: DHA, DHB ! Debye-Huckel parameters A and B
    double precision, intent(in),    dimension(NEQ,3)        :: DHP      ! per-species DH parameters C, D and E
  
    integer i
    DOUBLE PRECISION, dimension(NEQ) :: corr   ! the activity correction to the equilibrium constant
    integer,          dimension(NEQ) :: zstar  ! z* = sum p_i z_i^2 - z_s^2 
  
    do i=1,NEQ
        zstar(i) = -charge(NCP+i)**2 + sum(P(i,:)*charge(1:NCP)**2)
    end do

    corr = 0.0 
    CTRL = 0
    do i=1,NPT
        if (i .eq. 1) then
            call solve_point(C(1,:), T(1,:), Bref+corr, P, CTRL, .true., NCP, NEQ)
        else
            C(i,:) = C(i-1,:)
            call solve_point(C(i,:), T(i,:), Bref+corr, P, CTRL, .false., NCP, NEQ)
        end if
        if (CTRL > 50) exit  ! error unable to continue
        call outer(corr, C(i,:), charge, Bref, Iref, Iback(i), zstar, DHA, DHB, DHP, T(i,:), P, CTRL, NCP, NEQ)
        if (CTRL > 50) exit  ! error unable to continue
        call nr(C(i,:), T(i,:), Bref+corr, P, ctrl, NCP, NEQ)
    end do
end subroutine solve_vis


subroutine solve_solids(C, CS, T, B, P, KS, Q, CTRL, NPT, NCP, NEQ, NSL)
  ! Solve titration equilibrium when solids are present. The first initial value must be provided
  ! in the first row of C. The following points take as starting value the value of 
  ! the last iteration.
  implicit none
  integer, intent(in) :: NPT                                      ! number of titration points
  integer, intent(in) :: NCP                                      ! number of components
  integer, intent(in) :: NEQ                                      ! number of equilibria
  integer, intent(in) :: NSL                                      ! number of solids
  integer, intent(inout) :: CTRL                                  ! control value. 0 on no error, 1 otherwise
  integer, dimension(NEQ, NCP), intent(in) :: P                   ! The soluble stoichiometry array
  integer, dimension(NSL, NCP), intent(in) :: Q                   ! The solids stoichiometry array
  double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C   ! free concentation initial estimates (in) and refined values (out)
  double precision, dimension(NPT, NSL), intent(inout) :: CS      ! solid concentration initial estimates (in) and refined values (out)
  double precision, dimension(NPT, NCP), intent(in) :: T          ! analytical concentration values
  double precision, dimension(NEQ), intent(in) :: B               ! the equilibrium constants
  double precision, dimension(NSL), intent(in) :: KS              ! the solubility constant
  !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
  !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
  !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)

  integer :: i

  CTRL = 0
  do i=1,NPT
    if (i .eq. 1) then
        call solve_point_solid(C(1,:), CS(1,:), T(1,:), B, P, KS, Q, CTRL, .true., NEQ, NCP, NSL)
    else
        C(i,:) = C(i-1,:)
        CS(i,:) = CS(i-1,:)
        call solve_point_solid(C(i,:), CS(i,:), T(i,:), B, P, KS, Q, CTRL, .false., NEQ, NCP, NSL)
    end if
    if (CTRL > 50) exit  ! error unable to continue
  end do
end subroutine solve_solids


subroutine solve_solids_vis(C, CS, T, Bref, Iref, Iback, P, KS, Q, charge, DHA, DHB, DHP, CTRL, NPT, NCP, NEQ, NSL)
    ! Solve titration equilibrium when solids are present. The first initial value must be provided
    ! in the first row of C. The following points take as starting value the value of 
    ! the last iteration.
    implicit none
    integer, intent(in) :: NPT                                      ! number of titration points
    integer, intent(in) :: NCP                                      ! number of components
    integer, intent(in) :: NEQ                                      ! number of equilibria
    integer, intent(in) :: NSL                                      ! number of solids
    integer, intent(inout) :: CTRL                                  ! control value. 0 on no error, 1 otherwise
    integer, dimension(NEQ, NCP), intent(in) :: P                   ! The soluble stoichiometry array
    integer, dimension(NSL, NCP), intent(in) :: Q                   ! The solids stoichiometry array
    double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C   ! free concentation initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NSL), intent(inout) :: CS      ! solid concentration initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NCP), intent(in) :: T          ! analytical concentration values
    double precision, dimension(NEQ), intent(in) :: Bref            ! the equilibrium constants at I=Iref
    double precision, dimension(NSL), intent(in) :: KS              ! the solubility constant
    DOUBLE PRECISION, INTENT(IN)                             :: Iref     ! reference ionic strength
    DOUBLE PRECISION, INTENT(IN),    dimension(NPT)          :: Iback    ! background ionic strength
    integer,          INTENT(IN),    DIMENSION(NCP+NEQ)      :: charge   ! charge for every species
    DOUBLE PRECISION, INTENT(IN)                             :: DHA, DHB ! Debye-Huckel parameters A and B
    double precision, intent(in),    dimension(NEQ,3)        :: DHP      ! per-species DH parameters C, D and E
    !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
    !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
    !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
    !f2py integer, intent(hide), depend(q) :: nsl = shape(q,0)

    integer :: i
    DOUBLE PRECISION, dimension(NEQ) :: corr   ! the activity correction to the equilibrium constant
    integer,          dimension(NEQ) :: zstar  ! z* = sum p_i z_i^2 - z_s^2 
    do i=1,NEQ
        zstar(i) = -charge(NCP+i)**2 + sum(P(i,:)*charge(1:NCP)**2)
    end do

    corr = 0.0 
    CTRL = 0
    do i=1,NPT
      if (i .eq. 1) then
          call solve_point_solid(C(1,:), CS(1,:), T(1,:), Bref, P, KS, Q, CTRL, .true., NEQ, NCP, NSL)
      else
          C(i,:) = C(i-1,:)
          CS(i,:) = CS(i-1,:)
          call solve_point_solid(C(i,:), CS(i,:), T(i,:), Bref, P, KS, Q, CTRL, .false., NEQ, NCP, NSL)
          call outer_solid(corr, C(i,:), charge, Bref, Iref, Iback(i), zstar, DHA, DHB, DHP, T(i,:), P, CTRL, NCP, NEQ)
      end if
      if (CTRL > 50) exit  ! error unable to continue
    end do
end subroutine solve_solids_vis


subroutine solve_distribution(C, T, B, P, CTRL, REF, NPT, NCP, NEQ)
    ! Solve the equilibrium when one of the components has a fixed concentration, therefore said
    ! concentration is a given, not an unknown. The index of the fixed component is given in the
    ! variable REF, and the concentration is in C(:,REF). The contents of T(:,REF) is replaced
    ! with the corresponding value.
    implicit none
    integer, intent(in) :: NPT                                     ! number of titration points
    integer, intent(in) :: NCP                                     ! number of components
    integer, intent(in) :: NEQ                                     ! number of equilibria
    integer, intent(in) :: REF                                     ! component of reference
    integer, intent(inout) :: CTRL                                 ! control value. 0 on no error, 1 otherwise
    integer, dimension(NEQ, NCP), intent(in) :: P                  ! The stoichiometry array
    double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C  ! initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NCP), intent(inout) :: T      ! analytical concentration values
    double precision, dimension(NEQ), intent(in) :: B              ! the equilibrium constants
    !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
    !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
    !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
  
    double precision :: Cref
    double precision, dimension(NEQ)        :: fB  ! the frozen equilibrium constants
    integer,          dimension(NEQ, NCP-1) :: fP  ! The frozen stoichiometry array
    double precision, dimension(NEQ+NCP-1)  :: fC  ! frozen initial estimates (in) and refined values (out)
    double precision, dimension(NCP-1)      :: fT  ! frozen analytical concentration values
    integer i, j
    CTRL = 0
  
    call freeze_stoich(P, fP, REF, NCP, NEQ)

    do i=1,NPT
        Cref = C(i,REF)
        call freeze_component(Cref, B, fB, T(i,:), fT, P, REF, NCP, NEQ)
        if (i .eq. 1) then     ! first point
            do j=1,NCP         ! copy initial guess to frozenC
                if (j < REF) fC(j) = C(1,j)
                if (j > REF) fC(j-1) = C(1,j)
            end do
            call solve_point(fC, fT, fB, fP, CTRL, .true., NCP-1, NEQ)
        else
            call solve_point(fC, fT, fB, fP, CTRL, .false., NCP-1, NEQ)
        end if
        if (CTRL > 50) exit  ! error unable to continue
        do j=1,NCP+NEQ
            if (j < REF) C(i,j) = fC(j)
            if (j > REF) C(i,j) = fC(j-1)
        end do
        T(i,REF) = Cref + sum(P(:,REF)*C(i,NCP+1:NCP+NEQ))
    end do
end subroutine solve_distribution


subroutine solve_distribution_solids(C, CS, T, B, P, KS, Q, CTRL, REF, NPT, NEQ, NCP, NSL)
    ! Perform Newton-Raphson iterations per point when the beta matrix
    ! is 1D and solids are present. The first initial value must be provided
    ! in the first row of C. The following points take as starting value the value of 
    ! the last iteration.

    implicit none
    integer, intent(in) :: NPT                                    ! number of titration points
    integer, intent(in) :: NCP                                    ! number of components
    integer, intent(in) :: NEQ                                    ! number of equilibria
    integer, intent(in) :: NSL                                    ! number of solids
    integer, intent(in) :: REF                                    ! component of reference
    integer, intent(inout) :: CTRL                                ! control value. 0 on no error, 1 otherwise
    integer, dimension(NEQ, NCP), intent(in) :: P                 ! The soluble stoichiometry array
    integer, dimension(NSL, NCP), intent(in) :: Q                 ! The solids stoichiometry array
    double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C ! free concentation initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NSL), intent(inout) :: CS    ! solid concentration initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NCP), intent(inout) :: T     ! analytical concentration values
    double precision, dimension(NEQ), intent(in) :: B             ! the equilibrium constants
    double precision, dimension(NSL), intent(in) :: KS            ! the solubility constant
    !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
    !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
    !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
    !f2py integer, intent(hide), depend(q) :: nsl = shape(q,0)
  
    double precision :: Cref
    double precision, dimension(NEQ)        :: fB  ! the frozen equilibrium constants
    double precision, dimension(NSL)        :: fKS ! the frozen equilibrium constants
    integer,          dimension(NSL, NCP-1) :: fQ  ! the frozen equilibrium constants
    integer,          dimension(NEQ, NCP-1) :: fP  ! The frozen stoichiometry array
    double precision, dimension(NEQ+NCP-1)  :: fC  ! frozen initial estimates (in) and refined values (out)
    double precision, dimension(NCP-1)      :: fT  ! frozen analytical concentration values
    integer i, j

    CTRL = 0
    call freeze_stoich(P, fP, REF, NCP, NEQ)
    call freeze_stoich(Q, fQ, REF, NCP, NSL)
  
    ! solve point-wise
    point_loop: do i=1,NPT
        Cref = C(i,REF)
        ! Freeze component
        call freeze_component(Cref, B, fB, T(i,:), fT, P, REF, NCP, NEQ)
        call freeze_component_solid(Cref, Q, KS, fKS, REF, NCP, NSL)
  
        ! Solve point
        if (i == 1) then
            do j=1,NCP     ! copy first point to frozenC
                if (j < REF) fC(j) = C(1,j)
                if (j > REF) fC(j-1) = C(1,j)
            end do
            call solve_point_solid(fC, CS(1,:), fT, fB, fP, fKS, fQ, CTRL, .true., NEQ, NCP-1, NSL)
        else
            CS(i,:) = CS(i-1,:)
            call solve_point_solid(fC, CS(i,:), fT, fB, fP, fKS, fQ, CTRL, .false., NEQ, NCP-1, NSL)
        end if
        if (CTRL > 50) exit point_loop ! error unable to continue
        do j=1,NCP+NEQ
          if (j < REF) C(i,j) = fC(j)
          if (j > REF) C(i,j) = fC(j-1)
        end do
        T(i,REF) = Cref + sum(P(:,REF)*C(i,NCP+1:NCP+NEQ)) + sum(Q(:,REF)*CS(i,:))
    end do point_loop
end subroutine solve_distribution_solids

! ====================================================================================================
! END PUBLIC PYTHON MODULE INTERFACE
! ====================================================================================================


subroutine freeze_stoich(P, fP, REF, NCP, NEQ)
    ! freeze stoichiometry array, either soluble or solid
    integer, intent(in) :: NCP                         ! number of components
    integer, intent(in) :: NEQ                         ! number of equilibria
    integer, intent(in) :: REF                         ! component of reference
    integer, INTENT(OUT), dimension(NEQ, NCP-1) :: fP  ! The stoichiometry array
    integer, intent(in),  dimension(NEQ, NCP)   :: P   ! The stoichiometry array
    do i=1,NCP                         ! freeze stoichiometry array
        if (i<REF) fP(:,i) = P(:,i)
        if (i>REF) fP(:,i-1) = P(:,i)
    end do
end subroutine freeze_stoich


subroutine freeze_component(Cref, B, fB, T, fT, P, REF, NCP, NEQ)
    ! given a reference component, freeze all equilibrium parameters relative to it.
    ! The reference index is given in variable REF; the concentration to be frozen
    ! must be provided in Cref.
    integer, intent(in) :: NCP                                ! number of components
    integer, intent(in) :: NEQ                                ! number of equilibria
    integer, intent(in) :: REF                                ! component of reference
    double precision, intent(in) :: Cref                      ! the reference component concentration
    integer,          intent(in),  dimension(NEQ, NCP) :: P   ! The stoichiometry array
    double precision, intent(in),  dimension(NCP)      :: T   ! analytical concentration values
    double precision, intent(in),  dimension(NEQ)      :: B   ! the equilibrium constants
    double precision, intent(out), dimension(NEQ)      :: fB  ! the frozen equilibrium constants
    double precision, intent(out), dimension(NCP-1)    :: ft  ! analytical concentration values

    integer j

    do j=1,NEQ
        fB(j) = B(j)*Cref**P(j,REF)
    end do
    do j=1,NCP
        if (j<REF) fT(j) = T(j)
        if (j>REF) fT(j-1) = T(j)
    end do
end subroutine freeze_component


subroutine freeze_component_solid(Cref, Q, KS, fKS, REF, NCP, NSL)
    ! given a reference component, freeze all equilibrium parameters relative to it.
    ! The reference index is given in variable REF; the concentration to be frozen
    ! must be provided in Cref.
    integer,          intent(in)                       :: NSL  ! number of solids
    integer,          intent(in)                       :: REF  ! component of reference
    integer,          intent(in)                       :: NCP  ! number of components
    double precision, intent(in)                       :: Cref ! concentration of the reference component
    double precision, INTENT(IN),  dimension(NSL)      :: KS   ! solubility constants
    double precision, INTENT(OUT), dimension(NSL)      :: fKS  ! frozen solubility constants
    integer,          intent(in),  dimension(NSL, NCP) :: Q    ! The solids stoichiometry array

    INTEGER j

    do j=1,NSL
        fKS(j) = KS(j)/(Cref**Q(j, REF))
    end do
end subroutine freeze_component_solid


subroutine solve_point(C, T, B, P, CTRL, is_first, NCP, NEQ)
    ! Perform Newton-Raphson iterations per point when the beta matrix
    ! is 1D. The first initial value must be provided in the first row
    ! of C. The following points take as starting value the value of 
    ! the last iteration.
    implicit none
    integer, intent(in)    :: NCP        ! number of components
    integer, intent(in)    :: NEQ        ! number of equilibria
    integer, intent(inout) :: CTRL       ! control value. 0 on no error, 1 otherwise
    logical, intent(in)    :: is_first   ! whether it is first point or not
    integer,          intent(in),    dimension(NEQ, NCP) :: P  ! The stoichiometry array
    double precision, intent(inout), dimension(NEQ+NCP)  :: C  ! initial estimates (in) and refined values (out)
    double precision, intent(in),    dimension(NCP)      :: T  ! analytical concentration values
    double precision, intent(in),    dimension(NEQ)      :: B  ! the equilibrium constants
    !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
    !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
  
    DOUBLE PRECISION, parameter :: PCF_TOLERANCE = 0.25
    DOUBLE PRECISION, parameter :: PCF_TOLERANCE_FIRST = 1D-9
  
    CTRL = 0
    if (is_first) then
      call pcf(C, T, B, P, PCF_TOLERANCE_FIRST, NCP, NEQ)
    else
      call pcf(C, T, B, P, PCF_TOLERANCE, NCP, NEQ)
    end if
    call nr(C, T, B, P, CTRL, NCP, NEQ)
end subroutine solve_point


subroutine solve_point_solid(C, CS, T, B, P, KS, Q, CTRL, is_first, NEQ, NCP, NSL)
    ! Perform Newton-Raphson iterations per point when the beta matrix
    ! is 1D. The first initial value must be provided in the first row
    ! of C. The following points take as starting value the value of 
    ! the last iteration.
    implicit none
    integer,          intent(in)    :: NCP                      ! number of components
    integer,          intent(in)    :: NEQ                      ! number of equilibria
    integer,          intent(in)    :: NSL                      ! number of solids
    integer,          intent(inout) :: CTRL                     ! control value. 0 on no error, 1 otherwise
    logical,          intent(in)    :: is_first                 ! whether it is first point or not
    integer,          intent(in),    dimension(NEQ, NCP)  :: P  ! The stoichiometry array
    double precision, intent(inout), dimension(NEQ+NCP)   :: C  ! initial estimates (in) and refined values (out)
    double precision, intent(inout), dimension(NSL)       :: CS ! solid concentration initial estimates (in) and refined values (out)
    double precision, intent(in),    dimension(NCP)       :: T  ! analytical concentration values
    double precision, intent(in),    dimension(NEQ)       :: B  ! the equilibrium constants
    integer,          intent(in),    dimension(NSL, NCP)  :: Q  ! the equilibrium constants
    double precision, intent(in),    dimension(NSL)       :: KS ! the precipitation constants
  
    DOUBLE PRECISION, PARAMETER :: si_tolerance = 1e-6
    double PRECISION, dimension(NSL) :: si
    double PRECISION, dimension(NCP+NSL) :: F
    logical, dimension(NSL) :: do_solid
    integer :: i, its, NX
    DOUBLE PRECISION, parameter :: PCF_TOLERANCE = 0.25
    DOUBLE PRECISION, parameter :: PCF_TOLERANCE_FIRST = 1D-9
    integer, PARAMETER :: MAXITS=100
    do_solid(:) = .true.

    if (is_first) then
        call pcf_solids(C, CS, T, B, KS, P, Q, PCF_TOLERANCE_FIRST, NCP, NEQ, NSL)
    else
        call pcf_solids(C, CS, T, B, KS, P, Q, PCF_TOLERANCE, NCP, NEQ, NSL)
    end if
    do its=1,MAXITS
      call saturation_product(C, Q, NCP, NSL, si)
      SI = SI/KS
      do i=1,NSL
        if (si(i) < 1.0 + si_tolerance) then
            if (.not. do_solid(i)) then
                CS(i) = 0.0     ! if there is no solid, make CS=0.0
              end if
          do_solid(i) = .false.  ! remove from model
          F(NCP+i) = 0.0
          cycle
        end if
        if (si(i) > 1.0 + si_tolerance) then
          do_solid(i) = .true.
          if (CS(i) .eq. 0.0) then
              CS(i) = 1e-8
          end if
          F(NCP+i) = 1-SI(i)
          cycle
        end if
      end do
      NX = NCP + count(do_solid)
      call nr_solids_log(C, CS, T, B, P, KS, Q, ctrl, do_solid, NCP, NEQ, NSL, NX)
      call mb_solids(F, C, P, T, Q, CS, NCP, NEQ, NSL)
      if (.not. any(do_solid) .and. maxval(abs(F)) < 1e-20) exit    ! no solids to do
    end do
end subroutine solve_point_solid
