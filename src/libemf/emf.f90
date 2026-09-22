! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>
!
!----------------------------------------------------------------------------------------------------
! Routines for fitting equilibrium constants and other parameters from potentiometric data.
!
! HOW TO USE
! 1. call init_module with the parameters indicated
! 2. load the titrations one by one by calling load_titration
! 3. set the weighting scheme by calling one of the set weights routine
! 4. call fit_emf
! 5. copy the result from the public variables
! 6. call finalize to deallocate and cleanup
!----------------------------------------------------------------------------------------------------

module emf
    implicit none

    private

    ! variables for constraints
    type value_holder
        double precision, pointer :: origin_value   ! points to value in xvalues
        integer :: origin_index                     ! index of xvalues
        real, pointer :: destination_value          ! points to value in tmp_*
        integer :: destination_index1               ! index of EQ, CP or NT
        integer :: destination_index2               ! index of EQ, CP or NT
        integer :: id
        integer :: jtype                            ! type of variable (beta, init, etc.)
        double precision :: ratio                   ! 1.0 for free and reference; other for the rest
    end type value_holder

    type(value_holder), allocatable, dimension(:)         :: variables
    double precision,   allocatable, dimension(:), target :: xvalues
    double precision,   allocatable, dimension(:,:)       :: W
    double precision,   allocatable, dimension(:)         :: weights

    ! public variables for after fitting
    double precision, allocatable, public, protected, dimension(:,:) :: C           ! the free concentrations
    double precision, allocatable, public, protected, dimension(:,:) :: covariance  ! the covatiance matrix
    real, allocatable, target, public, protected, dimension(:)   :: electr_emf0   ! standard potential in mV
    real, allocatable, target, public, protected, dimension(:)   :: electr_slope  ! nernstian slope in unit
    real, allocatable, target, public, protected, dimension(:,:) :: titr_init     ! initial amount in mmol
    real, allocatable, target, public, protected, dimension(:,:) :: titr_buret    ! buret concentration in mmol/mL
    real, allocatable, target, public, protected, dimension(:)   :: lbeta         ! formation constants in log units
    real, allocatable, target, public, protected, dimension(:,:) :: junctionp     ! junction potential in mV.mL/mmmol

    ! private variables for internal use
    real,    allocatable,  dimension(:)   :: V            ! titre in mL
    real,    allocatable,  dimension(:)   :: V0           ! initial amount in mL
    real,    allocatable,  dimension(:)   :: E            ! experimental emf in mV
    integer, allocatable,  dimension(:,:) :: junctionc    ! indices of the components for the junction potential
    integer, allocatable,  dimension(:,:) :: junctionf    ! refinement flags for the junction potential
    integer, allocatable,  dimension(:)   :: eactive      ! indices of electroactive component
    integer, allocatable,  dimension(:)   :: flag_beta    ! refinement flags for formation constants
    integer, allocatable,  dimension(:,:) :: flag_init
    integer, allocatable,  dimension(:,:) :: flag_buret
    integer, allocatable,  dimension(:)   :: flag_emf0
    integer, allocatable,  dimension(:)   :: flag_slope
    real,    allocatable,  dimension(:)   :: err_emf      ! the error for potential, in mV, per titration
    real,    allocatable,  dimension(:)   :: err_vol      ! the error for volume, in mL, per titration
    integer, allocatable,  dimension(:)   :: titr_size    ! the size in experimental points per titration
    integer, allocatable,  dimension(:,:) :: P            ! the stoichiometry array
    real,    allocatable,  dimension(:,:) :: T            ! analytical concentration
    logical, allocatable,  dimension(:)   :: update_titr  ! flag to know which titration must be updated

    ! variables for handling customized refinement
    real, target, allocatable, dimension(:)   :: tmp_lbeta
    real, target, allocatable, dimension(:,:) :: tmp_init
    real, target, allocatable, dimension(:,:) :: tmp_buret
    real, target, allocatable, dimension(:)   :: tmp_emf0
    real, target, allocatable, dimension(:)   :: tmp_slope
    real, target, allocatable, dimension(:,:) :: tmp_junction  ! junction potential in mV.mL/mmmol
    real, pointer,             dimension(:)   :: use_lbeta
    real, pointer,             dimension(:,:) :: use_init
    real, pointer,             dimension(:,:) :: use_buret
    real, pointer,             dimension(:)   :: use_emf0
    real, pointer,             dimension(:)   :: use_slope
    real, pointer,             dimension(:,:) :: use_junction

    integer, public, protected :: NTT   ! number of titrations
    integer, public, protected :: NPT   ! number of titration points for all titrations
    integer, public, protected :: NCP   ! number of components
    integer, public, protected :: NEQ   ! number of equilibria
    character(len=40), public, protected :: emf_diagnostic
    real :: nernst
    integer :: tmp_loaded_titrs, tmp_space_left

    ! internal constants
    integer, parameter :: RF_CONSTANT = 0
    integer, parameter :: RF_REFINE = 1
    integer, parameter :: J_LBETA = 1
    integer, parameter :: J_TINIT = 2
    integer, parameter :: J_TBURE = 3
    integer, parameter :: J_EMF0  = 4
    integer, parameter :: J_SLOPE = 5
    integer, parameter :: J_JUNCT = 6
    double precision, parameter :: LG10=LOG(10.0)
    
    integer outu  ! unit for writing output

    PUBLIC :: init_module, load_titration, finalize, fit_emf, set_weight_gradient, set_weight_unit

contains

!----------------------------------------------------------------------------------------------------
! Initiate module. Allocates variables and some data. Call this one first.
!----------------------------------------------------------------------------------------------------
subroutine init_module(n_points, n_titrs, stoich, logbeta, flgbeta, n_junct, outn)
    integer, intent(in)                                 :: n_points  ! total number of experimental points
    integer, intent(in)                                 :: n_titrs   ! number of titrations to be loaded
    integer, intent(in), dimension(:,:)                 :: stoich    ! the stoichiometry array
    real,    intent(in), dimension(size(stoich, dim=1)) :: logbeta   ! the log10 of the formation constants
    integer, intent(in), dimension(size(stoich, dim=1)) :: flgbeta   ! the refinement flags for the formation constants
    integer, intent(in)                                 :: n_junct   ! number of junction components (usually zero)
    integer, intent(in)                                 :: outn

    emf_diagnostic = "no error"

    write (outn, '(A)') "Initiating module"

    if (n_titrs < 1) then
        emf_diagnostic = "n_titrs must be positive"
        return
    end if
    NTT = n_titrs

    if (n_points < 1) then
        emf_diagnostic = "n_points must be positive"
        return
    end if
    NPT = n_points

    NEQ = SIZE(stoich, dim=1)
    NCP = SIZE(stoich, dim=2)

    ALLOCATE(C(NPT,NCP+NEQ))
    C = 0.0
    ALLOCATE(V(NPT))
    V = 0.0
    ALLOCATE(V0(NTT))
    V0 = 0.0
    ALLOCATE(E(NPT))
    E = 0.0
    ALLOCATE(electr_emf0(NTT))
    electr_emf0 = 0.0
    ALLOCATE(electr_slope(NTT))
    electr_slope = 1.00
    ALLOCATE(eactive(NTT))
    eactive = 0
    ALLOCATE(titr_init(NTT,NCP)) 
    titr_init = 0.0
    ALLOCATE(titr_buret(NTT,NCP)) 
    titr_buret = 0.0
    ALLOCATE(lbeta(NEQ))
    ALLOCATE(flag_beta(NEQ))
    ALLOCATE(flag_init(NTT,NCP)) 
    flag_init = 0
    ALLOCATE(flag_buret(NTT,NCP)) 
    flag_buret = 0
    ALLOCATE(flag_emf0(NTT))
    flag_emf0 = 0
    ALLOCATE(flag_slope(NTT))
    flag_slope = 0
    ALLOCATE(err_emf(NTT))
    err_emf = 0.0
    ALLOCATE(err_vol(NTT))
    err_vol = 0.0
    ALLOCATE(titr_size(NTT))
    titr_size = 0
    ALLOCATE(P(NEQ, NCP))
    ALLOCATE(T(NPT,NCP))
    T = 0.0
    ALLOCATE(update_titr(NTT))
    update_titr = .false.
    ALLOCATE(W(NPT,NPT))
    W = 0.0
    ALLOCATE(weights(NPT))
    weights = 0.0
    allocate(junctionp(NTT,n_junct))
    allocate(junctionf(NTT,n_junct))
    allocate(junctionc(NTT,n_junct))

    P = stoich
    lbeta = REAL(logbeta * LG10)     ! convert from log10 to ln
    flag_beta = flgbeta

    tmp_loaded_titrs = 0
    tmp_space_left = NPT

    outu = outn
end subroutine init_module


!----------------------------------------------------------------------------------------------------
! Load single titration. The memory has been previously allocated by calling the init_module first.
!----------------------------------------------------------------------------------------------------
subroutine load_titration(titre, emf, titrparms, active, initV, emf0, nernst_slope, tflags, emf0_flag, slope_flag, &
        error_emf, error_vol, junctionpot, junctioncom, junctionflag)
    implicit none
    real,    intent(in), dimension(:)           :: titre       ! the titre in mL
    real,    intent(in), dimension(size(titre)) :: emf         ! experiental emf in mV
    real,    intent(in), dimension(NCP,2)       :: titrparms   ! titration parameters; dim(1)-> component;
                                                               ! dim(2)=1 -> init in mmol; dim(2)=2 -> buret c in mmol/mL
    real,    intent(in), dimension(:)   :: junctionpot  ! the junction potential in mV for each active species
    integer, intent(in), dimension(:)   :: junctioncom  ! the component index for the corresponding junction
    integer, intent(in), dimension(:)   :: junctionflag ! refinament flag for the juntion
    integer, intent(in), dimension(NCP,2)             :: tflags       ! refinement flags for the titration parameters
    real,    intent(in) :: initV                 ! initial volume in mL
    real,    intent(in) :: nernst_slope          ! the nernstian slope in x1
    real,    intent(in) :: emf0                  ! the standard potential in mL
    integer, intent(in) :: active                ! index of the electroactive component
    integer, intent(in) :: emf0_flag, slope_flag ! refinement flag for emf0 and slope
    real,    intent(in) :: error_emf             ! uncertainty for potential in mV
    real,    intent(in) :: error_vol             ! uncertainty for volume in mL

    integer n, bottom

    if (tmp_loaded_titrs .eq. NTT) then
        emf_diagnostic = "no room for more titrations"
        return
    end if

    if (all(shape(titrparms) .ne. shape(tflags))) then
        emf_diagnostic = "titrparms has wrong shape"
        return
    end if

    n = SIZE(titre, dim=1)
    if (n > tmp_space_left) then
        emf_diagnostic = "no room for more experimental points"
        return
    end if

    if ((active < 1) .or. (active > NCP)) then
        emf_diagnostic = "electroactive index out of bounds"
        return
    end if

    tmp_space_left = tmp_space_left - n
    tmp_loaded_titrs = tmp_loaded_titrs + 1
    if (tmp_loaded_titrs .eq. 1) then
        bottom = 1
    else
        bottom = 1+sum(titr_size(1:tmp_loaded_titrs-1))
    end if
    titr_size(tmp_loaded_titrs) = n
    V(bottom:bottom+n-1) = titre
    E(bottom:bottom+n-1) = emf
    V0(tmp_loaded_titrs) = initV
    eactive(tmp_loaded_titrs) = active

    electr_emf0(tmp_loaded_titrs) = emf0
    electr_slope(tmp_loaded_titrs) = nernst_slope

    titr_init(tmp_loaded_titrs,:) = titrparms(:, 1)
    titr_buret(tmp_loaded_titrs,:) = titrparms(:, 2)
    flag_init(tmp_loaded_titrs, :) = tflags(:,1)
    flag_buret(tmp_loaded_titrs, :) = tflags(:,2)
    flag_slope(tmp_loaded_titrs) = slope_flag
    flag_emf0(tmp_loaded_titrs) = emf0_flag
    err_emf(tmp_loaded_titrs) = error_emf
    err_vol(tmp_loaded_titrs) = error_vol

    junctionp(tmp_loaded_titrs,:) = junctionpot
    junctionc(tmp_loaded_titrs,:) = junctioncom
    junctionf(tmp_loaded_titrs,:) = junctionflag

    write (outu, '(I4,A)') n, " points loaded"
    
end subroutine load_titration


!----------------------------------------------------------------------------------------------------
! Select unit values as weight. The weight matrix will be a unit matrix.
!----------------------------------------------------------------------------------------------------
subroutine set_weight_unit()
    integer i
    W = 0.0
    do i=1,NPT
        W(i,i) = 1.0
    end do
    weights = 1.0
    write (outu, '(A)') "weighting scheme set to UNIT"
end subroutine set_weight_unit


!----------------------------------------------------------------------------------------------------
! Select weight matrix as proportional to the gradient.
! W_ij = 1/(sigma_emf**2 * gradient(i) + sigma_vol^2)  if i=j
! W_ij = 0  if i /= j
!----------------------------------------------------------------------------------------------------
subroutine set_weight_gradient()
    use misc, only: derivative_1d
    real, allocatable, dimension(:) :: dy
    integer i, j, blow, bhigh


    do i=1,NTT
        call titr_bounds(i, blow, bhigh)
    end do

    allocate(dy(NPT))
    do i=1,NTT
        call titr_bounds(i, blow, bhigh)
        call derivative_1d(titr_size(i), V(blow:bhigh), E(blow:bhigh), dy(blow:bhigh))
        do j=blow,bhigh
            W(j,j) = 1.0/(err_emf(i)**2 + dy(j)**2 * err_vol(i)**2)
            weights(j) = W(j,j)
        end do
    end do
    deallocate(dy)
    write (outu, '(A)') "weighting scheme set to AUTO"
end subroutine set_weight_gradient


!----------------------------------------------------------------------------------------------------
! Cleanup routine and deallocate memory assets.
!----------------------------------------------------------------------------------------------------
subroutine finalize()
    DEALLOCATE(C)
    DEALLOCATE(V)
    DEALLOCATE(V0)
    DEALLOCATE(E)
    DEALLOCATE(electr_emf0)
    DEALLOCATE(eactive)
    DEALLOCATE(electr_slope)
    DEALLOCATE(titr_init) 
    DEALLOCATE(titr_buret) 
    DEALLOCATE(lbeta)
    DEALLOCATE(flag_beta)
    DEALLOCATE(flag_init) 
    DEALLOCATE(flag_buret) 
    DEALLOCATE(flag_emf0) 
    DEALLOCATE(flag_slope) 
    deallocate(err_emf)
    deallocate(err_vol)
    DEallocate(titr_size)
    DEallocate(P)
    DEALLOCATE(T)
    DEALLOCATE(update_titr)
    DEALLOCATE(W)
    DEALLOCATE(weights)
    deallocate(junctionp)
    deallocate(junctionf)
    deallocate(junctionc)
    if (allocated(variables)) deallocate(variables)
    if (allocated(covariance)) deallocate(covariance)
end subroutine finalize


subroutine fit_emf(maxits, temp)
    use misc, only: calc_analc
    implicit none
    real,   intent(in) :: temp   ! in degrees Celsius
    integer, intent(in) :: maxits

    double precision :: cov
    real, parameter :: nernst_base = 0.08617333262 ! mV/K
    integer NX, NV
    integer i, j, blow, bhigh
    double precision, parameter :: TDAMP = 0.01
    double precision, parameter :: TMIN = 1D-10
    real, PARAMETER :: TOKELVIN = 273.15
    !character(len=8) :: cflag

    ! handle temperature
    nernst = nernst_base*(temp+TOKELVIN)  ! mV

    ! create constraints
    call count_dof(NX, NV)

    write (outu, '(2(A,I3),A)') "Fitting ", NX, " degrees of freedom and ", NV, " variables"

    allocate(xvalues(NX))
    allocate(variables(NV))
    call create_variables()

    ! Initialise T
    do i=1,NTT
        if (sum(flag_init(i,:))+sum(flag_buret(i,:)) > 0) update_titr(i) = .true.
        call titr_bounds(i, blow, bhigh)
        call calc_analc(V(blow:bhigh), V0(i), use_init(i,:), use_buret(i,:), T(blow:bhigh,:))
    end do

    ! Initialise C
    do i=1,NPT
        do j=1,NCP
            if (T(i,j) > 0.0) then
                C(i,j) = TDAMP*T(i,j)
            else
                C(i,j) = TMIN
            end if
        end do
    end do

    call lm(maxits)

    ! post processing
    ! update values
    if (allocated(tmp_lbeta)) then
        lbeta = tmp_lbeta / REAL(LG10)
        DEALLOCATE(tmp_lbeta)
    end if

    if (allocated(tmp_init)) then
        titr_init = tmp_init
        deallocate(tmp_init)
    end if

    if (allocated(tmp_buret)) then
        titr_buret = tmp_buret
        deallocate(tmp_buret)
    end if

    if (allocated(tmp_emf0)) then
        electr_emf0 = tmp_emf0
        deallocate(tmp_emf0)
    end if

    if (allocated(tmp_slope)) then
        electr_slope = tmp_slope
        deallocate(tmp_slope)
    end if

    if (allocated(tmp_junction)) then
        junctionp = tmp_junction
        deallocate(tmp_junction)
    end if

    ! print values
    write (outu, "(/,A)") "REFINED VARIABLES"
    write (outu, "(A)") "  #  PARM       VALUE       ERROR"
    do i=1,size(variables)
        cov = variables(i)%ratio*SQRT(covariance(variables(i)%origin_index,variables(i)%origin_index))
        select case (variables(i)%jtype)
            case (1)   ! log beta
                write (outu, "(I3,A6,2F12.4,*(I4))") i, "LOGB", variables(i)%origin_value/LG10, cov/LG10, &
                    P(i,:)
            case (2)   ! initial amount
                write (outu, "(I3,A6,2E12.4)") i, "init", variables(i)%origin_value, cov
            case (3)   ! buret concentration
                write (outu, "(I3,A6,2F12.4)") i, "buret", variables(i)%origin_value, cov
            case (4)   ! emf0
                write (outu, "(I3,A6,2F12.4)") i, "emf0", variables(i)%origin_value, cov
            case (5)   ! slope
                write (outu, "(I3,A6,2F12.4)") i, "slope", variables(i)%origin_value, cov
            case (6)   ! slope
                write (outu, "(I3,A6,2F12.4)") i, "junct", variables(i)%origin_value, cov
        end select
    end do

    ! deallocation
    DEALLOCATE(xvalues, variables)
end subroutine fit_emf


!----------------------------------------------------------------------------------------------------
! Run the Levenberg-Marquardt loop.
!----------------------------------------------------------------------------------------------------
subroutine lm(maxits)
    ! routine for Levenberg-Marquardt iterations.
    use misc, only: fit_sigma, lesolve
    implicit none
    integer, intent(in) :: maxits

    integer :: NX ! number of variables to optimize
    integer iter  ! iteration number

    real, parameter :: DAMPING_UPPER   = 1e10
    real, parameter :: DAMPING_LOWER   = 1e-12
    real, parameter :: DAMPING_UPF     = 8.0
    real, parameter :: DAMPING_LOWF    = 3.0
    real, parameter :: SIGMA_THRESHOLD = 1e-2
    real, parameter :: GRAD_THRESHOLD  = 1e-2
    real, parameter :: STEP_THRESHOLD  = 1e-5
    !double precision, PARAMETER :: RHO_THRESH = 1D-4
    double precision, allocatable, dimension(:,:) :: jac
    double precision, allocatable, dimension(:)   :: dx
    double precision, allocatable, dimension(:)   :: fvec, trial_fvec, F
    double precision, allocatable, dimension(:)   :: gradient
    double precision, allocatable, dimension(:,:) :: M, dampedM
    double precision :: gradient_norm, sumsq, damping, trial_sumsq
    double precision :: sigma
    integer i
    integer CTRL

    NX = size(xvalues)
    ALLOCATE(jac(NPT, NX))
    ALLOCATE(dx(NX))
    ALLOCATE(gradient(NX))
    ALLOCATE(fvec(NPT))
    ALLOCATE(trial_fvec(NPT))
    ALLOCATE(M(NX,NX))
    ALLOCATE(dampedM(NX,NX))
    ALLOCATE(F(NX))

    dx = 0.0
    jac = 0.0
    fvec = 0.0
    damping = 1D+2
    sumsq = HUGE(1D1)
    trial_sumsq = 0.0
    trial_fvec = 0.0

    iter = 0
    main_loop: do while (iter < maxits)
        call try_step(dx)

        call update_equil()
        call calc_resid(fvec)
        trial_sumsq = SUM(weights*fvec**2)

        if (iter .eq. 0) then
            sumsq = trial_sumsq
            call jacobian_lm(jac)

            M = MATMUL(MATMUL(TRANSPOSE(jac), W), jac)
            F = -MATMUL(MATMUL(TRANSPOSE(jac), W), fvec)
            dampedM = M
            do i=1,NX
                dampedM(i,i) = (1.0 + damping) * M(i,i)
            end do
            call lesolve(dampedM, F, dx, NX, CTRL)
        end if

        if (CTRL .ne. 0) then
            emf_diagnostic="linear system could not be solved"
            return
        end if

        if (iter .eq. 0) then
            iter = iter + 1
            cycle
        end if

        !rho = calc_rho((chisq - trial_chisq), dx, gradient, M, damping)
        !if (abs(rho) > RHO_THRESH) then                     ! sometimes rho can be negative very close to the solution
        if (trial_sumsq < sumsq) then
            ! step accepted
            write(outu, '(/,A,I4,/)') "Iteration", iter
            call report_variables(dx)

            xvalues = xvalues + dx
            sumsq = trial_sumsq

            call jacobian_lm(jac)

            M = MATMUL(MATMUL(TRANSPOSE(jac), W), jac)
            F = -MATMUL(MATMUL(TRANSPOSE(jac), W), fvec)
            gradient = 2*MATMUL(MATMUL(TRANSPOSE(jac),W),fvec)
            gradient_norm = NORM2(gradient)
            sigma = fit_sigma(fvec, weights, NPT, NX)

            write(outu, '(A,E14.6)') "   sumsq", sumsq

            ! Adaptive damping (very effective)
            !if (rho > 0.75) then
            !    damping = max(damping / 2*DAMPING_LOWF, DAMPING_LOWER)
            !else
            !    damping = max(damping / DAMPING_LOWF, DAMPING_LOWER)
            !end if
            damping = max(damping / DAMPING_LOWF, DAMPING_LOWER)

            ! Exit check
            ! Sigma check
            write(outu, '(A,E14.6)', advance='no') "   sigma", sigma
            if (sigma < SIGMA_THRESHOLD) then
                write(outu, '(A)') " YES"
                exit
            else
                write(outu, '(A)') " NO"
            end if

            ! Gradient norm check
            write(outu, '(A,E14.6)', advance='no') "   gnorm", gradient_norm
            if (gradient_norm < GRAD_THRESHOLD) then
                write(outu, '(A)') " YES"
                exit
            else
                write(outu, '(A)') " NO"
            end if

            ! step length check
            i = MAXLOC(abs(dx), dim=1)
            write(outu, '(A,E14.6)', advance='no') "   step ", dx(i)
            if (abs(dx(i)) < STEP_THRESHOLD) then
                write(outu, '(A)') " YES"
                exit
            else
                write(outu, '(A)') " NO"
            end if

            iter = iter + 1
        else
            ! step rejected
            damping = damping*DAMPING_UPF
            if (damping > DAMPING_UPPER) then
                emf_diagnostic="LM loop rejected too many steps"
                exit
            end if
            call lesolve(dampedM, F, dx, NX, CTRL)
            write (outu, "(A,E14.4)") "damping = ", damping
            !write (outu, "(2(A,E14.4))") "trial sumsq=", trial_sumsq, " sumsq=", sumsq
        end if

        dampedM = M
        do i=1,NX
            dampedM(i,i) = (1.0 + damping) * M(i,i)
        end do
        call lesolve(dampedM, F, dx, NX, CTRL)
    end do main_loop

    call calc_covariance(M)

    DEALLOCATE(jac, dx, fvec, M)
end subroutine lm

!----------------------------------------------------------------------------------------------------
! Update the free concentrations.
!----------------------------------------------------------------------------------------------------
subroutine update_equil()
    use misc, only: calc_analc
    integer i, blow, bhigh, CTRL
    external eq_000
    do i=1,NTT
        call titr_bounds(i, blow, bhigh)
        if (update_titr(i)) call calc_analc(V(blow:bhigh), V0(i), use_init(i,:), use_buret(i,:), T(blow:bhigh,:))
        call eq_000(C(blow:bhigh,:), DBLE(T(blow:bhigh,:)), EXP(DBLE(use_lbeta)), P, CTRL, titr_size(i), NCP, NEQ)
    end do
end subroutine update_equil


!----------------------------------------------------------------------------------------------------
! Update the free variables in each step.
!----------------------------------------------------------------------------------------------------
subroutine try_step(step)
    double precision, intent(in), dimension(:) :: step

    double precision :: istep
    integer i

    ! update tmp values
    do i=1,size(variables)
        istep = step(variables(i)%origin_index)
        variables(i)%destination_value = real(variables(i)%ratio*(variables(i)%origin_value + istep))
    end do
end subroutine try_step


!----------------------------------------------------------------------------------------------------
! Compute residuals.
!----------------------------------------------------------------------------------------------------
subroutine calc_resid(fvec)
    implicit none
    double precision, intent(out), DIMENSION(:) :: fvec

    integer i, j
    !integer lb, ub
    integer hact
    integer curr_titr
    integer curr_jcomp
    double precision calc_emf
    double precision junction

    curr_titr = 0

    do i=1,NPT
        if (i > sum(titr_size(1:curr_titr))) then             
            curr_titr = curr_titr + 1
            hact = eactive(curr_titr)
        end if
        junction=0.0 
        do j=1,size(use_junction,dim=2)
            curr_jcomp = junctionc(curr_titr,j)
            junction=junction+use_junction(curr_titr,j)*C(i,curr_titr)
        end do
        calc_emf = use_emf0(curr_titr)  &
                 + use_slope(curr_titr) * nernst * LOG(C(i,hact))  & 
                 + junction
        fvec(i) = E(i) - calc_emf
    end do

    ! l_titr: do i=1,NTT
    !     hact = eactive(i)
    !     call titr_bounds(i, lb, ub)
    !     fvec(lb:ub) = DBLE(                                 &
    !         E(lb:ub) - (                                    &
    !             use_emf0(i) +                               &
    !             use_slope(i) * nernst * LOG(C(lb:ub,hact))  &
    !             + sum(use_junction(i)*C(lb:ub,)
    !         )                                               &
    !     )
    !     if (i > sum(titr_size(1:curr_titr))) then
    !         curr_titr = curr_titr + 1
    !     end if
    ! end do l_titr

end subroutine calc_resid


!------------------------------------------------------------------------------------------------------------------------
! Compute the jacobian for the Levenberg-Marquardt iteration.
! The rows are the experimental points; the columns are the variables to be optimized in the
! following order: (1) equilibrium constants; (2) initial amounts; (3) buret concentrations;
! (4) standard potentials; (5) nernstian slopes; (6) constraints
subroutine jacobian_lm(jac)
    use misc, only: msysolve
    double precision, intent(out), dimension(:,:) :: jac

    double precision, DIMENSION(NCP, NCP) :: Am      ! the A-matrix
    double precision, allocatable, dimension(:,:) :: bvec
    integer i, j, k, l, CTRL
    integer curr_titr
    integer aux, eac
    integer bsize
    logical do_beta, do_init, do_buret, do_amatrix
    integer, DIMENSION(6,size(variables)) :: filter_jtype
    double precision tmpd
    integer tmpi

    jac = 0.0
    bsize = 0
    do_beta = allocated(tmp_lbeta)
    do_init = allocated(tmp_init)
    do_buret = allocated(tmp_buret)
    do_amatrix = do_beta .or. do_init .or. do_buret

    filter_jtype = 0
    do i=1,6
        do j=1,size(variables)
            if (variables(j)%jtype .eq. i) filter_jtype(i,j) = j
        end do
    end do

    bsize = 0
    do i=1,size(variables)
        if (variables(i)%jtype < 4) bsize = bsize + 1
    end do
    allocate(bvec(NCP,bsize))

    curr_titr = 1       ! track which titration are we currently in

    do i=1,NPT
        if (do_amatrix) then
            call amatrix(Am, C(i,:), P, NCP, NEQ)

            do j=1,size(variables)
                if (variables(j)%jtype > 4) exit
                aux = variables(j)%destination_index1
                select case (variables(j)%jtype)
                    case (J_LBETA)
                        do k=1,NCP
                            bvec(k, j) = -P(aux,k)*C(i,aux+NCP)
                        end do
                    case (J_TINIT)
                        if (variables(j)%destination_index2 .ne. curr_titr) then
                            bvec(:,j) = 0.0
                            cycle
                        end if
                        do k=1,NCP
                            if (aux .eq. k) then
                                bvec(k, j) = 1.0/(V(i)+V0(curr_titr))
                            else
                                bvec(k, j) = 0.0
                            end if
                        end do
                    case (J_TBURE)
                        if (variables(j)%destination_index2 .ne. curr_titr) then
                            bvec(:,j) = 0.0
                            cycle
                        end if
                        do k=1,NCP
                            if (aux .eq. k) then
                                bvec(k, j) = V(i)/(V(i)+V0(curr_titr))
                            else
                                bvec(k, j) = 0.0
                            end if
                        end do
                end select
            end do
            call msysolve(Am, bvec, CTRL)
            if (CTRL .eq. 99) then
                emf_diagnostic="Error computing jacobian"
                return
            end if
        end if

        do j=1,size(variables)
            eac = eactive(curr_titr)
            aux = variables(j)%origin_index
            select case (variables(j)%jtype)
                case (J_LBETA)
                    jac(i, aux) = jac(i, aux) - variables(j)%ratio*use_slope(curr_titr)*nernst*bvec(eac,j)
                    do k=1,size(junctionp,dim=2)
                        tmpi = junctionc(curr_titr,k)
                        if (tmpi > NCP) then
                            tmpd = DOT_PRODUCT(P(tmpi-NCP,:),bvec(:,j))
                            if (j .eq. tmpi) tmpd = tmpd + 1.0
                        else
                            tmpd = use_junction(curr_titr,k)*C(i,tmpi)*bvec(tmpi,j)
                        end if
                        jac(i, aux) = jac(i, aux) + tmpd
                    end do
                case (J_TINIT:J_TBURE)
                    if (curr_titr .eq. variables(j)%destination_index2) then
                        tmpd = - variables(j)%ratio*use_slope(curr_titr)*nernst*bvec(eac,j)
                        jac(i, aux) = jac(i, aux) + tmpd
                    end if
                case (J_EMF0)
                    if (curr_titr .eq. variables(j)%destination_index1) then
                        jac(i, aux) = jac(i, aux) - variables(j)%ratio
                    end if
                case (J_SLOPE)
                    if (curr_titr .eq. variables(j)%destination_index1) then
                        jac(i, aux) = jac(i, aux) - variables(j)%ratio*nernst
                    end if
                case (J_JUNCT)
                    if (curr_titr .eq. variables(j)%destination_index2) then
                        tmpi = junctionc(curr_titr,variables(j)%destination_index1)
                        tmpd = - variables(j)%ratio*C(i,tmpi)
                        jac(i, aux) = jac(i, aux) + tmpd
                    end if
            end select
        end do

        if (i > sum(titr_size(1:curr_titr))) then
            curr_titr = curr_titr + 1
        end if
    end do

    deallocate(bvec)
end subroutine jacobian_lm


subroutine create_variables()
    implicit none

    integer i, vcount, xcount
    real, dimension(2:9) :: vref
    integer, dimension(2:9) :: iref

    vcount = 0
    xcount = 0
    vref = 0.0
    iref = 0

    ! I. count formation constants
    if (sum(flag_beta) > 0) then        ! something to refine here
        allocate(tmp_lbeta(NEQ))
        tmp_lbeta = lbeta
        call aux_create_variables(tmp_lbeta, flag_beta, xcount, vcount, iref, vref, J_LBETA, 0)
        use_lbeta => tmp_lbeta
    else
        use_lbeta => lbeta
    end if
    
    ! II. count initial amounts
    if (sum(flag_init) > 0) then        ! something to refine here
        allocate(tmp_init(NTT,NCP))
        tmp_init = titr_init
        use_init => tmp_init
        do i=1,NTT
            call aux_create_variables(tmp_init(i,:), flag_init(i,:), xcount, vcount, iref, vref, J_TINIT, i)
        end do
    else
        use_init => titr_init
    end if

    ! III. count buret concentrations
    if (sum(flag_buret) > 0) then        ! something to refine here
        allocate(tmp_buret(NTT,NCP))
        tmp_buret = titr_buret
        use_buret => tmp_buret
        do i=1,NTT
            call aux_create_variables(tmp_buret(i,:), flag_buret(i,:), xcount, vcount, iref, vref, J_TBURE, i)
        end do
    else
        use_buret => titr_buret
    end if

    ! IV: count standard potentials
    if (sum(flag_emf0) > 0) then        ! something to refine here
        allocate(tmp_emf0(NTT))
        tmp_emf0 = electr_emf0
        use_emf0 => tmp_emf0
        call aux_create_variables(tmp_emf0, flag_emf0, xcount, vcount, iref, vref, J_EMF0, 0)
    else
        use_emf0 => electr_emf0
    end if

    ! V. count nernstian slopes
    if (sum(flag_slope) > 0) then        ! something to refine here
        allocate(tmp_slope(NTT))
        tmp_slope = electr_slope
        use_slope => tmp_slope
        call aux_create_variables(tmp_slope, flag_slope, xcount, vcount, iref, vref, J_SLOPE, 0)
    else
        use_slope => electr_slope
    end if
    
    ! VI. count junction potentials
    if (sum(junctionf) > 0) then        ! something to refine here
        allocate(tmp_junction(NTT,size(junctionp,dim=2)))
        tmp_junction = junctionp
        use_junction => tmp_junction
        do i=1,NTT
            call aux_create_variables(tmp_junction(i,:), junctionf(i,:), xcount, vcount, iref, vref, J_JUNCT, i)
        end do
    else
        use_junction => junctionp
    end if
end subroutine create_variables


subroutine aux_create_variables(vector, flags, xcount, vcount, iref, vref, jtype, index2)
    real, intent(in), target, DIMENSION(:) :: vector
    integer, intent(in), DIMENSION(:) :: flags
    integer, intent(inout) :: xcount, vcount
    real, intent(inout), dimension(2:9) :: vref
    integer, intent(inout), dimension(2:9) :: iref
    integer, intent(in) :: jtype
    integer, intent(in) :: index2

    integer i
    do i=1,size(flags)
        if (flags(i) == RF_CONSTANT) cycle
        if (flags(i) == RF_REFINE) then
            vcount = vcount + 1
            xcount = xcount + 1
            xvalues(xcount) = vector(i)
            variables(vcount)%origin_value => xvalues(xcount)
            variables(vcount)%origin_index = xcount
            variables(vcount)%destination_value => vector(i)
            variables(vcount)%ratio = 1.0
        end if
        if (flags(i) > 1) then
            vcount = vcount + 1
            if (iref(flags(i)) .eq. 0) then   ! constraint, first time
                xcount = xcount + 1
                xvalues(xcount) = vector(i)
                variables(vcount)%origin_value => xvalues(xcount)
                variables(vcount)%destination_value => vector(i)
                variables(vcount)%ratio = 1.0
                variables(vcount)%origin_index = xcount

                iref(flags(i)) = xcount
                vref(flags(i)) = vector(i)
            else                                  ! constraint already set
                variables(vcount)%origin_value => xvalues(iref(flags(i)))
                variables(vcount)%destination_value => vector(i)
                variables(vcount)%ratio = vector(i)/vref(flags(i))
                variables(vcount)%origin_index = xcount
            end if
        end if
        variables(vcount)%id = flags(i)
        variables(vcount)%destination_index1 = i
        variables(vcount)%destination_index2 = index2
        variables(vcount)%jtype = jtype
    end do
end subroutine aux_create_variables


subroutine count_dof(NX, NV)
    ! count degrees of freedom
    implicit none
    integer, intent(out) :: NX ! number of free variables
    integer, intent(out) :: NV ! total number of variables

    logical constr(2:9)
    integer i
    NX = 0
    NV = 0
    constr = .false.

    ! I. count formation constants
    call sub_count_dof(flag_beta, NX, NV, constr)

    ! II. count initial amounts
    do i=1,size(flag_init, dim=1)
        call sub_count_dof(flag_init(i,:), NX, NV, constr)
    end do

    ! III. count buret concentrations
    do i=1,size(flag_buret, dim=1)
        call sub_count_dof(flag_buret(i,:), NX, NV, constr)
    end do

    ! IV: count standard potentials
    call sub_count_dof(flag_emf0, NX, NV, constr)

    ! V. count nernstian slopes
    call sub_count_dof(flag_slope, NX, NV, constr)

    ! VI. count junction potentials
    do i=1,size(junctionf, dim=1)
        call sub_count_dof(junctionf(i,:), NX, NV, constr)
    end do
end subroutine count_dof


subroutine sub_count_dof(vec, NX, NV, constr)
    ! auxiliary routine to count degrees of freedom
    integer, intent(out) :: NX ! number of variables to optimize
    integer, intent(out) :: NV ! number of constrained variables
    logical, intent(inout) :: constr(:)
    integer, intent(in) :: vec(:)

    integer i
    do i=1,size(vec)
        if (vec(i) .eq. RF_REFINE) then
            NX = NX + 1
            NV = NV + 1
        else if (vec(i) > RF_REFINE) then
            NV = NV + 1
            if (.not. constr(vec(i))) then
                NX = NX + 1
                constr(vec(i)) = .true.
            end if
        end if
    end do
end subroutine sub_count_dof


subroutine titr_bounds(i, low, high)
    integer, intent(in) :: i
    integer, intent(out) :: low, high
    low = 1+sum(titr_size(1:(i-1)))  ! this should work also for i==1
    high = low + titr_size(i) - 1
end subroutine titr_bounds


function calc_chisq(resid, npars) result(chisq)
    implicit none
    double precision, intent(in), dimension(:) :: resid
    integer, intent(in) :: npars
    double precision :: chisq

    chisq = sum(resid**2)/(NPT-npars)
end function calc_chisq


function calc_rho(descent, step, gradient, M, damping) result(rho)
    double precision, INTENT(IN) :: descent
    double precision, INTENT(IN), DIMENSION(:) :: step
    double precision, INTENT(IN), DIMENSION(:) :: gradient
    double precision, INTENT(IN), DIMENSION(:,:) :: M
    double precision, INTENT(IN) :: damping

    double precision :: rho
    double precision :: predicted
    double precision, DIMENSION(size(M,dim=1)) :: D

    integer i

    D = [ (M(i,i), i=1,size(M,dim=1)) ]
    predicted = -DOT_PRODUCT(step, gradient) &
                -0.5*DOT_PRODUCT(step, damping*(D*step))
    rho = descent/predicted
end function calc_rho


! Given the Hessian matrix, calculate the covariance
subroutine calc_covariance(H)
    implicit none
    double precision, intent(in), dimension(:,:) :: H

    integer n, info, i, j

    n = size(H, dim=1)
    allocate(covariance(n,n))

    covariance = H

    ! 1. LU factorization
    call DPOTRF('U', n, covariance, n, info)
    if (info .ne. 0) then
        emf_diagnostic="COV could not be calculated"
        return
    end if

    ! 2. Invert using the factors
    call DPOTRI('U', n, covariance, n, info)
    if (info .ne. 0) then
        emf_diagnostic="COV could not be calculated"
        return
    end if

    ! 3. Copy the lower triangular matrix
    do i=1,n
        do j=i+1,n
            covariance(j,i) = covariance(i,j)
        end do
    end do
end subroutine calc_covariance


subroutine report_variables(step)
    double precision, intent(in), dimension(:) :: step
    double precision :: value_
    double precision :: step_val
    integer i

    write(outu, '(A)') "   #    previous        step         new"
    do i=1,size(variables)
        value_ = variables(i)%ratio*variables(i)%origin_value
        step_val = variables(i)%ratio*step(variables(i)%origin_index)
        select case (variables(i)%jtype)
            case (J_LBETA)
                write(outu, '(I4,3F12.4)') i, value_/LG10, step_val/LG10, (value_+step_val)/LG10
            case (J_TINIT)
                write(outu, '(I4,3E12.4)') i, value_, step_val, (value_+step_val)
            case default
                write(outu, '(I4,3F12.4)') i, value_, step_val, (value_+step_val)
        end select
    end do
end subroutine report_variables

end module emf
