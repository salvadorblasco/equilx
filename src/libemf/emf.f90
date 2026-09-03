! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>
!
!----------------------------------------------------------------------------------------------------
! Routines for fitting equilibrium constants and other parameters from potentiometric data.
!

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
    type(value_holder), allocatable, dimension(:) :: variables
    double precision, allocatable, target, dimension(:) :: xvalues
    DOUBLE PRECISION, allocatable, DIMENSION(:,:) :: W

    double precision, allocatable, public, protected, dimension(:,:)  :: C    ! the free concentrations
    real,    allocatable,         dimension(:)   :: V              ! titre in mL
    real,    allocatable,         dimension(:)   :: V0             ! initial amount in mL
    real,    allocatable,         dimension(:)   :: E              ! experimental emf in mV
    real,    allocatable, target, dimension(:)   :: E0             ! standard potential in mV
    real,    allocatable, target, dimension(:)   :: slope          ! nernstian slope in unit
    real,    allocatable, target, dimension(:,:) :: titr_init      ! initial amount in mmol
    real,    allocatable, target, dimension(:,:) :: titr_buret     ! buret concentration in mmol/mL
    real,    allocatable, target, dimension(:)   :: lbeta          ! formation constants in log units
    integer, allocatable,         dimension(:)   :: eactive        ! indices of electroactive component
    integer, allocatable,         dimension(:)   :: flag_beta      ! refinement flags for formation constants
    integer, allocatable,         dimension(:,:) :: flag_init
    integer, allocatable,         dimension(:,:) :: flag_buret
    integer, allocatable,         dimension(:)   :: flag_emf0
    integer, allocatable,         dimension(:)   :: flag_slope
    real,    allocatable,         dimension(:)   :: err_emf        ! the error for potential, in mV, per titration
    real,    allocatable,         dimension(:)   :: err_vol        ! the error for volume, in mL, per titration
    integer, allocatable,         dimension(:)   :: titr_size      ! the size in experimental points per titration
    integer, allocatable,         dimension(:,:) :: P              ! the stoichiometry array
    real,    allocatable, target, dimension(:)   :: tmp_lbeta
    real,    allocatable, target, dimension(:,:) :: tmp_init
    real,    allocatable, target, dimension(:,:) :: tmp_buret
    real,    allocatable, target, dimension(:)   :: tmp_emf0
    real,    allocatable, target, dimension(:)   :: tmp_slope
    real,    pointer,             dimension(:)   :: use_lbeta
    real,    pointer,             dimension(:,:) :: use_init
    real,    pointer,             dimension(:,:) :: use_buret
    real,    pointer,             dimension(:)   :: use_emf0
    real,    pointer,             dimension(:)   :: use_slope
    real,    allocatable,         dimension(:,:) :: T
    logical, allocatable,         dimension(:)   :: update_titr    ! flag to know which titration must be updated
    integer, public, protected :: NTT   ! number of titrations
    integer, public, protected :: NPT   ! number of titration points for all titrations
    integer, public, protected :: NCP   ! number of components
    integer, public, protected :: NEQ   ! number of equilibria
    character(len=40), public, protected :: emf_diagnostic
    real :: nernst
    integer :: tmp_loaded_titrs, tmp_space_left

    integer, parameter :: RF_CONSTANT = 0
    integer, parameter :: RF_REFINE = 1
    integer, parameter :: J_LBETA = 1
    integer, parameter :: J_TINIT = 2
    integer, parameter :: J_TBURE = 3
    integer, parameter :: J_EMF0  = 4
    integer, parameter :: J_SLOPE = 5

    PUBLIC :: init_module, load_titration, finalize, fit_emf, set_weight_gradient, set_weight_unit

contains

!----------------------------------------------------------------------------------------------------
! Initiate module. Allocates variables and some data. Call this one first.
!----------------------------------------------------------------------------------------------------
subroutine init_module(n_points, n_titrs, stoich, logbeta, flgbeta)
    integer, intent(in)                                 :: n_points  ! total number of experimental points
    integer, intent(in)                                 :: n_titrs   ! number of titrations to be loaded
    integer, intent(in), dimension(:,:)                 :: stoich    ! the stoichiometry array
    real,    intent(in), dimension(size(stoich, dim=1)) :: logbeta   ! the log10 of the formation constants
    integer, intent(in), dimension(size(stoich, dim=1)) :: flgbeta   ! the refinement flags for the formation constants

    emf_diagnostic = "no error"

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
    ALLOCATE(E0(NTT))
    E0 = 0.0
    ALLOCATE(slope(NTT))
    slope = 1.00
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

    P = stoich
    lbeta = logbeta * LOG(10.0)     ! convert from log10 to ln
    flag_beta = flgbeta

    tmp_loaded_titrs = 0
    tmp_space_left = NPT
end subroutine init_module


!----------------------------------------------------------------------------------------------------
! Load single titration. The memory has been previously allocated by calling the init_module first.
!----------------------------------------------------------------------------------------------------
subroutine load_titration(titre, emf, titrparms, active, initV, emf0, nernst_slope, tflags, emf0_flag, slope_flag, &
        error_emf, error_vol)
    implicit none
    real,    intent(in), dimension(:)           :: titre       ! the titre in mL
    real,    intent(in), dimension(size(titre)) :: emf         ! experiemntal emf in mV
    real,    intent(in), dimension(NCP,2)       :: titrparms
    integer, intent(in), dimension(:,:)         :: tflags
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

    E0(tmp_loaded_titrs) = emf0
    slope(tmp_loaded_titrs) = nernst_slope

    titr_init(tmp_loaded_titrs,:) = titrparms(:, 1)
    titr_buret(tmp_loaded_titrs,:) = titrparms(:, 2)
    flag_init(tmp_loaded_titrs, :) = tflags(:,1)
    flag_buret(tmp_loaded_titrs, :) = tflags(:,2)
    flag_slope(tmp_loaded_titrs) = slope_flag
    flag_emf0(tmp_loaded_titrs) = emf0_flag
    err_emf(tmp_loaded_titrs) = error_emf
    err_vol(tmp_loaded_titrs) = error_vol
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
end subroutine set_weight_unit


!----------------------------------------------------------------------------------------------------
! Select weight matrix as proportional to the gradient.
! W_ij = 1/(sigma_emf**2 * gradient(i) + sigma_vol^2)  if i=j
! W_ij = 0  if i /= j
!----------------------------------------------------------------------------------------------------
subroutine set_weight_gradient()
    real, allocatable, dimension(:) :: dy
    integer i, j, blow, bhigh
    allocate(dy(NPT))
    do i=1,NTT
        call titr_bounds(i, blow, bhigh)
        call derivative_1d(titr_size(i), V(blow:bhigh), E(blow:bhigh), dy(blow:bhigh))
        do j=blow,bhigh
            W(j,j) = 1.0/(err_emf(i)**2 + dy(j)**2 * err_vol(i)**2)
        end do
    end do
    deallocate(dy)
end subroutine set_weight_gradient


!----------------------------------------------------------------------------------------------------
! Cleanup routine and deallocate memory assets.
!----------------------------------------------------------------------------------------------------
subroutine finalize()
    DEALLOCATE(C)
    DEALLOCATE(V)
    DEALLOCATE(V0)
    DEALLOCATE(E)
    DEALLOCATE(E0)
    DEALLOCATE(eactive)
    DEALLOCATE(slope)
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
    if (allocated(variables)) deallocate(variables)
end subroutine finalize


subroutine fit_emf(temp)
    implicit none
    real,   intent(in) :: temp

    real, parameter :: nernst_base = 0.08617333262 ! mV/K
    integer NX, NV

    ! handle temperature
    nernst = nernst_base*temp  ! mV

    ! create constraints
    call count_dof(NX, NV)
    allocate(xvalues(NX))
    allocate(variables(NV))
    call create_variables()

    call lm()

    ! TODO post processing

    ! deallocation
    DEALLOCATE(xvalues, variables)
end subroutine fit_emf


!----------------------------------------------------------------------------------------------------
! Run the Levenberg-Marquardt loop.
!----------------------------------------------------------------------------------------------------
subroutine lm()
    ! routine for Levenberg-Marquardt iterations.
    implicit none

    integer :: NX ! number of variables to optimize
    integer iter  ! iteration number
    integer, parameter :: MAX_ITERS = 100
    real,    parameter :: DAMPING_UPPER = 1e10
    real,    parameter :: DAMPING_LOWER = 1e-12
    real,    parameter :: DAMPING_UPF = 8.0
    real,    parameter :: DAMPING_LOWF = 3.0
    double precision, allocatable, dimension(:,:) :: jac
    double precision, allocatable, dimension(:)   :: dx
    double precision, allocatable, dimension(:)   :: fvec, trial_fvec
    double precision, allocatable, dimension(:)   :: gradient
    double precision, allocatable, dimension(:,:) :: M, dampedM
    double precision :: gradient_norm, chisq, damping
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

    dx = 0.0
    jac = 0.0
    fvec = 0.0
    damping = 1D2

    main_loop: do iter=0,MAX_ITERS
        call try_step(dx)

        call update_equil()
        call calc_resid(fvec)
        call jacobian_lm(jac)

        M = MATMUL(MATMUL(TRANSPOSE(jac), W), jac)
        gradient = 2*MATMUL(MATMUL(TRANSPOSE(jac),W),fvec)
        gradient_norm = NORM2(gradient)
        chisq = DOT_PRODUCT(fvec, MATMUL(W,fvec))

        dampedM = M
        do i=1,NX
            dampedM(i,i) = (1.0 + damping) * M(i,i)
        end do

        if (iter .ne. 0) then
            call lesolve(dampedM, -fvec, dx, NX, CTRL)
            if (CTRL .eq. 99) then
                emf_diagnostic = "jacobian could not be inverted"
                return
            end if
        end if

    end do main_loop

    DEALLOCATE(jac, dx, fvec, M)
end subroutine lm

!----------------------------------------------------------------------------------------------------
! Update the free concentrations.
!----------------------------------------------------------------------------------------------------
subroutine update_equil()
    integer i, blow, bhigh, CTRL
    do i=1,NTT
        if (update_titr(i)) then
            blow = 1+sum(titr_size(1:(i-1)))  ! this should work also for i==1
            bhigh = blow + titr_size(i)
            call calc_analc(V(blow:bhigh), V0(i), use_init(i,:), use_buret(i,:), T(blow:bhigh,:))
            call eq_000(C(blow:bhigh,:), T(blow:bhigh,:), EXP(use_lbeta), P, CTRL, titr_size(i), NCP, NEQ)
        end if
    end do
end subroutine update_equil


subroutine try_step(step)
    double precision, intent(in), dimension(:) :: step

    integer i

    ! update tmp values
    do i=1,size(variables)
        variables(i)%destination_value = real(variables(i)%ratio*(variables(i)%origin_value + step(i)))
    end do
end subroutine try_step


subroutine calc_resid(fvec)
    implicit none
    double precision, intent(out), DIMENSION(:) :: fvec

    integer i
    integer lb, ub
    integer hact

    l_titr: do i=1,NTT
        hact = eactive(i)
        call titr_bounds(i, lb, ub)
        fvec(lb:ub) = E(lb:ub) - (use_emf0(i) - use_slope(i) * nernst * LOG(C(lb:ub,hact)))
    end do l_titr

end subroutine calc_resid


subroutine jacobian_lm(jac)
    ! Compute the jacobian for the Levenberg-Marquardt iteration.
    ! The rows are the experimental points; the columns are the variables to be optimized in the
    ! following order: (1) equilibrium constants; (2) initial amounts; (3) buret concentrations;
    ! (4) standard potentials; (5) nernstian slopes; (6) constraints
    double precision, intent(out), dimension(:,:) :: jac

    double precision, DIMENSION(NCP, NCP) :: Am      ! the A-matrix
    double precision, dimension(NCP, NEQ) :: dlcdlb  ! the result
    double precision, allocatable, dimension(:,:) :: bvec
    integer i, j, k, CTRL
    integer curr_titr
    integer aux, eac
    integer bsize
    logical do_beta, do_init, do_buret, do_amatrix
    integer, DIMENSION(5,size(variables)) :: filter_jtype

    jac = 0.0
    do_beta = allocated(tmp_lbeta)
    if (do_beta) bsize = bsize + NEQ
    do_init = allocated(tmp_init)
    do_buret = allocated(tmp_buret)
    do_amatrix = do_beta .or. do_init .or. do_buret

    filter_jtype = 0
    do i=1,5
        do j=1,size(variables)
            if (variables(j)%jtype .eq. i) filter_jtype(i,j) = j
        end do
    end do

    bsize = 0
    do i=1,size(variables)
        if (variables(i)%jtype < 3) bsize = bsize + 1
    end do
    allocate(bvec(NCP,bsize))

    curr_titr = 1       ! track which titration are we currently in

    do i=1,NPT
        if (do_amatrix) then
            call amatrix(Am, C(i,:), P, NCP, NEQ)
            bsize = 0

            do j=1,size(variables)
                if (j > 4) exit
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
                case (J_LBETA:J_TBURE)
                    jac(i, aux) = jac(i, aux) - variables(j)%ratio*use_slope(curr_titr)*nernst*bvec(eac,j)
                case (J_EMF0)
                    if (curr_titr .eq. variables(j)%destination_index1) then
                        jac(i, aux) = jac(i, aux) - variables(j)%ratio
                    end if
                case (J_SLOPE)
                    if (curr_titr .eq. variables(j)%destination_index1) then
                        jac(i, aux) = jac(i, aux) - variables(j)%ratio*nernst
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
        use_lbeta => tmp_lbeta
        call aux_create_variables(tmp_lbeta, flag_beta, xcount, vcount, iref, vref, J_LBETA, 0)
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
        tmp_emf0 = E0
        use_emf0 => tmp_emf0
        call aux_create_variables(tmp_emf0, flag_emf0, xcount, vcount, iref, vref, J_EMF0, 0)
    else
        use_emf0 => E0
    end if

    ! V. count nernstian slopes
    if (sum(flag_slope) > 0) then        ! something to refine here
        allocate(tmp_slope(NTT))
        tmp_slope = slope
        call aux_create_variables(tmp_slope, flag_slope, xcount, vcount, iref, vref, J_SLOPE, 0)
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
            variables(vcount)%origin_value => xvalues(xcount)
            variables(vcount)%origin_index = xcount
            variables(vcount)%destination_value => vector(i)
            variables(vcount)%ratio = 1.0
        end if
        if (flags(i) > 1) then
            vcount = vcount + 1
            if (iref(flags(i)) .eq. 0) then   ! constraint, first time
                xcount = xcount + 1
                variables(vcount)%origin_value => xvalues(xcount)
                variables(vcount)%destination_value => vector(i)
                variables(vcount)%ratio = 1.0
                variables(vcount)%origin_index = xcount

                iref(flags(i)) = xcount
                vref(flags(i)) = vector(i)
            else                                  ! constraint already set
                variables(vcount)%origin_value => xvalues(iref(flags(i)))
                variables(vcount)%destination_value => vector(i)
                variables(vcount)%ratio = vref(flags(i))/vector(i)
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
            if (constr(i)) then
                NX = NX + 1
                constr(i) = .true.
            end if
        end if
    end do
end subroutine sub_count_dof


subroutine calc_analc(titre, init_vol, init, buret, analc)
    ! Calculate the analytical concentration, given the titration paramers, 
    real, intent(in), dimension(:)  :: titre
    real, intent(in), dimension(:)  :: init
    real, intent(in), dimension(:)  :: buret
    real, intent(out), dimension(:,:)  :: analc
    real, INTENT(IN) :: init_vol

    integer i, j
    do i=1,size(titre)
        do j=1,NCP
            analc(i,j) = (init(j) + buret(j)*titre(i)) / (init_vol+titre(i))
        end do
    end do
end subroutine calc_analc


!------------------------------------------------------------------------------
! Compute numerical first derivative dY/dX of a 1-D array
! using second-order finite differences.
!
! Interior points  : central difference (O(h^2))
! Left boundary    : forward three-point stencil (O(h^2))
! Right boundary   : backward three-point stencil (O(h^2))
!
! The grid X need not be uniform.
!------------------------------------------------------------------------------
subroutine derivative_1d(n, x, y, dy)
  implicit none
  integer,          intent(in)  :: n          ! number of points (must be >= 3)
  real,     intent(in)  :: x(n)       ! abscissae (strictly increasing)
  real,     intent(in)  :: y(n)       ! function values
  real,     intent(out) :: dy(n)      ! derivative dy/dx

  integer :: i
  real(kind=8) :: h1, h2, p1, p2, p3

  if (n < 3) then
     write(*,*) 'derivative_1d: n must be >= 3'
     stop
  end if

  ! ---- left boundary (forward 2nd-order) ------------------------------------
  h1 = x(2) - x(1)
  h2 = x(3) - x(2)
  p1  = -(2.0d0*h1 + h2)/(h1*(h1+h2))
  p2  =  (h1 + h2)/(h1*h2)
  p3  = -h1/(h2*(h1+h2))
  dy(1) = real(p1*y(1) + p2*y(2) + p3*y(3))

  ! ---- interior points (central 2nd-order) ----------------------------------
  do i = 2, n-1
     h1 = x(i)   - x(i-1)
     h2 = x(i+1) - x(i)
     p1  = -h2 / (h1*(h1+h2))
     p2  =  (h2 - h1) / (h1*h2)
     p3  =  h1 / (h2*(h1+h2))
     dy(i) = real(p1*y(i-1) + p2*y(i) + p3*y(i+1))
  end do

  ! ---- right boundary (backward 2nd-order) ----------------------------------
  h1 = x(n-1) - x(n-2)
  h2 = x(n)   - x(n-1)
  p1  =  h2 / (h1*(h1+h2))
  p2  = -(h1 + h2)/(h1*h2)
  p3  =  (2.0d0*h2 + h1)/(h2*(h1+h2))
  dy(n) = real(p1*y(n-2) + p2*y(n-1) + p3*y(n))

end subroutine derivative_1d


subroutine titr_bounds(i, low, high)
    integer, intent(in) :: i
    integer, intent(out) :: low, high
    low = 1+sum(titr_size(1:(i-1)))  ! this should work also for i==1
    high = low + titr_size(i)
end subroutine titr_bounds


subroutine msysolve(A, b, CTRL)
    ! Solve lineal equations using LAPACK when A is symmetric. This is always the case when
    ! solids are not present. This routine requires
    ! the -llapack flag at compilation and, probably also -lblas in some systems.
    implicit none
    double precision, dimension(:,:), intent(in) :: A   ! The nxn array of coefficients
    double precision, dimension(:,:), intent(in) :: b   ! the n array of independent terms
    integer, intent(inout) :: CTRL                      ! control value. 0 on no error.
  
    ! For DGELSY 
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: nrhs, lda, ldb, lwork, info
    double precision, allocatable, dimension(:) :: work
    integer n

    ! External procedures defined in LAPACK
    external dsysv

    n = size(A, dim=1)
    nrhs = size(b, dim=2)
    lda = n
    ldb = n
    lwork = -1              ! workspace query first

    lwork = -1
    allocate(work(1))
    call dsysv("U", n, 1, A, lda, ipiv, b, ldb, work, lwork, info)
    lwork = nint(work(1))
    deallocate(work)

    allocate(work(lwork))
    call dsysv("U", n, 1, A, lda, ipiv, b, ldb, work, lwork, info)
    deallocate(work)

    if (info /= 0) then
        CTRL = 99    ! unable to invert matrix
        return
    end if
end subroutine msysolve

end module emf
