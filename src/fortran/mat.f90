! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>
!
!----------------------------------------------------------------------------------------------------
! Routines for fitting equilibrium constants and other parameters from experiemntal data in which
! a linearity dependence can be assumed: Y = AX; where Y is the experimental data (absorbance,
! chemical shift, emission or whatever) and X is the free concentrations array.
!
! HOW TO USE
! 1. call init_module with the parameters indicated
! 2. load the titrations one by one by calling load_titration
! 3. Optional additions
!   3a. set the weighting scheme by calling one of the set weights routine
!   3b. include variable ionic strength is aplicable with include_ionic
! 4. call fit_emf
! 5. copy the result from the public variables
! 6. call finalize to deallocate and cleanup
!----------------------------------------------------------------------------------------------------


module mat
    implicit none (type, external)

    private

    integer, public, parameter :: RF_CONSTANT = 0
    integer, public, parameter :: RF_REFINE = 1
    integer, public, parameter :: AC_TRANSPARENT = 0
    integer, public, parameter :: AC_ACTIVE = 1
    integer, public, parameter :: LE_GENERAL = 0    ! linear elimination allows positive and negative values
    integer, public, parameter :: LE_POSITIVE = 1   ! linear elimination allows only positive values

    integer :: NCP        ! number of components
    integer :: NEQ        ! number of equilibria
    integer :: NPT        ! number of titration points
    integer :: NHD        ! number of horizontal data
    integer :: NTT        ! number of titrations
    double precision,  allocatable :: analc(:,:)        ! Analytical concentration, in mmol/mL, dim NPT x NCP
    real,              allocatable :: absorbance(:,:)   ! NPT x NWL
    double precision,  allocatable, target :: extc(:,:)  ! NSP x NWL
    double precision,  allocatable :: conc(:,:)         ! NPT, NSP
    character(len=12), allocatable :: cmpnames(:)
    integer,           allocatable :: stoich(:,:)
    double precision,  allocatable, target :: lbeta(:)
    integer,           allocatable :: flag_lbeta(:)
    integer,           allocatable :: active(:)
    double precision,  allocatable :: covariance(:,:)  ! the covariance matrix
    double precision,  allocatable :: correlation(:,:)  ! the correlation matrix
    double precision :: weight_params(2)
    integer, allocatable :: titr_size(:)    ! the size in experimental points per titration

    ! variables for constraints
    type value_holder
        double precision, pointer :: origin_value => null()  ! points to value in xvalues
        integer :: origin_index                              ! index of xvalues
        double precision, pointer :: destination_value => null()         ! points to value in tmp_*
        integer :: destination_index1                        ! index of EQ, CP or NT
        integer :: destination_index2                        ! index of EQ, CP or NT
        integer :: id
        integer :: jtype                                     ! type of variable (beta, epsilon)
        double precision :: ratio                            ! 1.0 for free and reference; other for the rest
    end type value_holder

    type(value_holder), allocatable, dimension(:) :: variables
    double precision,   allocatable, target, dimension(:) :: xvalues
    double precision, allocatable, dimension(:) :: weights
    character(len=45), public :: msg
    integer :: outu
    double precision, parameter :: LG10=LOG(10.0)

    integer :: LOGLEVEL  ! set verbosity level; 0 -> mute output; 1 -> minimal output; 2 -> verbose output; 3->pedantic

    integer, parameter :: MAXITS = 20
    logical :: reweight_flag
    integer :: linear_elimination_type
    integer :: current_titration

    real, public :: stat_ssr   ! sum of squared residuals
    real, public :: stat_rms   ! root-mean-square of residuals
    integer, public :: stat_dof    ! degrees of freedom
    integer, public :: stat_nobs    ! number of observations

    real, allocatable, public :: out_lbeta(:)
    real, allocatable, public :: out_errlbeta(:)

    public init_module, finalize, lm, load_titration, fit
contains

!----------------------------------------------------------------------------------------------------
! Initialize module, allocate varables
!----------------------------------------------------------------------------------------------------
subroutine init_module(n_points, n_horiz, n_titrs, model, logbeta, flgbeta, activity, letype, verbosity, outn)
    use, intrinsic :: iso_fortran_env, only: output_unit
    integer, intent(in) :: n_points                  ! total number of experimental points
    integer, intent(in) :: n_titrs                   ! number of titrations to be loaded
    integer, intent(in) :: n_horiz                   ! number of horizontal data
    integer, intent(in) :: model(:,:)                ! the stoichiometry array
    real,    intent(in) :: logbeta(size(model, 1))   ! the log10 of the formation constants
    integer, intent(in) :: flgbeta(size(model, 1))   ! the refinement flags for the formation constants
    integer, intent(in) :: activity(:)
    integer, intent(in) :: verbosity                 ! level of verbosity
    integer, intent(in) :: letype                    ! type of linear elimination
    integer, intent(in), optional :: outn            ! the number for output unit

    integer :: i

    LOGLEVEL = verbosity

    if (present(outn)) then
        outu = outn
    else
        outu = output_unit
    end if
    if (LOGLEVEL>1) write (outu, "(A,/)") "Initiating module"

    if (letype == 0 .or. letype == 1) then
        linear_elimination_type = letype
    else
        msg = "letype parameter not allowed"
        return
    end if

    NPT = n_points
    NTT = n_titrs
    NHD = n_horiz
    NEQ = size(model, dim=1)
    NCP = size(model, dim=2)

    if (LOGLEVEL==3) write (outu, "(A,I4)") " .. NPT=", NPT
    if (LOGLEVEL==3) write (outu, "(A,I4)") " .. NTT=", NTT
    if (LOGLEVEL==3) write (outu, "(A,I4)") " .. NHD=", NHD
    if (LOGLEVEL==3) write (outu, "(A,I4)") " .. NEQ=", NEQ
    if (LOGLEVEL==3) write (outu, "(A,I4)") " .. NCP=", NCP

    allocate(weights(NPT*NHD))     ! variable allocation
    allocate(conc(NPT,NCP+NEQ))
    allocate(titr_size(NTT))
    allocate(analc(NPT, NCP))
    allocate(absorbance(NPT,NHD))
    allocate(lbeta(NEQ))
    allocate(flag_lbeta(NEQ))
    allocate(extc(NCP+NEQ, NHD))
    allocate(active(NCP+NEQ))
    allocate(stoich(NEQ,NCP))

    weights = 1.0               ! start with unit weights
    reweight_flag = .false.

    conc = 0.0
    extc = 0.0
    stoich = model

    titr_size = 0
    analc = 0.0
    absorbance = 0.0

    active = activity

    flag_lbeta = flgbeta

    if (LOGLEVEL==3) then
        write (outu, "(A)") " .. model to refine is"
        do i=1,NEQ
            write (outu, "(5X,I4)", advance='no') i
            write (outu, "(*(I4))", advance='no') model(i,:)
            write (outu, "(F12.4)", advance='no') logbeta(i)
            if (flgbeta(i) == 1) then
                write (outu, "(A)") " REFINE"
            else
                write (outu, "(/)")
            end if
        end do
    end if


    lbeta = logbeta * LG10
    current_titration = 0
end subroutine init_module


subroutine load_titration(titdata, titanalc)
    real, intent(in) :: titdata(:,:)
    real, intent(in) :: titanalc(:,:)

    integer :: low, high, dsize

    current_titration = current_titration + 1
    dsize = size(titdata, dim=1)
    titr_size(current_titration) = dsize

    call titr_bounds(current_titration, low, high)

    if (size(titdata, dim=2) /= NHD) then
        msg="ERROR: titdata has wrong dimension"
        if (LOGLEVEL>1) then
            write (outu, "(A,2(A,I4))") msg, " got", size(titdata, dim=2), " expected", NHD
        end if
        return
    end if

    absorbance(low:high,:) = titdata
    analc(low:high,:) = titanalc

    if (LOGLEVEL>1) write (outu, "(A)") " .. titration loaded"
end subroutine load_titration


subroutine set_weights_const(w)
    real, intent(in) :: w
    weights = w
end subroutine set_weights_const


subroutine fit()
    integer :: i, j

    do i=1,NCP
        conc(:,i) = analc(:,i) / 10.0
    end do
    where (conc < 0.0) conc=1D-12

    do i=1,NCP+NEQ
        if (active(i) < 0) then
            j = -active(i)
            extc(j,:) = absorbance(j,:) / analc(j,:)
        end if
    end do

    call lm()

    stat_nobs = NPT*NHD
    stat_rms = REAL(SQRT(stat_ssr/stat_nobs))

    allocate(out_lbeta(NEQ), out_errlbeta(NEQ))

    out_lbeta = REAL(lbeta / LG10)
    j = 0
    do i = 1,NEQ
        if (flag_lbeta(i) == 0) then
            out_errlbeta(i) = 0.0
        else
            j = j + 1
            out_errlbeta(i) = REAL(sqrt(covariance(j,j)))
        end if
    end do

end subroutine fit


!----------------------------------------------------------------------------------------------------
! Clear all variables and deallocate memory.
!----------------------------------------------------------------------------------------------------
subroutine finalize()
    if (LOGLEVEL>1) write (outu, "(/,A)") "Memory cleanup"

    if (allocated(analc))        deallocate(analc)
    if (allocated(cmpnames))     deallocate(cmpnames)
    if (allocated(extc))         deallocate(extc)
    if (allocated(active))       deallocate(active)
    if (allocated(weights))      deallocate(weights)
    if (allocated(conc))         deallocate(conc)
    if (allocated(titr_size))    deallocate(titr_size)
    if (allocated(absorbance))   deallocate(absorbance)
    if (allocated(lbeta))        deallocate(lbeta)
    if (allocated(flag_lbeta))   deallocate(flag_lbeta)
    if (allocated(stoich))       deallocate(stoich)
    if (allocated(out_lbeta))    deallocate(out_lbeta)
    if (allocated(out_errlbeta)) deallocate(out_errlbeta)
    if (allocated(variables))    deallocate(variables)
    if (allocated(xvalues))      deallocate(xvalues)
    if (allocated(covariance))   deallocate(covariance)
    if (allocated(correlation))  deallocate(correlation)

    if (LOGLEVEL>1) write (outu,"(A)") "Module terminated"
end subroutine finalize


!----------------------------------------------------------------------------------------------------
! Performe Levenberg-Marquardt iterations.
!----------------------------------------------------------------------------------------------------
subroutine lm()
    use stat, only: fit_sigma, calc_covariance, calc_correlation
    use linalg, only: las_cond

    double precision, allocatable :: jac(:,:)
    double precision, allocatable :: dx(:)
    double precision, allocatable :: fvec(:), trial_fvec(:)
    double precision, allocatable :: M(:,:), dampedM(:,:), diagM(:)
    double precision, allocatable :: resid(:,:)
    double precision, allocatable :: gradient(:)
    double precision :: gradient_norm, sumsq, damping, trial_sumsq
    double precision :: predicted, actual, rho, cond
    double precision :: relative_step
    double precision :: relative_grad
    double precision :: relative_sumsq
    double precision, parameter :: THRESHOLD_STEP = 1d-4
    double precision, parameter :: THRESHOLD_SSQ = 1d-3
    double precision, parameter :: THRESHOLD_GRAD = 1d-4
    double precision, parameter :: MAX_STEP = 1.0
    real, parameter :: DAMPING_UPPER = 1e10
    real, parameter :: DAMPING_LOWER = 1e-12
    logical :: exit_condition

    integer :: i, CTRL
    integer :: NX
    integer :: iter
    integer :: n_active

    n_active = sum(active)

    call create_variables(NX)

    if (LOGLEVEL>1) then
        write (outu, "(A)") " .. Allocating variables"
        write (outu, "(3(A,I5),/)") " .. NPT=", NPT, ", NWL=", NHD, ", NX=", NX
    end if

    allocate(resid(NPT,NHD))

    allocate(jac(NPT*NHD,NX))
    allocate(fvec(NPT*NHD))
    allocate(trial_fvec(NPT*NHD))

    allocate(M(NX,NX), dampedM(NX,NX), gradient(NX), diagM(NX))
    allocate(dx(NX))

    dx = 0.0
    jac = 0.0
    fvec = 0.0
    damping = 1D-2
    sumsq = HUGE(1D1)
    trial_sumsq = 0.0
    trial_fvec = 0.0

    if (LOGLEVEL>1) write(outu, "(A)") "Starting Levenberg-Marquardt iterations"

    call update_equil()                           ! update free concentrations
    call solve_extc_nnls(ctrl)
    call calc_resid(resid)                        ! calculate residuals
    if (reweight_flag) call reweight(resid, absorbance, .false.)
    resid = absorbance - matmul(conc, extc)
    fvec = reshape(resid, [ NPT*NHD ])
    sumsq = SUM(weights*fvec**2)
    call take_marquardt(damping, jac, M, dampedM)
    gradient = -MATMUL(TRANSPOSE(jac), weights*fvec)

    if (LOGLEVEL>1) write(outu,"(A,E14.6)") " ..  initial sumsq = ", sumsq

    iter = 1
    main_loop: do while (iter < maxits)
        call las_cond(dampedM, gradient, dx, cond, ctrl)
        where (abs(dx) > MAX_STEP) dx = SIGN(MAX_STEP, dx)     ! limit step

        call try_step(dx)                             ! update parameters to new step
        call update_equil()                           ! update free concentrations
        !call solve_extc()                             ! update extinction coefficients
        call solve_extc_nnls(ctrl)
        call calc_resid(resid)                        ! calculate residuals

        trial_fvec = reshape(resid, [ NPT*NHD ])
        trial_sumsq = SUM(weights*trial_fvec**2)

        ! after take_marquardt, store the diagonal once
        do i = 1, NX
            diagM(i) = M(i,i)
        end do

        ! \ldots later, after the solve that produced dx
        predicted = DOT_PRODUCT(gradient, dx) - 0.5d0 * DOT_PRODUCT(dx, MATMUL(M, dx))
        actual     = sumsq - trial_sumsq          ! positive when the step is good

        if (predicted > 100.0d0 * EPSILON(1.0d0) * (1.0d0 + abs(sumsq))) then
            rho = actual / predicted
        else if (predicted > 0.0 .and. actual > 0.0) then
            rho = 0.5
        else
            rho = -1.0d0          ! model is not decreasing → force rejection
        end if

        if (rho > 1D-3) then
            ! step accepted

            if (LOGLEVEL>1) then
                write(outu, "(/,A,I4,/)") "Iteration", iter
                write(outu, "(A)") "   #    previous        step         new"
                do i=1,NX
                    write(outu, "(I4,3F12.4)") i, xvalues(i)/LG10, dx(i)/LG10, (xvalues(i)+dx(i))/LG10
                end do
            end if

            xvalues = xvalues + dx

            if (reweight_flag .and. iter>10) call reweight(resid, absorbance, .true.)
            call update_equil()                           ! update free concentrations
            call solve_extc_nnls(ctrl)
            call calc_resid(resid)                        ! calculate residuals
            fvec = reshape(resid, [ NPT*NHD ])
            call take_marquardt(damping, jac, M, dampedM)
            gradient = -MATMUL(TRANSPOSE(jac), weights*fvec)
            sumsq = SUM(weights*fvec**2)
            gradient_norm = NORM2(gradient)

            if (LOGLEVEL>2) then
                write(outu, "(/,A8,E14.6)") "   sumsq", sumsq
                write(outu, "(A8,2E14.4)")  " weights", weight_params
                write(outu, "(A8,E14.6)")   " damping", damping
                write(outu, "(A8,E14.6,2(A,E12.4),A)")   "     rho", rho, " [actual=", actual, ", predicted=", predicted, "]"
                write(outu, "(A8,E14.6)")   "    cond", cond
            end if
            if (LOGLEVEL==1) write (outu, "(A)", advance='no') '.'

            ! Exit check
            relative_step = NORM2(dx) / (1.0d0 + NORM2(xvalues))
            relative_grad = MAXVAL(ABS(gradient)) / (1.0d0 + ABS(sumsq))
            relative_sumsq = ABS(actual) / (1.0d0 + sumsq)

            exit_condition = (relative_sumsq < THRESHOLD_SSQ) .and. (relative_step < THRESHOLD_STEP)

            if (LOGLEVEL > 1) then
                write(outu, "(A8,E14.6)", advance="no") " descent", relative_sumsq
                if (relative_sumsq < THRESHOLD_SSQ) then
                    write(outu, "(A)") " YES"
                else
                    write(outu, "(A)") " NO"
                end if

                write(outu, "(A,E14.6)", advance="no") "  rstep ", relative_step
                if (relative_step < THRESHOLD_STEP) then
                    write(outu, "(A)") " YES"
                else
                    write(outu, "(A)") " NO"
                end if

                write(outu, "(A,E14.6)", advance="no") "  rgradn", relative_grad
                if (relative_grad < THRESHOLD_GRAD) then
                    write(outu, "(A)") " YES"
                else
                    write(outu, "(A)") " NO"
                end if
            end if

            if (rho > 0.75d0) then
                damping = max(damping / 5.0d0, DAMPING_LOWER)   ! strong decrease
            else if (rho > 0.25d0) then
                damping = max(damping / 2.0d0, DAMPING_LOWER)   ! mild decrease
            else
                damping = max(damping / 1.1d0, DAMPING_LOWER)   ! mild decrease
            end if

            if (exit_condition) exit main_loop

            iter = iter + 1
        else

            damping = damping * 2.0d0   ! increase
            if (damping >= DAMPING_UPPER * 0.99d0) then
                msg = "LM stalled at maximum damping"
                if (LOGLEVEL>1) write (outu, "(A)") msg
                exit main_loop
            end if

            dampedM = M
            do i=1,NX
                dampedM(i,i) = (1.0 + damping) * M(i,i)
            end do
        end if
    end do main_loop

    if (iter == MAXITS) then
        msg = "iteration limit reached"
    else
        msg = "normal exit"
    end if

    allocate(covariance(size(M,1), size(M,2)))
    allocate(correlation(size(M,1), size(M,2)))
    call calc_covariance(M, covariance, msg)
    covariance = covariance * sumsq / (NPT*NHD-NX)
    call calc_correlation(covariance, correlation)
    stat_ssr = REAL(sumsq)
    stat_dof = NPT*NHD - NX

    deallocate(jac, dx, fvec, M)
end subroutine lm


subroutine take_marquardt(damping, jac, M, dampedM)
    use misc, only: lesolve
    double precision, intent(in)  :: damping
    double precision, intent(out) :: jac(:,:)
    double precision, intent(out) :: M(:,:), dampedM(:,:)

    integer :: i, j, NX
    NX = size(jac, dim=2)

    call calc_jacobian(jac)
    M = 0.0
    do i = 1,NX
        do j = 1,NX
            M(i,j) = M(i,j) + SUM(jac(:,i)*weights*jac(:,j))
        end do
    end do
    dampedM = M
    do i=1,NX
        dampedM(i,i) = (1.0 + damping) * M(i,i)
    end do
end subroutine take_marquardt


subroutine update_equil()
    use eq, only: eq_000
    integer :: CTRL

    call eq_000(conc, analc, EXP(lbeta), stoich, CTRL, NPT, NCP, NEQ)
end subroutine update_equil


subroutine try_step(step)
    double precision, intent(in), dimension(:) :: step

    integer :: i

    ! update tmp values
    do i=1,size(variables)
        variables(i)%destination_value = real(variables(i)%ratio*(variables(i)%origin_value + step(i)))
    end do
end subroutine try_step


subroutine create_variables(NX)
    integer, intent(out) :: NX
    integer :: i, xcount

    NX = COUNT(flag_lbeta==1)

    allocate(variables(NX))
    allocate(xvalues(NX))
    xcount = 0

    do i = 1,size(lbeta)
        if (flag_lbeta(i) == 1) then
            xcount = xcount + 1
            xvalues(xcount) = lbeta(i)
            variables(xcount)%origin_value => xvalues(xcount)
            variables(xcount)%origin_index = xcount
            variables(xcount)%destination_value => lbeta(i)
            variables(xcount)%destination_index1 = i
            variables(xcount)%jtype = 0
            variables(xcount)%ratio = 1.0
        end if
    end do
end subroutine create_variables


subroutine calc_resid(resid)
    double precision, intent(out), dimension(:,:)   :: resid

    resid = absorbance - MATMUL(conc, extc)
end subroutine calc_resid


subroutine calc_jacobian(jac)
    use misc, only: msysolve
    use err, only: calc_dlogcdlogb, calc_extdlogcdlogb

    double precision, intent(out), dimension(:,:) :: jac

    double precision, dimension(NCP+NEQ,NEQ,NPT) :: dc
    double precision, allocatable :: tmp_vec(:)
    double precision :: G(NPT,NHD)
    double precision, allocatable :: tmp(:)

    integer :: k, wl, pt, sp
    integer :: indexsp, nfree

    nfree = COUNT(active == AC_ACTIVE)
    allocate(tmp_vec(nfree))

    l_point: do k=1,NPT
        call calc_dlogcdlogb(conc(k,:), stoich, NCP, NEQ, dc(1:NCP,:,k))
        call calc_extdlogcdlogb(dc(:,:,k), stoich, NCP, NEQ)
    end do l_point

    ! convert dlogC/dlogB in dC/dlogB
    do k=1,NEQ
        dc(:,k,:) = dc(:,k,:) * transpose(conc)
    end do

    allocate(tmp(nfree))
    do k = 1, size(jac, dim=2)
       indexsp = variables(k)%destination_index1
       G = 0.0
       do wl = 1, NHD
          do pt = 1, NPT
             do sp = 1, NCP+NEQ
                G(pt,wl) = G(pt,wl) + dc(sp,indexsp,pt) * extc(sp,wl)
             end do
          end do
       end do

       jac(:,k) = -reshape(G, [NPT*NHD])

    end do
end subroutine calc_jacobian


subroutine reweight(resids, aobs, do_fit)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    double precision, intent(in) :: resids(:,:)
    real, intent(in) :: aobs(:,:)
    logical, intent(in) :: do_fit

    double precision, allocatable :: y(:)  ! residuals squared
    double precision, allocatable :: A2(:) ! absorbance squared
    double precision, allocatable :: w(:)  ! secondary weights
    double precision, allocatable :: X(:,:)  ! design matrix

    integer, parameter :: MAX_ITS = 10
    integer :: i, lwork
    integer :: n  ! number of data points
    double precision  :: work_query(1)
    double precision, allocatable :: work(:), bvec(:)
    integer :: iinfo
    double precision :: sigma2
    double precision, parameter :: TOL = 1D-5
    double precision, parameter :: MIN_A2 = 1D-12

    external dgels

    n = NHD*NPT
    allocate(A2(n))
    A2 = reshape(aobs, [n])**2

    if (do_fit) then
        allocate(y(n),  w(n), X(n,2), bvec(n))

        ! initial guess
        !weight_params = [1D-6, 0.0_dp]
        y = reshape(resids, [n])**2

        X(:,1) = 1.0_dp
        X(:,2) = A2
        bvec = y

        !----- query optimal work-space size for DGELS ----------------------------
        call dgels("N", n, 2, 1, X, n, bvec, n, work_query, -1, iinfo)
        lwork = int(work_query(1))
        allocate(work(lwork))
        call dgels("N", n, 2, 1, X, n, bvec, n, work, lwork, iinfo)

        weight_params(1) = SQRT(max(bvec(1), MIN_A2))
        weight_params(2) = SQRT(max(bvec(2), 0.0))
    end if

    ! recalculate weights
    do i=1,n
        sigma2 = weight_params(1)**2 + weight_params(2)**2 * A2(i)
        sigma2 = MAX(sigma2, 1.0D-30)   ! avoid division by zero
        weights(i) = 1.0D0 / sigma2
    end do
    deallocate(A2)
    if (do_fit) deallocate(y,  w, X, bvec, work)
end subroutine reweight


subroutine solve_extc()
    integer:: info
    double precision, allocatable, dimension(:,:) :: Acorr
    double precision, allocatable, dimension(:,:) :: Cfree
    double precision, allocatable, dimension(:,:) :: B
    integer :: nfree
    integer :: lda, ldb, lwork
    double precision, allocatable, dimension(:) :: work
    double precision, allocatable :: tau(:)

    integer :: i, j

    external dgeqrf
    external dormqr
    external dtrtrs
    external dorgqr

    nfree = COUNT(active == AC_ACTIVE)   ! number of species whose extinction coefficient is to be refined
    if (nfree == 0) then                     ! nothing to refine
       info = 0
       return
    end if

    allocate(Acorr(NPT, NHD))
    Acorr = absorbance
    l_species: do i=1,(NCP+NEQ)
        if (active(i) < 0) then
            do j=1,NHD
                Acorr(:,j) = Acorr(:,j) - conc(:,i) * extc(i,j)
            end do
        end if
    end do l_species

    allocate(Cfree(NPT,nfree))
    j = 0
    do i=1,NCP+NEQ
        if (active(i) == AC_ACTIVE) then
            j = j + 1
            Cfree(:,j) = conc(:,i)
        end if
    end do

    lda = NPT
    ldb = MAX(NPT, nfree)
    allocate(B(ldb, NHD))
    B(1:NPT,:) = Acorr
    if (NPT < nfree) B(NPT+1:ldb,:) = 0.0

    call workspace_queries(lwork, nfree, Cfree, B, lda, ldb)
    allocate(tau(nfree))
    allocate(work(lwork))

    ! QR factorisation
    call dgeqrf(NPT, nfree, Cfree, lda, tau, work, lwork, info)

    ! Qᵀ * B
    call dormqr("L","T", NPT, NHD, nfree, Cfree, lda, tau, B, ldb, work, lwork, info)

    ! R \ (Qᵀ B)
    call dtrtrs("U","N","N", nfree, NHD, Cfree, lda, B, ldb, info)

    ! (optional) form explicit thin Q for the Jacobian projection
    call dorgqr(NPT, nfree, nfree, Cfree, lda, tau, work, lwork, info)
    ! now Cfree contains Q

    if (info==0) then
        j = 0
        do i=1,NCP+NEQ
            if (active(i) == AC_ACTIVE) then
                j = j + 1
                extc(i,:) = B(j,:)
            end if
        end do
    end if
    where (extc < 0.0) extc = 0.0  ! unallow negative extinction coefficients

    !Q = Cfree
    deallocate(Cfree,B,work,Acorr)
end subroutine solve_extc


!==============================================================================
! Non-negative least squares for extinction coefficients
!
!   absorbance(NPT,NWL)  ≈  conc(NPT,NSP) * extc(NSP,NWL)
!   subject to  extc >= 0
!
! Uses Lawson-Hanson active-set + LAPACK DGELSY
!==============================================================================
subroutine solve_extc_nnls(info)
  use, intrinsic :: iso_fortran_env, only: dp => real64

  integer,  intent(out)   :: info          ! 0 = success, <0 = illegal arg, >0 = rank issues
  external dgelsy

  ! ---- local ---------------------------------------------------------------
  integer  :: jwl, i, j, k, rank, lwork, nfree, nbound, iter
  integer  :: maxit
  integer, allocatable :: idx(:), jpvt(:), free_idx(:), bound_idx(:)
  real(dp), allocatable :: A(:,:), b(:), x(:), work(:), dual(:), resid(:)
  real(dp) :: rcond, normb, tol, max_dual
  real(dp), parameter :: eps = epsilon(1.0_dp)
  logical  :: converged
  integer :: nsp

  nsp = NCP + NEQ

  info = 0
  if (npt < 1 .or. nsp < 1 .or. NHD < 1) then
     info = -1
     return
  end if

  maxit  = 3 * nsp          ! more than enough for spectroscopy
  rcond  = 1.0e-12_dp       ! relative rank tolerance for DGELSY
  tol    = 10.0_dp * eps    ! feasibility tolerance (relative)

  ! workspace for one wavelength
  allocate(A(npt, nsp), b(npt), x(nsp), dual(nsp), resid(npt))
  allocate(idx(nsp), jpvt(nsp), free_idx(nsp), bound_idx(nsp))

  ! query optimal lwork for DGELSY (m=npt, n=nsp, nrhs=1)
  lwork = -1
  allocate(work(1))
  call dgelsy(npt, nsp, 1, A, npt, b, npt, jpvt, rcond, rank, work, lwork, info)
  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))

  ! ---- main loop over wavelengths ----------------------------------------
  do jwl = 1, NHD

     ! copy data for this wavelength
     A = conc
     b = absorbance(:, jwl)
     normb = max(norm2(b), eps)
     tol   = 10.0_dp * eps * normb

     ! --- Lawson-Hanson active-set ---------------------------------------
     ! start with all variables free
     nfree  = nsp
     nbound = 0
     do i = 1, nsp
        free_idx(i) = i
     end do
     x = 0.0_dp

     converged = .false.
     do iter = 1, maxit

        ! 1. Build free sub-matrix and solve unconstrained LS
        if (nfree == 0) then
           x = 0.0_dp
           exit
        end if

        do j = 1, nfree
           A(:, j) = conc(:, free_idx(j))   ! re-copy (A was overwritten by DGELSY)
        end do
        b = absorbance(:, jwl)

        jpvt = 0
        call dgelsy(npt, nfree, 1, A, npt, b, npt, jpvt, rcond, rank, &
                    work, lwork, info)
        if (info < 0) return                 ! illegal argument

        ! solution is now in b(1:nfree)
        x = 0.0_dp
        do j = 1, nfree
           x(free_idx(j)) = b(j)
        end do

        ! 2. Feasibility: any negative free variable?
        if (any(x(free_idx(1:nfree)) < -tol)) then
           ! move the most negative variable to the bound set
           k = minloc(x(free_idx(1:nfree)), dim=1)
           nbound = nbound + 1
           bound_idx(nbound) = free_idx(k)
           ! compact free list
           free_idx(k:nfree-1) = free_idx(k+1:nfree)
           nfree = nfree - 1
           cycle
        end if

        ! 3. Compute residual and dual variables for bound set
        !    resid = b - C * x
        resid = absorbance(:, jwl)
        do j = 1, nsp
           if (x(j) > 0.0_dp) resid = resid - x(j) * conc(:, j)
        end do

        if (nbound == 0) then
           converged = .true.
           exit
        end if

        ! dual_i = c_i^T * resid   (for bound variables)
        dual = 0.0_dp
        do j = 1, nbound
           i = bound_idx(j)
           dual(i) = dot_product(conc(:, i), resid)
        end do

        ! 4. KKT test: all duals of bound variables >= -tol ?
        max_dual = maxval(dual(bound_idx(1:nbound)))
        if (max_dual <= tol) then          ! all duals non-negative (within tol)
           converged = .true.
           exit
        end if

        ! 5. Release the bound variable with the largest dual
        k = maxloc(dual(bound_idx(1:nbound)), dim=1)
        nfree = nfree + 1
        free_idx(nfree) = bound_idx(k)
        ! compact bound list
        bound_idx(k:nbound-1) = bound_idx(k+1:nbound)
        nbound = nbound - 1

     end do  ! active-set iteration

     if (.not. converged) then
        ! should almost never happen; keep the last feasible point
        info = max(info, 1)
     end if

     ! store result (force tiny negatives to zero)
     do i = 1, nsp
        extc(i, jwl) = max(x(i), 0.0_dp)
     end do

  end do  ! wavelength loop

  deallocate(A, b, x, dual, resid, idx, jpvt, free_idx, bound_idx, work)
end subroutine solve_extc_nnls


! Queries for the appropriate workspace
subroutine workspace_queries(lwork, nfree, Cfree, B, lda, ldb)
    integer, intent(out) :: lwork
    integer, intent(in) :: nfree, lda, ldb
    double precision, intent(in), dimension(:,:) :: Cfree, B

    integer :: lwork_geqrf, lwork_ormqr, lwork_orgqr, info
    double precision, allocatable :: work(:)
    double precision, allocatable :: tau(:)

    external dgeqrf
    external dormqr
    external dorgqr

    allocate(tau(nfree))

    ! 1. query for DGEQRF
    lwork = -1
    allocate(work(1))
    call dgeqrf(NPT, nfree, Cfree, lda, tau, work, lwork, info)
    lwork_geqrf = nint(work(1))
    deallocate(work)

    ! 2. query for DORMQR (side='L', trans='T')
    lwork = -1
    allocate(work(1))
    call dormqr("L","T", NPT, NHD, nfree, Cfree, lda, tau, B, ldb, work, lwork, info)
    lwork_ormqr = nint(work(1))
    deallocate(work)

    ! 3. query for DORGQR (only needed if you want the explicit Q later)
    lwork = -1
    allocate(work(1))
    call dorgqr(NPT, nfree, nfree, Cfree, lda, tau, work, lwork, info)
    lwork_orgqr = nint(work(1))
    deallocate(work, tau)

    ! final work-array size
    lwork = max(lwork_geqrf, lwork_ormqr, lwork_orgqr)
end subroutine workspace_queries


subroutine titr_bounds(i, low, high)
    integer, intent(in) :: i
    integer, intent(out) :: low, high
    low = 1+sum(titr_size(1:(i-1)))  ! this should work also for i==1
    high = low + titr_size(i) - 1
end subroutine titr_bounds

end module mat
