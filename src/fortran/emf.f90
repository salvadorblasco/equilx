! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>

module emf
    implicit none
    private
    double precision, allocatable, dimension(:,:)  :: C    ! the free concentrations
    real, allocatable, dimension(:)  :: V                  ! titre in mL
    real, allocatable, dimension(:)  :: V0                 ! initial amount in mL
    real, allocatable, dimension(:)  :: E                  ! experimental emf in mV
    real, allocatable, dimension(:)  :: E0                 ! standard potential in mV
    real, allocatable, dimension(:)  :: slope              ! nernstian slope in unit
    real, allocatable, dimension(:,:)  :: titr_init        ! initial amount in mmol
    real, allocatable, dimension(:,:)  :: titr_buret       ! buret concentration in mmol/mL
    real, ALLOCATABLE, DIMENSION(:) :: lbeta               ! formation constants in log10 units
    integer, allocatable, dimension(:)  :: flag_beta
    integer, allocatable, dimension(:,:)  :: flag_init
    integer, allocatable, dimension(:,:)  :: flag_buret
    integer, allocatable, dimension(:)  :: flag_emf0
    integer, allocatable, dimension(:)  :: flag_slope
    integer, public, protected :: NTT   ! number of titrations
    integer, public, protected :: NPT   ! number of titration points for all titrations
    integer, public, protected :: NCP   ! number of components
    integer, public, protected :: NEQ   ! number of equilibria
    real :: nernst
    integer, ALLOCATABLE, DIMENSION(:) :: titr_size 
    integer, allocatable, dimension(:,:) :: P
    integer, public, protected :: diagnostic  ! 0 no error; 1 size of data zero or negative; 2 P and Q incompatible;
                           ! 3  insucient memory allocated, increase n_points;
                           ! 4 no room for more titrations; 5 invalid weighting scheme
                           ! 99 -> linear system could not be solved
    integer :: tmp_loaded_titrs, tmp_space_left

    real, allocatable, DIMENSION(:) :: tmp_lbeta
    real, allocatable, DIMENSION(:,:) :: tmp_init
    real, allocatable, DIMENSION(:,:) :: tmp_buret
    real, allocatable, DIMENSION(:) :: tmp_emf0
    real, allocatable, DIMENSION(:) :: tmp_slope
    real, allocatable, DIMENSION(:,:) :: T
    logical, ALLOCATABLE, DIMENSION(:) :: update_titr  ! flag to know which titration must be updated

    integer, parameter :: RF_CONSTANT = 0
    integer, parameter :: RF_REFINE = 1

    PUBLIC :: init_module, load_titration, finalize, fit_emf

contains

subroutine init_module(n_points, n_titrs, stoich, logbeta, flgbeta)
    integer, INTENT(IN) :: n_points
    integer, INTENT(IN) :: n_titrs
    integer, intent(in), dimension(:,:) :: stoich
    double precision, INTENT(IN), DIMENSION(:) :: logbeta
    integer, intent(in), dimension(:) :: flgbeta

    diagnostic = 0

    if (n_titrs < 1) then
        diagnostic = 1
        return
    end if
    NTT = n_titrs

    if (n_points < 1) then
        diagnostic = 1
        return
    end if
    NPT = n_points

    NEQ = SIZE(stoich, dim=1)
    NCP = SIZE(stoich, dim=2)

    ALLOCATE(E(NPT))
    ALLOCATE(V(NPT))
    ALLOCATE(V0(NTT))
    ALLOCATE(E0(NTT))
    ALLOCATE(slope(NTT))
    ALLOCATE(titr_init(NTT,NCP)) 
    ALLOCATE(titr_buret(NTT,NCP)) 
    ALLOCATE(flag_init(NTT,NCP)) 
    ALLOCATE(flag_buret(NTT,NCP)) 
    allocate(titr_size(NTT))
    allocate(P(NEQ, NCP))
    ALLOCATE(lbeta(NEQ))
    ALLOCATE(flag_beta(NEQ))
    allocate(T(NPT,NCP))
    ALLOCATE(update_titr(NTT))
    update_titr = .false.

    P = stoich
    lbeta = logbeta * LOG(10.0)     ! convert from log10 to ln
    flag_beta = flgbeta

    tmp_loaded_titrs = 0
    tmp_space_left = NPT
end subroutine init_module


subroutine load_titration(titre, emf, titrparms, initV, emf0, nernst_slope, tflags, emf0_flag, slope_flag)
    implicit none
    real, intent(in), dimension(:) :: titre
    real, intent(in), dimension(size(titre)) :: emf
    real, intent(in), dimension(NCP,2) :: titrparms
    real, INTENT(IN) :: initV
    real, intent(in) :: emf0, nernst_slope
    integer, intent(in), dimension(:,:) :: tflags
    integer, INTENT(IN) :: emf0_flag, slope_flag

    integer n, space_left, bottom

    if (tmp_loaded_titrs .eq. NTT) then
        diagnostic = 4
        return
    end if

    if (all(shape(titrparms) .ne. shape(tflags))) then
        diagnostic = 5
        return
    end if

    n = SIZE(titre, dim=1)
    if (n > tmp_space_left) then
        diagnostic = 3
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

    E0(tmp_loaded_titrs) = emf0
    slope(tmp_loaded_titrs) = nernst_slope

    titr_init(tmp_loaded_titrs,:) = titrparms(:, 1)
    titr_buret(tmp_loaded_titrs,:) = titrparms(:, 2)
    flag_init(tmp_loaded_titrs, :) = tflags(:,1)
    flag_buret(tmp_loaded_titrs, :) = tflags(:,2)
end subroutine load_titration


subroutine finalize()
    DEALLOCATE(E)
    DEALLOCATE(V)
    DEALLOCATE(V0)
    DEALLOCATE(E0)
    DEALLOCATE(slope)
    DEALLOCATE(titr_init) 
    DEALLOCATE(titr_buret) 
    DEALLOCATE(flag_init) 
    DEALLOCATE(flag_buret) 
    DEallocate(titr_size)
    DEallocate(P)
    DEALLOCATE(lbeta)
    DEALLOCATE(flag_beta)
    DEALLOCATE(T)
    DEALLOCATE(update_titr)
end subroutine finalize


subroutine fit_emf(temp, weight_scheme)
    implicit none
    real,   intent(in) :: temp
    integer, intent(in) :: weight_scheme

    real, parameter :: nernst_base = 0.08617333262 !mV/K
    DOUBLE PRECISION, DIMENSION(NPT,NPT) :: W
    integer i

    ! handle temperature
    nernst = nernst_base*temp  ! mV

    ! handle weights
    if (weight_scheme .eq. 0) then      ! unit weights
        W = 0.0
        do i=1,NPT
            W(i,i) = 1.0
        end do
    else if (weight_scheme .eq. 1) then      ! unit weights
    else
        diagnostic = 6
    end if

    ! TODO
end subroutine fit_emf


subroutine lm(W)
    ! routine for Levenberg-Marquardt iterations.
    implicit none

    double precision, intent(in), dimension(:,:) :: W  ! the weights matrix

    integer :: NX ! number of variables to optimize
    integer iter  ! iteration numnber
    integer, parameter :: MAX_ITERS = 100
    real, parameter :: DAMPING_UPPER = 1e10
    real, parameter :: DAMPING_LOWER = 1e-12
    real, parameter :: DAMPING_UPF = 8.0
    real, parameter :: DAMPING_LOWF = 3.0
    real :: lm_damping
    double precision, allocatable, dimension(:,:) :: jac
    double precision, allocatable, dimension(:) :: dx
    double precision, allocatable, dimension(:) :: fvec, trial_fvec
    double precision, allocatable, dimension(:) :: gradient
    double precision, allocatable, dimension(:,:) :: M, dampedM
    double precision :: gradient_norm, chisq, damping
    integer i
    integer CTRL

    dx = 0.0
    jac = 0.0
    fvec = 0.0
    damping = 1D2

    call count_dof(NX)

    allocate(jac(NPT, NX))
    ALLOCATE(dx(NX))
    ALLOCATE(gradient(NX))
    ALLOCATE(fvec(NPT))
    ALLOCATE(trial_fvec(NPT))
    ALLOCATE(M(NX,NX))
    ALLOCATE(dampedM(NX,NX))

    dx = 0.0

    main_loop: do iter=0,MAX_ITERS
        if (iter .eq. 0) then
            call try_step(dx, 0)
        else
            call try_step(dx, 1)
        end if

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
                diagnostic = 99
                return
            end if
        end if

    end do main_loop

    DEALLOCATE(jac, dx, fvec, M)
end subroutine lm


subroutine update_equil()
    ! update the free concentrations
    integer i, blow, bhigh, CTRL
    do i=1,NTT
        if (update_titr(i)) then
            blow = 1+sum(titr_size(1:(i-1)))  ! this should work also for i==1
            bhigh = blow + titr_size(i)
            call calc_analc(V(blow:bhigh), V0(i), tmp_init(i,:), tmp_buret(i,:), T(blow:bhigh,:))
            call eq_000(C(blow:bhigh,:), T(blow:bhigh,:), 10.0**tmp_lbeta, P, CTRL, titr_size(i), NCP, NEQ)
        end if
    end do
end subroutine update_equil


subroutine try_step(step, task)
    double precision, intent(in), dimension(:) :: step
    integer, intent(in) :: task

    double precision, save, allocatable, DIMENSION(:) :: xvalues
    integer, save, allocatable, dimension(:) :: constraint_loc

    real, dimension(10) :: rconstr

    integer i, k

    if (task .eq. 0) then ! initiation
        allocate(xvalues(size(step)))
        xvalues = 0.0

        allocate(constraint_loc(size(step)))
        constraint_loc = 0


        rconstr = 0.0
        k = 0

        ! I. count formation constants
        if (any(flag_beta  .ne. RF_CONSTANT)) then
            allocate(tmp_lbeta(NEQ))
            tmp_lbeta = lbeta
            call aux_fillx(lbeta, flag_beta, xvalues, k, constraint_loc)
        end if

        ! II. count initial amounts
        if (any(flag_init  .ne. RF_CONSTANT)) then
            allocate(tmp_init(NTT, NCP))
            tmp_init = titr_init
            do i=1,size(flag_init, dim=1)
                call aux_fillx(titr_init(i,:), flag_init(i,:), xvalues, k, constraint_loc)
                if (any(flag_init(i,:)  .ne. RF_CONSTANT)) update_titr(i) = .true.
            end do
        end if

        ! III. count buret concentrations
        if (any(flag_buret .ne. RF_CONSTANT)) then
            allocate(tmp_buret(NTT, NCP))
            tmp_buret = titr_buret
            do i=1,size(flag_buret, dim=1)
                call aux_fillx(titr_buret(i,:), flag_buret(i,:), xvalues, k, constraint_loc)
                if (any(flag_buret(i,:)  .ne. RF_CONSTANT)) update_titr(i) = .true.
            end do
        end if

        ! IV: count standard potentials
        if (any(flag_emf0  .ne. RF_CONSTANT)) then
            allocate(tmp_emf0(NTT))
            tmp_emf0 = E0
            call aux_fillx(E0, flag_emf0, xvalues, k, constraint_loc)
        end if

        ! V. count nernstian slopes
        if (any(flag_slope .ne. RF_CONSTANT)) then
            allocate(tmp_slope(NTT))
            tmp_slope = slope
            call aux_fillx(slope, flag_slope, xvalues, k, constraint_loc)
        end if
    end if

    if (task .eq. 1) then ! trial step
        k = 0
        rconstr = 0.0
        if (allocated(tmp_lbeta)) then
            tmp_lbeta = lbeta
            call aux_takestep(flag_beta, tmp_lbeta, lbeta, xvalues+step, rconstr, k)
        end if
        if (allocated(tmp_init))   then
            do i=1,size(flag_init, dim=1)
                call aux_takestep(flag_init(i,:), tmp_init(i,:), titr_init(i,:), xvalues+step, rconstr, k)
            end do
        end if
        if (allocated(tmp_buret)) then
            do i=1,size(flag_buret, dim=1)
                call aux_takestep(flag_buret(i,:), tmp_buret(i,:), titr_buret(i,:), xvalues+step, rconstr, k)
            end do
        end if
        if (allocated(tmp_emf0)) then
            tmp_emf0 = E0
            call aux_takestep(flag_emf0, tmp_emf0, E0, xvalues+step, rconstr, k)
        end if
        if (allocated(tmp_slope)) then 
            tmp_slope = slope
            call aux_takestep(flag_slope, tmp_slope, slope, xvalues+step, rconstr, k)
        end if
    end if

    if (task .eq. 2) then ! cleanup
        deallocate(xvalues, constraint_loc)
        if (allocated(tmp_lbeta)) deallocate(tmp_lbeta)
        if (allocated(tmp_init))  deallocate(tmp_init)
        if (allocated(tmp_buret)) deallocate(tmp_buret)
        if (allocated(tmp_emf0))  deallocate(tmp_emf0)
        if (allocated(tmp_slope)) deallocate(tmp_slope)
    end if
end subroutine try_step


subroutine aux_takestep(flag, tmpv, refv, xval, ratio, counter)
    integer, INTENT(IN), DIMENSION(:) :: flag
    real, INTENT(INout), DIMENSION(size(flag)) :: tmpv
    real, INTENT(IN), DIMENSION(size(flag)) :: refv
    double precision, INTENT(IN), DIMENSION(size(flag)) :: xval
    real, INTENT(INOUT), DIMENSION(:) :: ratio
    integer, INTENT(INout) :: counter

    integer i

    do i=1,size(tmpv)
        if (flag(i) .eq. RF_REFINE) then
            counter = counter + 1
            tmpv(i) = xval(counter)
        else if (flag(i) > RF_REFINE) then
            if (ratio(flag(i)-1) .eq. 0.0) then
                counter = counter + 1
                tmpv(i) = xval(counter)
                ratio(flag(i)-1) = xval(counter)/refv(i)
            else
                tmpv(i) = ratio(flag(i)-1)*xval(counter)
            end if
        end if
    end do
end subroutine aux_takestep


subroutine aux_fillx(vec, flag, xv, counter, locv)
    real, intent(in), dimension(:) :: vec
    double precision, INTENT(INOUT), DIMENSION(:) :: xv
    integer, INTENT(INOUT), DIMENSION(:) :: locv
    integer, INTENT(IN), DIMENSION(:) :: flag
    integer, intent(inout) :: counter

    integer i

    do i=1,size(vec)
        if (flag(i) .eq. RF_REFINE) then
            counter = counter + 1
            xv(counter) = vec(i)
        else if (flag(i) > RF_REFINE) then
            if (locv(flag(i)-1) .eq. 0) then
                counter = counter + 1
                locv(flag(i)-1) = counter
            end if
        end if
    end do
end subroutine aux_fillx


subroutine jacobian_lm(jac)
    ! Compute the jacobian for the Levenberg-Marquardt iteration.
    ! The rows are the experimental points; the columns are the variables to be optimized in the
    ! following order: (1) equilibrium constants; (2) initial amounts; (3) buret concentrations;
    ! (4) standard potentials; (5) nernstian slopes; (6) constraints
    double precision, intent(out), dimension(:,:) :: jac

    double precision, dimension(NCP, NEQ) :: dlcdlb  ! the result
    !call calc_dlogcdlogb(C, P, NCP, NEQ, dlcdlb)
end subroutine jacobian_lm


subroutine count_dof(NX)
    ! count degrees of freedom
    implicit none
    integer, intent(out) :: NX ! number of variables to optimize

    logical constr(9)
    integer i, j
    NX = 0
    constr = .false.

    ! I. count formation constants
    call sub_count_dof(flag_beta, NX, constr)

    ! II. count initial amounts
    do i=1,size(flag_init, dim=1)
        call sub_count_dof(flag_init(i,:), NX, constr)
    end do

    ! III. count buret concentrations
    do i=1,size(flag_buret, dim=1)
        call sub_count_dof(flag_buret(i,:), NX, constr)
    end do

    ! IV: count standard potentials
    call sub_count_dof(flag_emf0, NX, constr)

    ! V. count nernstian slopes
    call sub_count_dof(flag_slope, NX, constr)
end subroutine count_dof


subroutine sub_count_dof(vec, NX, constr)
    ! auxiliary routine to count degrees of freedom
    integer, intent(inout) :: NX ! number of variables to optimize
    logical, intent(inout) :: constr(:)
    integer, intent(in) :: vec(:)

    integer i
    do i=1,size(vec)
        if (vec(i) .eq. RF_REFINE) then
            NX = NX + 1
        else if (vec(i) > RF_REFINE) then
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
  dy(1) = p1*y(1) + p2*y(2) + p3*y(3)

  ! ---- interior points (central 2nd-order) ----------------------------------
  do i = 2, n-1
     h1 = x(i)   - x(i-1)
     h2 = x(i+1) - x(i)
     p1  = -h2 / (h1*(h1+h2))
     p2  =  (h2 - h1) / (h1*h2)
     p3  =  h1 / (h2*(h1+h2))
     dy(i) = p1*y(i-1) + p2*y(i) + p3*y(i+1)
  end do

  ! ---- right boundary (backward 2nd-order) ----------------------------------
  h1 = x(n-1) - x(n-2)
  h2 = x(n)   - x(n-1)
  p1  =  h2 / (h1*(h1+h2))
  p2  = -(h1 + h2)/(h1*h2)
  p3  =  (2.0d0*h2 + h1)/(h2*(h1+h2))
  dy(n) = p1*y(n-2) + p2*y(n-1) + p3*y(n)

end subroutine derivative_1d

end module emf
