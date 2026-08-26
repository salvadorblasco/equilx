! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>

module emf
    implicit none
    private
    double precision, allocatable, dimension(:,:)  :: C
    real, allocatable, dimension(:)  :: V
    real, allocatable, dimension(:)  :: V0
    real, allocatable, dimension(:)  :: E
    real, allocatable, dimension(:)  :: E0
    real, allocatable, dimension(:)  :: slope
    real, allocatable, dimension(:,:)  :: titr_init
    real, allocatable, dimension(:,:)  :: titr_buret
    double precision, ALLOCATABLE, DIMENSION(:) :: lbeta
    integer, allocatable, dimension(:)  :: flag_beta
    integer, allocatable, dimension(:,:)  :: flag_init
    integer, allocatable, dimension(:,:)  :: flag_buret
    integer, allocatable, dimension(:)  :: flag_emf0
    integer, allocatable, dimension(:)  :: flag_slope
    integer, public, protected :: NTT   ! number of titrations
    integer, public, protected :: NPT   ! number of titration points for all titrations
    integer, public, protected :: NCP   ! number of components
    integer, public, protected :: NEQ   ! number of equilibria
    integer :: weight_scheme
    integer, ALLOCATABLE, DIMENSION(:) :: titr_size 
    integer, allocatable, dimension(:,:) :: P
    integer, public, protected :: diagnostic  ! 0 no error; 1 size of data zero or negative; 2 P and Q incompatible;
                           ! 3  insucient memory allocated, increase n_points;
                           ! 4 no room for more titrations
    integer :: tmp_loaded_titrs, tmp_space_left

    integer, parameter :: RF_CONSTANT = 0
    integer, parameter :: RF_REFINE = 1

    PUBLIC :: init_module, load_titration, finalize, fit_emf

contains

subroutine init_module(n_points, n_titrs, wsch, stoich, logbeta, flgbeta)
    integer, INTENT(IN) :: n_points
    integer, INTENT(IN) :: n_titrs
    INTEGER, INTENT(IN) :: wsch
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

    P = stoich
    lbeta = logbeta
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
end subroutine finalize


subroutine fit_emf()
    implicit none
end subroutine fit_emf


subroutine lm(W)
    ! routine for Levenberg-Marquardt iterations.
    implicit none

    double precision, intent(in), dimension(:,:) :: W  ! the weights matrix

    integer :: NX ! number of variables to optimize
    integer iter  ! iteration numnber
    integer, parameter :: MAX_ITERS = 100
    double precision :: lm_damping
    double precision, allocatable, dimension(:,:) :: jac
    double precision, allocatable, dimension(:) :: dx
    double precision, allocatable, dimension(:) :: fvec

    dx = 0.0
    jac = 0.0
    fvec = 0.0

    main_loop: do iter=1,MAX_ITERS
        call jacobian_lm(jac)
        call calc_resid(fvec)
    end do main_loop
end subroutine lm


subroutine calc_resid(fvec)
    double precision, dimension(:) :: fvec
end subroutine calc_resid


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

end module emf
