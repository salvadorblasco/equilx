! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>

module emf
    implicit none
    private
    double precision, allocatable, dimension(:,:)  :: C
    double precision, allocatable, dimension(:)  :: V
    double precision, allocatable, dimension(:)  :: E
    double precision, allocatable, dimension(:)  :: E0
    double precision, allocatable, dimension(:)  :: slope
    double precision, allocatable, dimension(:,:)  :: titr_init
    double precision, allocatable, dimension(:,:)  :: titr_buret
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
    integer, ALLOCATABLE, DIMENSION(:) :: cuts 
    integer, allocatable, dimension(:,:) :: P
    integer, public, protected :: diagnostic  ! 0 no error; 1 size of data zero or negative; 2 P and Q incompatible;
                           ! 3  insucient memory allocated, increase n_points;
                           ! 4 no room for more titrations
    integer :: tmp_loaded_titrs, tmp_space_left
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
    ALLOCATE(E0(NTT))
    ALLOCATE(slope(NTT))
    ALLOCATE(titr_init(NTT,NCP)) 
    ALLOCATE(titr_buret(NTT,NCP)) 
    ALLOCATE(flag_init(NTT,NCP)) 
    ALLOCATE(flag_buret(NTT,NCP)) 
    allocate(cuts(NTT))
    allocate(P(NEQ, NCP))
    ALLOCATE(lbeta(NEQ))
    ALLOCATE(flag_beta(NEQ))

    P = stoich
    lbeta = logbeta
    flag_beta = flgbeta

    tmp_loaded_titrs = 0
    tmp_space_left = NPT
end subroutine init_module


subroutine load_titration(curve, titrparms, emf0, nernst_slope, tflags, emf0_flag, slope_flag)
    implicit none
    double precision, intent(in), dimension(:,:) :: curve
    double precision, intent(in), dimension(:,:) :: titrparms
    double precision, intent(in) :: emf0, nernst_slope
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

    n = SIZE(curve, dim=1)
    if (n < tmp_space_left) then
        diagnostic = 3
        return
    end if

    tmp_space_left = tmp_space_left - n
    
    tmp_loaded_titrs = tmp_loaded_titrs + 1
    if (tmp_loaded_titrs .eq. 1) then
        bottom = 1
    else
        bottom = cuts(tmp_loaded_titrs)
    end if
    V(bottom:bottom+n) = curve(:,1)
    E(bottom:bottom+n) = curve(:,2)

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
    DEALLOCATE(E0)
    DEALLOCATE(slope)
    DEALLOCATE(titr_init) 
    DEALLOCATE(titr_buret) 
    DEALLOCATE(flag_init) 
    DEALLOCATE(flag_buret) 
    DEallocate(cuts)
    DEallocate(P)
    DEALLOCATE(lbeta)
    DEALLOCATE(flag_beta)
end subroutine finalize


subroutine fit_emf()
    implicit none
end subroutine fit_emf


subroutine lm(W, NX, NV)
    ! routine for Levenberg-Marquardt iterations.
    implicit none

    integer, intent(in) :: NX ! number of experimental points
    integer, intent(in) :: NV ! number of variables to optimize
    double precision, intent(in), dimension(NX,NX) :: W  ! the weights matrix

    integer iter  ! iteration numnber
    integer, parameter :: MAX_ITERS = 100
    double precision :: lm_damping

    main_loop: do iter=1,MAX_ITERS

    end do main_loop
end subroutine lm


subroutine jacobian_lm()
    ! Compute the jacobian for the Levenberg-Marquardt iteration.
    ! The rows are the experimental points; the columns are the variables to be optimized in the
    ! following order: (1) equilibrium constants; (2) initial amounts; (3) buret concentrations;
    ! (4) standard potentials; (5) nernstian slopes; (6) constraints
    double precision, dimension(NCP, NEQ) :: dlcdlb  ! the result
    call calc_dlogcdlogb(C, P, NCP, NEQ, dlcdlb)
end subroutine jacobian_lm

end module emf
