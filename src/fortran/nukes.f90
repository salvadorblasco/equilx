program nukes
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none (type, external)

    integer :: NCP, NEQ, NWL, NPT
    real, allocatable :: absorbance(:,:)
    real, allocatable :: analc(:,:)
    real, allocatable :: wavelength(:)
    real, allocatable :: lbeta(:)
    integer, allocatable :: global_stoich(:,:), active(:), flag_lbeta(:)
    character(len=12), allocatable :: cmpnames(:)
    integer, allocatable :: refindices(:)

    real, allocatable :: models_ssr(:)
    real, allocatable :: models_rms(:)
    real, allocatable :: models_aic(:)
    real, allocatable :: models_bic(:)

    integer, allocatable :: models_flags(:,:)
    integer :: num_models

    integer, parameter :: outu = output_unit    ! output unit to write to

    write(outu, "(A)") "     N   N U  U K  K EEEE  SSS"
    write(outu, "(A)") "     NN  N U  U K K  E    S   "
    write(outu, "(A)") "     N N N U  U KK   EE    SS "
    write(outu, "(A)") "     N  NN U  U K K  E       S"
    write(outu, "(A)") "     N   N  UU  K  K EEEE SSS "
    write(outu, "(A,/)") "by Salvador Blasco <salvador.blasco@uv.es>"
    write(outu, "(A,/)") "Starting program"

    call read_data()
    call generate_combinations(flag_lbeta, models_flags)

    num_models = size(models_flags, dim=2)
    allocate(models_ssr(num_models))
    allocate(models_rms(num_models))
    allocate(models_aic(num_models))
    allocate(models_bic(num_models))

    call cycle_models()
    call print_final()

    if (allocated(models_ssr)) deallocate(models_ssr)
    if (allocated(models_rms)) deallocate(models_rms)
    if (allocated(models_aic)) deallocate(models_aic)
    if (allocated(models_bic)) deallocate(models_bic)
contains

subroutine cycle_models()
    integer :: i
    do i=1,size(models_flags, dim=2)
        write (outu, "(A,I2,/)") "Trying model", i
        call fit_try_model(models_flags(:,i))
    end do
end subroutine cycle_models


subroutine fit_try_model(mflags)
    use mat, only: init_module, finalize, load_titration, fit, LE_POSITIVE
    integer, intent(in) :: mflags(:)

    integer :: cNEQ
    real, allocatable :: try_lbeta(:)
    integer, allocatable :: try_model(:,:), try_betaflags(:)
    integer :: i, j

    cNEQ = sum(mflags)
    allocate(try_model(cNEQ, NCP), try_lbeta(cNEQ), try_betaflags(cNEQ))

    j = 0
    do i=1,NEQ
        if (mflags(i) == 1) then
            j = j + 1
            try_model(j,:) = global_stoich(i,:)
            try_lbeta(j) = lbeta(i)
            select case (flag_lbeta(i))
                case (0)
                    try_betaflags(j) = 0
                case (1:2)
                    try_betaflags(j) = 1
                case default
                    try_betaflags(j) = 0
            end select
        end if
    end do

    call init_module(NPT, NWL, 1, try_model, try_lbeta, try_betaflags, refindices, LE_POSITIVE, 3)
    call load_titration(absorbance, analc)
    call fit()

    call report_model(try_model)

    call finalize()

    deallocate(try_model, try_lbeta, try_betaflags)

    write (outu, "(/,A)") "NUKES ended normally"
end subroutine fit_try_model


subroutine report_model(model)
    use mat, only: msg, out_lbeta, out_errlbeta, stat_ssr, stat_dof, stat_rms, stat_nobs
    integer, intent(in) :: model(:,:)

    real :: aic, bic
    integer :: i, nparms

    write (outu, "(/,A)") " .. output message: " // msg

    write (outu, "(A3, *(A12))", advance='no') "N", cmpnames(1:NCP)
    write (outu, "(2A8)") "logB", "error"
    do i=1,size(model, dim=1)
        write (outu, "(I3,*(I8))", advance='no') i, model(i,:)
        write (outu, "(2F12.4)") out_lbeta(i), out_errlbeta(i)
    end do

    nparms = stat_nobs - stat_dof
    aic = stat_nobs*LOG(stat_ssr/stat_nobs) + 2*nparms
    bic = stat_nobs*LOG(stat_ssr/stat_nobs) + nparms*LOG(REAL(stat_nobs))
    write (outu, "(5X,A,E14.4)") "SSR=", stat_ssr
    write (outu, "(5X,A,E14.4)") "RMS=", stat_rms
    write (outu, "(5X,A,E14.4)") "RXS=", stat_ssr/stat_dof
    write (outu, "(5X,A,E14.4)") "AIC=", aic
    write (outu, "(5X,A,E14.4)") "BIC=", bic
    write (outu, "(5X,A,I14)")   "DOF=", stat_dof
end subroutine report_model


subroutine print_final()
    integer :: i, best

    best = minloc(models_bic, dim=1)
    write (outu, "(A,/)") "Final report"
    write (outu, "(A6,4A14)") "Model", "SSR", "RMS", "BIC", "AIC"
    do i = 1,num_models
        write (outu, "(I6, 4E14.4)", advance='no') i, models_ssr(i), models_rms(i), models_bic(i), models_aic(i)
        if (i == best) then
            write (outu, "(A)") " best"
        else
            write (outu, "(/)")
        end if
    end do

end subroutine print_final

!------------------------------------------------------------------
!  Generate combinatorial pattern for a all nested sub-models.
!
!  mask -> array of length NEQ
!  On exit COMBOS(i,j) = 1 if equilibrium i is active in model j,
!                      = 0 otherwise.
!------------------------------------------------------------------
subroutine generate_combinations(mask, COMBOS)
  integer, intent(in) :: mask(:)  ! the mask for each combination; 0=always present; 1=turn on/off
  integer, allocatable, intent(out) :: COMBOS(:, :)  ! dim1 is the equilibria and dim2 is the combination number
    ! it contain 0 and 1 dependin on whether the constant is turned on or off.
  integer :: NCM

  integer :: j, i, bit, bmask, col

  !NEQ = size(mask)
  NCM = 2**count(mask /= 0)-1
  allocate(COMBOS(NEQ, NCM))

  ! build bit mask
  bmask = 0
  do i = 1, NEQ
     if (MASK(i) < 2) bmask = ior(bmask, ishft(1, i-1))
  end do

  ! enumerate and filter
  col = 0
  do j = 1, 2**NEQ - 1
    if (iand(j, bmask) == bmask) then
     col = col + 1
     do i = 1, NEQ
        bit = ishft(1, i-1)          ! 2**(i-1)
        if (iand(j, bit) /= 0) then
           COMBOS(i,col) = 1
        else
           COMBOS(i,col) = 0
        end if
      end do
    end if
  end do
end subroutine generate_combinations

!------------------------------------------------------------------
! read data from stdin
!------------------------------------------------------------------
subroutine read_data()
    use misc, only: build_species_names
    use mat, only: AC_TRANSPARENT, AC_ACTIVE
    !integer, intent(out) :: NCP, NEQ, NWL
    !real, allocatable, intent(out) :: lbeta(:), absorbance(:,:), wavelength(:)
    !integer, allocatable, intent(out) :: stoich(:,:), active(:), flag_lbeta(:)
    !character(len=*), allocatable, intent(out) :: cmpnames(:)

    integer :: ios
    character(len=64) :: title
    integer :: i, j
    real, allocatable :: tmpdata(:,:)
    real, allocatable :: tmpanalc(:,:), tmpbeta(:)
    integer, parameter :: MAX_DATA = 2000, MAX_EQ = 100
    integer :: tmp_npt, nline
    integer, allocatable :: tmpp(:,:), tmpactive(:), tmpkey(:), tmpint(:)
    character(len=256) :: line
    character(len=8), allocatable :: tmpnames(:)

    ! read title
    read(*, "(A)", iostat=ios) title
    nline = 1
    if (ios /= 0) stop "ERROR reading header"
    print "(A)", "Reading data file: " // title

    read(*, *, iostat=ios) NCP, tmp_npt
    nline = nline + 1
    if (ios /= 0) then
        write (outu, "(A,I4)") "ERROR reading on line ", nline
        stop
    end if
    print "(2(A,I4))", " .. components=", NCP, ", points=",tmp_NPT

    allocate(tmpanalc(tmp_NPT,NCP))
    allocate(tmpnames(NCP))
    allocate(tmpp(NCP, MAX_EQ))
    allocate(tmpbeta(MAX_EQ), tmpkey(MAX_EQ), tmpactive(MAX_EQ))

    do i=1,NCP
        read(*, "(A,I4)", iostat=ios) tmpnames(i), tmpactive(i)
        nline = nline + 1
        if (ios /= 0) then
            write (outu, "(A,I4)") "ERROR reading components on line ", nline
            stop
        end if
    end do

    ! read model
    allocate(tmpint(99))
    NEQ = 0
    do
        read(*, "(A)", iostat=ios) line
        nline = nline + 1
        if (ios /= 0) then
            write (outu, "(A,I4)") "ERROR model on line ", nline
            stop
        end if
        if (ios /= 0) stop "Error reading model"

        ! Remove leading and trailing blanks
        line = adjustl(trim(line))

        ! Stop if the line is blank
        if (len_trim(line) == 0) exit
        NEQ = NEQ + 1
        read(line, "(F8.4,*(I5))", iostat=ios) tmpbeta(NEQ), tmpint(1:NCP+2)
        if (ios /= 0) then
            write (outu, "(A,I4)") "ERROR model on line ", nline
            stop
        end if
        tmpp(:,NEQ) = tmpint(1:NCP)
        tmpkey(NEQ) = tmpint(NCP+1)
        tmpactive(NCP+NEQ) = tmpint(NCP+2)
    end do
    deallocate(tmpint)

    allocate(global_stoich(NEQ,NCP),lbeta(NEQ), active(NCP+NEQ), flag_lbeta(NEQ))

    global_stoich = transpose(tmpp(:,1:NEQ))
    lbeta = tmpbeta(1:NEQ)
    active = tmpactive(1:(NCP+NEQ))
    flag_lbeta = tmpkey(1:NEQ)
    deallocate(tmpp,tmpbeta,tmpactive,tmpkey)

    ! prepare components names
    allocate(cmpnames(NCP+NEQ))
    cmpnames(1:NCP) = tmpnames
    call build_species_names(global_stoich, cmpnames, NCP, NEQ)

    write (outu, "(A)") " .. components found"
    do i=1,NCP+NEQ
        write (outu, "(I5,2X,A12)", advance="no") i, cmpnames(i)
        if (active(i) < 0) then
            write (outu, "(A,I3,A)") "reference [point ", -active(i), "]"
        else if (active(i) == AC_TRANSPARENT) then
            write (outu, "(A)") "transparent"
        else if (active(i) == AC_ACTIVE) then
            write (outu, "(A)") "active"
        end if
    end do

    write (outu, "(A, I0, A)") " .. found ", NEQ, " equilibria"
    write (outu, "(A)") " .. reading analytical concentrations (values assumed to be mol/L)"
    do i=1,NCP
        read(*, *, iostat=ios) tmpanalc(:,i)
        if (ios /= 0) stop "Error reading components concentations"
    end do

    allocate(tmpdata(tmp_NPT+1,MAX_DATA))
    i = 1
    do
        read(*, *, iostat=ios) tmpdata(:,i)
        if (ios < 0) exit
        if (ios > 0) stop "Error reading data file"
        i = i + 1
    end do
    NWL = i-1

    allocate(wavelength(NWL))
    allocate(refindices(tmp_npt))

    refindices = 0
    NPT = tmp_NPT
    do i=1,NCP+NEQ
        if (active(i) < 0) then
            if (refindices(-active(i)) /= 0) stop "Duplicate index"
            refindices(-active(i)) = i
        end if
    end do

    allocate(absorbance(NPT,NWL))
    wavelength = tmpdata(1,:)
    allocate(analc(NPT,NCP))
    analc = 0.0
    j = 0
    do i=1,tmp_npt
        j = j + 1
        absorbance(j,1:NWL) = tmpdata(1+i,1:NWL)
        analc(j,:) = tmpanalc(i,:)
    end do

    write(outu,"(A,I5,A)") " .. read ", NWL, " wavelengths"

    deallocate(tmpdata, tmpanalc)
    write(outu,"(A,/)") " .. data loaded correctly"
end subroutine read_data

end program nukes

