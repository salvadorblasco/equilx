! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>

subroutine solve(C, T, B, P, CS, KS, Q, Iref, Iback, charge, VIS, VIP, freeze, CTRL, NPT, NEQ, NCP, NSL, NIP)
    implicit none
    integer, intent(in) :: NPT                                    ! number of titration points
    integer, intent(in) :: NCP                                    ! number of components
    integer, intent(in) :: NEQ                                    ! number of equilibria
    integer, intent(in) :: NSL                                    ! number of solids
    integer, intent(in) :: NIP                                    ! number of variable ionic strength parameters
    integer, intent(in) :: freeze                                    ! component of reference
    integer, intent(inout) :: CTRL                                ! control value. 0 on no error, 1 otherwise
    integer, dimension(NEQ, NCP), intent(in) :: P                 ! The soluble stoichiometry array
    integer, dimension(NSL, NCP), intent(in) :: Q                 ! The solids stoichiometry array
    double precision, dimension(NPT, NEQ+NCP), intent(inout) :: C ! free concentation initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NSL), intent(inout) :: CS    ! solid concentration initial estimates (in) and refined values (out)
    double precision, dimension(NPT, NCP), intent(inout) :: T     ! analytical concentration values
    double precision, dimension(NEQ), intent(in) :: B             ! the equilibrium constants
    double precision, dimension(NSL), intent(in) :: KS            ! the solubility constant
    DOUBLE PRECISION, INTENT(IN)                             :: Iref     ! reference ionic strength
    DOUBLE PRECISION, INTENT(IN),    dimension(NPT)          :: Iback    ! background ionic strength
    integer,          INTENT(IN),    DIMENSION(NCP+NEQ)      :: charge   ! charge for every species
    integer, intent(in) :: VIS ! variable ionic strength scheme. 0 -> constant ionic strength
    double precision, intent(in),    dimension(NIP)        :: VIP      ! parameters for variable ionic strength handling
    !f2py integer, intent(hide), depend(t) :: npt = shape(t,0)
    !f2py integer, intent(hide), depend(p) :: ncp = shape(p,1)
    !f2py integer, intent(hide), depend(p) :: neq = shape(p,0)
    !f2py integer, intent(hide), depend(q) :: nsl = shape(q,0)
    !f2py integer, intent(hide), depend(vip) :: nip = shape(vip,0)
    
    logical solids, variable_ionic, frozen_comp
    DOUBLE PRECISION                             :: DHA, DHB ! Debye-Huckel parameters A and B
    double precision,    dimension(NEQ,3)        :: DHP      ! per-species DH parameters C, D and E

    solids = NSL .ne. 0
    variable_ionic = VIS > 0
    frozen_comp = freeze > 0

    if (variable_ionic) then
        if (size(VIP,1) .ne. (3*NEQ+2)) then
            CTRL = 96
            return
        end if
        DHA = VIP(1)
        DHB = VIP(2)
        DHP = reshape(VIP(3:),shape=[NEQ, 3])
    end if

    if (.not. solids .and. .not. variable_ionic .and. .not. frozen_comp) then
        call eq_000(C, T, B, P, CTRL, NPT, NCP, NEQ)
    else if (solids .and. .not. variable_ionic .and. .not. frozen_comp) then
        call eq_0s0(C, T, B, P, CS, KS, Q, CTRL, NPT, NCP, NEQ, NSL)
    else if (.not. solids .and. variable_ionic .and. .not. frozen_comp) then
        call eq_00i(C, T, B, P, Iref, Iback, charge, DHA, DHB, DHP, CTRL, NPT, NCP, NEQ)
    else if (.not. solids .and. .not. variable_ionic .and. frozen_comp) then
        call eq_f00(C, T, B, P, freeze, CTRL, NPT, NCP, NEQ)
    else if (solids .and. .not. variable_ionic .and. frozen_comp) then
        call eq_fs0(C, T, B, P, CS, KS, Q, freeze, CTRL, NPT, NEQ, NCP, NSL)
    else if (solids .and. variable_ionic .and. .not. frozen_comp) then
        call eq_0si(C, T, B, P, CS, KS, Q, Iref, Iback, charge, DHA, DHB, DHP, CTRL, NPT, NCP, NEQ, NSL)
    else if (.not. solids .and. variable_ionic .and. frozen_comp) then
        call eq_f0i(C, T, B, P, Iref, Iback, charge, DHA, DHB, DHP, freeze, CTRL, NPT, NCP, NEQ)
    else if (solids .and. variable_ionic .and. frozen_comp) then
        call eq_fsi(C, T, B, P, CS, KS, Q, Iref, Iback, charge, DHA, DHB, DHP, freeze, CTRL, NPT, NEQ, NCP, NSL)
    end if
end subroutine solve

