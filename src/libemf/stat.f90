module stat
    implicit none (type, external)
    public
contains

pure function calc_chisq(resids, weights) result(chsq)
    double precision, intent(in), dimension(:) :: resids
    double precision, intent(in), dimension(:) :: weights

    double precision :: chsq

    chsq = sum(weights*resids**2)
end function calc_chisq


! Calculate the fit's sigma value for a given set of residuals and weights.
function fit_sigma(residuals, weights, npoints, nparams) result(sigma)
    double precision, intent(in), dimension(npoints) :: residuals
    double precision, intent(in), dimension(npoints) :: weights
    integer, intent(in) :: npoints
    integer, intent(in) :: nparams
    double precision :: sigma

    sigma = SQRT(sum(weights*residuals**2)/(sum(weights)*(npoints-nparams)))
end function fit_sigma



end module stat

