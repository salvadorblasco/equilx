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


! Given the Hessian matrix, calculate the covariance
subroutine calc_covariance(H, cov, msg)
    double precision, intent(in), dimension(:,:) :: H
    double precision, intent(out), dimension(:,:) :: cov
    character(len=45), intent(out) :: msg

    integer :: n, info, i, j

    external dpotrf
    external dpotri

    n = size(H, dim=1)

    cov = H

    ! 1. LU factorization
    call DPOTRF('U', n, cov, n, info)
    if (info /= 0) then
        msg="COV could not be calculated"
        return
    end if

    ! 2. Invert using the factors
    call DPOTRI('U', n, cov, n, info)
    if (info /= 0) then
        msg="COV could not be calculated"
        return
    end if

    ! 3. Copy the lower triangular matrix
    do i=1,n
        do j=i+1,n
            cov(j,i) = cov(i,j)
        end do
    end do
end subroutine calc_covariance


subroutine calc_correlation(cov, cor)
    double precision, intent(in) :: cov(:,:)
    double precision, intent(out) :: cor(:,:)

    integer :: i, j
    double precision :: di, dj

    do i=1,size(cov, dim=1)
        di = sqrt(cov(i,i))
        do j=1,size(cov, dim=2)
            dj = sqrt(cov(j,j))
            cor(i,j) = cov(i,j) / (di*dj)
        end do
    end do
end subroutine calc_correlation

end module stat

