! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
! by Salvador Blasco <salvador.blasco@uv.es>
!
! Miscellaneous routines

module misc

contains

! Calculate the fit's sigma value for a given set of residuals and weights.
function fit_sigma(residuals, weights, npoints, nparams) result(sigma)
    double precision, intent(in), dimension(npoints) :: residuals
    double precision, intent(in), dimension(npoints) :: weights
    integer, intent(in) :: npoints
    integer, intent(in) :: nparams
    double precision :: sigma

    sigma = SQRT(sum(weights*residuals**2)/(sum(weights)*(npoints-nparams)))
end function fit_sigma


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
    call dsysv("U", n, nrhs, A, lda, ipiv, b, ldb, work, lwork, info)
    lwork = nint(work(1))
    deallocate(work)

    allocate(work(lwork))
    call dsysv("U", n, nrhs, A, lda, ipiv, b, ldb, work, lwork, info)
    deallocate(work)

    if (info /= 0) then
        CTRL = 99    ! unable to invert matrix
        return
    end if
end subroutine msysolve


subroutine lesolve(A, b, x, n, CTRL)
  ! Solve system of lineal equations using LAPACK. This routine requires
  ! the -llapack flag at compilation and, probably also -lblas in some systems.
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: A   ! The nxn array of coefficients
  double precision, dimension(n),   intent(in) :: b   ! the n array of independent terms
  double precision, dimension(n),   intent(out) :: x  ! the n array where the result is returned
  integer, intent(inout) :: CTRL                      ! control value. 0 on no error.

  integer, dimension(n) :: ipiv   ! pivot indices
  integer :: info, nrhs, lda, ldb

  ! External procedures defined in LAPACK
  external dgesv

  ! Store b in x to prevent it from being overwritten by LAPACK
  x = b

  ! Set up solve by LAPACK routine DGESV
  nrhs = 1 ! number of right hand sides in b
  lda = n  ! leading dimension of a
  ldb = n  ! leading dimension of b

  call dgesv(n, nrhs, A, lda, ipiv, x, ldb, info)

  CTRL = info
end subroutine lesolve


subroutine calc_analc(titre, init_vol, init, buret, analc)
    ! Calculate the analytical concentration, given the titration paramers, 
    real, intent(in), dimension(:)  :: titre
    real, intent(in), dimension(:)  :: init
    real, intent(in), dimension(:)  :: buret
    real, intent(out), dimension(:,:)  :: analc
    real, INTENT(IN) :: init_vol

    integer i, j

    do i=1,size(titre)
        do j=1,size(init)
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

end module misc




!def covariance(J, W):
!    """Compute covariance matrix.
!
!    Parameters:
!        J (:class:`numpy.ndarray`): the jacobian
!        W (:class:`numpy.ndarray`): the weights matrix
!    Returns:
!        :class:`numpy.ndarray`: an (*p*, *p*)-sized array representing
!            the covariance matrix.
!    """
!    aux2 = J.T @ np.diag(np.diag(W)) @ J
!    return np.linalg.pinv(aux2)
!
!
!def fitting_errors(covar):
!    """Return the standard errors of the fitted parameters.
!
!    Extracts the square root of the diagonal elements of the covariance
!    matrix, which equals the standard deviation (1 σ) of each parameter.
!
!    Parameters
!    ----------
!    covar : numpy.ndarray
!        Covariance matrix of shape ``(p, p)``, as returned by
!        :func:`covariance`.
!
!    Returns
!    -------
!    numpy.ndarray
!        1-D array of length *p* with the standard errors.
!    """
!    return np.sqrt(np.diag(covar))
!
!
!def correlation_matrix(covar):
!    """Compute the Pearson correlation matrix from a covariance matrix.
!
!    Normalises *covar* by the outer product of its diagonal elements:
!
!    .. math::
!
!        R_{ij} = \\frac{\\text{Cov}_{ij}}{\\sqrt{\\text{Cov}_{ii}\\,\\text{Cov}_{jj}}}
!
!    Parameters
!    ----------
!    covar : numpy.ndarray
!        Covariance matrix of shape ``(p, p)``.
!
!    Returns
!    -------
!    numpy.ndarray
!        Correlation matrix of shape ``(p, p)`` with values in ``[-1, 1]``.
!    """
!    D = np.diag(covar)
!    nD = len(D)
!    return covar/np.sqrt(np.dot(D.reshape((nD, 1)), D.reshape((1, nD))))
!
!
!def extrapoly(x0, X, Y):
!    r"""Polynomial extrapolation.
!
!    Given a list of *x* and *y* points, and an *x0* point, this routine
!    calculates the polynomials that goes through the given points and
!    evaluates them for *x0*. It evaluates *m* extrapolations at once out of
!    *n* polynomials of degree *g* -1
!
!    Parameters:
!        x0 (:class:`numpy.ndarray`): An (m, n)-array with the values for the
!            polynomials to be evaluated.
!        X (:class:`numpy.ndarray`): An (n, g)-array with the *x* values of the
!            data.
!        Y (:class:`numpy.ndarray`): An (n, g)-array with the *y* values of the
!            data. Axis-0 must be of the same length than X.
!
!    Returns:
!        :class:`numpy.ndarray`: An (m, n)-array with calculated *y* values
!             that produce the evaluation of the polynomials at *x0*
!    """
!    def polyexp(array, axis, n):
!        """Given any array and one axis of it, return the same array with one
!        extra dimmension in which the new dimmension is the array raised to
!        0..n
!
!        >>> a = np.array([[1, 1], [2, 2]])
!        >>> polyexp(a, -1, 2)
!        [[[1, 1], [1, 1]],
!         [[1, 1], [2, 2]],
!         [[1, 1], [4, 4]]]
!        """
!        a = np.expand_dims(array, axis)
!        s = [1] * a.ndim
!        s[axis] = n+1
!        b = np.reshape(np.arange(n+1), tuple(s))
!        return np.power(a, b)
!
!    g, n = Y.shape
!    Q = polyexp(X.T, -1, g-1)
!    assert Q.shape == (n, g, g)
!    A = np.linalg.solve(Q, Y.T)
!    assert A.shape == (n, g)
!    m = x0.shape[1]
!    assert x0.shape[0] == n
!    xexp = polyexp(x0.T, -1, g-1)
!    assert xexp.shape == (m, n, g)
!    r = np.sum(xexp * A[np.newaxis, ...], axis=-1)
!    return np.squeeze(r)
!
!
!def quadratic_extrapolation(a, b, c):
!    """Perform 3-point extrapolation.
!
!    Given three ordered points of an evenly spaced curve, return the
!    quadratic extrapolation for estimating the following point.
!    """
!    return a - 3*b + 3*c
!
!
!def nearest(a, b):
!    """Find nearest element in a sorted list of numbers.
!
!    Given a two list of indices, returns a sorted list of the indices of
!    **a** which are closer to those of **b**
!    """
!    return np.argsort(np.abs(a[:, None]-b[None, :]), axis=0)
!
!
!def sample_size_change(data, new_size):
!    """Resample a free concentration array.
!
!    Parameters:
!        data (:class:`numpy.ndarray`): The data to be resized. Dimmension 0
!            must be the sample number. Lenght of dimmension is assumed to
!            be the size of the old data.
!        new_size (int): The length of the new data.
!    Returns:
!        :class:`numpy.ndarray`: The resized array along dimmension 0.
!    """
!    from scipy.interpolate import interp1d
!    old_size = data.shape[0]
!    f_interp = interp1d(np.arange(old_size), data, axis=0, assume_sorted=True)
!    new_x = np.linspace(0, old_size-1, new_size)
!    return f_interp(new_x)
!
!
!def weighting_slope(x, y, error_x, error_y):
!    r"""Calculate weighting scheme acording to Gans et al.
!
!    .. math:: w = s_E^2 + \left(\frac{\partial E}{\partial V}\right)^2  s_V^2
!
!    Parameters:
!        x (iterable): volume in mL
!        y (iterable): emf in mV
!        error_x (iterable): error of volume in mL
!        error_y (iterable): error of emf in mV
!
!    Returns:
!        :class:`numpy.ndarray`: An array containing the calculated weights
!    """
!    dydx = np.gradient(x, y)
!    yield from (1/(error_y**2 + d**2 * error_x**2) for d in dydx)
!
!
!def m_matrix(jacobian, weights):
!    """Compute the normal-equations matrix :math:`M = J^T W J`.
!
!    Parameters
!    ----------
!    jacobian : numpy.ndarray
!        Jacobian matrix of shape ``(n, p)``.
!    weights : numpy.ndarray
!        Square weights matrix of shape ``(n, n)`` (typically diagonal).
!
!    Returns
!    -------
!    numpy.ndarray
!        Normal-equations matrix of shape ``(p, p)``.
!    """
!    return np.dot(np.dot(jacobian.T, weights), jacobian)
!
!
!def error_params(jacobian, weights):
!    """Compute the variance of each parameter from the Jacobian and weights.
!
!    Calculates ``diag(inv(J^T W J))`` which represents the variances of the
!    fitted parameters under the weighted least-squares model.
!
!    Parameters
!    ----------
!    jacobian : numpy.ndarray
!        Jacobian matrix of shape ``(n, p)``.
!    weights : numpy.ndarray
!        Square weights matrix of shape ``(n, n)`` (typically diagonal).
!
!    Returns
!    -------
!    numpy.ndarray
!        1-D array of length *p* containing the variance of each parameter.
!    """
!    M = m_matrix(jacobian, weights)
!    return np.diag(np.linalg.inv(M))
!
!
!def fit_sigma(residuals: np.ndarray, weights: np.ndarray, npoints: int, nparams: int) -> float:
!    """Calculate the fit's sigma value for a given set of residuals and weights."""
!    return np.sqrt(np.sum(weights*residuals**2))/(npoints-nparams)
