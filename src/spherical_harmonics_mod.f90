module spherical_harmonics_mod
  use my_prec, only: mp
  implicit none

  real(mp), parameter :: pi = acos(-1.0_mp)

contains

  integer pure function n_sph_coeff(lmax)
    integer, intent(in) :: lmax
    n_sph_coeff = (lmax + 1) * (lmax + 1)
  end function n_sph_coeff

  subroutine real_sph_basis(lmax, lat_deg, lon_deg, basis)
    integer, intent(in) :: lmax
    real(mp), intent(in) :: lat_deg, lon_deg
    real(mp), intent(out) :: basis(:)
    real(mp), allocatable :: p(:, :)
    real(mp) :: lat_rad, lon_rad, x, norm
    integer :: l, m, k

    if (size(basis) /= n_sph_coeff(lmax)) error stop 'Wrong basis size'

    allocate(p(0:lmax, 0:lmax))
    lat_rad = lat_deg * pi / 180.0_mp
    lon_rad = lon_deg * pi / 180.0_mp
    x = sin(lat_rad)

    call associated_legendre(lmax, x, p)

    k = 0
    do l = 0, lmax
      k = k + 1
      norm = sph_norm(l, 0)
      basis(k) = norm * p(l, 0)
      do m = 1, l
        norm = sph_norm(l, m)
        k = k + 1
        basis(k) = norm * p(l, m) * cos(real(m, mp) * lon_rad)
        k = k + 1
        basis(k) = norm * p(l, m) * sin(real(m, mp) * lon_rad)
      end do
    end do
  end subroutine real_sph_basis

  function sph_norm(l, m) result(norm)
    integer, intent(in) :: l, m
    real(mp) :: norm, log_ratio

    log_ratio = log_factorial(l - m) - log_factorial(l + m)
    norm = sqrt(real(2 * l + 1, mp) * exp(log_ratio))
    if (m > 0) norm = norm * sqrt(2.0_mp)
  end function sph_norm

  function log_factorial(n) result(value)
    integer, intent(in) :: n
    real(mp) :: value
    integer :: i

    value = 0.0_mp
    do i = 2, n
      value = value + log(real(i, mp))
    end do
  end function log_factorial

  subroutine associated_legendre(lmax, x, p)
    integer, intent(in) :: lmax
    real(mp), intent(in) :: x
    real(mp), intent(out) :: p(0:, 0:)
    real(mp) :: somx2, fact
    integer :: l, m

    p = 0.0_mp
    p(0, 0) = 1.0_mp
    if (lmax == 0) return

    somx2 = sqrt(max(0.0_mp, (1.0_mp - x) * (1.0_mp + x)))
    fact = 1.0_mp
    do m = 1, lmax
      p(m, m) = -fact * somx2 * p(m - 1, m - 1)
      fact = fact + 2.0_mp
    end do

    do m = 0, lmax - 1
      p(m + 1, m) = x * real(2 * m + 1, mp) * p(m, m)
    end do

    do m = 0, lmax
      do l = m + 2, lmax
        p(l, m) = (real(2 * l - 1, mp) * x * p(l - 1, m) - &
                   real(l + m - 1, mp) * p(l - 2, m)) / real(l - m, mp)
      end do
    end do
  end subroutine associated_legendre

  function evaluate_sph(lmax, coeff, lat_deg, lon_deg) result(value)
    integer, intent(in) :: lmax
    real(mp), intent(in) :: coeff(:), lat_deg, lon_deg
    real(mp) :: value
    real(mp), allocatable :: basis(:)

    allocate(basis(n_sph_coeff(lmax)))
    call real_sph_basis(lmax, lat_deg, lon_deg, basis)
    value = sum(coeff * basis)
  end function evaluate_sph

end module spherical_harmonics_mod
