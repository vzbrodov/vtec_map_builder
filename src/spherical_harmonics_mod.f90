module spherical_harmonics_mod
  ! Вещественные сферические гармоники для представления глобальной карты VTEC.
  ! Для каждого l используются:
  !   m = 0:       C_l_0 * P_l0(sin(lat))
  !   m = 1..l:    C_l_m * P_lm(sin(lat))*cos(m*lon)
  !             +  S_l_m * P_lm(sin(lat))*sin(m*lon)
  ! Здесь используется нормированный базис, чтобы МНК был численно устойчивее.
  use my_prec, only: mp
  implicit none

  real(mp), parameter :: pi = acos(-1.0_mp)

contains

  integer pure function n_sph_coeff(lmax)
    ! Число вещественных коэффициентов до степени lmax включительно: (lmax+1)^2.
    integer, intent(in) :: lmax
    n_sph_coeff = (lmax + 1) * (lmax + 1)
  end function n_sph_coeff

  subroutine real_sph_basis(lmax, lat_deg, lon_deg, basis)
    ! Возвращает строку матрицы дизайна A для одной точки (lat, lon).
    ! basis(:) затем умножается на вектор коэффициентов, чтобы получить VTEC.
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
    ! Аргумент присоединенных полиномов Лежандра: sin(latitude).
    ! Это эквивалентно cos(colatitude).
    x = sin(lat_rad)

    call associated_legendre(lmax, x, p)

    ! Заполняем вещественный базис в том же порядке, в каком потом пишутся коэффициенты.
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
    ! Нормировка сферических гармоник без множителя 1/sqrt(4*pi).
    ! Она уменьшает разброс масштабов между разными l,m и делает МНК устойчивее.
    integer, intent(in) :: l, m
    real(mp) :: norm, log_ratio

    log_ratio = log_factorial(l - m) - log_factorial(l + m)
    norm = sqrt(real(2 * l + 1, mp) * exp(log_ratio))
    if (m > 0) norm = norm * sqrt(2.0_mp)
  end function sph_norm

  function log_factorial(n) result(value)
    ! log(n!) нужен для устойчивого вычисления нормировочного множителя.
    ! Через логарифмы не переполняемся на больших l.
    integer, intent(in) :: n
    real(mp) :: value
    integer :: i

    value = 0.0_mp
    do i = 2, n
      value = value + log(real(i, mp))
    end do
  end function log_factorial

  subroutine associated_legendre(lmax, x, p)
    ! Рекуррентное вычисление присоединенных полиномов Лежандра P_lm(x).
    ! Используется стандартная ненормированная рекурсия, нормировка применяется выше.
    integer, intent(in) :: lmax
    real(mp), intent(in) :: x
    real(mp), intent(out) :: p(0:, 0:)
    real(mp) :: somx2, fact
    integer :: l, m

    p = 0.0_mp
    ! Начальное значение P_00.
    p(0, 0) = 1.0_mp
    if (lmax == 0) return

    somx2 = sqrt(max(0.0_mp, (1.0_mp - x) * (1.0_mp + x)))
    fact = 1.0_mp
    ! Сначала диагональные элементы P_mm.
    do m = 1, lmax
      p(m, m) = -fact * somx2 * p(m - 1, m - 1)
      fact = fact + 2.0_mp
    end do

    ! Затем элементы P_(m+1),m.
    do m = 0, lmax - 1
      p(m + 1, m) = x * real(2 * m + 1, mp) * p(m, m)
    end do

    ! Остальная таблица строится по l.
    do m = 0, lmax
      do l = m + 2, lmax
        p(l, m) = (real(2 * l - 1, mp) * x * p(l - 1, m) - &
                   real(l + m - 1, mp) * p(l - 2, m)) / real(l - m, mp)
      end do
    end do
  end subroutine associated_legendre

  function evaluate_sph(lmax, coeff, lat_deg, lon_deg) result(value)
    ! Вычисляет VTEC по готовым коэффициентам в одной точке сетки.
    integer, intent(in) :: lmax
    real(mp), intent(in) :: coeff(:), lat_deg, lon_deg
    real(mp) :: value
    real(mp), allocatable :: basis(:)

    allocate(basis(n_sph_coeff(lmax)))
    call real_sph_basis(lmax, lat_deg, lon_deg, basis)
    value = sum(coeff * basis)
  end function evaluate_sph

end module spherical_harmonics_mod
