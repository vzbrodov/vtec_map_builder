module vtec_model_mod
  ! Высокоуровневые процедуры моделирования VTEC:
  ! - подгонка коэффициентов сферических гармоник для каждой эпохи;
  ! - регуляризация высоких степеней;
  ! - расчет карты на регулярной широтно-долготной сетке.
  use my_prec, only: mp
  use station_data_mod, only: observation_set
  use spherical_harmonics_mod, only: n_sph_coeff, real_sph_basis, evaluate_sph, pi
  use least_squares_mod, only: weighted_least_squares
  implicit none

contains

  subroutine fit_vtec_series(obs, lmax, ridge, weight_radius_deg, coeffs, residual_rms, nused)
    ! Для каждой эпохи независимо строит МНК-систему по всем доступным станциям.
    ! На выходе coeffs(:,it) - коэффициенты карты для времени obs%times(it).
    !
    ! ridge управляет гладкостью: чем больше, тем сильнее подавляются высокие степени.
    ! weight_radius_deg > 0 включает плотностные веса станций; 0 оставляет равные веса.
    type(observation_set), intent(in) :: obs
    integer, intent(in) :: lmax
    real(mp), intent(in) :: ridge, weight_radius_deg
    real(mp), allocatable, intent(out) :: coeffs(:, :), residual_rms(:)
    integer, allocatable, intent(out) :: nused(:)
    real(mp), allocatable :: a(:, :), y(:), w(:), coeff(:), basis(:), station_weights(:), ridge_diag(:)
    integer :: it, is, iobs, ncoef
    logical :: ok

    ncoef = n_sph_coeff(lmax)
    allocate(coeffs(ncoef, obs%ntimes), residual_rms(obs%ntimes), nused(obs%ntimes))
    allocate(a(obs%nstations, ncoef), y(obs%nstations), w(obs%nstations), coeff(ncoef), basis(ncoef))
    allocate(station_weights(obs%nstations))
    allocate(ridge_diag(ncoef))
    coeffs = 0.0_mp
    residual_rms = -1.0_mp
    nused = 0
    ! Веса станций не зависят от времени, поэтому считаются один раз.
    call compute_density_weights(obs, weight_radius_deg, station_weights)
    ! Диагональ регуляризации тоже одинакова для всех эпох.
    call make_degree_damping(lmax, ridge, ridge_diag)

    do it = 1, obs%ntimes
      ! Собираем матрицу A и вектор y только из станций, у которых есть VTEC в эту эпоху.
      iobs = 0
      do is = 1, obs%nstations
        if (.not. obs%has_value(is, it)) cycle
        iobs = iobs + 1
        call real_sph_basis(lmax, obs%stations(is)%lat_deg, obs%stations(is)%lon_deg, basis)
        a(iobs, :) = basis
        y(iobs) = obs%vtec(is, it)
        w(iobs) = station_weights(is)
      end do

      nused(it) = iobs
      ! Если наблюдений меньше коэффициентов, обычная МНК-задача недоопределена.
      if (iobs < ncoef) cycle

      call weighted_least_squares(a(:iobs, :), y(:iobs), w(:iobs), ridge_diag, coeff, ok)
      if (.not. ok) cycle
      coeffs(:, it) = coeff
      residual_rms(it) = compute_rms(a(:iobs, :), y(:iobs), coeff)
    end do
  end subroutine fit_vtec_series

  subroutine make_degree_damping(lmax, ridge, ridge_diag)
    ! Формирует регуляризацию по степени l.
    ! Штраф растет как [l(l+1)]^2, поэтому низкие глобальные формы почти свободны,
    ! а мелкомасштабные осцилляции высоких l подавляются сильнее.
    integer, intent(in) :: lmax
    real(mp), intent(in) :: ridge
    real(mp), intent(out) :: ridge_diag(:)
    integer :: l, m, k
    real(mp) :: degree_scale

    k = 0
    do l = 0, lmax
      degree_scale = real(l * (l + 1), mp)
      degree_scale = degree_scale * degree_scale

      k = k + 1
      ridge_diag(k) = ridge * degree_scale
      do m = 1, l
        k = k + 1
        ridge_diag(k) = ridge * degree_scale
        k = k + 1
        ridge_diag(k) = ridge * degree_scale
      end do
    end do
  end subroutine make_degree_damping

  subroutine compute_density_weights(obs, radius_deg, weights)
    ! Опциональные веса по плотности сети.
    ! Если radius_deg <= 0, все веса равны 1.
    ! Если включить, станции в плотном кластере получают меньший вес.
    ! В текущих экспериментах основной режим - radius_deg = 0, без разреживания сети.
    type(observation_set), intent(in) :: obs
    real(mp), intent(in) :: radius_deg
    real(mp), intent(out) :: weights(:)
    real(mp) :: radius_rad, density, total_weight
    integer :: i, j

    weights = 1.0_mp
    if (radius_deg <= 0.0_mp) return

    radius_rad = radius_deg * pi / 180.0_mp
    do i = 1, obs%nstations
      density = 0.0_mp
      do j = 1, obs%nstations
        density = density + exp(-0.5_mp * (angular_distance_rad( &
          obs%stations(i)%lat_deg, obs%stations(i)%lon_deg, &
          obs%stations(j)%lat_deg, obs%stations(j)%lon_deg) / radius_rad) ** 2)
      end do
      weights(i) = 1.0_mp / max(density, tiny(1.0_mp))
    end do

    ! Нормируем веса так, чтобы средний вес оставался около 1.
    total_weight = sum(weights)
    if (total_weight > 0.0_mp) weights = weights * real(obs%nstations, mp) / total_weight
  end subroutine compute_density_weights

  function angular_distance_rad(lat1_deg, lon1_deg, lat2_deg, lon2_deg) result(distance)
    ! Сферическое угловое расстояние между двумя точками в радианах.
    real(mp), intent(in) :: lat1_deg, lon1_deg, lat2_deg, lon2_deg
    real(mp) :: distance
    real(mp) :: lat1, lat2, dlon, c

    lat1 = lat1_deg * pi / 180.0_mp
    lat2 = lat2_deg * pi / 180.0_mp
    dlon = (lon1_deg - lon2_deg) * pi / 180.0_mp
    c = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon)
    distance = acos(max(-1.0_mp, min(1.0_mp, c)))
  end function angular_distance_rad

  function compute_rms(a, y, coeff) result(rms)
    ! RMS невязки на станциях для одной эпохи.
    ! Это диагностическая величина, а не единственный критерий качества карты.
    real(mp), intent(in) :: a(:, :), y(:), coeff(:)
    real(mp) :: rms
    real(mp) :: residual
    integer :: i

    rms = 0.0_mp
    do i = 1, size(y)
      residual = y(i) - sum(a(i, :) * coeff)
      rms = rms + residual * residual
    end do
    rms = sqrt(rms / real(max(1, size(y)), mp))
  end function compute_rms

  subroutine make_latlon_grid(lat_step, lon_step, lat_grid, lon_grid)
    ! Создает регулярную сетку широт [-90,90] и долгот [-180,180].
    real(mp), intent(in) :: lat_step, lon_step
    real(mp), allocatable, intent(out) :: lat_grid(:), lon_grid(:)
    integer :: nlat, nlon, i

    nlat = int(180.0_mp / lat_step) + 1
    nlon = int(360.0_mp / lon_step) + 1
    allocate(lat_grid(nlat), lon_grid(nlon))

    do i = 1, nlat
      lat_grid(i) = -90.0_mp + real(i - 1, mp) * lat_step
    end do
    do i = 1, nlon
      lon_grid(i) = -180.0_mp + real(i - 1, mp) * lon_step
    end do
  end subroutine make_latlon_grid

  subroutine evaluate_maps(lmax, coeffs, lat_grid, lon_grid, maps)
    ! Вычисляет VTEC на всей регулярной сетке для всех эпох.
    ! maps(lon_index, lat_index, time_index) затем пишется в CSV.
    integer, intent(in) :: lmax
    real(mp), intent(in) :: coeffs(:, :), lat_grid(:), lon_grid(:)
    real(mp), allocatable, intent(out) :: maps(:, :, :)
    integer :: it, ilat, ilon

    allocate(maps(size(lon_grid), size(lat_grid), size(coeffs, 2)))
    do it = 1, size(coeffs, 2)
      do ilat = 1, size(lat_grid)
        do ilon = 1, size(lon_grid)
          maps(ilon, ilat, it) = evaluate_sph(lmax, coeffs(:, it), lat_grid(ilat), lon_grid(ilon))
        end do
      end do
    end do
  end subroutine evaluate_maps

end module vtec_model_mod
