program build_vtec_maps
  ! Главная программа намеренно оставлена тонкой:
  ! 1) прочитать параметры запуска;
  ! 2) загрузить наблюдения;
  ! 3) найти коэффициенты сферических гармоник;
  ! 4) посчитать карту на сетке;
  ! 5) записать CSV-результаты.
  use my_prec, only: mp
  use station_data_mod, only: observation_set
  use data_io_mod, only: load_observations, write_coefficients, write_grid
  use vtec_model_mod, only: fit_vtec_series, make_latlon_grid, evaluate_maps
  implicit none

  type(observation_set) :: obs
  real(mp), allocatable :: coeffs(:, :), residual_rms(:), lat_grid(:), lon_grid(:), maps(:, :, :)
  integer, allocatable :: nused(:)
  character(len=256) :: data_dir, coords_file, output_dir
  integer :: doy, year, lmax
  real(mp) :: ridge, lat_step, lon_step, weight_radius_deg

  ! Параметры можно передать через командную строку; если их нет, берутся значения ниже.
  call read_options(data_dir, coords_file, output_dir, doy, year, lmax, ridge, lat_step, lon_step, weight_radius_deg)
  call execute_command_line('mkdir -p ' // trim(output_dir))

  ! Загружаем VTEC и координаты станций.
  call load_observations(trim(data_dir), doy, year, trim(coords_file), obs)
  if (obs%nstations == 0 .or. obs%ntimes == 0) error stop 'No VTEC observations were loaded'

  ! Основной расчет: коэффициенты по эпохам и значения карты на регулярной сетке.
  call fit_vtec_series(obs, lmax, ridge, weight_radius_deg, coeffs, residual_rms, nused)
  call make_latlon_grid(lat_step, lon_step, lat_grid, lon_grid)
  call evaluate_maps(lmax, coeffs, lat_grid, lon_grid, maps)

  ! Результаты пишутся в CSV, чтобы их легко читать Python-скриптами.
  call write_coefficients(trim(output_dir) // '/coefficients.csv', obs%times, lmax, coeffs)
  call write_grid(trim(output_dir) // '/vtec_grid.csv', obs%times, lat_grid, lon_grid, maps)
  call write_diagnostics(trim(output_dir) // '/fit_diagnostics.csv', obs%times, nused, residual_rms)

  write(*, '(A,I0,A,I0,A)') 'Loaded ', obs%nstations, ' stations and ', obs%ntimes, ' epochs.'
  write(*, '(A,I0,A,I0,A)') 'lmax=', lmax, ', coefficients per epoch=', size(coeffs, 1), '.'
  write(*, '(A,F6.2,A)') 'density weight radius=', weight_radius_deg, ' deg.'
  write(*, '(A)') 'Wrote:'
  write(*, '(A)') '  ' // trim(output_dir) // '/coefficients.csv'
  write(*, '(A)') '  ' // trim(output_dir) // '/vtec_grid.csv'
  write(*, '(A)') '  ' // trim(output_dir) // '/fit_diagnostics.csv'

contains

  subroutine read_options(data_dir, coords_file, output_dir, doy, year, lmax, ridge, lat_step, lon_step, weight_radius_deg)
    ! Читает до 10 позиционных аргументов:
    ! data_dir coords_file output_dir doy year lmax lat_step lon_step ridge weight_radius_deg.
    ! Это простой интерфейс без внешних зависимостей и конфигурационных файлов.
    character(len=*), intent(out) :: data_dir, coords_file, output_dir
    integer, intent(out) :: doy, year, lmax
    real(mp), intent(out) :: ridge, lat_step, lon_step, weight_radius_deg
    character(len=256) :: arg

    ! Значения по умолчанию удобны для быстрого тестового запуска.
    data_dir = '2023'
    coords_file = 'station_coordinates_from_tec.csv'
    output_dir = 'output'
    doy = 320
    year = 2023
    lmax = 6
    ridge = 1.0e-8_mp
    lat_step = 5.0_mp
    lon_step = 5.0_mp
    weight_radius_deg = 0.0_mp

    call get_command_argument(1, arg)
    if (len_trim(arg) > 0) data_dir = trim(arg)
    call get_command_argument(2, arg)
    if (len_trim(arg) > 0) coords_file = trim(arg)
    call get_command_argument(3, arg)
    if (len_trim(arg) > 0) output_dir = trim(arg)
    call get_command_argument(4, arg)
    if (len_trim(arg) > 0) read(arg, *) doy
    call get_command_argument(5, arg)
    if (len_trim(arg) > 0) read(arg, *) year
    call get_command_argument(6, arg)
    if (len_trim(arg) > 0) read(arg, *) lmax
    call get_command_argument(7, arg)
    if (len_trim(arg) > 0) read(arg, *) lat_step
    call get_command_argument(8, arg)
    if (len_trim(arg) > 0) read(arg, *) lon_step
    call get_command_argument(9, arg)
    if (len_trim(arg) > 0) read(arg, *) ridge
    call get_command_argument(10, arg)
    if (len_trim(arg) > 0) read(arg, *) weight_radius_deg
  end subroutine read_options

  subroutine write_diagnostics(filename, times, nused, residual_rms)
    ! Диагностика по эпохам: сколько станций использовано и RMS невязки.
    ! По этому файлу удобно сравнивать варианты lmax/ridge.
    character(len=*), intent(in) :: filename
    real(mp), intent(in) :: times(:), residual_rms(:)
    integer, intent(in) :: nused(:)
    integer :: unit, i

    open(newunit=unit, file=filename, status='replace', action='write')
    write(unit, '(A)') 'time_ut,n_observations,residual_rms'
    do i = 1, size(times)
      write(unit, '(F8.3,",",I0,",",ES16.8)') times(i), nused(i), residual_rms(i)
    end do
    close(unit)
  end subroutine write_diagnostics

end program build_vtec_maps
