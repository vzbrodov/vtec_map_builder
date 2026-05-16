module data_io_mod
  use my_prec, only: mp
  use station_data_mod, only: station_meta, observation_set
  use string_utils, only: lower_string, upper_string, replace_char
  implicit none

contains

  subroutine read_station_coordinates(filename, stations)
    character(len=*), intent(in) :: filename
    type(station_meta), allocatable, intent(out) :: stations(:)
    character(len=512) :: line, clean
    character(len=32) :: code
    integer :: unit, ios, n, i
    real(mp) :: lat, lon, h

    n = 0
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Cannot open station coordinate file'
    read(unit, '(A)', iostat=ios) line
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) n = n + 1
    end do
    rewind(unit)

    allocate(stations(n))
    read(unit, '(A)', iostat=ios) line
    i = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      clean = replace_char(line, ',', ' ')
      read(clean, *, iostat=ios) code, lat, lon, h
      if (ios /= 0) cycle
      i = i + 1
      stations(i)%code = upper_string(adjustl(code))
      stations(i)%lat_deg = lat
      stations(i)%lon_deg = lon
      stations(i)%height_m = h
    end do
    close(unit)

    if (i < n) stations = stations(:i)
  end subroutine read_station_coordinates

  subroutine load_observations(data_dir, doy, year, coords_file, obs)
    character(len=*), intent(in) :: data_dir, coords_file
    integer, intent(in) :: doy, year
    type(observation_set), intent(out) :: obs
    type(station_meta), allocatable :: all_stations(:), used(:)
    real(mp), allocatable :: tmp_times(:), tmp_vtec(:), vtec_work(:, :)
    logical, allocatable :: mask_work(:, :)
    character(len=512) :: path
    logical :: exists
    integer :: i, nused, ntimes, nrows

    call read_station_coordinates(coords_file, all_stations)
    allocate(used(size(all_stations)))
    nused = 0
    ntimes = 0

    do i = 1, size(all_stations)
      path = station_data_path(data_dir, all_stations(i)%code, doy, year)
      inquire(file=trim(path), exist=exists)
      if (.not. exists) cycle

      call read_station_vtec(trim(path), tmp_times, tmp_vtec)
      nrows = size(tmp_times)
      if (nrows == 0) cycle

      if (ntimes == 0) then
        ntimes = nrows
        allocate(obs%times(ntimes))
        obs%times = tmp_times
        allocate(vtec_work(size(all_stations), ntimes))
        allocate(mask_work(size(all_stations), ntimes))
        vtec_work = 0.0_mp
        mask_work = .false.
      end if

      nused = nused + 1
      used(nused) = all_stations(i)
      vtec_work(nused, 1:min(ntimes, nrows)) = tmp_vtec(1:min(ntimes, nrows))
      mask_work(nused, 1:min(ntimes, nrows)) = .true.
    end do

    obs%nstations = nused
    obs%ntimes = ntimes
    if (nused == 0 .or. ntimes == 0 .or. .not. allocated(vtec_work)) then
      allocate(obs%stations(0))
      allocate(obs%times(0))
      allocate(obs%vtec(0, 0))
      allocate(obs%has_value(0, 0))
      return
    end if

    allocate(obs%stations(nused))
    allocate(obs%vtec(nused, ntimes))
    allocate(obs%has_value(nused, ntimes))
    obs%stations = used(:nused)
    obs%vtec = vtec_work(:nused, :)
    obs%has_value = mask_work(:nused, :)
  end subroutine load_observations

  function station_data_path(data_dir, station_code, doy, year) result(path)
    character(len=*), intent(in) :: data_dir, station_code
    integer, intent(in) :: doy, year
    character(len=512) :: path
    character(len=16) :: code_lower

    code_lower = lower_string(trim(station_code))
    write(path, '(A,"/",A,"/",A,"_",I3.3,"_",I4,".dat")') &
      trim(data_dir), trim(code_lower), trim(code_lower), doy, year
  end function station_data_path

  subroutine read_station_vtec(filename, times, vtec)
    character(len=*), intent(in) :: filename
    real(mp), allocatable, intent(out) :: times(:), vtec(:)
    character(len=512) :: line
    integer :: unit, ios, n, i
    real(mp) :: t, val

    n = 0
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      allocate(times(0), vtec(0))
      return
    end if
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) t, val
      if (ios == 0) n = n + 1
    end do
    rewind(unit)

    allocate(times(n), vtec(n))
    i = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) t, val
      if (ios /= 0) cycle
      i = i + 1
      times(i) = t
      vtec(i) = val
    end do
    close(unit)
  end subroutine read_station_vtec

  subroutine write_coefficients(filename, times, lmax, coeffs)
    character(len=*), intent(in) :: filename
    real(mp), intent(in) :: times(:), coeffs(:, :)
    integer, intent(in) :: lmax
    integer :: unit, it, ic
    character(len=32) :: name

    open(newunit=unit, file=filename, status='replace', action='write')
    write(unit, '(A)', advance='no') 'time_ut'
    do ic = 1, size(coeffs, 1)
      name = coeff_name(ic, lmax)
      write(unit, '(",",A)', advance='no') trim(name)
    end do
    write(unit, *)
    do it = 1, size(times)
      write(unit, '(F10.4)', advance='no') times(it)
      do ic = 1, size(coeffs, 1)
        write(unit, '(",",ES18.9)', advance='no') coeffs(ic, it)
      end do
      write(unit, *)
    end do
    close(unit)
  end subroutine write_coefficients

  function coeff_name(index, lmax) result(name)
    integer, intent(in) :: index, lmax
    character(len=32) :: name
    integer :: l, m, k

    k = 0
    do l = 0, lmax
      k = k + 1
      if (k == index) then
        write(name, '("C_",I0,"_0")') l
        return
      end if
      do m = 1, l
        k = k + 1
        if (k == index) then
          write(name, '("C_",I0,"_",I0)') l, m
          return
        end if
        k = k + 1
        if (k == index) then
          write(name, '("S_",I0,"_",I0)') l, m
          return
        end if
      end do
    end do
    name = 'unknown'
  end function coeff_name

  subroutine write_grid(filename, times, lat_grid, lon_grid, maps)
    character(len=*), intent(in) :: filename
    real(mp), intent(in) :: times(:), lat_grid(:), lon_grid(:), maps(:, :, :)
    integer :: unit, it, ilat, ilon

    open(newunit=unit, file=filename, status='replace', action='write')
    write(unit, '(A)') 'time_ut,latitude,longitude,vtec'
    do it = 1, size(times)
      do ilat = 1, size(lat_grid)
        do ilon = 1, size(lon_grid)
          write(unit, '(F8.3,",",F9.3,",",F9.3,",",ES16.8)') &
            times(it), lat_grid(ilat), lon_grid(ilon), maps(ilon, ilat, it)
        end do
      end do
    end do
    close(unit)
  end subroutine write_grid

end module data_io_mod
