module station_data_mod
  use my_prec, only: mp
  implicit none

  type :: station_meta
    character(len=16) :: code = ''
    real(mp) :: lat_deg = 0.0_mp
    real(mp) :: lon_deg = 0.0_mp
    real(mp) :: height_m = 0.0_mp
  end type station_meta

  type :: observation_set
    integer :: nstations = 0
    integer :: ntimes = 0
    type(station_meta), allocatable :: stations(:)
    real(mp), allocatable :: times(:)
    real(mp), allocatable :: vtec(:, :)
    logical, allocatable :: has_value(:, :)
  end type observation_set

end module station_data_mod
