module station_data_mod
  ! Общие типы данных для всего проекта: описание станции и массив VTEC-наблюдений.
  ! Модуль не содержит вычислений, только структуры, которыми обмениваются остальные блоки.
  use my_prec, only: mp
  implicit none

  type :: station_meta
    ! Метаданные GNSS-приемника из station_coordinates_from_tec.csv.
    character(len=16) :: code = ''
    real(mp) :: lat_deg = 0.0_mp
    real(mp) :: lon_deg = 0.0_mp
    real(mp) :: height_m = 0.0_mp
  end type station_meta

  type :: observation_set
    ! Набор наблюдений после чтения файлов TayAbsTec.
    ! vtec(station,time) хранит VTEC, has_value показывает, есть ли значение.
    integer :: nstations = 0
    integer :: ntimes = 0
    type(station_meta), allocatable :: stations(:)
    real(mp), allocatable :: times(:)
    real(mp), allocatable :: vtec(:, :)
    logical, allocatable :: has_value(:, :)
  end type observation_set

end module station_data_mod
