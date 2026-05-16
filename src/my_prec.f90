module my_prec
  ! Единое место, где задается точность вещественных расчетов.
  ! Чтобы перейти, например, с real(8) на real(4), достаточно изменить mp здесь.
  implicit none

  integer, parameter :: mp = 8
end module my_prec
