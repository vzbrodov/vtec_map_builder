module least_squares_mod
  use my_prec, only: mp
  implicit none

contains

  subroutine weighted_least_squares(a, y, w, ridge_diag, coeff, ok)
    real(mp), intent(in) :: a(:, :), y(:), w(:), ridge_diag(:)
    real(mp), intent(out) :: coeff(:)
    logical, intent(out) :: ok
    real(mp), allocatable :: normal(:, :), rhs(:)
    integer :: i, j, k, nobs, ncoef

    nobs = size(y)
    ncoef = size(coeff)
    allocate(normal(ncoef, ncoef), rhs(ncoef))
    normal = 0.0_mp
    rhs = 0.0_mp

    do i = 1, nobs
      do j = 1, ncoef
        rhs(j) = rhs(j) + w(i) * a(i, j) * y(i)
        do k = 1, ncoef
          normal(j, k) = normal(j, k) + w(i) * a(i, j) * a(i, k)
        end do
      end do
    end do

    do j = 1, ncoef
      normal(j, j) = normal(j, j) + ridge_diag(j)
    end do

    call solve_linear_system(normal, rhs, coeff, ok)
  end subroutine weighted_least_squares

  subroutine solve_linear_system(a, b, x, ok)
    real(mp), intent(inout) :: a(:, :)
    real(mp), intent(inout) :: b(:)
    real(mp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(mp) :: factor, pivot_abs, tmp
    integer :: n, i, j, k, pivot

    n = size(b)
    ok = .true.

    do k = 1, n - 1
      pivot = k
      pivot_abs = abs(a(k, k))
      do i = k + 1, n
        if (abs(a(i, k)) > pivot_abs) then
          pivot_abs = abs(a(i, k))
          pivot = i
        end if
      end do
      if (pivot_abs <= tiny(1.0_mp)) then
        ok = .false.
        x = 0.0_mp
        return
      end if
      if (pivot /= k) call swap_rows(a, b, k, pivot)
      do i = k + 1, n
        factor = a(i, k) / a(k, k)
        a(i, k) = 0.0_mp
        do j = k + 1, n
          a(i, j) = a(i, j) - factor * a(k, j)
        end do
        b(i) = b(i) - factor * b(k)
      end do
    end do

    if (abs(a(n, n)) <= tiny(1.0_mp)) then
      ok = .false.
      x = 0.0_mp
      return
    end if

    do i = n, 1, -1
      tmp = b(i)
      do j = i + 1, n
        tmp = tmp - a(i, j) * x(j)
      end do
      x(i) = tmp / a(i, i)
    end do
  end subroutine solve_linear_system

  subroutine swap_rows(a, b, i, j)
    real(mp), intent(inout) :: a(:, :), b(:)
    integer, intent(in) :: i, j
    real(mp) :: tmp
    real(mp), allocatable :: row(:)

    allocate(row(size(a, 2)))
    row = a(i, :)
    a(i, :) = a(j, :)
    a(j, :) = row
    tmp = b(i)
    b(i) = b(j)
    b(j) = tmp
  end subroutine swap_rows

end module least_squares_mod
