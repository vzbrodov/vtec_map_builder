module string_utils
  ! Небольшие строковые функции, чтобы не размазывать обработку имен станций
  ! и CSV-строк по модулям ввода-вывода.
  implicit none

contains

  pure function lower_string(text) result(out)
    ! Возвращает строку в нижнем регистре. Используется для путей вида 2023/01kr.
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, c

    out = text
    do i = 1, len(text)
      c = iachar(out(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c + 32)
    end do
  end function lower_string

  pure function upper_string(text) result(out)
    ! Возвращает строку в верхнем регистре. Используется для сравнения кодов станций.
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, c

    out = text
    do i = 1, len(text)
      c = iachar(out(i:i))
      if (c >= iachar('a') .and. c <= iachar('z')) out(i:i) = achar(c - 32)
    end do
  end function upper_string

  pure function replace_char(text, old_char, new_char) result(out)
    ! Простая замена символов. Нужна для чтения CSV через обычный list-directed read.
    character(len=*), intent(in) :: text
    character(len=1), intent(in) :: old_char, new_char
    character(len=len(text)) :: out
    integer :: i

    out = text
    do i = 1, len(out)
      if (out(i:i) == old_char) out(i:i) = new_char
    end do
  end function replace_char

end module string_utils
