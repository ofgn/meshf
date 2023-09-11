module error_handling
   implicit none
   type :: error_type
      integer :: code
      character(len=128) :: procedure
      character(len=512) :: message
   end type error_type
contains
   subroutine handle_error(err)
      type(error_type), intent(in) :: err
      write (*, "(A, I0)") "Error Code: ", err%code
      write (*, "(2A)") "Procedure: ", trim(err%procedure)
      write (*, "(2A)") "Message: ", trim(err%message)
      stop
   end subroutine handle_error
end module error_handling