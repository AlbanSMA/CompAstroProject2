program test
    implicit none

    abstract interface
        subroutine myproc(a)
            integer :: a
        end subroutine
    end interface

    procedure(myproc), pointer :: f

end program