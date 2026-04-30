module interfaces
    implicit none

    abstract interface
      subroutine myproc(temp, den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
        real(8), intent(in)  :: temp, den
        real(8), intent(out) :: fr, dfrdt, dfrdd, rr, drrdt, drrdd
      end subroutine myproc
    end interface

    type :: proc_ptr
      procedure(myproc), pointer, nopass :: ptr => null()
    end type proc_ptr
end module interfaces