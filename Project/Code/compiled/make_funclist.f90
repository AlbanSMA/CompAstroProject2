module interfaces
    implicit none
abstract interface
    subroutine myproc(T, rho,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
        real, intent(in) :: T, rho
        real, intent(out) :: fr,dfrdt,dfrdd,rr,drrdt,drrdd
    end subroutine myproc
end interface

  type :: proc_ptr
    procedure(myproc), pointer, nopass :: ptr => null()
  end type proc_ptr
end module interfaces