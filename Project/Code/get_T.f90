subroutine get_T(time, T)
    real(8), intent(in) :: time
    real(8), intent(out) :: T
    T = (4e10/sqrt(time))
end subroutine