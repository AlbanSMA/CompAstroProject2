subroutine get_T(time, T)
    real, intent(in) :: time
    real, intent(out) :: T
    T = (1.3e10/sqrt(time))!/1e9
end subroutine