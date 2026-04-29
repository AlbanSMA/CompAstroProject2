subroutine get_rho(time, rho)
    real, intent(in) :: time
    real, intent(out) :: rho
    exp_rate = sqrt(time)/(1.3e10)
    rho = rho/exp_rate
end subroutine
