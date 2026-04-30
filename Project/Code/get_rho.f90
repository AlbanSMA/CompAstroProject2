subroutine get_rho(time, rho, rho_0)
    real(8), intent(in) :: time, rho_0
    real(8), intent(out) :: rho
    real(8) :: pi, cstt
    pi = 3.14159265
    cstt = 1.7e10

    exp_rate = sqrt(time)/cstt
    rho = (3*rho_0)/(4*pi*(exp_rate**3))
end subroutine
