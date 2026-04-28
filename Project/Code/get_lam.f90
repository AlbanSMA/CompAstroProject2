subroutine get_lam(T9, lam, dt)
    implicit none
    real, intent(in) :: T9, dt
    real, intent(out) :: lam(7)
    real :: lam1, lam2, lam3, lam4, lam5, lam6, lam7
    real :: rho_b, h
    h = 1e-6
    rho_b = h * (T9**3)

    lam1 = 2.5e4*rho_b
    lam2 = 3.9e8*rho_b * (T9**(-2./3)) * exp(-4.26 * (T9**(-1./3))) * (1 + 0.0979 * (T9**(1./3)) + 0.642 * (T9**(2./3)) + 0.440 * T9)
    lam3 = 2.6e9*rho_b * (T9**(-3./2)) * exp(-2.99 * (T9**(-1)))
    lam4 = 1.38e9*rho_b * (T9**(-3./2)) * exp(-0.745 * (T9**(-1)))
    lam5 = 2.23e3*rho_b * (T9**(-2./3)) * exp(-3.72 * (T9**(-1./3))) * (1 + 0.112 * (T9**(1./3)) + 3.38 * (T9**(2./3)) + 2.65 * T9)
    lam6 = 3.9e8*rho_b * (T9**(-2./3)) * exp(-4.26 * (T9**(-1./3))) * (1 + 0.0979 * (T9**(1./3)) + 0.642 * (T9**(2./3)) + 0.440 * T9)
    lam7 = 7.06e8*rho_b
    lam = dt*(/lam1, lam2, lam3, lam4, lam5, lam6, lam7/)
    return
end subroutine