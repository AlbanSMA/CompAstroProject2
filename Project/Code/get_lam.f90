subroutine get_lam(rho, lam, temp)
    real(8), intent(in) :: rho, temp
    real(8), intent(out) :: lam(8)
    real(8) :: lam0, lam1, lam2, lam3, lam4, lam5, lam6, lam7
    real(8) :: T9

    T9 = temp/1e9

    lam0 = 2.5e4
    lam1 = 3.9e8 * (T9**(-2./3)) * exp(-4.26*(T9**(-1./3))) * (1 + 0.0979 * (T9**(1./3)) + 0.642 * (T9**(2./3)) + 0.440 * T9)
    lam2 = 2.6e9 * (T9**(-3./2)) * exp(-2.99*(T9**((-1))))
    lam3 = 1.38e9 * (T9**(-3./2)) * exp(-0.745*(T9**(-1)))
    lam4 = 2.23e3 * (T9**(-2./3)) * exp(-3.72*(T9**(-1./3))) * (1 + 0.112*(T9**(1./3)) + 3.38*(T9**(2./3)) + 2.65*T9)
    lam5 = 3.9e8 * (T9**(-2./3)) * exp(-4.26*(T9**(-1./3))) * (1 + 0.0979 * (T9**(1./3)) + 0.642 * (T9**(2./3)) + 0.440 * T9)
    lam6 = 7.06e8
    lam7 = 1.0/888.0

    lam = rho*(/lam0, lam1, lam2, lam3, lam4, lam5, lam6, lam7/)
end subroutine