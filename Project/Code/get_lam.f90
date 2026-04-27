subroutine get_lam(rho_b, T9, lam, dt)
    implicit none
    real, intent(in) :: rho_b, T9, dt
    real, intent(out) :: lam(1,7)
    real :: lam1, lam2, lam3, lam4, lam5, lam6, lam7

    lam1 = dt*(2.5*(10**4)*rho_b)
    lam2 = dt*(3.9*(10**8)*rho_b*(T9**(-2.0/3))*exp(-4.26*(T9**(-1.0/3)))*(1 + 0.0979*(T9**(1.0/3)) + 0.642*(T9**(2.0/3)) + 0.440*T9))
    lam3 = dt*(2.60*(10**9)*rho_b*(T9**(-3.0/2))*exp(-2.99*(T9**(-1))))
    lam4 = dt*(1.38*(10**9)*rho_b*(T9**(-3.0/2))*exp(-0.745*(T9**(-1))))
    lam5 = dt*(2.23*(10**3)*rho_b*(T9**(-2.0/3))*exp(-3.72*(T9**(-1.0/3)))*(1 + 0.112*(T9**(1.0/3)) + 3.38*(T9**(2.0/3)) + 2.65*T9))
    lam6 = dt*(3.9*(10**8)*rho_b*(T9**(-2.0/3))*exp(-4.26*(T9**(-1.0/3)))*(1 + 0.0979*(T9**(1.0/3)) + 0.642*(T9**(2.0/3)) + 0.440*T9))
    lam7 = dt*(7.06*(10**8)*rho_b)

    lam(1,:) = (/lam1, lam2, lam3, lam4, lam5, lam6, lam7/)
    return
end subroutine