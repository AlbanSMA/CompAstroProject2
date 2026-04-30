subroutine make_A(ns, lam, A)
    real(8), intent(in) :: ns(7), lam(8)
    real(8), intent(out) :: A(7)
    real(8) :: dn0, dn1, dn2, dn3, dn4, dn5, dn6
    
    dn0 = -lam(1)*ns(1)*ns(2) + lam(2)*(ns(3)**2) + lam(3)*ns(6)*ns(3) + lam(7)*ns(6)*ns(2) - lam(5)*ns(3)*ns(1) + lam(8)*ns(2)
    dn1 = -lam(1)*ns(1)*ns(2) + lam(4)*ns(5)*ns(3) + lam(6)*(ns(3)**2) - lam(7)*ns(6)*ns(2) - lam(8)*ns(2)
    dn2 = lam(1)*ns(1)*ns(2) - lam(2)*(ns(3)**2) - lam(3)*ns(6)*ns(3) - lam(4)*ns(5)*ns(3) - lam(5)*ns(3)*ns(1) - lam(6)*(ns(3)**2)
    dn3 = lam(1)*ns(1)*ns(2) + lam(5)*ns(3)*ns(1)
    dn4 = lam(2)*(ns(3)**2) - lam(4)*ns(5)*ns(3) + lam(7)*ns(6)*ns(2)
    dn5 = -lam(3)*ns(6)*ns(3) + lam(6)*(ns(3)**2) + lam(5)*ns(3)*ns(1) - lam(7)*ns(6)*ns(2)
    dn6 = lam(3)*ns(6)*ns(3) + lam(4)*ns(5)*ns(3)
    A = (/dn0, dn1, dn2, dn3, dn4, dn5, dn6/)
    return
end subroutine