subroutine make_A(ns, lam, A)
    real, intent(in) :: ns(7), lam(7)
    real, intent(out) :: A(1,7)
    real :: dn0, dn1, dn2, dn3, dn4, dn5, dn6
    
    dn0 = -lam(1)*ns(1)*ns(2) + lam(2)*ns(3)**2 + lam(3)*ns(6)*ns(3) + lam(7)*ns(6)*ns(2) - lam(5)*ns(3)*ns(1)
    dn1 = -lam(1)*ns(1)*ns(2) + lam(4)*ns(5)*ns(3) + lam(6)*ns(3)**2 - lam(7)*ns(6)*ns(2)
    dn2 = lam(1)*ns(1)*ns(2) - lam(2)*ns(3)**2 - lam(3)*ns(6)*ns(3) - lam(4)*ns(5)*ns(3) - lam(5)*ns(3)*ns(1) - lam(6)*ns(3)**2
    dn3 = lam(1)*ns(1)*ns(2) + lam(5)*ns(3)*ns(1)
    dn4 = lam(2)*ns(3)**2 - lam(4)*ns(5)*ns(3) + lam(7)*ns(6)*ns(2)
    dn5 = -lam(3)*ns(6)*ns(3) + lam(6)*ns(3)**2 + lam(5)*ns(3)*ns(1) - lam(7)*ns(6)*ns(2)
    dn6 = lam(3)*ns(6)*ns(3) + lam(4)*ns(5)*ns(3)

    A = reshape((/dn0, dn1, dn2, dn3, dn4, dn5, dn6/), shape(A))
    return
end subroutine