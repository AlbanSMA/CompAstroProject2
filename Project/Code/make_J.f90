subroutine make_J(ns, lam, J)
    real, intent(in) :: ns(7), lam(7)
    real, intent(out) :: J(7,7)
    ! Declare the variables
    real :: J00, J01, J02, J03, J04, J05, J06
    real :: J10, J11, J12, J13, J14, J15, J16
    real :: J20, J21, J22, J23, J24, J25, J26
    real :: J30, J31, J32, J33, J34, J35, J36
    real :: J40, J41, J42, J43, J44, J45, J46
    real :: J50, J51, J52, J53, J54, J55, J56
    real :: J60, J61, J62, J63, J64, J65, J66

    ! First row
    J00 = -lam(1)*ns(2) - lam(5)*ns(3)
    J01 = -lam(1)*ns(1) + lam(7)*ns(6)
    J02 = 2*lam(2)*ns(3) + lam(3)*ns(6) - lam(5)*ns(1)
    J03 = 0
    J04 = 0
    J05 = lam(3)*ns(3) + lam(7)*ns(2)
    J06 = 0

    J(:,1) = (/J00, J01, J02, J03, J04, J05, J06/)

    ! Second row
    J10 = -lam(1)*ns(2)
    J11 = -lam(1)*ns(1) - lam(7)*ns(6)
    J12 = lam(4)*ns(5) + 2*lam(6)*ns(3)
    J13 = 0
    J14 = lam(4)*ns(3)
    J15 = -lam(7)*ns(2)
    J16 = 0

    J(:,2) = (/J10, J11, J12, J13, J14, J15, J16/)

    ! Third row
    J20 = lam(1)*ns(2) - lam(5)*ns(3)
    J21 = lam(1)*ns(1)
    J22 = -2*lam(2)*ns(3) - lam(3)*ns(6) - lam(4)*ns(5) - lam(5)*ns(1) - 2*lam(6)*ns(3)
    J23 = 0
    J24 = -lam(4)*ns(3)
    J25 = -lam(3)*ns(3)
    J26 = 0

    J(:,3) = (/J20, J21, J22, J23, J24, J25, J26/)

    ! Fourth row
    J30 = lam(1)*ns(2) + lam(5)*ns(3)
    J31 = lam(1)*ns(1)
    J32 = lam(5)*ns(1)
    J33 = 0
    J34 = 0
    J35 = 0
    J36 = 0

    J(:,4) = (/J30, J31, J32, J33, J34, J35, J36/)

    ! Fifth row
    J40 = 0
    J41 = lam(7)*ns(6)
    J42 = 2*lam(2)*ns(3) - lam(4)*ns(5)
    J43 = 0
    J44 = -lam(4)*ns(3)
    J45 = lam(7)*ns(2)
    J46 = 0

    J(:,5) = (/J40, J41, J42, J43, J44, J45, J46/)

    ! Sixth row
    J50 = lam(5)*ns(3)
    J51 = -lam(7)*ns(6)
    J52 = -lam(3)*ns(6) + 2*lam(6)*ns(3) + lam(5)*ns(1)
    J53 = 0
    J54 = 0
    J55 = -lam(3)*ns(3) - lam(7)*ns(2)
    J56 = 0

    J(:,6) = (/J50, J51, J52, J53, J54, J55, J56/)

    ! Seventh row
    J60 = 0
    J61 = 0
    J62 = lam(3)*ns(6) + lam(4)*ns(5)
    J63 = 0
    J64 = lam(4)*ns(3)
    J65 = lam(3)*ns(3)
    J66 = 0

    J(:,7) = (/J60, J61, J62, J63, J64, J65, J66/)

    J = reshape(J, shape(J))
    return
end subroutine