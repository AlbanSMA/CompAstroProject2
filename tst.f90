program tst
    implicit none
    integer, parameter :: dp = selected_real_kind(30, 1000)
    real(kind=dp) :: X


    X = exp(296.6_dp)
    print *, X
end program