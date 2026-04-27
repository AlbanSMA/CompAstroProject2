program main
    implicit none
    ! Declare the variables
        ! matrices A, J : matrix and Jacobian of dn/dt
        ! matrices d, ns : n_i+1 - n_i and list of n_i

        ! IPIV, INFO : for solver
        ! N, NRHS : dim 2 of J, dim 1 of A
        ! LDA, LDB : dim 1 of J, dim 2 of A
        ! lam : list of reaction rates
        ! steps : number of steps want to run for
        ! J*d = -A


        ! n0, n1, n2, n3, n4, n5, n6 : initial values of n
        ! dt : size of a step, in seconds
        ! T9 : list of temperatures over 10^9
        ! rho_b : density of baryons
        ! lst_nsdt : array to store ns after every step
        ! lst_dt : array to store t for saving the data

        ! OUT_FILE and fu : file and opener to save the data

    ! Matrices initialisation
    real :: J(7, 7), J_copy(7, 7)
    real :: A(1,7), A_copy(1,7), d(1,7), ns(1,7), new_ns(1,7), lam(1,7)

    ! Parameters for the matrix and number of steps
    integer :: IPIV(7)
    integer, parameter :: N = 7, NRHS = 1, LDA = 7, LDB = 7, steps = 3000
    integer :: INFO, i

    ! Variables : number densities, dt, Temperature and baryon density
    real :: n0, n1, n2, n3, n4, n5, n6
    real :: dt, T9(1, steps)
    real :: rho_b

    ! Get storing list and steps number
    real lst_nsdt(1, 7, steps)
    real lst_dt(1, steps)

    ! Make data file to store the results and plotting file
    character(len=*), parameter :: OUT_FILE = "data.txt"
    integer :: fu

    ! Initial values of ns
    n0 = 2.0*(10.0**5)
    n1 = 3.0*(10.0**5)
    n2 = 0
    n3 = 0
    n4 = 0
    n5 = 0
    n6 = 0
    ns = reshape((/n0, n1, n2, n3, n4, n5, n6/), shape(ns))
    print *, shape(ns)

    ! Values of rho_b, dt, and T9
    rho_b = 3*(10.0**(-31))
    dt = (0.1+100)/steps
    call linspace((10.0**11)/(10.0**9), (10.0**9)/(10.0**9), T9)
    print *, shape(T9), dt

    ! Initialise for the loop
    d = 0
    new_ns = ns+d
    do i=1, steps
        ! Get the reaction rates
        call get_lam(rho_b, T9(:,i), lam, dt)
        print *, shape(lam)

        ! Get A and its Jacobian
        call make_A(new_ns, lam, A)
        A_copy = -A
        print *, shape(A), shape(A_copy)
        print *, A_copy
    
        call make_J(new_ns, lam, J)
        J_copy = J
        print *, shape(J), shape(J_copy)
        
        ! Solve the matrix
        call dgesv(N, NRHS, J_copy, LDA, IPIV, A_copy, LDB, INFO)
        print *, shape(J_copy), shape(A_copy)

        ! get d = ns_(i+1) - ns from the result, which is A_copy
        d = A_copy
        print*, d
        new_ns = ns + d  ! new_ns is ns + d

        ! store the results in a list for later
        lst_nsdt(:,:,i) = new_ns
        lst_dt(:,i) = lst_dt(:,i-1)+dt
        ns = new_ns
        print *, ns
    end do

    ! Now save this data to a text file
    open(action="write", file=OUT_FILE, newunit=fu, status="replace")
    do i=1, steps
        write(fu,*) lst_dt(:,i), lst_nsdt(:,:,i)
    end do
    close(fu)

end program
