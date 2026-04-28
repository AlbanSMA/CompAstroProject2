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

    ! Force quadruple precision for large numbers
    !integer, parameter :: dp = selected_real_kind(20, 900)


    ! Matrices initialisation
    real :: J(7, 7), J_copy(7, 7)
    real :: A(7), A_copy(7), d(7), ns(7), new_ns(7), lam(7)

    ! Parameters for the matrix
    integer :: IPIV(7)
    integer, parameter :: N = 7, NRHS = 0, LDA = 7, LDB = 7
    integer :: INFO, i

    ! Variables : 
    ! number of steps (=resolution)
    ! number densities, dt, Temperature
    ! top and bottom temperature
    ! smst, div and tstep are to define arrays of dt and T9 later
    integer, parameter :: steps = 1000000
    real :: n0, n1, n2, n3, n4, n5, n6
    real :: dt(steps), T9(steps)
    real :: tstart, tstop

    ! Get storing list and steps number
    real :: lst_nsdt(7, steps)
    real :: lst_dt(steps)

    ! Make data file to store the results
    character(len=*), parameter :: OUT_FILE = "data.txt"
    integer :: fu


!-----------------------------------------------------------------------------

    ! Initial values of ns
    n0 = 2e3
    n1 = 3e3
    n2 = 0
    n3 = 0
    n4 = 0
    n5 = 0
    n6 = 0
    ns = reshape((/n0, n1, n2, n3, n4, n5, n6/), shape(ns))


!-----------------------------------------------------------------------------

    ! Values of t, dt and T9, by creating range(start, stop, step) functions
        ! t and dt
    tstart = 1e-1
    tstop = 1.5e1
        
    dt = (tstop-tstart)/steps

            ! list of ts
    lst_dt(1) = tstart
    do i=2, steps
        lst_dt(i) = lst_dt(i-1)+dt(i)
    end do

    T9 = ((1.3e10)/sqrt(lst_dt))/1e9

!-----------------------------------------------------------------------------

    ! Initialise for the loop
    d = 0.
    new_ns = ns+d
    lst_nsdt(:,1) = ns

    do i=1, steps-1
        ! Get the reaction rates
        call get_lam(T9(i), lam, dt)
        lam = lam*dt(i)
        !print*, lam

        ! Get A and its Jacobian
        call make_A(new_ns, lam, A)
        A_copy = A
        !print*,A
    
        !call make_J(new_ns, lam, J)
        J_copy = J
        
        ! Solve the matrix
        call dgesv(N, NRHS, J_copy, LDA, IPIV, A_copy, LDB, INFO)

        ! get d = ns_(i+1) - ns from the result, which is A_copy
        d = A_copy
        print*,"d",d
        print*,"ns",ns

        ! new_ns is ns + d
        new_ns = ns + d
        !print*, new_ns

        ! store the results in a list for later
        lst_nsdt(:,i+1) = new_ns

        ! reinitialise
        ns = new_ns
    end do

!-----------------------------------------------------------------------------
    ! Now save this data to a text file
    open(action="write", file=OUT_FILE, newunit=fu, status="replace")
    do i=1, steps
        write(fu,*) lst_dt(i), lst_nsdt(:,i)
    end do
    close(fu)

end program
