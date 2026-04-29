program main
    use interfaces
    use rates
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
    real :: J(7, 7), Id(7,7), A_prime(7,7), A_min(7)
    real :: A(7), d(7), ns(7), new_ns(7), lam(7)!, X(7)
    
    ! Parameters for the matrix
    integer :: IPIV(7)
    integer, parameter :: N = 7, NRHS = 0, LDA = 7, LDB = 7
    integer :: INFO, i, k

    ! Variables : 
    ! number of steps (=resolution)
    ! number densities, dt, Temperature
    ! top and bottom temperature
    ! smst, div and tstep are to define arrays of dt and T9 later
    integer, parameter :: steps = 5000
    real :: n0, n1, n2, n3, n4, n5, n6
    !real :: x0, x1, x2, x3, x4, x5, x6
    real :: tstart, tstop, time, T, dt, rho!, NA
    real :: fr, dfrdt, dfrdd, rr, drrdd
    type(proc_ptr) :: func(7)

    real, external :: rate_png


    ! Get storing list and steps number
    real :: lst_nsdt(7, steps)
    real :: lst_dt(steps)

    ! Make data file to store the results
    character(len=*), parameter :: OUT_FILE = "data.txt"
    integer :: fu


!-----------------------------------------------------------------------------
    ! Initial values of ns and density
    n0 = 0.3
    n1 = 0.7
    n2 = 0
    n3 = 0
    n4 = 0
    n5 = 0
    n6 = 0
    ns = (/n0, n1, n2, n3, n4, n5, n6/)

    rho = 1.

    ! Make a list of the rate functions from ratelib
    func(1)%ptr => rate_png
    func(2)%ptr => rate_ddp
    func(3)%ptr => rate_he3dp
    func(4)%ptr => rate_tdn
    func(5)%ptr => rate_dpg
    func(6)%ptr => rate_ddn
    func(7)%ptr => rate_tpn

    ! Id matrix:
    do i=1,7
        do k=1,7
        if (i /= k) then
            Id(i,k) = 0
        end if
        if (i == k) then
            Id(i,k) = 1
        end if
        end do
    end do

!-----------------------------------------------------------------------------

    ! Range of t and first dt
        ! t and dt
    tstart = 1e-1
    tstop = 1e2
        
        ! get first dt
    dt = 1e-5

!-----------------------------------------------------------------------------
    ! Initialise for the loop
    lst_nsdt(:,1) = ns
    new_ns = ns
    A_prime = 1e5

    ! Loop
    iloop : do i=1, steps-1
        ! Adaptative time steps and Newton-Raphson
        whileloop : do while (sum(A_prime) > 1.)
            dt = 0.9*dt
            time = lst_dt(i) + dt
            call get_T(time, T)
            print*, T

            call get_rho(time, rho)
            print*, rho

            ! Get the reaction rates
            do k=1, 7
                call func(k)(T, rho, fr, dfrdt, dfrdd, rr, drr1dd)
                lam(k) = fr
            end do
            !call get_lam(T9, lam, rho)
            !print*,lam

            ! Get A and its Jacobian
            call make_A(new_ns, lam, A)
            !A = (X/(rho*NA))*A
            call make_J(new_ns, lam, J)
            !J = (X/(rho*NA))*J

            !print*,A

            ! A'(Y_i) = I-dt*J(Y_i)
            A_prime = abs(Id - dt*J)
            !print*, sum(A_prime), sum(dt*J)

        end do whileloop

        ! Want -A, not A
        print *, i
        A_min = -A

        ! Solve the matrix
        call dgesv(N, NRHS, J, LDA, IPIV, A_min, LDB, INFO)
        if (INFO == 0) then
            print*, "Unsolvable matrix"
            stop
        end if


        ! get d = ns_(i+1) - ns from the result, which is A
        d = A
        print*, "here"

        ! new_ns is ns + d
        new_ns = ns + d

        ! store the results in a list for later
        lst_nsdt(:,i+1) = new_ns
        lst_dt(i+1) = lst_dt(i)+dt


        ! If tstop is reached, then exit
        if (lst_dt(i) >= tstop) then
            exit
            print*,"exit"
        end if

        ! reinitialise
        ns = new_ns
    end do iloop

!-----------------------------------------------------------------------------
    ! Now save this data to a text file
    open(action="write", file=OUT_FILE, newunit=fu, status="replace")
    do i=1, steps
        write(fu,*) lst_dt(i), lst_nsdt(:,i)
    end do
    close(fu)

end program
