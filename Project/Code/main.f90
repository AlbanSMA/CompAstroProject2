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
        ! A'(Y_i)*d = -A(Y_i)
        ! A'(Y_i) = abs(I - dt*J(Y_i))


        ! n0, n1, n2, n3, n4, n5, n6 : initial values of n
        ! dt : size of a step, in seconds
        ! temp : temperature
        ! rho, rho_0 : total density
        ! lst_nsdt : array to store ns after every step
        ! lst_dt : array to store t for saving the data

        ! OUT_FILE and fu : file and opener to save the data


    ! Matrices initialisation, density and reaction rates
    real(8) :: J(7, 7), Id(7,7), A_prime(7,7), A_min(7)
    real(8) :: A(7), d(7), ns(7), new_ns(7), lam(8)
    
    ! Parameters for the matrix
    integer :: IPIV(7)
    integer, parameter :: N = 7, NRHS = 0, LDA = 7, LDB = 7, steps=100000
    integer :: INFO, i, k, size

    ! Variables : 
    ! densities, time range, temperature
    ! total density, constant density outside of loop
    ! variables for rates : fr, dfrdt etc...
    ! save funcs from ratelib to list 
    real(8) :: n0, n1, n2, n3, n4, n5, n6
    real(8) :: tstart, tstop, time, temp, dt, rho, rho_0, NA, den
    !real(8) :: fr, dfrdt, dfrdd, rr, drrdt, drrdd
    type(proc_ptr) :: func(8)


    ! Get storing list and steps number
    real(8) :: lst_nsdt(7,steps)
    real(8) :: lst_dt(steps)

    ! Make data file to store the results
    character(len=*), parameter :: OUT_FILE = "data.txt"
    integer :: fu

!-----------------------------------------------------------------------------
    ! Initial values of ns and density
    rho = 2e-31
    NA = 6.02e23
    den = 1/(rho*NA)

    n0 = 0.85*den
    n1 = 0.15*den
    n2 = 0
    n3 = 0
    n4 = 0
    n5 = 0
    n6 = 0
    ns = (/n0, n1, n2, n3, n4, n5, n6/)

    ! Make a list of the rate functions from ratelib
    func(1)%ptr => rate_png
    func(2)%ptr => rate_ddp
    func(3)%ptr => rate_he3dp
    func(4)%ptr => rate_tdn
    func(5)%ptr => rate_dpg
    func(6)%ptr => rate_ddn
    func(7)%ptr => rate_tpn
    func(8)%ptr => rate_weaknp

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
    tstart = 1.
    tstop = 6e5

    dt = 1e-3

!-----------------------------------------------------------------------------
    ! Initialise for the loop
    lst_nsdt(:,1) = ns
    lst_dt(1) = tstart
    new_ns = ns
    A_prime = 1e5

    rho_0 = rho
    time = tstart

    i = 2
    ! Loop
    iloop : do while (time < tstop)

        ! Reset A_prime and dt
        A_prime = 1e5
        ! Adaptative time steps
        whileloop : do while (sum(A_prime) > 8.)
            dt = 0.98*dt
            time = lst_dt(i-1) + dt

            call get_T(time, temp)

            ! Get the reaction rates
            call get_lam(rho, lam, temp)
            lam = lam !*1e13 -> added to accelerate the nucleosynthesis

            ! Get A and its Jacobian
            call make_A(new_ns, lam, A)
            call make_J(new_ns, lam, J)

            ! A'(Y_i) = abs(I-dt*J(Y_i))
            A_prime = abs(Id - dt*J)

            ! Stop the iteration if dt is too small
            if (dt<((tstop-tstart)/steps)) then
                dt = ((tstop-tstart)/steps)
                exit
            end if
        end do whileloop

        print*, i
        print*, sum(A_prime)

        A_min = A

        ! Solve the matrix
        call dgesv(N, NRHS, A_prime, LDA, IPIV, A_min, LDB, INFO)
        if (INFO /= 0) then
            print*, "Unsolvable matrix"
            stop
        end if

        ! get d = ns_(i+1) - ns from the result, which is A
        d = A_min
        print*, d

        ! new_ns is ns + d
        new_ns = ns + d
        print*, new_ns

        ! store the results in a list for later
        lst_nsdt(:,i) = new_ns
        lst_dt(i) = time

        ! reinitialise
        ns = new_ns
        rho = rho_0
        i = i + 1

        if (i >= steps) then
            exit
        end if
    end do iloop

    size = i
    lst_nsdt = lst_nsdt(:, 1:size)

!-----------------------------------------------------------------------------
    ! Now save this data to a text file
    open(action="write", file=OUT_FILE, newunit=fu, status="replace")
    do i=1, steps/2
        write(fu,*) lst_dt(i), lst_nsdt(:,i)
    end do
    close(fu)

end program
