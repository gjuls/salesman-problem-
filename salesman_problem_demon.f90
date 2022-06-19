module functions
    implicit none
    public :: place_cities, exchange, factor
    integer, public :: N
    integer, public, parameter :: double = selected_real_kind(13)

    real(kind = double), dimension(:,:), allocatable :: distances

contains

    !initialises cities coordinates and distances
    subroutine place_cities(initial_distance, X, Y)
        integer :: sizer, i, ii, iii
        integer, dimension(:), allocatable :: seed
        real(kind = double) :: L
        real(kind = double), intent(out) :: initial_distance
        real(kind = double), dimension(:), intent(out), allocatable :: X, Y
        open(unit = 1, file="place_cities.dat", status="replace", action="write")

        print*, "Insert number of cities :"
        read*, N
        
        L = sqrt(real(N))

        call random_seed(sizer)
        allocate(seed(sizer))
        call random_seed(put = seed)
        allocate(X(N))
        allocate(Y(N))
        allocate(distances(N, N))
        X = 0
        Y = 0
        distances = 0
        !coordinates generation in a LxL square
        call random_number(X)
        X = X*L
        call random_number(Y)
        Y = Y*L
        !distances calculation
        do i = 1, N
            !print the coordinates
            write(unit = 1, fmt = *) i, X(i), Y(i)
            do ii = 1, N
            !if (i > ii)   VOLENDO SFRUTTARE LA SIMMETRIA
            distances(i, ii) = sqrt( (X(ii) - X(i))**2 + (Y(ii) - Y(i) )**2 )
            end do
            !final point == initial one
            if (i == N) write(unit = 1, fmt = *) 1, X(1), Y(1)
        end do
        !calculate initial distance
        initial_distance = 0
        do iii = 1, N-1
            initial_distance = initial_distance + distances(iii, iii + 1)
        end do
        initial_distance = initial_distance + distances(N, 1)
        print*, "Initial distance :", initial_distance
        !call system('gnuplot -p place_cities.gpl')
        close(1)
        
    end subroutine place_cities


    subroutine exchange(travel, annealing_distance)
        integer, dimension(N) :: travel
        real(kind = double), dimension(2) :: rnd
        integer :: index1, index2, tmp1, tmp2, ii
        real(kind = double), intent(out) :: annealing_distance
        annealing_distance = 0
        call random_number(rnd)
        
        index1 = 1 + FLOOR(N*rnd(1))
        index2 = 1 + FLOOR(N*rnd(2))
        tmp1 = travel(index1)
        tmp2 = travel(index2)
        travel(index1) = tmp2
        travel(index2) = tmp1

        do ii = 1, N-1
            annealing_distance = annealing_distance + distances(travel(ii), travel(ii+1))
        end do
        annealing_distance = annealing_distance + distances(travel(N), travel(1))

    end subroutine exchange


    recursive function factor(n) result (factorial_result)
        integer, intent (in) :: n
        integer :: factorial_result
 
         if (n <= 0) then
         factorial_result = 1
            else
         factorial_result = n*factor(n-1)
         end if
    end function factor

end module functions


program salesman_program_brute
    use functions
    real(kind = double), dimension(:), allocatable ::  X, Y
    integer ::  ii, iii, nmcs, accept, maxbin
    integer, dimension(:), allocatable :: travel, histo
    integer, dimension(:), allocatable :: travel_minimum
    real(kind = double) :: initial_distance, i_distance, ii_distance
    real(kind = double) :: start, finish, dE, delta
    real(kind = double) :: sum_energy, sum_demon, energy, absolute_minimum, demon

    open(unit = 8, file = "demon_minimum_perT.dat", status="replace", action="write")
    open(unit = 9, file = "histo_demon.dat", status = 'replace', action = 'write')
    open(unit = 10, file = "initial_config.dat", status = 'replace', action = 'write')
    open(unit = 11, file = "final_config.dat", status = 'replace', action = 'write')
    open(unit = 13, file = "guo_demon.dat", status = 'replace', action = 'write')
    

    !generate initial condition
    call place_cities(initial_distance, X, Y)

    do j = 1, N
        write(unit = 10, fmt = *) j, X(j), Y(j)
        if (j == N) write(unit = 10, fmt = *) N, X(1), Y(1)
    end do
    
    !histogram
    maxbin = 50*N
    delta = 0.1
    allocate(histo(0 : + maxbin))
    histo = 0

    !initialisations
    accept = 0
    nmcs = N*100000                                                                             ! <------------
    sum_energy = 0
    sum_demon = 0
    demon = sqrt(real(N)) / 4   !from book init.

    !print*, 'Initial value of demon :', demon

    allocate(travel(N)) 
    allocate(travel_minimum(N))
  
    !initial order of visit - ordine di generazione
    travel = 0
    do ii = 1, N
        travel(ii) = ii
    end do    
    
    i_distance = initial_distance
    ii_distance = 0
    energy = initial_distance

    demon = sqrt(real(N))/4

    do while(demon > 0.5)
        absolute_minimum = 1000

        call cpu_time(start)


        do iii = 1, nmcs
            sum_demon = 0
            sum_energy = 0
            call exchange(travel, ii_distance)

            dE =  (ii_distance - i_distance) 
            !cerco e salvo un nuovo minimo
            if ( dE <=  0.0 ) then
                demon = demon + abs(dE)

                accept = accept + 1
                energy = energy  - abs(dE)

                sum_energy = sum_energy + energy
                sum_demon = sum_demon + demon

            else if (dE > 0 .and. demon >= 0) then
                demon = demon - abs(dE)
                energy = energy + abs(dE)

                accept = accept +1

                sum_energy = sum_energy + energy
                sum_demon = sum_demon + demon

                !check that the istantaneous sum is constant = demon_0 + d_0
                write(12,*) energy + demon
            else 
            !no increase, no new config

            end if

            if (ii_distance < absolute_minimum) then
                absolute_minimum = ii_distance
                travel_minimum = travel
                !for a single run - guo
                write(13, *) iii, absolute_minimum

            end if

            !prepare new config
            i_distance = ii_distance

        end do

        call cpu_time(finish)

        !minimum for each T
        write(8,*) demon, absolute_minimum

        ibin = nint(absolute_minimum/delta)
        if (ibin < maxbin) histo(ibin) = histo(ibin) + 1

        demon = demon - 0.1
        
    end do

    print*, 'Absolut Minimum :', absolute_minimum
    !print the final configuration
    do j = 1, N
        write(unit = 11, fmt = *) travel_minimum(j), X(travel_minimum(j)), Y(travel_minimum(j))
        if (j == N) write(unit = 11, fmt = *) travel_minimum(1), X(travel_minimum(1)), Y(travel_minimum(1))
    end do

    !file histogram
    do ibin = 1 , maxbin
        write(9 , *) ibin*delta, histo(ibin)
    end do

    !check results
    !print*, sum_energy/nmcs/( (   (sqrt(real(N))/4 - demon)/0.1) * nmcs ) , sum_demon/nmcs/(((sqrt(real(N))/4 - demon)/0.1)*nmcs),&
    !                 accept/nmcs/( ((sqrt(real(N))/4 - demon)/0.1)*nmcs )

    close(9)
    close(10)
    close(11)

    call system('gnuplot -p histo_demon.gpl')
    call system('gnuplot -p final_config.gpl')
    call system('gnuplot -p initial_config.gpl')

    print*, "Time :", finish - start, "s"

end program 