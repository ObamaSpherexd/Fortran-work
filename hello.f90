program convection_solver_gnuplot
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307) ! 15 digits after point, up to ^307
    integer :: N, Nt, m, i
    integer :: print_step
    real(dp) :: l, c, a, p, r, h, tau, T
    real(dp), allocatable :: x(:), u(:), u_new(:) ! arrays of unknown length
    
    
   
    real(dp), allocatable :: ubig(:,:) 
    real(dp) :: t_new, pi_val

    ! === Parameters ===
    l = 2.0_dp
    c = 1.0_dp
    a = 0.7_dp
    p = 0.6_dp
    N = 100
    r = 0.8_dp
    
    
    print_step = 10 

    h = l / real(N, dp)
    tau = r * h / c
    T = 0.2_dp * l / c
    Nt = int(T / tau)
    if (Nt < 1) Nt = 1

    pi_val = acos(-1.0_dp)

    ! === Allocate arrays ===
    allocate(x(0:N)) ! arrays of known size
    allocate(u(0:N), u_new(0:N))
    allocate(ubig(0:N, 0:Nt))

    ! === Spatial grid ===
    do i = 0, N
        x(i) = real(i, dp) * h
    end do

    ! === Initial condition ===
    do i = 0, N
        u(i) = abs(sin(2.0_dp * pi_val * (x(i)/l)**4))**3
    end do
    ubig(:, 0) = u

    !openning a file for gnuplot data (x t u)

    open(unit=10, file='solution.dat', status='replace')
    
    
    write(10, *) "# x t u"

    ! t=0
    do i = 0, N
        write(10, '(3F15.6)') x(i), 0.0_dp, u(i)
    end do
    write(10, *) "" 

    ! initial state in terminal
    print *, "--- Time t = 0.0 ---"
    call print_console_slice(x, u, N, 5) ! space step 5

    ! === Time stepping (upwind scheme) ===
    do m = 0, Nt - 1
        t_new = real(m + 1, dp) * tau

        ! Left boundary
        u_new(0) = exp(-p * t_new**2)

        
        do i = 1, N
            u_new(i) = (r * u_new(i-1) + u(i) + tau * (x(i)/l) * (1.0_dp - x(i)/l) * exp(-a * t_new)) / (1.0_dp + r)
        end do

        u = u_new
        ubig(:, m+1) = u

        ! writing file for Gnuplot
        do i = 0, N
            write(10, '(3F15.6)') x(i), t_new, u(i)
        end do
        write(10, *) "" ! devider

        ! console log with steps
        if (mod(m + 1, print_step) == 0 .or. m == Nt - 1) then
            print *, "--- Time t =", t_new, "---"
            call print_console_slice(x, u, N, 5)
        end if
    end do

    close(10)

    ! === generating Gnuplot ===
    call create_gnuplot_script()

    print *, ' Data saved to solution.dat'
    print *, ' Gnuplot script created: plot_script.gnu'
    print *, ' To view 3D plot, run in terminal:'
    print *, '   gnuplot -persist plot_script.gnu'
    

contains

    ! subroutine for pickme style log
    subroutine print_console_slice(x_arr, u_arr, n, step_space)
        integer, intent(in) :: n, step_space
        real(dp), intent(in) :: x_arr(:), u_arr(:)
        integer :: k
        
        print '(A)', " Index |      X      |      U      "
        print '(A)', "-------|-------------|-------------"
        do k = 0, n, step_space
            print '(I5, " | ", F11.6, " | ", F11.6)', k, x_arr(k), u_arr(k)
        end do
        print *
    end subroutine print_console_slice

    !  Gnuplot script
    subroutine create_gnuplot_script()
        integer :: unit_gnu
        unit_gnu = 20
        open(unit=unit_gnu, file='plot_script.gnu', status='replace')
        
        write(unit_gnu, *) "set terminal wxt size 800,600 enhanced font 'Verdana,10'"
        write(unit_gnu, *) "set title 'Convection Solver: u(x,t)'"
        write(unit_gnu, *) "set xlabel 'x'"
        write(unit_gnu, *) "set ylabel 't'"
        write(unit_gnu, *) "set zlabel 'u'"
        write(unit_gnu, *) "set grid"
        write(unit_gnu, *) "set hidden3d" ! clearer image
        write(unit_gnu, *) "set dgrid3d 50,50 splines" ! evening the surface
        write(unit_gnu, *) "set style data lines"
        write(unit_gnu, *) ""
        write(unit_gnu, *) "# Plot as 3D surface"
        write(unit_gnu, *) " splot 'solution.dat' using 1:2:3 with lines title 'Solution'"

        
        close(unit_gnu)
    end subroutine create_gnuplot_script

end program convection_solver_gnuplot