program wave_equ
    implicit none
    ! precision
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! Parameters for initial and boundary conditions
    real(dp) :: c(4), n_val, p(2), v(2)
    real(dp), parameter :: pi = 3.141592653589793_dp
    integer, parameter :: N = 50, M = 30
    real(dp), parameter :: l = 1.0_dp, a = 1.0_dp
    
    ! Bounds for mu2(x)
    real(dp), parameter :: alpha = 0.3_dp, beta = 0.7_dp

    ! Grid and time variables
    !real(dp) :: h, tau, tau0, T0, T, courant, t
    real(dp)::h,t0,tau0,t,courant,tau
    integer :: k, kmax, j, frame_count

    real(dp), allocatable :: u_prev(:), u_now(:), u_next(:), x(:)

    ! --- RANDOM INITIALIZATION ---
    call random_number(c)
    call random_number(n_val) ! renamed to n_val to avoid confusion
    call random_number(p)
    call random_number(v)
    
    n_val = 1.0_dp + n_val * 2.0_dp
    c = 1.0_dp + c * 9.0_dp
    p = 0.2_dp + p * 1.8_dp
    v = 1.0_dp + v * 9.0_dp

    ! --- GRID PARAMS ---
    h = l / real(N, dp)
    T0 = l / a
    tau0 = h / a
    
    ! Courant number condition c = a*tau/h < 1. 
    ! tau = 0.5 * h/a => c = 0.5 (Stable)
    tau = 0.5_dp * tau0  
    T = 5.0_dp * T0
    courant = a * tau / h
    kmax = nint(T / tau)

    ! --- ALLOCATION ---
    allocate(x(0:N), u_prev(0:N), u_now(0:N), u_next(0:N))

    ! Fill spatial grid
    do j = 0, N
        x(j) = j * h
    end do

    ! --- INITIAL CONDITIONS (t=0) ---
    u_prev(0) = mu34(3, 0.0_dp)
    u_prev(N) = mu34(4, 0.0_dp)
    do j = 1, N - 1
        u_prev(j) = mu1(x(j))
    end do

    ! --- INITIAL CONDITIONS (t=tau) using Taylor expansion ---
    u_now(0) = mu34(3, tau)
    u_now(N) = mu34(4, tau)
    do j = 1, N - 1
        ! Formula: u(tau) = u(0) + tau*u_t(0) + 0.5*tau^2*u_tt(0)
        ! u_tt = a^2 * u_xx + f. Here f=0.
        u_now(j) = mu1(x(j)) + tau * mu2(x(j)) + &
                   0.5_dp * tau**2 * (a**2 * d2mu1_dx2(x(j)))
    end do

    ! --- SAVE INITIAL FRAME (t=0) ---
    frame_count = 0
    call save_frame(frame_count, 0.0_dp, x, u_prev, N)

    print *, '=== 1D Wave Equation Animation ==='
    print *, 'L =', l, 'a =', a
    print *, 'N =', N, 'Steps =', kmax
    print *, 'Saving frames to disk...'

    ! --- MAIN TIME CYCLE ---
    do k = 1, kmax
        t = k * tau

        ! Explicit scheme for internal points
        do j = 1, N - 1
            u_next(j) = 2.0_dp * (1.0_dp - courant**2) * u_now(j) - u_prev(j) + &
                        courant**2 * (u_now(j + 1) + u_now(j - 1))
        end do

        ! Boundary conditions at current time t
        u_next(0) = mu34(3, t)
        u_next(N) = mu34(4, t)

        ! Save frame every M steps
        if (mod(k, M) == 0) then
            frame_count = frame_count + 1
            call save_frame(frame_count, t, x, u_next, N)
        end if

        ! Shift arrays: u_prev <- u_now <- u_next
        u_prev = u_now
        u_now = u_next
    end do

    ! --- FINALIZE ANIMATION ---
    call create_animation_script(frame_count)

    print *, 'Calculation complete.'
    print *, 'Frames saved: ', frame_count + 1
    print *, 'To view animation, run: gnuplot plot_animation.gnu'
    print *, 'To create GIF (requires ImageMagick): bash animate_wave.sh'

    deallocate(x, u_prev, u_now, u_next)

contains

    ! --- USER FUNCTIONS ---

    ! Second derivative of mu1: d2/dx2 (c1 * sin(pi*n*x/l))
    ! = -c1 * (pi*n/l)^2 * sin(...)
    function d2mu1_dx2(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res
        res = -c(1) * (pi * n_val / l)**2 * sin(pi * n_val * x / l)
    end function d2mu1_dx2

    ! Initial displacement u(x,0)
    function mu1(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res
        res = c(1) * sin(pi * n_val * x / l)
    end function mu1

    ! Initial velocity u_t(x,0)
    function mu2(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res
        if (x >= alpha .and. x <= beta) then
            res = c(2)
        else
            res = 0.0_dp
        end if
    end function mu2

    ! Boundary conditions mu3(t) and mu4(t)
    ! Fixed: Array p has size 2, so for i=3 use p(1), for i=4 use p(2)
    function mu34(i, t) result(res)
        integer, intent(in) :: i
        real(dp), intent(in) :: t
        real(dp) :: res
        real(dp) :: p_idx
        
        ! Map i=3 -> index 1, i=4 -> index 2
        p_idx = real(i - 2, dp) 
        
        res = c(i) * (1.0_dp - exp(-p(i - 2) * t / T0))
    end function mu34

    ! --- ANIMATION SUBROUTINES ---

    subroutine save_frame(frame_num, t_val, x_arr, u_arr, N_nodes)
        integer, intent(in) :: frame_num, N_nodes
        real(dp), intent(in) :: t_val
        real(dp), intent(in) :: x_arr(0:N_nodes), u_arr(0:N_nodes)
        integer :: j
        character(len=50) :: filename
        
        write(filename, '("frame_", I3.3, ".dat")') frame_num
        
        open(unit=30, file=filename, status='replace')
        write(30, *) '# Frame ', frame_num, ' - Time t = ', t_val
        write(30, *) '# x u(x,t)'
        
        do j = 0, N_nodes
            write(30, '(2F15.6)') x_arr(j), u_arr(j)
        end do
        
        close(30)
    end subroutine save_frame

    subroutine create_animation_script(n_frames)
        integer, intent(in) :: n_frames
        integer :: unit_gnu, f
        real(dp) :: dt_print
        
        ! Calculate time step between frames for correct labels
        ! T = 5*T0, kmax steps, output every M steps
        ! n_frames = kmax / M
        ! dt_print = (M * tau)
        ! But we can just approximate or calculate strictly. 
        ! Here we assume linear progression.
        
        unit_gnu = 21
        open(unit=unit_gnu, file='plot_animation.gnu', status='replace')
        
        write(unit_gnu, *) '# Gnuplot script for 1D Wave Equation Animation'
        write(unit_gnu, *) 'set terminal pngcairo size 1200,600 enhanced font "Verdana,12"'
        write(unit_gnu, *) 'set output "anim_frame_%03d.png"'
        write(unit_gnu, *) ''
        write(unit_gnu, *) 'set title "1D Wave Equation - Time Evolution" font ",16"'
        write(unit_gnu, *) 'set xlabel "x" font ",13"'
        write(unit_gnu, *) 'set ylabel "u(x,t)" font ",13"'
        write(unit_gnu, *) ''
        write(unit_gnu, *) 'set grid lt 0 lw 0.5 lc rgb "#cccccc"'
        write(unit_gnu, *) 'set border lw 1.5'
        write(unit_gnu, *) 'set yrange [-2:2]' ! Adjust range as needed
        write(unit_gnu, *) ''
        
        write(unit_gnu, *) '# Animation loop'
        write(unit_gnu, '(A,I0,A)') 'do for [frame=0:', n_frames, '] {'
        
        write(unit_gnu, *) '    set output sprintf("anim_frame_%03d.png", frame)'
        
        ! Approximate time label (Assuming uniform dt between outputs)
        ! dt = M * tau. 
        ! Since we don't pass tau here easily, we use a generic label or simple index
        write(unit_gnu, *) '    set title sprintf("1D Wave Equation - Frame %d", frame) font ",16"'
        
        write(unit_gnu, *) '    filename = sprintf("frame_%03d.dat", frame)'
        write(unit_gnu, *) '    plot filename using 1:2 with lines lw 3 lc rgb "navy" title "u(x,t)"'
        write(unit_gnu, *) '}'
        write(unit_gnu, *) ''
        write(unit_gnu, *) 'print "Animation frames saved: anim_frame_*.png"'
        
        close(unit_gnu)

        ! Create shell script for GIF conversion
        open(unit=22, file='animate_wave.sh', status='replace')
        write(22, *) '#!/bin/bash'
        write(22, *) '# Script to create animated GIF from wave frames'
        write(22, *) 'echo "Creating animated GIF..."'
        write(22, '(A,I0,A)') 'convert -delay 10 -loop 0 anim_frame_000.png anim_frame_*.png wave_animation.gif'
        write(22, *) 'echo "Animation saved: wave_animation.gif"'
        close(22)
    end subroutine create_animation_script

end program wave_equ