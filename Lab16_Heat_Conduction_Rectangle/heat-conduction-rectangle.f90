program heat_rectangle
    implicit none

    ! precision
    integer,parameter :: dp=selected_real_kind(15,307)

    ! task params

    real(dp), parameter :: lx=1.0_dp
    real(dp), parameter :: ly=1.0_dp
    real(dp), parameter :: T_total=5.0_dp
    real(dp), parameter :: a=1.0_dp

    ! grid params
    integer :: N=50,M=50,kmax=500, out_step=50

    ! counters
    integer :: i,j,k,out_count

    ! arrays
    real(dp), allocatable :: u(:,:), u_new(:,:), x(:),y(:)

    ! grid steps
    real(dp):: hx,hy,tau,tau_max
    real(dp):: t_current,rx,ry

    ! random coefs for funnies
    real(dp) :: alpha(0:5),beta(0:5),mu(0:5),nu(0:5)
    real(dp) :: p(5),q(5),r(5),s(5)
    real(dp) :: gamma(1:5),c(1:5)

    call random_number(alpha)
    call random_number(beta)
    call random_number(mu)
    call random_number(nu)
    alpha=alpha*lx
    beta=beta*lx
    mu=mu*ly
    nu=nu*ly
    ! all coefs are random and correspond to the lengths


    call random_number(p)
    call random_number(q)
    call random_number(r)
    call random_number(s)
    p=p*3
    q=q*3
    r=r*3
    s=s*3


    call random_number(gamma)
    call random_number(c)
    gamma=0.1+9.9*gamma
    c=0.1+1.9*c
    hx=lx/real(N,dp)
    hy=ly/real(M,dp)

    tau_max=1.0_dp/(2.0_dp*a**2*(1.0_dp/hx**(2)+1.0_dp*hy**(2)))
    tau=0.5_dp*tau_max

    print *, '=== Heat Equation in Rectangle ==='
    print *, 'Lx =', lx, 'Ly =', ly, 'T =', T_total
    print *, 'N =', N, 'M =', M, 'hx =', hx, 'hy =', hy
    print *, 'tau_max =', tau_max, 'tau =', tau

    ! arrays
    allocate(u(0:N,0:M))
    allocate(u_new(0:N,0:M))
    allocate(x(0:N))
    allocate(y(0:M))

    ! fill coords
    do j = 0, N
        x(j)=real(j,dp)*hx       
    end do
    do i = 0, M
        y(i)=real(i,dp)*hy
    end do

    ! init condition t=0
    do j=0,N
        do i=0,M
            u(j,i)=u0(x(j),y(i))
        end do
    end do


    ! Console output header
    print *, ''
    print '(A10, A15, A15, A15, A15)', 'Time', 'u(0,0)', 'u(Lx/2,Ly/2)', 'u(Lx,0)', 'u(Lx,Ly)'
    print '(A10, A15, A15, A15, A15)', '----------', '----------', '----------', '----------', '----------'
    print '(F10.4, 4F15.6)', 0.0_dp, u(0, 0), u(N/2, M/2), u(N, 0), u(N, M)

    ! ==================================================================
    ! Time stepping (explicit scheme 16.4)
    ! ==================================================================
    t_current = 0.0_dp
    out_count = 0

    do k=1,kmax
        t_current=t_current+tau

        if (t_current>T_total) exit

        ! u_new - unify 4 borders
        do i = 0, M
            u_new(0,i)=psi(1,y(i),t_current)
            u_new(N,i)=psi(2,y(i),t_current)
        end do

        do j = 0, N
            u_new(j,0)=psi(3,x(j),t_current)
            u_new(j,M)=psi(4,x(j),t_current)
        end do

        ! 16.5 scheme

        ! coefs so more UI pleasant
        rx=a**2*tau/(hx**2)
        ry=a**2*tau/(hy**2)
        do j=1,N-1
            do i=1,M-1
                u_new(j,i)=u(j,i) &
                + rx*(u(j+1,i)-2.0_dp*u(j,i)+u(j-1,i)) &
                + ry*(u(j,i+1)-2.0_dp*u(j,i)+u(j,i-1)) &
                + tau*f(x(j),y(i),t_current)
            end do
        end do
        
        ! update u
        u=u_new

        ! Console output (less frequent)
        if (mod(k, out_step) == 0 .or. k == kmax) then
            out_count = out_count + 1
            print '(F10.4, 4F15.6)', t_current, u(0, 0), u(N/2, M/2), u(N, 0), u(N, M)
        end if
    end do
    
    ! ==================================================================
    ! OPTIMIZED: Write ONLY FINAL STATE for heatmap
    ! ==================================================================
    open(unit=10, file='heatmap_data.dat', status='replace')
    write(10, *) '# x y u(x,y) - heatmap data (final state t=', T_total, ')'
    write(10, *) '# Assignment #16 - Optimized output'
    
    ! Write final state with blank lines between scans for gnuplot pm3d
    do i = 0, M
        do j = 0, N
            write(10, '(3F15.6)') x(j), y(i), u(j, i)
        end do
        write(10, *) ''  ! separator for gnuplot
    end do
    close(10)

    ! Generate gnuplot script for HEATMAP
    call create_heatmap_script()

    print *, ''
    print *, '=== Results saved ==='
    print *, 'Heatmap data: heatmap_data.dat'
    print *, 'Gnuplot script: plot_heatmap.gnu'
    print *, 'To view: gnuplot -persist plot_heatmap.gnu'

    ! Deallocate
    deallocate(u, u_new, x, y)

contains

    function eta(z) result(res)
        real(dp), intent(in) :: z
        integer :: res
        if (z>=0.0_dp) then
            res=1
        else 
            res=0
        end if
    end function eta

    function u0(x,y) result(res)
        real(dp), intent(in) :: x,y
        real(dp) :: res
        ! res=c(1)*eta(x-alpha(0))*eta(beta(0)-x) &
        ! *eta(y-mu(0))*eta(nu(0)-y)
        res=1
        
    end function u0

    function psi(i,var,t) result(res)
        integer, intent(in) :: i
        real(dp), intent(in) :: var,t
        real(dp) :: res
        res=0
        ! if (i==1 .or. i==2) then
        !     res=c(i)*eta(var-mu(i))*eta(nu(i)-var) &
        !     *(1.0_dp-exp(-gamma(i)*t))
        ! else if (i==3 .or. i==4) then
        !     res=c(i)*eta(var-alpha(i))*eta(beta(i)-var) &
        !     *(1.0_dp-exp(-gamma(i)*t))
        ! end if
    end function psi

    function f(x,y,t) result(res)
        real(dp), intent(in) :: x,y,t
        real(dp) :: res
        !res=c(5)*u0(x,y)*(1.0_dp-exp(-gamma(5)))
        res=32*(x*(1-x)+y*(1-y))
    end function f
 

    ! ==================================================================
    ! HEATMAP gnuplot script generator
    ! ==================================================================

    ! Gnuplot script generator
    subroutine create_heatmap_script()
        integer :: unit_gnu
        unit_gnu = 20
        open(unit=unit_gnu, file='plot_heatmap.gnu', status='replace')
        
        write(unit_gnu, *) '# Gnuplot script for 2D heatmap'
        write(unit_gnu, *) 'set terminal pngcairo size 1200,1000 enhanced font "Verdana,11"'
        write(unit_gnu, *) 'set output "heatmap.png"'
        write(unit_gnu, *) ''
        write(unit_gnu, *) 'set title "Heat Equation in Rectangle\\nAssignment #16" font ",16"'
        write(unit_gnu, *) 'set xlabel "x" font ",13"'
        write(unit_gnu, *) 'set ylabel "y" font ",13"'
        write(unit_gnu, *) ''
        write(unit_gnu, *) 'set grid front lt 0 lw 0.5 lc rgb "#cccccc"'
        write(unit_gnu, *) 'set border lw 1.5'
        write(unit_gnu, *) ''
        write(unit_gnu, *) '# HEATMAP VIEW (top-down)'
        write(unit_gnu, *) 'set view map'
        write(unit_gnu, *) 'set pm3d map'
        write(unit_gnu, *) ''
        write(unit_gnu, *) '# Color palette (viridis-like)'
        write(unit_gnu, *) 'set palette defined (0 "#440154", 0.2 "#3b528b", 0.4 "#21918c", 0.6 "#5ec962", 0.8 "#fde725", 1 "#ffffff")'
        write(unit_gnu, *) ''
        write(unit_gnu, *) '# Colorbar'
        write(unit_gnu, *) 'set cbrange [*:*]'
        write(unit_gnu, *) 'set colorbox vertical user origin 0.92, 0.15 size 0.03, 0.7'
        write(unit_gnu, *) 'set cblabel "u(x,y)" font ",11" offset 2,0'
        write(unit_gnu, *) ''
        write(unit_gnu, *) '# Plot'
        write(unit_gnu, *) 'splot "heatmap_data.dat" using 1:2:3 with pm3d title "Temperature"'
        write(unit_gnu, *) ''
        write(unit_gnu, *) 'print "Heatmap saved to: heatmap.png"'
        
        close(unit_gnu)
    end subroutine create_heatmap_script
    
end program heat_rectangle