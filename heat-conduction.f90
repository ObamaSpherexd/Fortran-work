program heat_conduction
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    integer :: N, Nt, i, m
    integer :: save_step,print_step
    real(dp) :: l, c, a2, T,c,phi1,phi2,b,f0,beta1,beta2
    real(dp) :: sigma, theta, nu,mu
    real(dp) allocatable :: x(:), u(:), u_new(:)
    real(dp) allocatable :: xi_coef(:), eta_coef(:)

    ! parameters
    l=2.0_dp ! length
    a2=1.0_dp ! conduction coef
    T=0.5_dp ! modeling time
    N=50 ! spatial steps
    c=1.0_dp ! start coef
    phi1=0.0_dp ! 1 border
    phi2=1.0_dp ! 2 border
    b=1.0_dp ! damping factor
    f0=1.0_dp ! source amp
    beta1=0.3_dp ! 1 source area
    beta2=0.7_dp ! 2 source area

    save_step=5 ! write every 5 steps
    print_step=10 ! print every 10 steps


    mu=2*h**2/(tau*a2)
    sugma=2+mu
    nu=2*h**2/a2
    theta=2-nu
   ! g(i)=-2*h**2/a2*f(x(i),t_new+tau/2)-(u(i+1)-2-)
    
    ! GRID
    h=l/real(N,dp)
    ! tau*a2/h=lambda=0.5 for example
    lambda=0.5_dp
    tau=2.0_dp*lambda*h**2/a2

    Nt=int(T/tau) ! time steps
    tau=T/real(nt,dp) ! correcting tau to fit T
    lambda=a2*tau/(2.0_dp*h**2)
    ! memory allocation
    allocate (x(0,N))
    allocate (u(0:N), u_new(0:N))
    allocate (xi_coef(0:N), eta_coef(0:N))

    ! GRID INIT
    do i=1,N
        x(i)=real(i,dp)*h
    end do

    ! init condition
    ! u0=c(l-x)/l
    do i=1,N
        u(i)=c*(l-x(i))/l
    end do

    ! opening file for data
    open(unit=10, file='heat_solution.dat',status='replace')
    write(10, *) "# x t u"

    ! t=0 data
    do i=0,N
        write(10, '(3F15.6)') x(i), 0.0_dp, u(i)
    end do
    write(10, *) ""

    print *, ' Heat conduction solver'
    print *, 'N= ', N, "Nt= ", Nt, "tau= ", tau
    call print_console_slice(x,u,N,5)

    ! Time cycle
    do m=0,Nt-1
        t_new=real(m,dp)*tau
        t_next=t_new+tau

        ! border coefs
        ! phi1=p1(1-e^-bt)
        ! phi2=p2
        phi1=p1*(1.0_dp-exp(-b*t_new))
        phi2=p2

        ! forward approach
        xi_coef(0)=0.0_dp
        eta_coef(0)=phi1

        do i =1,N
            if (i==0) then
                kn_val=phi1
            else if (i==N) then 
                kn_val=phi2
            else
                f_val=0.0_dp
                if (0<=beta1<=x(i)<=beta2<=l) then
                    f_val=f0*(1.0_dp-exp(-x(i)*(t_new+0.5_dp*tau)))
                end if

                kn_val=lambda*u(i-1)+(1.0_dp-2.0_dp*lambda)
 ! UNFINISHED, CHECK WITH THE BOOK, THERE HAS TO BE A BETTER WAY



    
    
        
        if (mod(i+1,save_step)==0 .or. m==Nt-1) then
            do i=1,N
                write(10,'(3F15.6)') x(i) t_next, u(i)
            end do
            write(10, *) "" ! block devider for plot
        end if

    ! console output
        if (mod(m+1, print_step)==0 .or. m==Nt-1) then
            print *, '--- time t = ',t_next,'---'
            call print_console_slice(x,u,N,5)
        end if
    end do
    close(10)

    !creating gru script
    call create_gnuplot_script
    print *, ' data saved to heat_solution.dat'
    print *, ' plot created: plot_heat.gnu'
    print *, ' to view run: gnuplot -persist plot_heat_gnu'



contains
    subroutine print_console_slice(x_arr,  u_arr, n, step_space)
        integer, intent(in) :: n, step_space
        real, intent(in) ::  x_arr(:), u_arr(:)
        integer :: k

        print '(A)', " Index |      X      |      U      "
        print '(A)', "-------|-------------|-------------"
        do k = 0, n, step_space
            print '(I5, " | ", F11.6, " | ", F11.6)', k, x_arr(k), u_arr(k)
        end do
        print *
    end subroutine print_console_slice
        
    subroutine create_gnuplot_script()
        integer :: unit_gnu
        unit_gnu = 20
        open(unit=unit_gnu, file='plot_heat.gnu', status='replace')
        
        write(unit_gnu, *) "set terminal wxt size 800,600 enhanced font 'Verdana,10'"
        write(unit_gnu, *) "set title 'Heat Equation"
        write(unit_gnu, *) "set xlabel 'x'"
        write(unit_gnu, *) "set ylabel 't'"
        write(unit_gnu, *) "set zlabel 'u'"
        write(unit_gnu, *) "set grid"
        write(unit_gnu, *) "set pm3d depthorder" ! quick print
        write(unit_gnu, *) "set style data pm3d"
        write(unit_gnu, *) "set view 60, 30"
        write(unit_gnu, *) ""
        write(unit_gnu, *) "splot 'heat_solution.dat' using 1:2:3 title 'Solution'"

        close(unit_gnu)
    end subroutine create_gnuplot_script

    
        

end program heat_conduction