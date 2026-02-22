program heat_conduction
   implicit none

   ! precision
   integer, parameter :: dp = selected_real_kind(15, 307)

   ! task params
   real(dp), parameter :: l = 1.0_dp      ! stick len
   real(dp), parameter :: T_total = 1.0_dp ! model time
   real(dp), parameter :: a_coef = 1.0_dp ! heat coef

   ! variant params
   real(dp), parameter :: c = 1.0_dp      ! init coef
   real(dp), parameter :: p1 = 1.0_dp     ! left border init
   real(dp), parameter :: p2 = 1.0_dp     ! right border init
   real(dp), parameter :: b_param = 1.0_dp ! damp factor
   real(dp), parameter :: f0 = 1.0_dp     ! source amp
   real(dp), parameter :: beta1 = 0.2_dp  ! start source area
   real(dp), parameter :: beta2 = 0.8_dp  ! end source area

   ! grid params
   integer, parameter :: N = 50
   integer :: Nt, m, j
   integer :: print_step, save_step
   real(dp) :: h, tau, t_current
   real(dp) :: sigma, mu, nu, theta, ratio
   real(dp) :: pi_val

   ! arrays (dynamic)
   real(dp), allocatable :: x(:)
   real(dp), allocatable :: u(:)
   real(dp), allocatable :: u_new(:)
   real(dp), allocatable :: g(:)
   real(dp), allocatable :: xi(:)
   real(dp), allocatable :: eta(:)


   ! init params
   print_step = 10
   save_step = 5
   h = l / real(N, dp)

   ! Tau eval
   ratio = 0.5_dp
   tau = ratio * h**2 / a_coef**2
   Nt = int(T_total / tau)
   if (Nt < 1) Nt = 1
   tau = T_total / real(Nt, dp) ! fit tau

   ! Coefs
   mu = 2.0_dp * h**2 / (tau * a_coef**2)
   sigma = 2.0_dp + mu
   nu = 2.0_dp * h**2 / a_coef**2
   theta = 2.0_dp - mu

   pi_val = acos(-1.0_dp)


   allocate(x(0:N))
   allocate(u(1:N+1), u_new(1:N+1))
   allocate(g(1:N+1))
   allocate(xi(1:N+1))
   allocate(eta(1:N+1))

   ! Grid
   do j = 0, N
      x(j) = real(j, dp) * h
   end do

   ! init condition
   do j = 1, N+1
      u(j) = u0(x(j-1))
   end do

   ! open file for gnu
   open(unit=10, file='heat_solution.dat', status='replace')
   write(10, *) "# x t u"

   !  t=0
   do j = 0, N
      write(10, '(3F15.6)') x(j), 0.0_dp, u(j+1)
   end do
   write(10, *) "" ! Block separator for splot

   ! output init
   print *, "--- Time t = 0.0 ---"
   call print_console_slice(x, u, N, 5)

   ! time steps
   call cpu_time(t_current)
   t_current = 0.0_dp

   do m = 1, Nt
      t_current = real(m, dp) * tau

      ! border
      u_new(1) = phi1(t_current)
      u_new(N+1) = phi2(t_current)

      ! right part e
      do j = 2, N
         g(j) = -nu * f_source(x(j-1), t_current + tau/2.0_dp) &
            - (u(j+1) - theta * u(j) + u(j-1))
      end do

      ! forward

      xi(2) = 0.0_dp
      eta(2) = u_new(1)

      do j = 3, N
         xi(j) = 1.0_dp / (sigma - xi(j-1))
         eta(j) = (eta(j-1) - g(j-1)) * xi(j)
      end do

      ! backward
      do j = N, 2, -1
         u_new(j) = xi(j+1) * u_new(j+1) + eta(j+1)
      end do

      u = u_new

      ! write in file
      if (mod(m, save_step) == 0 .or. m == Nt) then
         do j = 0, N
            write(10, '(3F15.6)') x(j), t_current, u(j+1)
         end do
         write(10, *) "" ! separator
      end if

      ! console log
      if (mod(m, print_step) == 0 .or. m == Nt) then
         print *, "--- Time t =", t_current, "---"
         call print_console_slice(x, u, N, 5)
      end if
   end do

   close(10)

   ! gnu script
   call create_gnuplot_script()

   print *, ' Data saved to heat_solution.dat'
   print *, ' Gnuplot script created: plot_heat.gnu'
   print *, ' To view 3D plot, run in terminal:'
   print *, '   gnuplot -persist plot_heat.gnu'

contains

   subroutine print_console_slice(x_arr, u_arr, n, step_space)
      integer, intent(in) :: n, step_space
      real(dp), intent(in) :: x_arr(:), u_arr(:)
      integer :: k

      print '(A)', " Index |      X      |      U      "
      print '(A)', "-------|-------------|-------------"
      do k = 0, n, step_space
         print '(I5, " | ", F11.6, " | ", F11.6)', k, x_arr(k), u_arr(k+1)
      end do
      print *
   end subroutine print_console_slice


   subroutine create_gnuplot_script()
      integer :: unit_gnu
      unit_gnu = 20
      open(unit=unit_gnu, file='plot_heat.gnu', status='replace')


      write(unit_gnu, *) "set terminal wxt size 800,600 enhanced font 'Verdana,10'"
      write(unit_gnu, *) "set title 'Heat Conduction Solver: u(x,t)'"
      write(unit_gnu, *) "set xlabel 'x'"
      write(unit_gnu, *) "set ylabel 't'"
      write(unit_gnu, *) "set zlabel 'u'"
      write(unit_gnu, *) "#set grid"
      write(unit_gnu, *) "set pm3d depthorder" ! Better performance for surfaces
      write(unit_gnu, *) "set style data pm3d"
      write(unit_gnu, *) ""
      write(unit_gnu, *) "# Plot as 3D surface without interpolation"
      write(unit_gnu, *) "splot 'heat_solution.dat' using 1:2:3 title 'Solution'"

      close(unit_gnu)
   end subroutine create_gnuplot_script

   ! phys funcs
   function u0(x) result(res)
      real(dp), intent(in) :: x
      real(dp) :: res
      !res = c * (l - x) / l
      !res=x/l
      res=abs(cos(3.14_dp*2.0_dp*(x/(l**3))))**5
   end function u0

   function phi1(t) result(res)
      real(dp), intent(in) :: t
      real(dp) :: res
      !res = p1 * (1.0_dp - exp(-b_param * t))
      !res=1
      res=p1*(1.0_dp-exp(-b_param*t))
   end function phi1

   function phi2(t) result(res)
      real(dp), intent(in) :: t
      real(dp) :: res
      !res = p2
      !res=0
      res=p2*exp(-b_param*t)
   end function phi2

   function f_source(x, t) result(res)
      real(dp), intent(in) :: x, t
      real(dp) :: res
      !if (beta1 <= x .and. x <= beta2) then
      !   res = f0
      !else
      !   res = 0.0_dp
      !end if
      if (beta1<=x .and. x<=beta2) then
         res=f0*(1.0_dp-exp(-x*t))
      else
         res=0.0_dp
      end if
      !res=0
   end function f_source

end program heat_conduction
