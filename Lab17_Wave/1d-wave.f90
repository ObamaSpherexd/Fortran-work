program wave_equ
   implicit none
   ! precision
   integer, parameter :: dp = selected_real_kind(15, 307)

   ! Parameters for initial and boundary conditions
   real(dp) :: c(4), n_val, p(2), v(2)
   real(dp), parameter :: pi = 3.141592653589793_dp
   integer, parameter :: N = 50, M = 10
   real(dp), parameter :: l = 1.0_dp, a = 1.0_dp

   ! Bounds for mu2(x)
   real(dp), parameter :: alpha = 0.3_dp, beta = 0.7_dp

   ! Grid and time variables
   !real(dp) :: h, tau, tau0, T0, T, courant, t
   real(dp)::h,t0,tau0,t,courant,tau
   integer :: k, kmax, j, frame_count

   real(dp), allocatable :: u_prev(:), u_now(:), u_next(:), x(:)


   call random_number(c)
   call random_number(n_val)
   call random_number(p)
   call random_number(v)
   n_val=1.0_dp+n_val*2.0_dp
   c=1.0_dp+c*9.0_dp
   p=0.2_dp+p*1.8_dp
   v=1.0_dp+v*9.0_dp

   ! grid params

   h=l/real(N,dp)
   T0=l/a
   tau0=h/a

   tau=0.5_dp*tau0
   T=5.0_dp*T0

   courant=a*tau/h

   kmax=nint(T/tau)

   allocate(x(0:N), u_prev(0:N),u_now(0:N), u_next(0:N))

   !grid
   do j=0,N
      x(j)=j*h
   end do

   ! init conditions at t=0, u(x,0)=mu1(x)
   u_prev(0)=mu34(3,0.0_dp)
   u_prev(N)=mu34(4,0.0_dp)


   do j=1,N-1
      u_prev(j)=mu1(x(j))
   end do

   u_now(0)=mu34(3,tau)
   u_now(N)=mu34(4,tau)
   do j=1,N-1
      u_now(j)=mu1(x(j))+tau*mu2(x(j))+ &
         0.5_dp*tau**2 *(a**2*d2mu1_dx2(x(j)))
   end do


   ! --- SAVE INITIAL FRAME (t=0) ---
   frame_count = 0
   call save_frame(frame_count, 0.0_dp, x, u_prev, N)

   print *, '=== 1D Wave Equation Animation ==='
   print *, 'L =', l, 'a =', a
   print *, 'N =', N, 'Steps =', kmax
   print *, 'Saving frames to disk...'

   ! main time cycle
   ! Open 3D data file
   open(unit=40, file='wave_solution.dat', status='replace')
   write(40, *) "# x t u"

   ! t=0 data
   do j = 0, N
      write(40, '(3F15.6)') x(j), 0.0_dp, u_prev(j)
   end do
   write(40, *) "" ! Block separator for splot
   !$omp parallel do default(none)
   do k=1,kmax
      t=k*tau
      do j=1,N-1
         u_next(j)=2.0_dp*(1.0_dp-courant**2)*u_now(j)-u_prev(j)+ &
            courant**2*(u_now(j+1)+u_now(j-1))
      end do
      ! border on current step
      u_next(0)=mu34(3,t)
      u_next(N)=mu34(4,t)


      ! Save frame every M steps
      if (mod(k, M) == 0) then
         frame_count = frame_count + 1
         call save_frame(frame_count, t, x, u_next, N)

         ! Save 3D data
         do j = 0, N
            write(40, '(3F15.6)') x(j), t, u_next(j)
         end do
         write(40, *) "" ! separator
      end if

      u_prev=u_now
      u_now=u_next
   end do
   !$omp end parallel do

   close(40)

   ! --- FINALIZE ANIMATION ---
   call create_animation_script(frame_count)

   ! --- CREATE AND RENDER 3D PLOT ---
   call create_and_render_3d_plot()

   print *, 'Calculation complete.'
   print *, 'Frames saved: ', frame_count + 1
   print *, 'To view animation, run: gnuplot plot_animation_wave.gnu'
   print *, 'To create GIF (requires ImageMagick): bash animate_wave.sh'
   print *, '3D surface plot saved: wave_3d_surface.png'

   deallocate(x, u_prev, u_now, u_next)


contains
   function d2mu1_dx2(x) result(res)
      real(dp), intent(in) :: x
      real(dp) :: res
      !res=-(pi**2*c(1)*n_val**2*sin(pi*n*x/l))/l**2
      res=0
   end function d2mu1_dx2

   function mu1(x) result(res)
      real(dp), intent(in) :: x
      real(dp) :: res
      !res=c(1)*sin(pi*n*x/l)
      res=0
   end function mu1

   function mu2(x) result(res)
      real(dp), intent(in) :: x
      real(dp) :: res

      !if (alpha<=x .and. x<=beta) then
      !    res=c(2)
      !else
      !   res=0.0_dp
      !end if
      res=(pi*a/l)*sin(pi*x/l)
   end function mu2

   function mu34(i,t) result(res)
      real(dp), intent(in) :: t
      integer, intent(in) :: i
      real(dp) :: res
      res=0
      !res=c(i)*(1-exp((-p(i)*t)/T0))

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

   ! === ИСПРАВЛЕННЫЙ GNUPLOT SCRIPT ===
   subroutine create_animation_script(n_frames)
      integer, intent(in) :: n_frames
      integer :: unit_gnu

      unit_gnu = 21
      open(unit=unit_gnu, file='plot_animation_wave.gnu', status='replace')

      write(unit_gnu, *) '# Gnuplot script for 1D Wave Equation Animation'
      write(unit_gnu, *) 'set terminal pngcairo size 1200,600 enhanced font "Verdana,12"'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'set title "1D Wave Equation - Time Evolution" font ",16"'
      write(unit_gnu, *) 'set xlabel "x" font ",13"'
      write(unit_gnu, *) 'set ylabel "u(x,t)" font ",13"'
      write(unit_gnu, *) 'set grid lt 0 lw 0.5 lc rgb "#cccccc"'
      write(unit_gnu, *) 'set border lw 1.5'
      ! Автоматический масштаб по Y (убирает обрезание)
      write(unit_gnu, *) 'set yrange [-2:2]'
      write(unit_gnu, *) ''

      write(unit_gnu, *) '# Animation loop'
      write(unit_gnu, '(A,I0,A)') 'do for [frame=0:', n_frames, '] {'
      write(unit_gnu, *) '    set output sprintf("anim_frame_wave%03d.png", frame)'
      write(unit_gnu, *) '    set title sprintf("1D Wave Equation - Frame %d (t=%.2f)", frame, frame * 30 * 0.02) font ",16"'
      write(unit_gnu, *) '    filename = sprintf("frame_%03d.dat", frame)'
      write(unit_gnu, *) '    plot filename using 1:2 with lines lw 3 lc rgb "navy" title "u(x,t)"'
      write(unit_gnu, *) '}'
      write(unit_gnu, *) 'print "Animation frames saved: anim_frame_wave*.png"'
      close(unit_gnu)

      ! === ИСПРАВЛЕННЫЙ BASH SCRIPT ===
      open(unit=22, file='animate_wave.sh', status='replace')
      write(22, *) '#!/bin/bash'
      write(22, *) 'echo "Creating animated GIF..."'
      ! -delay 10 = 0.1 сек на кадр. Меняйте это число для регулировки скорости.
      write(22, '(A,I0,A)') 'convert -delay 15 -loop 0 anim_frame_wave000.png anim_frame_wave*.png wave_animation.gif'
      write(22, *) 'echo "Animation saved: wave_animation.gif"'
      close(22)
   end subroutine create_animation_script

   subroutine create_and_render_3d_plot()
      integer :: unit_gnu
      character(len=300) :: cmd

      ! Create gnuplot script for 3D surface
      unit_gnu = 23
      open(unit=unit_gnu, file='plot_3d_wave.gnu', status='replace')

      write(unit_gnu, *) '# Gnuplot script for 3D Wave Equation Surface'
      write(unit_gnu, *) 'set terminal pngcairo size 1200,800 enhanced font "Verdana,11"'
      write(unit_gnu, *) 'set output "wave_3d_surface.png"'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'set title "1D Wave Equation: u(x,t) Surface" font ",14"'
      write(unit_gnu, *) 'set xlabel "x" font ",12"'
      write(unit_gnu, *) 'set ylabel "t" font ",12"'
      write(unit_gnu, *) 'set zlabel "u(x,t)" font ",12"'
      write(unit_gnu, *) ''
      write(unit_gnu, *) '# 3D surface with pm3d coloring'
      write(unit_gnu, *) 'set pm3d depthorder'
      write(unit_gnu, *) 'set pm3d implicit'
      write(unit_gnu, *) 'set palette defined (-1 "blue", 0 "white", 1 "red")'
      write(unit_gnu, *) 'set view 60, 30, 1.0, 1.0'
      write(unit_gnu, *) 'set hidden3d'
      write(unit_gnu, *) 'set contour base'
      write(unit_gnu, *) 'set cntrparam levels 20'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'set key off'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'splot "wave_solution.dat" using 1:2:3 with pm3d title "Wave Solution"'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'print "3D plot saved: wave_3d_surface.png"'

      close(unit_gnu)

      ! Run gnuplot to render the 3D plot
      print *, 'Rendering 3D surface plot...'

      ! Try with LD_PRELOAD to avoid snap library conflicts
      cmd = 'LD_PRELOAD=/lib/x86_64-linux-gnu/libpthread.so.0 ' // &
         'gnuplot plot_3d_wave.gnu 2>/dev/null || ' // &
         'gnuplot plot_3d_wave.gnu 2>/dev/null'
      call system(trim(cmd))

   end subroutine create_and_render_3d_plot

end program wave_equ
