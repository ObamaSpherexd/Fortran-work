program poisson
   implicit none

   ! precision
   integer, parameter:: dp=selected_real_kind(15,307)

   ! measurements
   real(dp),parameter :: lx=2.0_dp
   real(dp),parameter :: ly=2.0_dp

   ! variant params
   real(dp),parameter :: c=10.0_dp
   real(dp),parameter :: p=2.0_dp
   real(dp),parameter :: q=1.9_dp
   real(dp),parameter:: r=1.2_dp
   real(dp),parameter:: s=1.0_dp

   real(dp)::omega,d_inv,max_res,f_val,residual

   ! grid params
   integer, parameter ::N=50
   integer,parameter::M=50

   ! limits
   integer,parameter::kmax=600
   real(dp),parameter::eps=1.0e-6_dp

   real(dp),allocatable:: u(:,:),x(:),y(:)

   integer::i,j,k

   real(dp):: hx,hy

   real(dp):: pi=3.14159265358979_dp

   real(dp):: hx_inv,hy_inv,d

   character(len=*), parameter :: output_dir = '/home/artem/Study/fortran_projects/my_first_fortran/Lab18_Poisson/'

   hx=lx/real(N,dp)
   hy=ly/real(M,dp)

   hx_inv=1.0_dp/(hx*hx)
   hy_inv=1.0_dp/(hy*hy)

   d=2.0_dp*(hx_inv+hy_inv)
   d_inv=1.0_dp/d

   omega=2.0_dp/(1.0_dp+sin(pi/real(N,dp)))

   print *, '=== Poisson Equation in Rectangle ==='
   print *, 'Assignment #18, Variant 2'
   print *, 'Domain: Lx =', lx, ' Ly =', ly
   print *, 'Grid: N =', N, ' M =', M, ' hx =', hx, ' hy =', hy
   print *, 'RHS: f(x,y) = c*(lx-x)^p * x^q * y^r * (ly-y)^s'
   print *, '  c =', c, ' p =', p, ' q =', q, ' r =', r, ' s =', s
   print *, 'BC: psi = 1 (everywhere)'
   print *, 'SOR method: omega =', omega
   print *, 'Tolerance: eps =', eps
   print *, 'Max iterations: kmax =', kmax
   print *, ''

   allocate(u(0:N,0:M))
   allocate(x(0:N))
   allocate(y(0:M))

   do j=0,N
      x(j)=hx*real(j,dp)
   end do
   do i=0,M
      y(i)=hy*real(i,dp)
   end do

   !bounds
   do i=0,M
      u(0,i)=1.0_dp
      u(N,i)=1.0_dp
   end do
   do j = 0, N
      u(j,0)=1.0_dp
      u(j,M)=1.0_dp
   end do

   ! init guess
   do j = 1, N-1
      do i=1,M-1
         u(j,i)=1.0_dp
      end do
   end do

   ! ================================================================
   ! SOR iterations (Successive Over-Relaxation)
   ! ================================================================
   print '(A)', 'Starting SOR iterations...'
   print '(A10, A20)', 'Iter', 'Max Residual'
   print '(A10, A20)', '----------', '--------------------'

   do k=1,kmax
      max_res=0.0_dp
      do j = 1, N-1
         do i = 1, M-1
            f_val=f(x(j),y(i))

            f_val=d_inv*( &
               hx_inv*(u(j+1,i)+u(j-1,i)) + &
               hy_inv*(u(j,i+1)+u(j,i-1)) + &
               f_val)
            residual=abs(f_val-u(j,i))
            u(j,i)=(1.0_dp-omega)*u(j,i)+omega*f_val

            if (residual>max_res) then
               max_res=residual
            end if
         end do
      end do

      ! print progress every 30 iterations
      if (mod(k, 30) == 0) then
         print '(I10, E20.6)', k, max_res
      end if

      ! check convergence
      if (max_res < eps) then
         print *,''
         print *, 'Converged after',k,'iterations'
         print *, 'Max residual:',max_res
         exit
      end if

      if (k == kmax) then
         print '(A)', 'Not converged, max iterations reached'
         print *, 'Max residual:',max_res
      end if
   end do


   ! ================================================================
   ! Output results
   ! ================================================================
   print *, ''
   print *, '=== Solution Summary ==='
   print '(A15, F12.6)', 'u(0,0)     = ', u(0, 0)
   print '(A15, F12.6)', 'u(Lx/2,Ly/2)= ', u(N / 2, M / 2)
   print '(A15, F12.6)', 'u(Lx,Ly)   = ', u(N, M)


   ! ================================================================
   ! Save solution for gnuplot
   ! ================================================================
   open(unit = 10, file = trim(output_dir)//'potential_data.dat', status = 'replace')
   write(10, *) '# x y u(x,y) - potential distribution'
   write(10, *) '# Assignment #18 - Poisson Equation, Variant 2'
   write(10, *) '# f(x,y) = c*(lx-x)^p * x^q * y^r * (ly-y)^s'
   write(10, *) '# BC: psi = 1 everywhere'
   write(10, *) '# Iterations:', k

   ! write data with blank lines between y-scans for gnuplot pm3d
   do i = 0, M
      do j = 0, N
         write(10, '(3F15.6)') x(j), y(i), u(j, i)
      end do
      write(10, *) ''  ! separator for gnuplot
   end do
   close(10)

   ! ================================================================
   ! Generate gnuplot script for potential distribution
   ! ================================================================
   call create_gnuplot_script(output_dir)

   ! also create contour plot script
   call create_contour_script(output_dir)

   print *, ''
   print *, '=== Results saved ==='
   print *, 'Potential data: potential_data.dat'
   print *, 'Gnuplot heatmap: plot_potential.gnu'
   print *, 'Gnuplot contour:  plot_contour.gnu'
   print *, 'To view heatmap:  gnuplot -persist plot_potential.gnu'
   print *, 'To view contour:  gnuplot -persist plot_contour.gnu'

   ! deallocate
   deallocate(u, x, y)


contains
   function f(x,y) result(res)
      real(dp), intent(in) :: x,y
      real(dp) :: res
      res=c*(lx-x)**p*x**q*y**r*(ly-y)**s
      !res=c*(sin(pi*(x/lx)**q))**p*(sin(pi*(x/ly)**s))**r
   end function f


   ! ================================================================
   ! Generate gnuplot script for heatmap (pm3d map view)
   ! ================================================================
   subroutine create_gnuplot_script(output_dir)
      character(len=*), intent(in) :: output_dir
      integer :: unit_gnu
      character(len=:), allocatable :: png_path, dat_path

      png_path = trim(output_dir)//'potential_heatmap.png'
      dat_path = trim(output_dir)//'potential_data.dat'

      unit_gnu = 20
      open(unit = unit_gnu, file = trim(output_dir)//'plot_potential.gnu', &
         status = 'replace')

      write(unit_gnu, *) '# Gnuplot script for 2D potential heatmap'
      write(unit_gnu, *) '# Assignment #18 - Poisson Equation, Variant 2'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'set terminal pngcairo size 1200,1000 enhanced '//&
         'font "Verdana,11"'
      write(unit_gnu, *) 'set output "'//trim(png_path)//'"'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'set title "Poisson Equation - Potential Distribution'//&
         '\\nAssignment #18, Variant 2" font ",16"'
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
      write(unit_gnu, *) 'set palette defined (0 "#440154", 0.2 "#3b528b", '//&
         '0.4 "#21918c", 0.6 "#5ec962", 0.8 "#fde725", 1 "#ffffff")'
      write(unit_gnu, *) ''
      write(unit_gnu, *) '# Colorbar'
      write(unit_gnu, *) 'set cbrange [*:*]'
      write(unit_gnu, *) 'set colorbox vertical user origin 0.92, 0.15 '//&
         'size 0.03, 0.7'
      write(unit_gnu, *) 'set cblabel "u(x,y)" font ",11" offset 2,0'
      write(unit_gnu, *) ''
      write(unit_gnu, *) '# Plot'
      write(unit_gnu, '(A)') 'splot "'//trim(dat_path)//'" using 1:2:3 '//&
         'with pm3d title "Potential"'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'print "Potential heatmap saved to: '//&
         'potential_heatmap.png"'

      close(unit_gnu)
   end subroutine create_gnuplot_script

   ! ================================================================
   ! Generate gnuplot script for contour plot
   ! ================================================================
   subroutine create_contour_script(output_dir)
      character(len=*), intent(in) :: output_dir
      integer :: unit_gnu
      character(len=:), allocatable :: png_path, dat_path

      png_path = trim(output_dir)//'potential_contour.png'
      dat_path = trim(output_dir)//'potential_data.dat'

      unit_gnu = 21
      open(unit = unit_gnu, file = trim(output_dir)//'plot_contour.gnu', &
         status = 'replace')

      write(unit_gnu, *) '# Gnuplot script for 2D potential with contours'
      write(unit_gnu, *) '# Assignment #18 - Poisson Equation, Variant 2'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'set terminal pngcairo size 1200,1000 enhanced '//&
         'font "Verdana,11"'
      write(unit_gnu, *) 'set output "'//trim(png_path)//'"'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'set title "Poisson Equation - Potential Contours'//&
         '\\nAssignment #18, Variant 2" font ",16"'
      write(unit_gnu, *) 'set xlabel "x" font ",13"'
      write(unit_gnu, *) 'set ylabel "y" font ",13"'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'set grid front lt 0 lw 0.5 lc rgb "#cccccc"'
      write(unit_gnu, *) 'set border lw 1.5'
      write(unit_gnu, *) ''
      write(unit_gnu, *) '# 3D surface with contours'
      write(unit_gnu, *) 'set view 60, 30'
      write(unit_gnu, *) 'set contour base'
      write(unit_gnu, *) 'set cntrparam bspline'
      write(unit_gnu, *) 'set cntrparam points 5'
      write(unit_gnu, *) 'set isosamples 50,50'
      write(unit_gnu, *) ''
      write(unit_gnu, *) '# Color palette'
      write(unit_gnu, *) 'set palette defined (0 "#440154", 0.2 "#3b528b", '//&
         '0.4 "#21918c", 0.6 "#5ec962", 0.8 "#fde725", 1 "#ffffff")'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'unset key'
      write(unit_gnu, *) ''
      write(unit_gnu, *) '# Plot'
      write(unit_gnu, '(A)') 'splot "'//trim(dat_path)//'" using 1:2:3 '//&
         'with pm3d title "Potential"'
      write(unit_gnu, *) ''
      write(unit_gnu, *) 'print "Potential contour saved to: '//&
         'potential_contour.png"'

      close(unit_gnu)
   end subroutine create_contour_script


end program poisson
