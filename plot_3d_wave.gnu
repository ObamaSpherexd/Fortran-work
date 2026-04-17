 # Gnuplot script for 3D Wave Equation Surface
 set terminal pngcairo size 1200,800 enhanced font "Verdana,11"
 set output "wave_3d_surface.png"
 
 set title "1D Wave Equation: u(x,t) Surface" font ",14"
 set xlabel "x" font ",12"
 set ylabel "t" font ",12"
 set zlabel "u(x,t)" font ",12"
 
 # 3D surface with pm3d coloring
 set pm3d depthorder
 set pm3d implicit
 set palette defined (-1 "blue", 0 "white", 1 "red")
 set view 60, 30, 1.0, 1.0
 set hidden3d
 set contour base
 set cntrparam levels 20
 
 set key off
 
 splot "wave_solution.dat" using 1:2:3 with pm3d title "Wave Solution"
 
 print "3D plot saved: wave_3d_surface.png"
