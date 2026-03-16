 # Gnuplot script for 2D heat equation
 set terminal pngcairo size 1200,800 enhanced font "Verdana,10"
 set output "heat_2d_surface.png"
 
 set title "Heat Equation in Rectangle\\nAssignment #16" font ",14"
 set xlabel "x" font ",12"
 set ylabel "y" font ",12"
 set zlabel "u(x,y,t)" font ",12"
 set grid
 set hidden3d
 
 # Color palette
 set palette defined (0 "#440154", 0.2 "#3b528b", 0.4 "#21918c", 0.6 "#5ec962", 0.8 "#fde725", 1 "#ffffff")
 set pm3d
 
 # 3D surface plot
 splot "heat_2d_data.dat" using 1:2:4 with pm3d title "u(x,y,t)"
 
 # Alternative: heatmap view (top-down)
 # set view map
 # set pm3d map
 # replot
