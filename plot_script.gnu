 set terminal wxt size 800,600 enhanced font 'Verdana,10'
 set title 'Convection Solver: u(x,t)'
 set xlabel 'x'
 set ylabel 't'
 set zlabel 'u'
 set grid
 set hidden3d
 set dgrid3d 50,50 splines
 set style data lines
 
 # Plot as 3D surface
  splot 'solution.dat' using 1:2:3 with lines title 'Solution'
