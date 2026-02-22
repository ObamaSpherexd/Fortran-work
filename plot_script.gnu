 set terminal wxt size 800,600 enhanced font 'Verdana,10'
 set title 'Convection Solver: u(x,t)'
 set xlabel 'x'
 set ylabel 't'
 set zlabel 'u'
 #set grid
 set pm3d depthorder
 set style data pm3d
 
 # Plot as 3D surface without interpolation
 splot 'solution.dat' using 1:2:3 title 'Solution'
