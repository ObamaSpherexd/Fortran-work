 # Heatmap script for heat conduction solution
 set terminal wxt size 1200,800 enhanced font 'Verdana,10'
 set title 'Heat Conduction: u(x,t)\\nVariant 2' font ',14'
 set xlabel 'x'
 set ylabel 't'
 set zlabel 'Temperature u'
 set grid
 set hidden3d
 set dgrid3d 50,50 splines
 
 splot 'solution.dat' using 1:2:3 with lines title 'Solution'
