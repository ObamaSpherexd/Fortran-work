 # Gnuplot script for 2D potential with contours
 # Assignment #18 - Poisson Equation, Variant 2
 
 set terminal pngcairo size 1200,1000 enhanced font "Verdana,11"
 set output "/home/artem/Study/fortran_projects/my_first_fortran/Lab18_Poisson/potential_contour.png"
 
 set title "Poisson Equation - Potential Contours\\nAssignment #18, Variant 2" font ",16"
 set xlabel "x" font ",13"
 set ylabel "y" font ",13"
 
 set grid front lt 0 lw 0.5 lc rgb "#cccccc"
 set border lw 1.5
 
 # 3D surface with contours
 set view 60, 30
 set contour base
 set cntrparam bspline
 set cntrparam points 5
 set isosamples 50,50
 
 # Color palette
 set palette defined (0 "#440154", 0.2 "#3b528b", 0.4 "#21918c", 0.6 "#5ec962", 0.8 "#fde725", 1 "#ffffff")
 
 unset key
 
 # Plot
splot "/home/artem/Study/fortran_projects/my_first_fortran/Lab18_Poisson/potential_data.dat" using 1:2:3 with pm3d title "Potential"
 
 print "Potential contour saved to: potential_contour.png"
