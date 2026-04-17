 # Gnuplot script for 2D potential heatmap
 # Assignment #18 - Poisson Equation, Variant 2
 
 set terminal pngcairo size 1200,1000 enhanced font "Verdana,11"
 set output "potential_heatmap.png"
 
 set title "Poisson Equation - Potential Distribution\\nAssignment #18, Variant 2" font ",16"
 set xlabel "x" font ",13"
 set ylabel "y" font ",13"
 
 set grid front lt 0 lw 0.5 lc rgb "#cccccc"
 set border lw 1.5
 
 # HEATMAP VIEW (top-down)
 set view map
 set pm3d map
 
 # Color palette (viridis-like)
 set palette defined (0 "#440154", 0.2 "#3b528b", 0.4 "#21918c", 0.6 "#5ec962", 0.8 "#fde725", 1 "#ffffff")
 
 # Colorbar
 set cbrange [*:*]
 set colorbox vertical user origin 0.92, 0.15 size 0.03, 0.7
 set cblabel "u(x,y)" font ",11" offset 2,0
 
 # Plot
 splot "potential_data.dat" using 1:2:3 with pm3d title "Potential"
 
 print "Potential heatmap saved to: potential_heatmap.png"
