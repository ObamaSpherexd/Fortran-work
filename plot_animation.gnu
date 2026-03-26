 # Gnuplot script for animated heatmap
 set terminal pngcairo size 1200,1000 enhanced font "Verdana,11"
 
 set title "Heat Equation in Rectangle - Time Evolution\\nAssignment #16" font ",16"
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
 
 # Animation loop
do for [frame=0:100] {
    filename = sprintf("frame_%3.3d.dat", frame)
    set output sprintf("anim_frame_%3.3d.png", frame)
     t_label = (frame == 0) ? "t = 0.00" : sprintf("t = %.2f", frame *    5.0000000000000003E-002 )
     set title sprintf("Heat Equation in Rectangle - Time Evolution\\nt = %.2f s", frame *    5.0000000000000003E-002 ) font ",16"
     splot filename using 1:2:3 with pm3d notitle
 }
 
 print "Animation frames saved: anim_frame_*.png"
