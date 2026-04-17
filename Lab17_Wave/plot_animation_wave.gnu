 # Gnuplot script for 1D Wave Equation Animation
 set terminal pngcairo size 1200,600 enhanced font "Verdana,12"
 
 set title "1D Wave Equation - Time Evolution" font ",16"
 set xlabel "x" font ",13"
 set ylabel "u(x,t)" font ",13"
 set grid lt 0 lw 0.5 lc rgb "#cccccc"
 set border lw 1.5
 set yrange [*:*]
 
 # Animation loop
do for [frame=0:50] {
     set output sprintf("anim_frame_wave%03d.png", frame)
     set title sprintf("1D Wave Equation - Frame %d (t=%.2f)", frame, frame * 30 * 0.02) font ",16"
     filename = sprintf("frame_%03d.dat", frame)
     plot filename using 1:2 with lines lw 3 lc rgb "navy" title "u(x,t)"
 }
 print "Animation frames saved: anim_frame_wave*.png"
