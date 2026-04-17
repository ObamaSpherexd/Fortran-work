set yrange [-0.3:0.3]
set grid
do for [i=0:60] {
   plot "wave_data.txt" index i using 1:2 with lines lw 2 title "Wave"
   pause 0.05
}
