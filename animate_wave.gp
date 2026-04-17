# Скрипт для анимации волны (метод прогонки)
set terminal gif animate delay 5 size 800,600 optimize
set output "wave_implicit.gif"
set xrange [0:  1.00]
set yrange [-0.3:0.3]
set xlabel "x"
set ylabel "u(x,t)"
set title "Волновое уравнение (неявная схема, прогонка)"
set grid
set key outside
do for [i=0:60] {
   t = i*.0500
   plot "wave_data.txt" index i using 1:2 with lines lw 2 lc rgb "blue" title sprintf("t = %.3f", t)
}
set output
