
# wave_animation.plt - Анимация решения волнового уравнения

set terminal gif animate delay 10 size 800,600
set output 'wave_animation.gif'

set title 'Одномерное волновое уравнение - Вариант 2'
set xlabel 'x'
set ylabel 'u(x,t)'
set grid
set xrange [0:1]
set yrange [-2:2]

# Count number of data files
n_files = system("ls wave_data_*.dat 2>/dev/null | wc -l")

if (n_files == 0) {
    print "No data files found! Run the Fortran program first."
    exit 1
}

# Animate through all data files
do for [i=0:n_files-1] {
    filename = sprintf('wave_data_%03d.dat', i)
    t_val = i * 20 * 0.5 * (1.0/50.0)
    plot filename using 1:2 with lines lw 3 lc rgb 'blue' title sprintf('t = %.3f', t_val)
}

set output
print 'Animation saved as wave_animation.gif'
