program plot_3d_surface
    implicit none
    integer, parameter :: nx = 50, ny = 50
    real :: x(nx), y(ny), z(nx, ny)
    real :: xmin, xmax, ymin, ymax
    integer :: i, j
    real :: dx, dy
    character(len=100) :: data_file, script_file

    ! Настройка диапазона
    xmin = -5.0; xmax = 5.0
    ymin = -5.0; ymax = 5.0
    
    dx = (xmax - xmin) / real(nx - 1)
    dy = (ymax - ymin) / real(ny - 1)

    ! 1. Генерация данных (пример: волна)
    do i = 1, nx
        x(i) = xmin + (i-1)*dx
        do j = 1, ny
            y(j) = ymin + (j-1)*dy
            ! Функция Z = sin(sqrt(x^2 + y^2))
            z(i,j) = sin(sqrt(x(i)**2 + y(j)**2))
        end do
    end do

    ! Имена файлов
    data_file = 'data.dat'
    script_file = 'plot_script.gnu'

    ! 2. Запись данных в файл
    ! Формат: X Y Z (по строкам для каждой точки сетки)
    open(unit=10, file=data_file, status='replace', action='write')
    do i = 1, nx
        do j = 1, ny
            write(10, '(3F12.6)') x(i), y(j), z(i,j)
        end do
        ! Пустая строка разделяет линии сетки для правильного 3D отображения в gnuplot
        write(10, *) 
    end do
    close(10)

    ! 3. Создание скрипта для Gnuplot
    open(unit=11, file=script_file, status='replace', action='write')
    write(11, *) 'set terminal wxt size 800,600'  ! Окно вывода (wxt для Linux)
    write(11, *) 'set title "3D Surface from Fortran"'
    write(11, *) 'set xlabel "X axis"'
    write(11, *) 'set ylabel "Y axis"'
    write(11, *) 'set zlabel "Z axis"'
    write(11, *) 'set dgrid3d'                   ! Интерполяция, если сетка не идеальна
    write(11, *) 'set surface'
    write(11, *) 'set hidden3d'                  ! Убрать скрытые линии
    write(11, *) 'splot "', trim(data_file), '" using 1:2:3 with lines'
    close(11)

    ! 4. Вызов Gnuplot из Fortran
    ! Команда: gnuplot -persist plot_script.gnu
    ! -persist оставляет окно открытым после завершения скрипта
    call system('gnuplot -persist ' // trim(script_file))

    print *, "График построен! Данные сохранены в ", trim(data_file)

end program plot_3d_surface