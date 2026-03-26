program wave_equ
    implicit none
    ! precision
    integer,parameter :: dp=selected_real_kind(15,307)

    real(dp) :: c(1:4),n, p(2),v(2)
    real(dp):: pi=3.14159265359

    call random_number(c)
    call random_number(n)
    call random_number(p)
    call random_number(v)
    n=1+n*2
    c=1+c*9
    p=0.2+p*1.8
    v=1+v*9

contains
    function mu1(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res
        res=c(1)*sin(pi*n*x/l)
    end function mu1

    function mu2(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res

        if (alpha<=x .and. x<=beta) then
            res=c(2)
        else 
            res=0
        end if
    end function mu2

    function mu34(i,x) result(res)
        real(dp), intent(in) :: x
        integer, intent(in) :: i
        real(dp) :: res

        res=c(i)*(1-exp((-p(i)*t)/T0))
        
    end function mu34

    
end program wave_equ