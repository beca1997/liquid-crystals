program maier_saupe_2d
    implicit none
    
    real                ::  x, u_z              !variables x and uz
    real                ::  dx, du              !x and u_z stepsize
    integer             ::  N                   !number of data points
    real                ::  x_low, x_upp        !lower and upper limit for x 
    real                ::  a, b, h             !limits of integration and stepsize
    real                ::  I                   !I(X)
    real                ::  num, den            !will store value of integrated num and den

    integer             ::  j, k                !just step variables for loops
    real, external      ::  den_int, num_int    !functions to integrate
    
    a=0.
    b=1.
    h=0.001
    x_low=0.
    x_upp=88.
    dx=0.001
    N=int((x_upp-x_low)/dx)
    x=x_low
    open(unit=21, file="maier3d.csv", status="replace")
    do j=1,N
        call trapez(num_int, a, b, h, x, num)
        call trapez(den_int, a, b, h, x, den)
        x=x+dx
        I=num/den
        write(21,*) x,",",I
    enddo
    close(21)
end program maier_saupe_2d

function    num_int(u_z, x)
    implicit none
    real, intent(in)    ::  x, u_z
    real                ::  num_int
    num_int=2.*(u_z**2-1./2.)*exp(x*u_z**2)
    return
end function

function    den_int(u_z, x)
    implicit none
    real, intent(in)    ::  x, u_z
    real                ::  den_int
    den_int=exp(x*u_z**2)
    return
end function


subroutine trapez(f, a, b, h, T, integral)
    implicit none
    real    ::  f
    real    ::  a, b 
    real    ::  h 
    real, intent(out)    ::  integral
    real    ::  x
    real    ::  area
    real    ::  T
    x=a
    integral=0.
    
    do while (x.le.b)
        area=(f(x, T)+f(x+h,T))/2*((x+h)-x)
        integral=integral+area
        x=x+h
    enddo
    return
end subroutine 