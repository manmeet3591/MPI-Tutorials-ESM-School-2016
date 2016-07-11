! Serial Trapezoid Rule Program                 
! [f(x0)/2 + f(xn)/2 + f(x1) + ... + f(xn-1)]*h 

program serial_trapezoid

  real :: integral ! definite integral result
  real,parameter :: a=0.0 ! left end point
  real,parameter :: b=1.0 ! right end point
  integer,parameter :: N=100000 ! number of subdivisions
  real :: h          !  base width of subdivision
  real :: x
  integer :: i
  real :: fa,fb,fx

  h = (b-a)/N
  call f(a,fa)
  call f(b,fb)      
  integral = (fa+fb)/2.0
  x = a

  do i=1,N-1
      x = x+h
      call f(x,fx)
      integral = integral + fx
  
  end do
  integral = integral*h

  print *,"WITH N=", N, " TRAPEZOIDS, INTEGRAL=", integral

end program serial_trapezoid

subroutine f(x,y)

real :: x,y
y = exp(x*x)

end subroutine f

