! Serial dot product program 

program serial_dot

integer,parameter :: N=2000

  real :: x(N)
  real :: y(N)
  real :: prod
  integer :: i

  do i=1,N-1
    x(i) = 0.01 * i
    y(i) = 0.03 * i
  end do

  call dotProduct(x,y,N,prod)
  print *,"dotProduct = ", prod

end program serial_dot

subroutine dotProduct(x, y, n, prod) 
  integer :: n
  real :: x(n), y(n)
  real :: prod
  integer :: i

  prod = 0.0
  do i=1,n-1
    prod = prod + x(i)*y(i)
  end do
  
end subroutine dotProduct
