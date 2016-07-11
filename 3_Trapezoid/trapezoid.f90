! MPI Trapezoid Rule Program                 */
! [f(x0)/2 + f(xn)/2 + f(x1) + ... + f(xn-1)]*h */
        
        program trapezoid

        use mpi
        implicit none

        real :: integral          ! definite integral result
        real, parameter ::  a=0.0 ! left end point
        real, parameter ::  b=1.0 ! right end point
        real, parameter ::  N=100000 ! number of subdivisions
        real :: h                 ! base width of subdivision
        real :: x,y,local_y1,local_y2, y1,y2
        integer :: i
        integer :: my_rank
        integer ::  numprocs
  
!  we will need some local variables 
        real :: local_a, local_b
        integer :: local_N
        real :: lcl_integral
        integer :: ierr
        integer :: dest = 0
        real ::  recv ! a variable to receive results 
        integer :: status(MPI_STATUS_SIZE)  
        h = (b-a)/N               ! we assume we use the same integration step on all processes 

!  MPI programming begins 
        call MPI_Init(ierr)


        call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

!  Find out what the local values are on each process   
        local_N = N/numprocs
        local_a = a + my_rank * local_N * h
        local_b = local_a + local_N * h

!  begins local integration */
        x = local_a
        call f(local_a,local_y1)
        call f(local_b,local_y2)
        lcl_integral = (local_y1 +local_y2)/2.0

        do i=1,local_N-1
           x = x+h
           call f(x,y)
           lcl_integral = lcl_integral + y
        end do

        lcl_integral = lcl_integral*h

!  send the local results to Process 0 

        if ( my_rank .ne. dest ) then

        call MPI_Send(lcl_integral, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, ierr)
        print *,"local integral in processor number ",my_rank," =",lcl_integral  
!  Process 0 receives and sums up the results 
        else 
        print *,"local integral in processor number ",my_rank," =",lcl_integral  
        integral = lcl_integral
        do i=1,numprocs-1 
                call MPI_Recv(recv, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, status, ierr)
                integral = integral + recv
        end do

              print *, "WITH N=", N, " TRAPEZOIDS, INTEGRAL=", integral
        
        end if

!  MPI programming ends 
        call MPI_Finalize(ierr)

        end program trapezoid        

        subroutine f(x,y)
        
        real :: x,y
        y = exp(x*x)  

        end subroutine f
