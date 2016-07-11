        program trapezoid
        
        use mpi
        implicit none
        integer,parameter :: ROOT = 0
        real :: integral ! definite integral result
        real :: a = 0.0 ! left end point
        real :: b = 1.0 ! right end point
        integer :: N ! number of subdivisions
        real :: h    ! base width of subdivision
        real :: x,fx
        integer :: i
        integer :: my_rank
        integer :: numprocs
        integer :: ierr
        integer :: status(MPI_STATUS_SIZE) 
 
! we will need some local variables 
        real :: local_a, local_b, flocal_a, flocal_b
        real :: local_N
        real :: lcl_integral
  
! MPI programming begins 
        call MPI_Init(ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
        
        if (my_rank .eq. ROOT) then

                print *,"Enter values of "
                print *," a = "
                read(*,*) a
                print *," b = "
                read(*,*) b
                print *," N = "
                read(*,*) N

        end if

        call MPI_Bcast(a, 1, MPI_REAL, ROOT, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(b, 1, MPI_REAL, ROOT, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(N, 1, MPI_INT,    ROOT, MPI_COMM_WORLD, ierr)

        h = (b-a)/N   ! we assume we use the same integration step on all processes 
  
! Find out what the local values are on each process */  
        local_N = N / numprocs
        local_a = a + my_rank * local_N * h
        local_b = local_a + local_N * h
  
! begins local integration 
        x = local_a
        call f(local_a,flocal_a)
        call f(local_b,flocal_b)
        lcl_integral = (flocal_a+flocal_b)/2.0

        do i=1,local_N-1

        x = local_a + i*h
        call f(x,fx)
        lcl_integral = lcl_integral + fx

        end do
        lcl_integral = lcl_integral*h
 
! Reduce and send result to ROOT 
        call MPI_Reduce(lcl_integral, integral, 1, MPI_REAL, MPI_SUM, ROOT, MPI_COMM_WORLD, ierr)
!       To send the updated value of integral to every processor
!       call MPI_Allreduce(lcl_integral, integral, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr) 

        print *,"Process ",my_rank," WITH N=", N, " TRAPEZOIDS, INTEGRAL=", integral

        call MPI_Finalize(ierr)

        end program trapezoid

subroutine f(x,y)

real :: x,y
y = exp(x*x)

end subroutine f
