!  Parallel dot product
        program dotproduct

        use mpi
        implicit none
        integer,parameter ::  N=2000
        integer :: i
        real ::  prod
        integer ::  my_rank
        integer :: num_procs
        integer :: ierr
        integer :: local_N
        real,allocatable,dimension(:) :: local_x, local_y 
        real :: local_prod
!--------------Definitions end-------------------------
        call MPI_Init(ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  
        local_N = N / num_procs ! assuming N is totally divisible by num_procs 
        allocate(local_x(1:local_N))
        allocate(local_y(1:local_N)) 
        do i=1,local_N 
                local_x(i) = 0.01 * (i + my_rank * local_N)
                local_y(i) = 0.03 * (i + my_rank * local_N)
        end do
        call dot(local_x,local_y,local_N,local_prod)
        call MPI_Reduce(local_prod, prod, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if (my_rank == 0) then 
        print *,"dotProduct = ", prod
        end if
        print *,"from my_rank = ", my_rank, " local_prod = ",local_prod

        call MPI_Finalize(ierr)

        end program dotproduct

subroutine dot(x, y, n, prod) 
  real,dimension(n) :: x,y
  integer :: i
  real :: prod 
        
  prod = 0.0        
  do i=1,n 
    prod = prod + x(i)*y(i)
  end do

end subroutine dot
