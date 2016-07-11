        program helloworld
        
        use mpi
        
        integer :: rank
        integer :: numtasks
        integer :: ierr

        call MPI_Init(ierr) ! Initializes MPI calls
  
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)   ! obtains the rank of current MPI process 
        call MPI_Comm_size(MPI_COMM_WORLD, numtasks, ierr) !  obtains the total number of MPI processes */

        if (rank == 0 ) then 
        print *,"hello world. I am MASTER"
        else 
        print *,"hello world from process", rank," of ", numtasks
        end if
        call MPI_Finalize(ierr) ! Finalizes MPI calls 

        end program helloworld

