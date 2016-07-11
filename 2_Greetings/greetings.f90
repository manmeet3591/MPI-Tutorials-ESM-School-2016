        program greetings
  
        use mpi
      
        integer,parameter :: ROOT = 0;
        integer :: my_rank
        integer :: recv_rank
        integer :: numtasks
        integer :: p
        integer :: ierr
        integer status(MPI_STATUS_SIZE)

        call MPI_Init(ierr)     ! Initializes MPI calls
        call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr) !  obtains the rank of current MPI process 
        call MPI_Comm_size(MPI_COMM_WORLD, numtasks, ierr) ! obtains the total number of MPI processes 

        if (my_rank .ne. ROOT ) then
                call MPI_Send(my_rank, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, ierr)
        else 
                do p = 1,numtasks-1
                        call MPI_Recv(recv_rank, 1, MPI_INT, p, 0, MPI_COMM_WORLD, status, ierr)
                        print *,"Greetings from Process ", recv_rank
                end do
        end if
        call MPI_Finalize(ierr) !        Finalizes MPI calls 
  
        end program greetings

