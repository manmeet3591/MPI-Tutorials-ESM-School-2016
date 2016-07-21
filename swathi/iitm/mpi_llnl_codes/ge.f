C                        Example 26 gauss_elimination.f
c
c
c  Objective            : To solve the system of Linear Equations
c                         using Gaussian Elimination without pivoting 
c                         on 'p' processors.
c
c  Input                : Read files (mdatgaus.inp) for Matrix A
c        	      	     and (vdatgaus.inp) for Vector b
c
c  Output               : The solution of matrix system of linear 
c                         equations Ax=b on processor 0.
c
c  Description          : Input matrix is stored in n by n format.
c                         Columnwise cyclic distribution of input
c                         matrix is used for partitioning of the 
c                         matrix.
c
c  Necessary conditions : Number of Processes should be less than
c                         or equal to 8.Matrix size for Matrix A 
c                         and vector size for vector b should be
c                         equally striped, that is Matrix size and 
c                         Vector Size should be divisible by Number
c                         of processes used.
c
c******************************************************************
c
c
c
      program main

      include 'mpif.h'

      PARAMETER ( NMAX = 100)
      integer MyRank,Numprocs,NoofRows,nlocal_proc,myrow, ierr

      integer  recv_send,send_recv_tag
      integer  idest, NoofCols

      integer i,j,k,STATUS(MPI_STATUS_SIZE)
      real Matrix_A,Input_B,Cols_MatA,b1,y,Buffer_Pivot,
     $	   X_buffer,x1

      dimension Matrix_A(NMAX,NMAX),Input_B(NMAX),Cols_MatA(NMAX),
     $		Buffer_Pivot(NMAX), X_buffer(NMAX)

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, MyRank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, Numprocs, ierr )

c
c     Read matrices A,B at process #0
c
      if ( MyRank .eq. 0) then
	open(unit=12, file = './data/mdatgaus.inp')
         read(12,*) NoofRows, NoofCols

c	  print *, NoofRows, NoofCols
         write(6,*) 'Input Matrix' 
         do i = 1,NoofRows
            read(12,*) (Matrix_A(i,j),j=1,NoofCols)
            write(6,75) (Matrix_A(i,j),j=1,NoofCols)
         enddo
75    format(8(2x,f8.3))

	  open(unit=13, file = './data/vdatgaus.inp')
         write(6,*) 'Input Vector' 
         read(13,*) NoofRows
         read(13,*) (Input_B(i), i =1,NoofRows)
         write(6,75) (Input_B(i), i = 1, NoofRows)

	 close(12)
	 close(13)
      endif

c     Broadcast the number of rows of Matrix A
c
      call MPI_BCAST(NoofRows,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(mod(NoofRows,Numprocs) .ne. 0 ) then
         if(MyRank .eq. 0) 
     $	    print*,"Matrix cannot be evenly striped among processes"
	 goto 30
      endif

c     Distribute A,B to other nodes
c     Cyclic rowwise distribution
      nlocal_proc = NoofRows/Numprocs

      if ( MyRank .eq. 0 ) then
         do i = 1,NoofRows
            do j = 1,NoofRows
               Cols_MatA(j) = Matrix_A(i,j)
            enddo
            b1 = Input_B(i)
            idest = mod(i,Numprocs)
            recv_send = (i/Numprocs+1)*2
            if ( idest .ne. 0) then
c
c           Node 0 send corresponding row of 
c           matrix A and element of B to other nodes
c
c
               call MPI_SEND( Cols_MatA,NoofRows, MPI_REAL,idest,
     $		    recv_send, MPI_COMM_WORLD,STATUS,ierr)
               send_recv_tag = recv_send + 1
c
               call MPI_SEND( b1,1, MPI_REAL,idest,send_recv_tag,
     $              MPI_COMM_WORLD,STATUS,ierr)
            endif
         enddo
      else
c
         do recv_send = 2, nlocal_proc*2,2

c           Receive information from node 0

            call MPI_RECV( Cols_MatA,NoofRows,
     $           MPI_REAL,0,recv_send,MPI_COMM_WORLD,STATUS,ierr)
            send_recv_tag = recv_send + 1

            call MPI_RECV( b1,1,
     $           MPI_REAL,0,send_recv_tag,MPI_COMM_WORLD,STATUS,ierr)
            myrow = MyRank + (recv_send/2-1)*Numprocs

            do k = 1, NoofRows
               Matrix_A(myrow,k) = Cols_MatA(k) 
            enddo

            Input_B(myrow) = b1
         enddo
      endif
c
c
c         Gaussian_elimination


      do k = 1, NoofRows

c          Decide whose turn it is to do normalize A(k,k)
c
         if (MyRank .eq. mod(k,Numprocs)) then
            do j = k+1, NoofRows
               Matrix_A(k,j) = Matrix_A(k,j) / Matrix_A(k,k)
               Buffer_Pivot(j) = Matrix_A(k,j)
            enddo
            Input_B(k) = Input_B(k)/Matrix_A(k,k)
            y = Input_B(k)
            Matrix_A(k,k) = 1.0
         endif 

c        Broadcast the coefficients Matrix_A and result Input_B
c
         call MPI_BCAST(y,1,MPI_REAL,mod(k,Numprocs),
     $        MPI_COMM_WORLD,ierr)
         call MPI_BCAST(Buffer_Pivot,NoofRows,MPI_REAL,mod(k,Numprocs),
     $        MPI_COMM_WORLD,ierr)

         do i = k+1,NoofRows
            if (MyRank .eq. mod(i,Numprocs)) then

c              update the coefficient Matrix_A and result Input_B

               do j = k+1,NoofRows
                  Matrix_A(i,j) = Matrix_A(i,j) - 
     $				  Matrix_A(i,k)*Buffer_Pivot(j)
               enddo
               Input_B(i) = Input_B(i) - Matrix_A(i,k)*y
               Matrix_A(i,k) = 0.0
            endif
         enddo
      enddo

c
c     Triangle Solver
      do k = NoofRows, 1, -1

c         Decide whose turn to solve the X_buffer(k) value

         if (MyRank .eq. mod(k,Numprocs)) then
            X_buffer(k) = Input_B(k)/Matrix_A(k,k)
            x1 = X_buffer(k)
         endif

c        Broadcast this X_buffer(k) value
c        Update result Input_B and eliminate X_buffer(k) term

         if ( k .ne. 0)
     $        call MPI_BCAST(x1,1,MPI_REAL,mod(k,Numprocs),
     $        MPI_COMM_WORLD,ierr)
         do i = k, 1, -1
            if (MyRank .eq. mod(i,Numprocs)) then
               Input_B(i) = Input_B(i) - Matrix_A(i,k)*x1
            endif

c           Accumulate X_buffer(k) value in node 0

            if (MyRank .eq. 0) X_buffer(k) = x1
         enddo
      enddo

cc         Node 0 output the final result X_buffer.
c
      if (MyRank .eq. 0) then
         write(6,*)
         write(6,*) 'Results for solution of matrix system of equations'
         do i = 1,NoofRows
            write(6,97) i,X_buffer(i)
 97         format(1x,'x[',i1,'] = ',f5.2)
         enddo
      endif

 30   call MPI_FINALIZE(ierr)

      stop
      end

c	***************************************************
c

