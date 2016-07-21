********************
c dot product
c xTx
c
      program main

      include 'mpif.h'

      PARAMETER ( no_of_rows=8,no_of_cols=8)
      integer myrank,numprocs,no_of_rows, ierr, elem_per_proc

      integer  idest, ist, iend
      real :: dot_prod_myproc,dot_product
      real :: vec_prod_myproc,vec_product

      integer i,j,k,status(MPI_STATUS_SIZE)

      real x(no_of_cols),y(no_of_cols)
      real a(no_of_rows,no_of_cols)

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      x(1) = 1.0 ; x(2) = 2.0 ; x(3) = 3.0 ; x(4) = 4.0; x(5) = 5.0; x(6) = 6.0 ; x(7) = 7.0; x(8) = 8.0
      y = 2.0*x
      a = 1.0
      elem_per_proc = no_of_rows/numprocs
      if ( mod(no_of_rows,numprocs) .ne. 0 ) then
        write(*,*) 'No_of_rows is not a multiple of Numprocs,',no_of_rows, numprocs
        call MPI_Finalize(ierr)
      endif
c
      ist = myrank*elem_per_proc + 1
      iend = ist + elem_per_proc -1
      call dot_prod(ist,iend,no_of_rows,x(ist:iend),dot_prod_myproc)
      call vec_prod(ist,iend,no_of_rows,x(ist:iend),y(ist:iend),vec_prod_myproc)
      call MPI_REDUCE(dot_prod_myproc,dot_product,1,MPI_REAL,MPI_SUM,0,      
     & MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(vec_prod_myproc,vec_product,1,MPI_REAL,MPI_SUM,0,      
     & MPI_COMM_WORLD,ierr)
      if (myrank .eq. 0) then 
       write(*,*)'Dot product of', x,dot_product
       write(*,*)'Vec product of', x,y,vec_product
      endif
      call MPI_Finalize(ierr)
      stop
      end
      subroutine dot_prod(ist,iend,n,x,dprod)
      integer, intent(in) :: ist,iend
      real, dimension(ist:iend), intent(in)  :: x
      integer :: i
      real, intent(out) :: dprod
      dprod = 0.0
      do i = ist,iend
        dprod = dprod + x(i)*x(i)
      enddo
      return
      end

      subroutine vec_prod(ist,iend,n,x,y,vprod)
      real, dimension(ist:iend), intent(in)  :: x,y
      integer, intent(in) :: ist,iend,n
      integer :: i
      real, intent(out) :: vprod
      vprod = 0.0
      do i = ist,iend
        vprod = vprod + x(i)*y(i)
      enddo
      return
      end
  
      subroutine mat_vec(ist,iend,jst,jend,ni,nj,a,v,y)
      integer :: ist,iend,jst,jend,i,j,k
      real, dimension(ist:iend,jst:jend),intent(in) :: a
      real, dimension(jst:jend),intent(in) :: v
      real, dimension(ist:iend),intent(out) :: y
      y = 0.0
      do j = jst,jend
        do i = ist, iend
          y(i) = y(i) + a(i,j)*y(j)
        enddo
      enddo
      return
      end
      
       
