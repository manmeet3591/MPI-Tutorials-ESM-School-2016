
c For mat_vec product
c
      program main

      include 'mpif.h'

      PARAMETER ( no_of_rows=8,no_of_cols=8)
      integer myrank,numprocs,no_of_rows, no_of_cols, ierr

      integer  ist, iend
      real :: dot_prod_myproc,dot_product
      real :: vec_prod_myproc,vec_product

      integer i,j,k,status(MPI_STATUS_SIZE)
      integer, dimension(2) :: domain_layout=(/2,2/)

      real x(no_of_cols),y(no_of_rows),yall(no_of_rows)
      real a(no_of_rows,no_of_cols)

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      if ( numprocs .ne. domain_layout(1)*domain_layout(2) ) then
        write(*,*) 'Domain layout not consistent with numprocs'
        call MPI_Finalize(ierr)
      endif

      if ( mod(no_of_cols,domain_layout(1)) .ne. 0 ) then
        write(*,*) 'No_of cols not consistent with domain_layout(1) ,',no_of_cols,domain_layout(1)
        call MPI_Finalize(ierr)
      endif

      if ( mod(no_of_rows,domain_layout(2)) .ne. 0 ) then
        write(*,*) 'No_of rows not consistent with domain_layout(2) ,',no_of_rows,domain_layout(2)
        call MPI_Finalize(ierr)
      endif

      rows_per_proc = no_of_rows/domain_layout(2)
      cols_per_proc = no_of_cols/domain_layout(1)
      col_block_index = mod(myrank,domain_layout(1)) + 1
      row_block_index = (myrank - col_block_index +1)/domain_layout(1) + 1
c
      ist = (row_block_index - 1)*rows_per_proc + 1
      iend = ist +  rows_per_proc - 1

      jst = (col_block_index - 1)*cols_per_proc + 1
      jend = jst +  cols_per_proc - 1
C
      
C Put some data
      x(1) = 1.0 ; x(2) = 2.0 ; x(3) = 3.0 ; x(4) = 4.0; x(5) = 5.0; x(6) = 6.0 ; x(7) = 7.0; x(8) = 8.0
      do j = 1, no_of_cols
        do i = 1, no_of_rows
          a(i,j) = (i-1)*no_of_cols + j
        enddo
      enddo
C

      call mat_vec(ist,iend,jst,jend,no_of_rows,no_of_cols,a(ist:iend,jst:jend),x(jst:jend),y(ist:iend))
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call MPI_REDUCE(y,yall,no_of_rows,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if (myrank .eq. 0) then 
       write(*,*)'yall', yall
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
          y(i) = y(i) + a(i,j)*v(j)
        enddo
      enddo
      return
      end
      
       
