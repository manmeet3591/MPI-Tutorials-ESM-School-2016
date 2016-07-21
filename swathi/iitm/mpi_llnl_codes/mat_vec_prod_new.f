
c For mat_vec product
c
      program main

      include 'mpif.h'

      PARAMETER ( ni=8,nj=8)
      integer myrank,numprocs,ni, nj, ierr

      integer  ist, iend
      real :: dot_prod_myproc,dot_product
      real :: vec_prod_myproc,vec_product

      integer i,j,k,status(MPI_STATUS_SIZE)
      integer, dimension(2) :: domain_layout=(/2,2/)

      real x(nj),y(ni),yall(ni)
      real a(ni,nj)

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      if ( numprocs .ne. domain_layout(1)*domain_layout(2) ) then
        write(*,*) 'Domain layout not consistent with numprocs'
        call MPI_Finalize(ierr)
      endif

      if ( mod(nj,domain_layout(1)) .ne. 0 ) then
        write(*,*) 'No_of cols not consistent with domain_layout(1) ,',nj,domain_layout(1)
        call MPI_Finalize(ierr)
      endif

      if ( mod(ni,domain_layout(2)) .ne. 0 ) then
        write(*,*) 'No_of rows not consistent with domain_layout(2) ,',ni,domain_layout(2)
        call MPI_Finalize(ierr)
      endif


      
C Put some data
      x(1) = 1.0 ; x(2) = 2.0 ; x(3) = 3.0 ; x(4) = 4.0; x(5) = 5.0; x(6) = 6.0 ; x(7) = 7.0; x(8) = 8.0
      do j = 1, nj
        do i = 1, ni
          a(i,j) = (i-1)*nj + j
        enddo
      enddo
C
      
      call mat_vec_mpp(ni,nj,a,x,y,domain_layout,myrank)

      if (myrank .eq. 0) then 
       write(*,*)'y', y
      endif
      call MPI_Finalize(ierr)
      stop
      end
C
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
  
      subroutine mat_vec_proc(ist,iend,jst,jend,ni,nj,a,v,y)
C This is with in a processor
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
   
      subroutine mat_vec_mpp(ni,nj,a,x,y,domain_layout,myrank)
      real  a(ni,nj),x(nj),y(ni)
      integer :: ni,nj,i,j,k,ierr
      integer,  dimension(2) :: domain_layout
      integer :: rows_per_proc,cols_per_proc,col_block_index,row_block_index
      integer :: ist,iend,jst,jend
C
      rows_per_proc = ni/domain_layout(2)
      cols_per_proc = nj/domain_layout(1)
      col_block_index = mod(myrank,domain_layout(1)) + 1
      row_block_index = (myrank - col_block_index +1)/domain_layout(1) + 1
c
      ist = (row_block_index - 1)*rows_per_proc + 1
      iend = ist +  rows_per_proc - 1

      jst = (col_block_index - 1)*cols_per_proc + 1
      jend = jst +  cols_per_proc - 1
C
      call mat_vec_proc(ist,iend,jst,jend,ni,nj,a(ist:iend,jst:jend),x(jst:jend),y(ist:iend))
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call MPI_REDUCE(y,y,ni,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      return
      end


     

     
      
       
