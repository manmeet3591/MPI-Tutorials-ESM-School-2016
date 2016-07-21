! Conjugate gradient; f90 version.
! Notation is as in the Kincaid & Cheney pg 211
!
! Parallelisation only for matrix vectr product. Vector-vector product and 
! matrix add, subtract etc are done in serial fully in each proc (Redundant computation).
!
      program conjugate_gradient
      use mpi

      PARAMETER ( n=1000)
      integer myrank,numprocs,n, ierr

      integer  ist, iend

      integer i,j,status(MPI_STATUS_SIZE)
      integer, dimension(2) :: domain_layout=(/2,2/)
      integer :: m, k  = 0

      real x(n),b(n),r(n),v(n),z(n),xe(n)
      real ax(n),av(n)
      real a(n,n)
      real :: c,t,d, delta = 100, vtv = 1.0e8
      double precision :: start_time, end_time, elapsed_time

!      data (a(1,i),i=1,n) /420.,210.,140.,105./
!      data (a(2,i),i=1,n) /210.,140.,105.,84./
!      data (a(3,i),i=1,n) /140.,105.,84.,70./
!      data (a(4,i),i=1,n) /105.,84.,70.,60./
!      data (b(i),i=1,n) /875.,539.,399.,319./
      data m/15/
!      data (x(i),i=1,n) /0.,0.,0.,0./

! Put in some data
      xe = 1.0
      x = 0.0
      do i = 1,n
          b(i) = 0.0
        do j = 1,n
          a(i,j) =  i*j + i**2 + j**2
          b(i) = b(i) + a(i,j)*xe(j)
        enddo
      enddo   

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      if ( numprocs .ne. domain_layout(1)*domain_layout(2) ) then
        write(*,*) 'Domain layout not consistent with numprocs'
        call MPI_Finalize(ierr)
      endif

      if ( mod(n,domain_layout(1)) .ne. 0 ) then
        write(*,*) 'Mat size  consistent with domain_layout(1) ,',n,domain_layout(1)
        call MPI_Finalize(ierr)
      endif

      if ( mod(n,domain_layout(2)) .ne. 0 ) then
        write(*,*) 'Mat size not consistent with domain_layout(2) ,',n,domain_layout(2)
        call MPI_Finalize(ierr)
      endif
      start_time = MPI_WTIME()


!
! Initial residual      
      call mat_vec_mpp(n,n,a,x,ax,domain_layout,myrank)
      r = b - ax
      v  = r
      call vec_prod(1,n,n,r,r,c)
      do while( k .le. m .and. vtv .ge. delta)
        k = k + 1 
        call vec_prod(1,n,n,v,v,vtv)   
        call mat_vec_mpp(n,n,a,v,z,domain_layout,myrank)
        call vec_prod(1,n,n,v,z,vtz)
        t = c/vtz
        x = x + t*v
        r = r - t*z
        call vec_prod(1,n,n,r,r,d)
        v = r + d*v/c
        c = d
        z = 0.0
        if (myrank .eq. 0) write(*,*)'k,vtv,x',k,vtv,x
      enddo
      end_time = MPI_WTIME()
      elapsed_time = end_time - start_time
      write(*,*)'st,end,elap',start_time,end_time,elapsed_time


      call MPI_Finalize(ierr)
      stop
      end
!
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
  
      subroutine mat_vec_proc(ist,iend,jst,jend,nrow,ncol,a,v,y)
! This is with in a processor
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
   
      subroutine mat_vec_mpp(nrow,ncol,a,x,y,domain_layout,myrank)
      use mpi
      real  a(nrow,ncol),x(ncol),y(nrow)
      integer :: nrow,ncol,i,j,k,ierr
      integer,  dimension(2) :: domain_layout
      integer :: rows_per_proc,cols_per_proc,col_block_index,row_block_index
      integer :: ist,iend,jst,jend
!
      rows_per_proc = nrow/domain_layout(1)
      cols_per_proc = ncol/domain_layout(2)
      col_block_index = mod(myrank,domain_layout(2)) + 1
      row_block_index = (myrank - col_block_index +1)/domain_layout(2) + 1
!
      ist = (row_block_index - 1)*rows_per_proc + 1
      iend = ist +  rows_per_proc - 1

      jst = (col_block_index - 1)*cols_per_proc + 1
      jend = jst +  cols_per_proc - 1
!
      y = 0.0
      call mat_vec_proc(ist,iend,jst,jend,nrow,ncol,a(ist:iend,jst:jend),x(jst:jend),y(ist:iend))
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(y,y,nrow,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
!      call MPI_BCAST(y,nrow,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      return
      end


     

     
      
       
