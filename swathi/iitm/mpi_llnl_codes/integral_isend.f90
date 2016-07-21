
      Program Example1_3
      use mpi
!#######################################################################
!#
!# This is an MPI example on parallel integration
!# It demonstrates the use of :
!#
!# * MPI_Init, MPI_Comm_rank, MPI_Comm_size, MPI_Finalize
!# * MPI_Recv, MPI_Isend, MPI_Wait
!# * MPI_ANY_SOURCE, MPI_ANY_TAG
!#
!# Dr. Kadin Tseng
!# Scientific Computing and Visualization
!# Boston University
!# 1998
!#
!#######################################################################
      implicit none
      integer n, p, i, j, proc, ierr, master
      real h, a, b, integral, pi
      integer req(1)

      integer myid, source, tag, status(MPI_STATUS_SIZE)
      real my_int, integral_sum
      real*8 st_time,end_time,elap_time

      data master/0/    ! processor 0 collects integral sums from other processors

!**Starts MPI processes ...

      call MPI_Init(ierr)                              ! starts MPI
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)   ! get current proc id
      call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)      ! get number of procs

      st_time = MPI_WTime()
      pi = acos(-1.0)   !  = 3.14159...
      a = 0.0           ! lower limit of integration
      b = pi/2.         ! upper limit of integration
      n = 5000000           ! number of increment within each process
      tag = 123         ! set the tag to identify this particular job
      h = (b-a)/n/p     ! length of increment

      my_int = integral(a,myid,h,n)
      write(*,*)'myid=',myid,',  my_int=',my_int

      if(myid .eq. master) then               ! the following is serial
        integral_sum = my_int
        do proc=1,p-1
          call MPI_Recv(my_int, 1, MPI_REAL, & 
            MPI_ANY_SOURCE, MPI_ANY_TAG, &  ! more efficient, less prone to deadlock
            MPI_COMM_WORLD, status, ierr) ! root receives my_int from proc
          integral_sum = integral_sum + my_int
        enddo
      else
        call MPI_Isend(my_int, 1, MPI_REAL, master, tag, &
          MPI_COMM_WORLD, req(1), ierr)     ! send my_int to master 
        call other_work('more work')  ! because of Isend, gets here immediately
        call MPI_Wait(req(1), status, ierr)  ! wait for nonblock send ...
      endif

!**integral_sums from all procs have been collected and summed ...

      if(myid .eq. master) then
        write(*,*)'The Integral =',integral_sum
      endif
      end_time = MPI_WTime()
      elap_time = end_time - st_time
      if (myid .eq. master) write(*,*)'elap_time=',elap_time
      call MPI_Finalize(ierr)                       ! let MPI finish up ...

      stop
      end

      subroutine other_work(header)
      implicit none
      integer i,j,nmax,l,m
      real k
      character*(*) header
      write(*,'(a)') header
      nmax = 100000
! Do some junk work
!      do l = 1,nmax
      k = 0.0
      do i = 1, nmax
        do j = 1,nmax
          k = k+ i +j
        enddo
      enddo
!      enddo
!
      write(*,*)'junk work done,nmax,k',nmax,k
      return
      end

      real function integral(a,i,h,n)
      implicit none
      integer n, i, j
      real h, h2, aij, a, fct, x

      integral = 0.0                ! initialize integral
      h2 = h/2.
      do j=0,n-1                    ! sum over all "j" integrals
        aij = a + (i*n +j)*h        ! lower limit of "j" integral
        integral = integral + fct(aij+h2)*h
      enddo

      return
      end

      real function fct(x)
      implicit none
      real x
      fct = cos(x)
      end


