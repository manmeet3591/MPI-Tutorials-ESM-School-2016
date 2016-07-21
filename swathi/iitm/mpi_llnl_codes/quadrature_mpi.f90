! A collection of subroutines to do 1-D numerical integration
  module num_integ_routines
  implicit none
  contains
  subroutine trap(func,a,b,c,n)
  real, external :: func
  real :: a,b,c,delx
  integer :: n,i
  delx = (b-a)/n
  c = func(a) + func(b)
  do i = 1, n-1
    c = c + 2.0* func(a + i*delx)
  enddo
  c = c*delx/2
  return
  end subroutine trap
!
  subroutine simpson(func,a,b,c,n)
  real, external :: func
  real :: a,b,c,delx
  integer :: n,i
  delx = (b-a)/n
  c = func(a) + func(b)
  do i = 1,n/2
    c = c + 4.0*func(a + (2*i-1)*delx)
  enddo
  do i = 2,n/2
    c = c + 2.0*func(a + (2*i -2)*delx)
  enddo
  c = c*delx/3
  return
  end subroutine simpson
!
  subroutine gauss_quad(func,a,b,c,n)
! For now a and b have to be -1 and 1; we will do change of base later.
  real, external :: func 
  real :: a,b,c
  integer :: n,i
  real, dimension(n) :: pt,wt
  if ( a .eq. -1.0 .and. b .eq. 1.0) then
  
    call  gauss_pt_wt(n,pt,wt) 
    c = 0.0
    do i = 1,n
      c = c + wt(i)*func(pt(i))
    enddo
  else 
    write(*,*) 'a and b must be -1 and 1 respectively'
    stop
  endif
  return
  end subroutine gauss_quad

  subroutine gauss_pt_wt(n,pt,wt)
  real, dimension(n) :: pt,wt
  integer :: n
  select case (n)
   case (1)
    pt(1) = 0.0 ; wt(1) = 2.0
   case (2)
    pt(1) = -1.0/sqrt(3.0) ; pt(2) = 1.0/sqrt(3.0)
    wt(2) = 1.0 ; wt(2) = 1.0
   case (3) 
    pt(1) =  -sqrt(0.6) ; pt(2) = 0.0 ; pt(3) = sqrt(0.6)
    wt(1) = 5.0/9.0 ; wt(2) = 8.0/9.0 ; wt(3) = 5.0/9.0
  end select
  end subroutine gauss_pt_wt 
   

  function func(x)
  real :: func,x
  func = x*x
  end function func
  end module num_integ_routines
!
  program num_integ_mpi
  use mpi
  use num_integ_routines
  real :: a=-1.0,b=1.0,c
  real :: a_proc,b_proc,c_proc
  integer :: n = 100,ng=3
  integer :: myrank,numprocs,ierr

  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
! Divide up the region a-b into numprocs subintervals
  delta = (b-a)/numprocs
  a_proc = a + myrank*delta
  b_proc = a_proc + delta


  call trap(func,a_proc,b_proc,c,n)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(c,c,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)  
  if (myrank .eq. 0) write(*,*)'c (trap)=',c

!  call simpson(func,a,b,c,n)
!  write(*,*)'c (simp)=',c
!
!  call gauss_quad(func,a,b,c,ng)
!  write(*,*)'c (gauss)=',c
  call MPI_Finalize(ierr) 
  stop
  end   
