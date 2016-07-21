! Swathi..Test program to solve u<sub>t = delsq(u) on a unit square
!
! with zero bc's and sin(pi*x)*sin(pi*y) initial condition
PROGRAM  heat_eqn
        use mpp_mod
        use mpp_domains_mod
        use mpp_io_mod


        implicit none

#ifdef use_netCDF
#include <netcdf.inc>
#endif


	integer :: ni, nj, halo = 1,n, unit, ntsteps
	real, dimension(:),allocatable :: xt0,yt0
	real, allocatable,dimension(:,:) :: u,unew
	real :: dtts,delx,dely, cfl
	integer :: ioun,npes,pe, iscomp,iecomp,jscomp,jecomp
	integer, dimension(2) :: domain_layout=(/1,1/)
	type(domain2D),  allocatable, dimension(:), target :: local_domains  ! domains with halos for all processors
	type(domain2D),  pointer                           :: local_domain   ! domain with halos on this processor
	type(domain2D) :: domain
	!
	type(domain1D) :: xdom,ydom
        type(axistype) :: x,y, t 
        type(fieldtype) :: f
	integer :: i,j,io_status
	namelist /heat_eqn_nml/ ni,nj,dtts,ntsteps
!
	call mpp_init()
	call mpp_io_init()
	npes = mpp_npes()
	pe = mpp_pe()


! Read some namelist inputs
	call mpp_open(ioun,'heat_eqn_nml',action=MPP_RDONLY,form=MPP_ASCII)
	read(ioun,heat_eqn_nml,iostat=io_status)
	write(stdout(),heat_eqn_nml)
! Make some checks
	if ( nj <= 0 .OR. nj <= 0 .OR. dtts <= 0 .OR. ntsteps <= 0) then
          write(stdout(),*)'nj,nj,dtts,ntsteps',nj,nj,dtts,ntsteps
	  call mpp_error(FATAL,'=>Error one of the above is <= 0')
        endif

! Next find the ni and nj's that will
        delx = 1.0/(ni+1)
	dely = 1.0/(nj+1)
! Do a cfl check
	cfl = 0.5*delx**2*dely**2/(delx**2+dely**2)
        write(stdout(),*)'cfl=',cfl
	if ( dtts > cfl)  call mpp_error(FATAL,'CFL violation')
	allocate (xt0(ni))
	allocate (yt0(nj))
	do i = 1,ni
	  xt0(i) = i*delx
	enddo
	  write(stdout(),*)'xt0',xt0
	do j = 1,nj
	  yt0(j) = j*dely
        enddo
	  write(stdout(),*)'yt0',yt0
!
	allocate (local_domains(0:npes-1))
	local_domain => local_domains(pe)
        call mpp_define_layout ((/1,ni,1,nj/),npes,domain_layout)
	if (pe ==0) then
       	  write(stdout(),*)'domain layout',domain_layout
	  write(stdout(),*)'ni,nj',ni,nj
        endif
	call mpp_define_domains( (/1,ni,1,nj/), domain_layout, local_domains(pe), &
		xhalo=halo, yhalo=halo )
	call mpp_get_compute_domain (local_domains(pe), iscomp, iecomp, jscomp, jecomp)
	write(stdout(),*)'pe,iscomp,iecomp,jscomp,jecomp',pe,iscomp,iecomp,jscomp,jecomp
	call mpp_get_domain_components( local_domains(pe), xdom, ydom )
	!
! Setup for netcdf distributed write
        if( pe.EQ.0 )print *, 'netCDF distributed write'
  	call mpp_open( unit, 'heat_gen_output', action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  	call mpp_write_meta( unit, x, 'X', 'km', 'X distance', domain=xdom, data=xt0 )
  	call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', domain=ydom, data=yt0 )
  	call mpp_write_meta( unit, t, 'T', 'sec', 'Time' )
  	call mpp_write_meta( unit, f, (/x,y,t/), 'TEMP', 'degrees', 'heat gen' )
  	call mpp_write( unit, x )
  	call mpp_write( unit, y )

	allocate (u(iscomp-halo:iecomp+halo,jscomp-halo:jecomp+halo))
	allocate (unew(iscomp-halo:iecomp+halo,jscomp-halo:jecomp+halo))
! Initial conditions
	do j = jscomp,jecomp
	  do i = iscomp, iecomp
	    u(i,j) = sin(3.1414*xt0(i))*sin(3.1414*yt0(j))
          enddo
        enddo
        call mpp_update_domains(u,local_domain)
!	write(stdout(),*)'init u',u
! The boundary conditions
        if (jscomp == 1) u(:,jscomp-1) = 0.0
	if (jecomp == nj) u(:,jecomp+1) = 0.0
	if (iscomp == 1) u(iscomp-1,:) = 0.0
	if (iecomp == ni) u(iecomp+1,:) = 0.0
!
	do n = 1,ntsteps 
	  do j = jscomp,jecomp
	    do i = iscomp,iecomp
	      unew(i,j) = u(i,j)*(1-2*dtts*(delx**2 + dely**2)/(delx**2*dely**2)) + & 
	                  dtts*((u(i-1,j) + u(i+1,j))/delx**2 + (u(i,j-1) + u(i,j+1))/dely**2)
            enddo
          enddo
          u(iscomp:iecomp,jscomp:jecomp) = unew(iscomp:iecomp,jscomp:jecomp)
	  call mpp_update_domains(u,local_domain)
          call mpp_write( unit, f , local_domain, u,  n*1. )
        enddo
        call mpp_close(unit)
!	write(stdout(),*)'iscomp,iecomp,n,u',iscomp,iecomp,n,u
	call mpp_exit()
END PROGRAM heat_eqn
