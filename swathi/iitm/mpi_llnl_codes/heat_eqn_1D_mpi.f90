!
      module heat_types
      type domain1D
        integer :: isc, iec, xhalo
        integer :: isd,ied
        integer :: nb_left, nb_right
        integer :: pe,npes
        integer :: ni
        integer :: ni_local
        real :: delx
        real,dimension(:),allocatable :: x_local
      end type domain1D
      end module heat_types

! For solving the 1D transient heat conduction problem.
!
      program heat_eqn_1D
      use mpi
      use heat_types
      include 'netcdf.inc'
      type(domain1D), allocatable, dimension(:),target  :: domains
      type(domain1D), pointer :: local_domain
      real, allocatable,  dimension(:,:) :: u
      real, allocatable,  dimension(:) :: usave,uexact


      integer :: this_pe,npes_tot, ierr
      integer :: ni_local,nb_left,nb_right
      integer :: xhalo = 1
      logical :: periodic = .false. 
      integer ::  isc, iec,n
      integer :: tau=0,taup1=1,nt_steps = 100
      integer :: ni = 20000 , io_int = 100000
      real :: del_t = 0.0000000001, del_t_cfl, time=0.0
      double precision :: st_time_comp,end_time_comp,st_time_comm,end_time_comm
      double precision :: st_time_io,end_time_io,st_time_wait,end_time_wait
      double precision :: st_time_all,end_time_all,elap_time_all
      double precision :: st_time_init,end_time_init,st_time_bc,end_time_bc
      double precision :: elap_time_init=0.0,elap_time_bc=0.0
      double precision :: tot_time_comp=0.0, tot_time_comm=0.0, tot_time_io = 0.0
      double precision :: tot_time_wait=0.0,tot_time_all, tot_time_allproc

        

      integer i,j,k,status(MPI_STATUS_SIZE)
      integer :: domain_layout=4

      real :: delx, pi = 3.14145
!
      integer ncid,ier,timeid,tempid,temp_exactid,ltime
      character (len=12) :: fname
      namelist /heat_eqn_nml/ ni,del_t,nt_steps,domain_layout,io_int
      

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, this_pe, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, npes_tot, ierr )
      st_time_all = MPI_WTime()
      open(10,file='heat_eqn_nml')
      read(10,nml=heat_eqn_nml)
      write(*,*)'ni,del_t,nt_steps,domain_layout',ni,del_t,nt_steps,domain_layout

      if ( mod(ni,domain_layout) .ne. 0 ) then
        write(*,*) 'ni not consistent with domain_layout ,',ni,domain_layout
        call MPI_Finalize(ierr)
        stop
      endif

! CFL Check
      delx = 1.0/(ni+1)
      del_t_cfl = 0.5*delx**2
      if ( del_t > del_t_cfl)  then 
        write(*,*) 'CFL violation, delt,del_t_cfl',del_t,del_t_cfl
         call MPI_Finalize(ierr)
         stop
      endif

!
      allocate (domains(0:npes_tot-1))
      local_domain => domains(this_pe)
!
      st_time_init = MPI_WTime()
      call define_domains(ni,delx,domain_layout,this_pe,domains,npes_tot,xhalo,periodic)
      call get_domain_components(local_domain,isc,iec,isd,ied, &
                 ni_local,ni,this_pe,nb_left,nb_right)

! 
      allocate(u(isd:ied,0:1))
      u(:,:) = 0.0
      allocate(usave(isc:iec),uexact(isc:iec))
      
!     Boundary conditions
! Not sure how this works for periodic condidtions and halo > 1.
      call boundary_cond(u,isd,ied,nb_left,nb_right,ni)
      
!     Initial conditions
      call init_cond(u,isc,iec,isd,ied,local_domain)
      end_time_init = MPI_WTime()
      elap_time_init = elap_time_init + end_time_init - st_time_init
! Allocate buffers for sendrecv

!
      st_time_wait = MPI_WTIME()     
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      end_time_wait = MPI_WTIME()
      tot_time_wait = tot_time_wait + end_time_wait - st_time_wait
      call update_domains(u(:,tau),local_domain)

!
!
      write (fname, "(a8,i4.4)") "heat.nc.", this_pe
      call nc_setup(local_domain,fname)
      ier=nf_open(fname,NF_WRITE,ncid)
      ier=nf_inq_varid(ncid,'Temp',tempid)
      ier=nf_inq_varid(ncid,'Temp_exact',temp_exactid)
      ier=nf_inq_varid(ncid,'time',timeid)
      ltime = 0

      
! Do the time stepping loop
      do n = 1,nt_steps
        tau = mod(n-1,2)
        taup1 = mod(n,2)
        time = time + del_t
        st_time_comp = MPI_WTIME() 
        call update_model(u,tau,taup1,isc,iec,isd,ied,delx,del_t)
        call exact_sol(uexact,isc,iec,isd,ied,local_domain,time)
        end_time_comp = MPI_WTIME()
        tot_time_comp = tot_time_comp + end_time_comp - st_time_comp
        st_time_wait = MPI_WTIME() 
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end_time_wait = MPI_WTIME()
        tot_time_wait = tot_time_wait + end_time_wait - st_time_wait
        st_time_comm = MPI_WTIME() 
        call update_domains(u(:,taup1),local_domain)
        end_time_comm = MPI_WTIME()
        tot_time_comm = tot_time_comm + end_time_comm - st_time_comm
        if (mod(n,io_int) .eq.0) then
          ltime=ltime+1
          do i = isc,iec
            usave(i) = u(i,taup1)
          enddo
          st_time_io = MPI_WTIME() 
          ier=NF_PUT_VARA_REAL(ncid,tempid,(/1,ltime/),(/ni_local,1/),usave(isc:iec))
          ier=NF_PUT_VARA_REAL(ncid,temp_exactid,(/1,ltime/),(/ni_local,1/),uexact(isc:iec))
          ier=NF_PUT_VARA_REAL(ncid,timeid,(/ltime/),(/1/),time)
          end_time_io = MPI_WTIME()
          tot_time_io = tot_time_io + end_time_io - st_time_io
        endif

      enddo
      tot_time_all = tot_time_comp + tot_time_comm + tot_time_io + tot_time_wait
      call MPI_ALLREDUCE(tot_time_all,tot_time_allproc,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
      end_time_all = MPI_WTime()
      elap_time_all =  end_time_all - st_time_all
      
      write(*,100)this_pe,tot_time_comp,tot_time_comm,tot_time_io,tot_time_wait,tot_time_all,tot_time_allproc,elap_time_all,elap_time_init
100   format(2x,'pe,tot_time_comp,comm,io,wait,all, allproc,elapall,elapinit',2x,i3,8(1pe10.3,2x))
      ier = nf_close(ncid)



      call MPI_Finalize(ierr)
      stop
      end
!
   
      subroutine define_domains(ni,delx,domain_layout,this_pe,domains,npes_tot,xhalo,periodic)
      use mpi
      use heat_types
      integer this_pe,npes_tot,xhalo,ierr
      type(domain1D) :: domains(0:npes_tot-1)
      integer :: ni
      integer :: domain_layout
      integer :: is_per_proc,is_block_index
      integer :: isc,iec,nb_left, nb_right
      integer :: isd,ied
      real :: delx
      real, allocatable, dimension(:) :: x_global  
      logical :: periodic
!
      is_per_proc = ni/domain_layout
!
      is_block_index = mod(this_pe,domain_layout) + 1
!
      isc = (is_block_index - 1)*is_per_proc + 1; isd = isc-xhalo
      iec = isc +  is_per_proc - 1; ied = iec+xhalo



! Find the neighbours
      if (isc .eq. 1 .and. periodic .eq. .false.) then 
        nb_left = MPI_PROC_NULL
      elseif (isc .eq. 1 .and. periodic .eq. .true.) then
        nb_left = this_pe + domain_layout - 1 
      else 
        nb_left = this_pe -1
      endif
!
      if ( iec .eq. ni .and. periodic .eq. .false.) then 
         nb_right = MPI_PROC_NULL
      elseif ( iec .eq. ni .and. periodic .eq. .true.) then
         nb_right = this_pe - domain_layout + 1 
      else
         nb_right = this_pe + 1
      endif
!
!
      allocate (x_global(0:ni+1))
      do i = 0,ni+1
        x_global(i) = delx*i
      enddo


! Tranfer all to domain type
      domains(this_pe)%isc = isc ; domains(this_pe)%iec = iec
      domains(this_pe)%isd = isd ; domains(this_pe)%ied = ied
      domains(this_pe)%nb_left = nb_left; domains(this_pe)%nb_right = nb_right
      domains(this_pe)%pe = this_pe ; domains(this_pe)%npes = npes_tot
      domains(this_pe)%xhalo = xhalo 
      domains(this_pe)%ni = ni 
      domains(this_pe)%ni_local = iec - isc + 1 
      domains(this_pe)%delx = delx 
      allocate (domains(this_pe)%x_local(isc:iec))
      domains(this_pe)%x_local(isc:iec) = x_global(isc:iec) 
!
      end subroutine define_domains

      subroutine update_domains(u,local_domain)
      use mpi
      use heat_types
      type(domain1D) :: local_domain
      real, dimension(local_domain%isd:local_domain%ied) :: u
      real, allocatable, dimension(:)  :: sl_buf,sr_buf,rl_buf,rr_buf
      real, allocatable, dimension(:)  :: sb_buf,st_buf,rb_buf,rt_buf 
      integer :: status(MPI_STATUS_SIZE)

      integer :: nlr ,xhalo,  tag=0
      integer :: isc,iec,isd,ied
      integer :: ni_local,ni,this_pe,nb_left,nb_right

      call get_domain_components(local_domain,isc,iec,isd,ied, &
                 ni_local,ni,this_pe,nb_left,nb_right)
      xhalo = local_domain%xhalo 
 
!
      allocate (sl_buf(xhalo),sr_buf(xhalo),rl_buf(xhalo),rr_buf(xhalo))
      sl_buf = 0.0; sr_buf = 0.0 
      rl_buf = 0.0; rr_buf = 0.0 
! Fill buffers
      sl_buf(1:xhalo) = u(isc:isc+xhalo-1)
      sr_buf(1:xhalo) = u(iec-xhalo+1:iec)
!
!


      call MPI_SENDRECV(sl_buf,xhalo,MPI_REAL,nb_left,tag,rr_buf, xhalo, MPI_REAL,nb_right, &
                      tag,MPI_COMM_WORLD, status, ierr)

      call MPI_SENDRECV(sr_buf,xhalo,MPI_REAL,nb_right,tag,rl_buf, xhalo, MPI_REAL,nb_left, &
                      tag,MPI_COMM_WORLD, status, ierr)
! Now put them back in the u array

      u(isc-xhalo:isc-1) =  rl_buf(1:xhalo)
      u(iec+1:iec+xhalo) =  rr_buf(1:xhalo)
!        
      end subroutine update_domains
!
      subroutine nc_setup(local_domain,fname)
      use heat_types
      type(domain1D) :: local_domain 
      integer :: ier
      character(len=12) ::fname
      integer :: ncid,xdimid,timedimid
      integer :: xid
      integer :: timeid,tempid,temp_exactid
      integer :: npes_tot
      integer :: isc,iec,isd,ied
      integer :: ni_local,ni,this_pe,nb_left,nb_right
      include 'netcdf.inc'
      call get_domain_components(local_domain,isc,iec,isd,ied, &
                 ni_local,ni,this_pe,nb_left,nb_right)

      npes_tot = local_domain%npes

!
      ier=nf_create(fname,nf_clobber,ncid)
      ier=NF_DEF_DIM(ncid,'x',ni_local,xdimid)
      ier=NF_DEF_DIM(ncid,'time',nf_unlimited,timedimid)
      ier=NF_DEF_VAR(ncid,'x',NF_REAL,1,(/xdimid/),xid)
      ier=NF_DEF_VAR(ncid,'time',NF_REAL,1,(/timedimid/),timeid)
      ier=NF_DEF_VAR(ncid,'Temp',NF_REAL,2,(/xdimid,timedimid/),tempid)
      ier=NF_DEF_VAR(ncid,'Temp_exact',NF_REAL,2,(/xdimid,timedimid/),temp_exactid)
      ier=NF_PUT_ATT_TEXT(ncid,xid,'long_name',1,'x')
      ier=NF_PUT_ATT_TEXT(ncid,xid,'units',1,'m')
      ier=NF_PUT_ATT_INT(ncid,xid,'domain_decomposition',NF_INT,4,(/1,ni,isc,iec/))
      ier=NF_PUT_ATT_TEXT(ncid,timeid,'long_name',10,'Julian Day')
      ier=NF_PUT_ATT_TEXT(ncid,timeid,'units',30,'days since 1994-01-01 00:00:00')
      ier=NF_PUT_ATT_TEXT(ncid,tempid,'units',6,'Kelvin')
!
! Global attributes for mppnccombine
      ier = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'filename',12,fname)
      ier = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NumFilesInSet',NF_INT,1,npes_tot)

      ier=NF_ENDDEF(ncid)
      ier=NF_PUT_VARA_REAL(ncid,xid,(/1/),(/ni_local/),local_domain%x_local(isc:iec))
      ier=nf_close(ncid)
      end subroutine nc_setup
!
      subroutine get_domain_components(local_domain,isc,iec,isd,ied, &
                 ni_local,ni,this_pe,nb_left,nb_right)
      use heat_types
      type(domain1D) :: local_domain
      integer :: isc,iec,isd,ied
      integer :: ni_local,ni,this_pe,nb_left,nb_right
      isc = local_domain%isc; iec = local_domain%iec 
      isd = local_domain%isd; ied = local_domain%ied
      ni_local = local_domain%ni_local 
      ni = local_domain%ni 
      nb_left = local_domain%nb_left; nb_right = local_domain%nb_right 
      this_pe = local_domain%pe
      return
      end
!
      subroutine update_model(u,tau,taup1,isc,iec,isd,ied,delx,del_t)
      integer :: tau,taup1,isc,iec,isd,ied
      real :: delx,del_t
      real :: u(isd:ied,0:1)
      do i = isc,iec
        u(i,taup1) = u(i,tau)*(1-2*del_t/delx**2) + &
                     del_t*(u(i-1,tau) + u(i+1,tau))/delx**2 
      enddo
      return
      end
!
      subroutine init_cond(u,isc,iec,isd,ied,local_domain)
      use heat_types
      type(domain1D) :: local_domain
      integer ::  isc,iec
      real :: u(isd:ied,0:1)
      real :: pi
      pi = 4.0*atan(1.0)
      do i = isc,iec
        u(i,0) = sin(pi*local_domain%x_local(i))
      enddo
      return
      end
!
      subroutine exact_sol(uexact,isc,iec,isd,ied,local_domain,t)
      use heat_types
      type(domain1D) :: local_domain
      integer ::  isc,iec
      real :: uexact(isc:iec)
      real :: pi,t
      pi = 4.0*atan(1.0)
      do i = isc,iec
        uexact(i) = sin(pi*local_domain%x_local(i))*exp(-pi*pi*t)
      enddo
      return
      end
!   
      subroutine boundary_cond(u,isd,ied,nb_left,nb_right,ni)
      use mpi
      integer :: isd,ied,nb_left,nb_right
      real :: u(isd:ied,0:1)
      if (nb_left .eq. MPI_PROC_NULL) u(0,:) = 0.0
      if (nb_right .eq. MPI_PROC_NULL) u(ni+1,:) = 0.0
      return
      end
      

  






      
!




      







     
      
       
