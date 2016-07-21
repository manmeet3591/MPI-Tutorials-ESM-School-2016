!
      module heat_types
      type domain2D
        integer :: isc, iec, jsc, jec,xhalo,yhalo
        integer :: isd,ied,jsd,jed
        integer :: nb_left, nb_right, nb_bot, nb_top
        integer :: pe,npes
        integer :: ni,nj
        integer :: ni_local, nj_local
        real :: delx,dely
        real,dimension(:),allocatable :: x_local,y_local
      end type domain2D
      end module heat_types

! For solving the 2D transient heat conduction problem.
!
      program heat_eqn_2D
      use mpi
      use heat_types
      include 'netcdf.inc'
      type(domain2D), allocatable, dimension(:),target  :: domains
      type(domain2D), pointer :: local_domain
      real, allocatable,  dimension(:,:,:) :: u
      real, allocatable,  dimension(:,:) :: usave


      integer :: this_pe,npes_tot, ierr
      integer :: ni_local,nj_local,nb_left,nb_right,nb_bot,nb_top
      integer :: xhalo = 1, yhalo = 1
      logical, dimension(2) :: periodic = (/.false.,.false./) 
      integer ::  isc, iec,n
      integer :: tau=0,taup1=1,nt_steps = 100
      integer :: ni = 20000 ,nj = 15000, io_int = 100000
      real :: del_t = 0.0000000001, del_t_cfl
      double precision :: st_time_comp,end_time_comp,st_time_comm,end_time_comm
      double precision :: st_time_io,end_time_io,st_time_wait,end_time_wait
      double precision :: st_time_all,end_time_all,elap_time_all
      double precision :: st_time_init,end_time_init,st_time_bc,end_time_bc
      double precision :: elap_time_init=0.0,elap_time_bc=0.0
      double precision :: tot_time_comp=0.0, tot_time_comm=0.0, tot_time_io = 0.0
      double precision :: tot_time_wait=0.0,tot_time_all, tot_time_allproc

        

      integer i,j,k,status(MPI_STATUS_SIZE)
      integer, dimension(2) :: domain_layout=(/4,3/)

      real :: delx,dely, pi = 3.14145
!
      integer ncid,ier,timeid,tempid,ltime
      character (len=12) :: fname
      namelist /heat_eqn_nml/ ni,nj,del_t,nt_steps,domain_layout,io_int
      

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, this_pe, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, npes_tot, ierr )
      st_time_all = MPI_WTime()
      open(10,file='heat_eqn_nml')
      read(10,nml=heat_eqn_nml)
      write(*,*)'ni,nj,del_t,nt_steps,domain_layout',ni,nj,del_t,nt_steps,domain_layout
      if ( npes_tot .ne. domain_layout(1)*domain_layout(2) ) then
        write(*,*) 'Domain layout not consistent with numprocs'
        call MPI_Finalize(ierr)
        stop
      endif

      if ( mod(ni,domain_layout(1)) .ne. 0 ) then
        write(*,*) 'ni not consistent with domain_layout(1) ,',ni,domain_layout(1)
        call MPI_Finalize(ierr)
        stop
      endif

      if ( mod(nj,domain_layout(2)) .ne. 0 ) then
        write(*,*) 'nj not consistent with domain_layout(2) ,',nj,domain_layout(2)
        call MPI_Finalize(ierr)
        stop
      endif
! CFL Check
      delx = 1.0/(ni+1); dely=1.0/(nj+1)
      del_t_cfl = 0.5*delx**2*dely**2/(delx**2+dely**2)
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
      call define_domains(ni,nj,delx,dely,domain_layout,this_pe,domains,npes_tot,xhalo,yhalo,periodic)
      call get_domain_components(local_domain,isc,iec,jsc,jec,isd,ied,jsd,jed, &
                 ni_local,nj_local,ni,nj,this_pe,nb_left,nb_right,nb_bot,nb_top)

! 
      allocate(u(isd:ied,jsd:jed,0:1))
      u(:,:,:) = 0.0
      allocate(usave(isc:iec,jsc:jec))
      
!     Boundary conditions
! Not sure how this works for periodic condidtions and halo > 1.
      call boundary_cond(u,isd,ied,jsd,jed,nb_left,nb_right,nb_bot,nb_top,ni,nj)
      
!     Initial conditions
      call init_cond(u,isc,iec,jsc,jec,isd,ied,jsd,jed,local_domain)
      end_time_init = MPI_WTime()
      elap_time_init = elap_time_init + end_time_init - st_time_init
! Allocate buffers for sendrecv

!
      st_time_wait = MPI_WTIME()     
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      end_time_wait = MPI_WTIME()
      tot_time_wait = tot_time_wait + end_time_wait - st_time_wait
      call update_domains(u(:,:,tau),local_domain)

!
!
      write (fname, "(a8,i4.4)") "heat.nc.", this_pe
      call nc_setup(local_domain,fname)
      ier=nf_open(fname,NF_WRITE,ncid)
      ier=nf_inq_varid(ncid,'Temp',tempid)
      ier=nf_inq_varid(ncid,'time',timeid)
      ltime = 0

      
! Do the time stepping loop
      do n = 1,nt_steps
        tau = mod(n-1,2)
        taup1 = mod(n,2)
        st_time_comp = MPI_WTIME() 
        call update_model(u,tau,taup1,isc,iec,jsc,jec,isd,ied,jsd,jed,delx,dely,del_t)
        end_time_comp = MPI_WTIME()
        tot_time_comp = tot_time_comp + end_time_comp - st_time_comp
        st_time_wait = MPI_WTIME() 
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end_time_wait = MPI_WTIME()
        tot_time_wait = tot_time_wait + end_time_wait - st_time_wait
        st_time_comm = MPI_WTIME() 
        call update_domains(u(:,:,taup1),local_domain)
        end_time_comm = MPI_WTIME()
        tot_time_comm = tot_time_comm + end_time_comm - st_time_comm
        if (mod(n,io_int) .eq.0) then
          ltime=ltime+1
          do j = jsc,jec
            do i = isc,iec
              usave(i,j) = u(i,j,taup1)
            enddo
          enddo
          st_time_io = MPI_WTIME() 
          ier=NF_PUT_VARA_REAL(ncid,tempid,(/1,1,ltime/),(/ni_local,nj_local,1/),usave(isc,jsc))
          ier=NF_PUT_VARA_REAL(ncid,timeid,(/ltime/),(/1/),n*1.0)
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
   
      subroutine define_domains(ni,nj,delx,dely,domain_layout,this_pe,domains,npes_tot,xhalo,yhalo,periodic)
      use mpi
      use heat_types
      integer this_pe,npes_tot,xhalo,yhalo,ierr
      type(domain2D) :: domains(0:npes_tot-1)
      integer :: ni,nj
      integer,  dimension(2) :: domain_layout
      integer :: is_per_proc,js_per_proc,js_block_index,is_block_index
      integer :: isc,iec,jsc,jec,nb_left, nb_right,nb_bot,nb_top
      integer :: isd,ied,jsd,jed
      real :: delx,dely
      real, allocatable, dimension(:) :: x_global,y_global  
      logical, dimension(2) :: periodic
!
      is_per_proc = ni/domain_layout(1)
      js_per_proc = nj/domain_layout(2)
!
      is_block_index = mod(this_pe,domain_layout(1)) + 1
      js_block_index = (this_pe - is_block_index +1)/domain_layout(1) + 1
!
      isc = (is_block_index - 1)*is_per_proc + 1; isd = isc-xhalo
      iec = isc +  is_per_proc - 1; ied = iec+xhalo

      jsc = (js_block_index - 1)*js_per_proc + 1; jsd = jsc-yhalo
      jec = jsc +  js_per_proc - 1; jed = jec + yhalo

! Find the neighbours
      if (isc .eq. 1 .and. periodic(1) .eq. .false.) then 
        nb_left = MPI_PROC_NULL
      elseif (isc .eq. 1 .and. periodic(1) .eq. .true.) then
        nb_left = this_pe + domain_layout(1) - 1 
      else 
        nb_left = this_pe -1
      endif
!
      if ( iec .eq. ni .and. periodic(1) .eq. .false.) then 
         nb_right = MPI_PROC_NULL
      elseif ( iec .eq. ni .and. periodic(1) .eq. .true.) then
         nb_right = this_pe - domain_layout(1) + 1 
      else
         nb_right = this_pe + 1
      endif
!
      if ( jsc .eq. 1 .and. periodic(2) .eq. .false.) then
        nb_bot = MPI_PROC_NULL
      elseif ( jsc .eq. 1 .and. periodic(2) .eq. .true.) then
        nb_bot = this_pe + npes_tot - domain_layout(1)
      else
        nb_bot = this_pe - domain_layout(1)
      endif
!
      if ( jec .eq. nj .and. periodic(2) .eq. .false.) then
        nb_top = MPI_PROC_NULL
      elseif ( jec .eq. nj .and. periodic(2) .eq. .true.) then
        nb_top = this_pe  - npes_tot + domain_layout(1)
      else
        nb_top = this_pe + domain_layout(1)
      endif
!
      allocate (x_global(0:ni+1),y_global(0:nj+1))
      do i = 0,ni+1
        x_global(i) = delx*i
      enddo
      do j = 0,nj+1
        y_global(j) = dely*j
      enddo

! Tranfer all to domain type
      domains(this_pe)%isc = isc ; domains(this_pe)%iec = iec
      domains(this_pe)%jsc = jsc ; domains(this_pe)%jec = jec
      domains(this_pe)%isd = isd ; domains(this_pe)%ied = ied
      domains(this_pe)%jsd = jsd ; domains(this_pe)%jed = jed
      domains(this_pe)%nb_left = nb_left; domains(this_pe)%nb_right = nb_right
      domains(this_pe)%nb_bot = nb_bot ; domains(this_pe)%nb_top = nb_top
      domains(this_pe)%pe = this_pe ; domains(this_pe)%npes = npes_tot
      domains(this_pe)%xhalo = xhalo ; domains(this_pe)%yhalo = yhalo
      domains(this_pe)%ni = ni ; domains(this_pe)%nj = nj
      domains(this_pe)%ni_local = iec - isc + 1 ; domains(this_pe)%nj_local = jec - jsc + 1
      domains(this_pe)%delx = delx ; domains(this_pe)%dely = dely
      allocate (domains(this_pe)%x_local(isc:iec),domains(this_pe)%y_local(jsc:jec))
      domains(this_pe)%x_local(isc:iec) = x_global(isc:iec) 
      domains(this_pe)%y_local(jsc:jec) = y_global(jsc:jec)
!
      end subroutine define_domains

      subroutine update_domains(u,local_domain)
      use mpi
      use heat_types
      type(domain2D) :: local_domain
      real, dimension(local_domain%isd:local_domain%ied,local_domain%jsd:local_domain%jed) :: u
      real, allocatable, dimension(:,:)  :: sl_buf,sr_buf,rl_buf,rr_buf
      real, allocatable, dimension(:,:)  :: sb_buf,st_buf,rb_buf,rt_buf 
      integer :: status(MPI_STATUS_SIZE)

      integer :: nlr , nbt, xhalo, yhalo, nj_local_tot, ni_local_tot, tag=0
      integer :: isc,iec,jsc,jec,isd,ied,jsd,jed
      integer :: ni_local,nj_local,ni,nj,this_pe,nb_left,nb_right,nb_bot,nb_top

      call get_domain_components(local_domain,isc,iec,jsc,jec,isd,ied,jsd,jed, &
                 ni_local,nj_local,ni,nj,this_pe,nb_left,nb_right,nb_bot,nb_top)
      xhalo = local_domain%xhalo ; yhalo = local_domain%yhalo

      nj_local_tot = nj_local*xhalo ; ni_local_tot = ni_local*yhalo
!
      allocate (sl_buf(xhalo,nj_local),sr_buf(xhalo,nj_local),rl_buf(xhalo,nj_local),rr_buf(xhalo,nj_local))
      allocate (sb_buf(ni_local,yhalo),st_buf(ni_local,yhalo),rb_buf(ni_local,yhalo),rt_buf(ni_local,yhalo))
      sl_buf = 0.0; sr_buf = 0.0 ; sb_buf = 0.0 ; st_buf = 0.0
      rl_buf = 0.0; rr_buf = 0.0 ; rb_buf = 0.0 ; rt_buf = 0.0
! Fill buffers
      sl_buf(1:xhalo,1:nj_local) = u(isc:isc+xhalo-1,jsc:jec)
      sr_buf(1:xhalo,1:nj_local) = u(iec-xhalo+1:iec,jsc:jec)
!
      sb_buf(1:ni_local,1:yhalo) = u(isc:iec,jsc:jsc+yhalo-1) 
      st_buf(1:ni_local,1:yhalo) = u(isc:iec,jec-yhalo+1:jec) 
!

      call MPI_SENDRECV(sb_buf,ni_local_tot,MPI_REAL,nb_bot,tag,rt_buf, ni_local_tot, MPI_REAL,nb_top, &
                      tag,MPI_COMM_WORLD, status, ierr)
      call MPI_SENDRECV(st_buf,ni_local_tot,MPI_REAL,nb_top,tag,rb_buf, ni_local_tot, MPI_REAL,nb_bot, & 
                      tag,MPI_COMM_WORLD, status, ierr)

      call MPI_SENDRECV(sl_buf,nj_local_tot,MPI_REAL,nb_left,tag,rr_buf, nj_local_tot, MPI_REAL,nb_right, &
                      tag,MPI_COMM_WORLD, status, ierr)

      call MPI_SENDRECV(sr_buf,nj_local_tot,MPI_REAL,nb_right,tag,rl_buf, nj_local_tot, MPI_REAL,nb_left, &
                      tag,MPI_COMM_WORLD, status, ierr)
! Now put them back in the u array

      u(isc-xhalo:isc-1,jsc:jec) =  rl_buf(1:xhalo,1:nj_local)
      u(iec+1:iec+xhalo,jsc:jec) =  rr_buf(1:xhalo,1:nj_local)
      u(isc:iec,jsc-xhalo:jsc-1) =  rb_buf(1:ni_local,1:yhalo)
      u(isc:iec,jec+1:jec+yhalo) =  rt_buf(1:ni_local,1:yhalo)
!        
      end subroutine update_domains
!
      subroutine nc_setup(local_domain,fname)
      use heat_types
      type(domain2D) :: local_domain 
      integer :: ier
      character(len=12) ::fname
      integer :: ncid,xdimid,ydimid,timedimid
      integer :: xid,yid
      integer :: timeid,tempid
      integer :: npes_tot
      integer :: isc,iec,jsc,jec,isd,ied,jsd,jed
      integer :: ni_local,nj_local,ni,nj,this_pe,nb_left,nb_right,nb_bot,nb_top
      include 'netcdf.inc'
      call get_domain_components(local_domain,isc,iec,jsc,jec,isd,ied,jsd,jed, &
                 ni_local,nj_local,ni,nj,this_pe,nb_left,nb_right,nb_bot,nb_top)

      npes_tot = local_domain%npes

!
      ier=nf_create(fname,nf_clobber,ncid)
      ier=NF_DEF_DIM(ncid,'x',ni_local,xdimid)
      ier=NF_DEF_DIM(ncid,'y',nj_local,ydimid)
      ier=NF_DEF_DIM(ncid,'time',nf_unlimited,timedimid)
      ier=NF_DEF_VAR(ncid,'x',NF_REAL,1,(/xdimid/),xid)
      ier=NF_DEF_VAR(ncid,'y',NF_REAL,1,(/ydimid/),yid)
      ier=NF_DEF_VAR(ncid,'time',NF_REAL,1,(/timedimid/),timeid)
      ier=NF_DEF_VAR(ncid,'Temp',NF_REAL,3,(/xdimid,ydimid,timedimid/),tempid)
      ier=NF_PUT_ATT_TEXT(ncid,xid,'long_name',1,'x')
      ier=NF_PUT_ATT_TEXT(ncid,xid,'units',1,'m')
      ier=NF_PUT_ATT_INT(ncid,xid,'domain_decomposition',NF_INT,4,(/1,ni,isc,iec/))
      ier=NF_PUT_ATT_TEXT(ncid,yid,'long_name',1,'y')
      ier=NF_PUT_ATT_TEXT(ncid,yid,'units',1,'m')
      ier=NF_PUT_ATT_INT(ncid,yid,'domain_decomposition',NF_INT,4,(/1,nj,jsc,jec/))
      ier=NF_PUT_ATT_TEXT(ncid,timeid,'long_name',10,'Julian Day')
      ier=NF_PUT_ATT_TEXT(ncid,timeid,'units',31,'days since 1994-01-01 00:00:00')
      ier=NF_PUT_ATT_TEXT(ncid,tempid,'units',6,'Kelvin')
!
! Global attributes for mppnccombine
      ier = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'filename',12,fname)
      ier = NF_PUT_ATT_INT(ncid,NF_GLOBAL,'NumFilesInSet',NF_INT,1,npes_tot)

      ier=NF_ENDDEF(ncid)
      ier=NF_PUT_VARA_REAL(ncid,xid,(/1/),(/ni_local/),local_domain%x_local(isc:iec))
      ier=NF_PUT_VARA_REAL(ncid,yid,(/1/),(/nj_local/),local_domain%y_local(jsc:jec))
      ier=nf_close(ncid)
      end subroutine nc_setup
!
      subroutine get_domain_components(local_domain,isc,iec,jsc,jec,isd,ied,jsd,jed, &
                 ni_local,nj_local,ni,nj,this_pe,nb_left,nb_right,nb_bot,nb_top)
      use heat_types
      type(domain2D) :: local_domain
      integer :: isc,iec,jsc,jec,isd,ied,jsd,jed
      integer :: ni_local,nj_local,ni,nj,this_pe,nb_left,nb_right,nb_bot,nb_top
      isc = local_domain%isc; iec = local_domain%iec 
      jsc = local_domain%jsc ; jec = local_domain%jec
      isd = local_domain%isd; ied = local_domain%ied
      jsd = local_domain%jsd ; jed = local_domain%jed
      ni_local = local_domain%ni_local; nj_local = local_domain%nj_local 
      ni = local_domain%ni ; nj = local_domain%nj
      nb_left = local_domain%nb_left; nb_right = local_domain%nb_right 
      nb_bot = local_domain%nb_bot ; nb_top = local_domain%nb_top
      this_pe = local_domain%pe
      return
      end
      subroutine update_model(u,tau,taup1,isc,iec,jsc,jec,isd,ied,jsd,jed,delx,dely,del_t)
      integer :: tau,taup1,isc,iec,jsc,jec,isd,ied,jsd,jed
      real :: delx,dely,del_t,delx_sq,dely_sq,factor
      real :: u(isd:ied,jsd:jed,0:1)
      delx_sq = delx*delx; dely_sq = dely*dely
      factor = 2*del_t*(delx_sq + dely_sq)/(delx_sq*dely_sq) 
      do j = jsc,jec
        do i = isc,iec
          u(i,j,taup1) = u(i,j,tau)*(1-factor) + &
                           del_t*((u(i-1,j,tau) + u(i+1,j,tau))/delx_sq + &
                                  (u(i,j-1,tau) + u(i,j+1,tau))/dely_sq)
        enddo
      enddo
      return
      end
      subroutine init_cond(u,isc,iec,jsc,jec,isd,ied,jsd,jed,local_domain)
      use heat_types
      type(domain2D) :: local_domain
      integer ::  isc,iec,jsc,jec
      real :: u(isd:ied,jsd:jed,0:1)
      real :: pi
      pi = 4.0*atan(1.0)
      do j = jsc,jec
        do i = isc,iec
          u(i,j,0) = sin(pi*local_domain%x_local(i))*sin(pi*local_domain%y_local(j))
        enddo
      enddo
      return
      end
!   
      subroutine boundary_cond(u,isd,ied,jsd,jed,nb_left,nb_right,nb_bot,nb_top,ni,nj)
      use mpi
      integer :: isd,ied,jsd,jed,nb_left,nb_right,nb_bot,nb_top
      real :: u(isd:ied,jsd:jed,0:1)
      if (nb_left .eq. MPI_PROC_NULL) u(0,:,:) = 0.0
      if (nb_right .eq. MPI_PROC_NULL) u(ni+1,:,:) = 0.0
      if (nb_bot .eq. MPI_PROC_NULL) u(:,0,:) = 0.0
      if (nb_top .eq. MPI_PROC_NULL) u(:,nj+1,:) = 0.0
      return
      end
      

  






      
!




      







     
      
       
