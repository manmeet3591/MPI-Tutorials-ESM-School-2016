module ocean_freesurf_mod
  !-----------------------------------------------------------------------
  !                   GNU General Public License                        
  !                                                                      
  ! This program is free software; you can redistribute it and/or modify it and  
  ! are expected to follow the terms of the GNU General Public License  
  ! as published by the Free Software Foundation; either version 2 of   
  ! the License, or (at your option) any later version.                 
  !                                                                      
  ! MOM is distributed in the hope that it will be useful, but WITHOUT    
  ! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
  ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
  ! License for more details.                                           
  !                                                                      
  ! For the full text of the GNU General Public License,                
  ! write to: Free Software Foundation, Inc.,                           
  !           675 Mass Ave, Cambridge, MA 02139, USA.                   
  ! or see:   http://www.gnu.org/licenses/gpl.html                      
  !-----------------------------------------------------------------------
  ! 
  ! <CONTACT EMAIL= "Ronald.Pacanowski@noaa.gov">Ronald Pacanowski </CONTACT>
  ! <CONTACT EMAIL= "Zhi.Liang@noaa.gov">Zhi Liang </CONTACT>

  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
  ! Time step free surface height
  !<OVERVIEW>
  !</OVERVIEW>
  !<DESCRIPTION>
  ! This module steps free surface height forward in time using a 
  ! leap-frog time stepping scheme.
  !</DESCRIPTION>

  use mpp_mod,                only: mpp_chksum, mpp_error, FATAL, stdout, stdlog
  use mpp_domains_mod,        only: mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod,        only: mpp_update_domains, domain2d
  use fms_mod,                only: open_namelist_file, check_nml_error, file_exist
  use fms_mod,                only: read_data, write_data, write_version_number, close_file
  use time_manager_mod,       only: time_type, increment_time
  use time_manager_mod,       only: operator( /= ), operator( == ), operator( > )
  use diag_manager_mod,       only: register_diag_field, send_data
  use ocean_types_mod,        only: ocean_grid_type, ocean_velocity_type
  use ocean_types_mod,        only: ocean_time_type, ocean_freesurf_type
  use ocean_horz_diffuse_mod, only: ocean_horz_diffuse

  implicit none
  private


  !--- public interface ----

  public :: ocean_freesurf_init, update_ocean_freesurf, ocean_freesurf_end

  !---namelist interface ----
  !<NAMELIST NAME="ocean_freesurf_nml">
  ! <DATA NAME="anomoly_amplitude" TYPE="real" DEFAULT="0.2">
  !   free surface height anomoly (meter).
  ! </DATA>
  ! <DATA NAME="anomoly_longitude" TYPE="real" DEFAULT="10.0">
  !   longitude (degrees) of free surface height anomoly center.
  ! </DATA>
  ! <DATA NAME="anomoly_latitude" TYPE="real" DEFAULT="0.0">
  !   latitude (degrees) of free surface height anomoly center.
  ! </DATA>  
  ! <DATA NAME="holding_time" TYPE="real" DEFAULT="0.0">
  !   time (days) to maintain the initial free elevation
  !   after which it is allowed to evolve.  
  ! </DATA>
  ! <DATA NAME="anomoly_scale" TYPE="real" DEFAULT="1.0">
  !  lateral gaussian scale (degrees) for free surface height anomoly
  ! </DATA>
  ! <DATA NAME="newtonian_timescale" TYPE="real" DEFAULT="100.0">
  !   time scale (days) for newtonian damping of free surface elevation to zero
  ! </DATA>
  ! <DATA NAME="horz_diffuse_on" TYPE="logical" DEFAULT="FALSE">
  !   Set true to add lateral diffusion.
  ! </DATA>
  ! <DATA NAME="freshwater_flux_on" TYPE="logical" DEFAULT="FALSE">
  !   Set true to add a fresh water flux (m/s).
  ! </DATA> 
  !</NAMELIST>

  real              :: holding_time        = 0.0
  real              :: anomoly_amplitude   = 0.2
  real              :: anomoly_longitude   = 10.0
  real              :: anomoly_latitude    = 0.0
  real              :: anomoly_scale       = 1.0
  real              :: newtonian_timescale = 100.0
  logical           :: horz_diffuse_on     = .false.
  logical           :: freshwater_flux_on  = .false.
  namelist /ocean_freesurf_nml/ anomoly_amplitude, anomoly_longitude, anomoly_latitude, anomoly_scale, &
                                newtonian_timescale, holding_time, horz_diffuse_on, freshwater_flux_on

  real,                    parameter :: secpday = 86400.0
  type(time_type)                    :: Time_free_surf
  type(domain2d),            pointer :: Domain
  type(ocean_grid_type),     pointer :: Grid
  type(ocean_time_type),     pointer :: Time
  logical                            :: nonlinear
  integer                            :: isc, iec, jsc, jec
  integer                            :: id_eta
  real, dimension(:,:),  allocatable :: wflux      ! fresh water flux at the tracer cell center (m/s)
  logical                            :: module_is_initialized = .FALSE.
  real                               :: newt_coeff

  !--- version information variables ----
  character(len=128) :: version = '$Id: ocean_freesurf.f90,v 1.1.2.1.2.3 2005/03/01 15:41:31 z1l Exp $'
  character(len=128) :: tagname = '$Name:  $'


  contains

  !###################################################################
    ! <SUBROUTINE NAME="ocean_freesurf_init">
    !
    ! <DESCRIPTION>
    ! Initialization routine
    ! </DESCRIPTION>
  subroutine ocean_freesurf_init(Domain_in, Grid_in, Time_in, nonlinear_in, Freesurf)
    type(domain2d), target,       intent(in) :: Domain_in
    type(ocean_grid_type), target,intent(in) :: Grid_in
    type(ocean_time_type), target,intent(in) :: Time_in
    logical,                      intent(in) :: nonlinear_in
    type(ocean_freesurf_type), intent(inout) :: Freesurf

    character(len=32) :: filename = 'INPUT/ocean_freesurf.res.nc'
    character(len=32) :: filename_wflux = 'INPUT/wflux.nc'
    integer           :: i, j, tau, taum1, taup1, isd, ied, jsd, jed
    integer           :: ioun, io_status, ierr
    real              :: rade, emax, argy, arg


    if ( module_is_initialized ) then 
       call mpp_error(FATAL, '==>Error from ocean_freesurf_mod (ocean_freesurf_init): module already initialized')
    endif

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

    ! provide for namelist over-ride of defaults 
    ioun = open_namelist_file()
    read  (ioun, ocean_freesurf_nml,iostat=io_status)
    write (stdout(),'(/)')
    write (stdout(), ocean_freesurf_nml)  
    write (stdlog(), ocean_freesurf_nml)
    ierr = check_nml_error(io_status, 'ocean_freesurf_nml')
    call close_file(ioun)

    Domain    => Domain_in
    Grid      => Grid_in
    Time      => Time_in
    nonlinear = nonlinear_in

    taum1     = Time%taum1
    tau       = Time%tau
    taup1     = Time%taup1   

    Time_free_surf = increment_time(Time%Time_init, nint(holding_time*secpday), 0) 


    call mpp_get_compute_domain(Domain,isc,iec,jsc,jec)        
    call mpp_get_data_domain   (Domain,isd,ied,jsd,jed)

    allocate(Freesurf%eta(isd:ied,jsd:jed,3))
    allocate(Freesurf%ud(isd:ied,jsd:jed), Freesurf%vd(isd:ied,jsd:jed) )

    !--- initialize to trivial value
    Freesurf%ud  = 0.0
    Freesurf%vd  = 0.0
    Freesurf%eta = 0.0

    !--- initial condition ---------------------------------------------
    ! The freesurface elevation is zero with an optional gaussian anomoly added

    if(Time%Model_time == Time%Time_init) then 
       write(stdout(),*)'==> Note: The shallow water model will start from initial conditions.' 
       if(holding_time > 0) write(stdout(),*) '=> Warning: The initial free surface anomoly will be maintained for ', &
            holding_time,' days'
       rade = 1.0/(anomoly_scale**2)
       emax = 0.0

       do j = jsc, jec
          argy = (Grid%yt(j) - anomoly_latitude)**2
          do i = isc, iec
             arg = rade*((Grid%xt(i) - anomoly_longitude)**2 + argy)
             arg = min(arg, 15.0)
             Freesurf%eta(i,j,:) = anomoly_amplitude*exp(-arg)*Grid%tmask(i,j)
             emax = max(emax, Freesurf%eta(i,j,1))
          enddo
       enddo

       if(emax > 0) write(stdout(),*)'=>NOTE: maximum free surface anomoly is ', emax
    else
       ! --- read restart file if exists, otherwise it will crash.
       if(file_exist(filename) ) then
          call read_data(filename,'eta',Freesurf%eta(isc:iec,jsc:jec,tau), Domain,timelevel=1) 
          call read_data(filename,'eta',Freesurf%eta(isc:iec,jsc:jec,taup1), Domain,timelevel=2) 
       else
          call mpp_error(FATAL,'==> Error from ocean_freesurf_mod: Expected file INPUT/ocean_freesurf.res.nc does not exist')
       endif
    endif

    call mpp_update_domains(Freesurf%eta, Domain)

    !--- boundary condition.
    if (freshwater_flux_on) then
        allocate(wflux(isc:iec,jsc:jec))
        wflux = 0.0
       if(file_exist(filename_wflux) ) then
          write(stdout(),*)'=>NOTE: wflux is being read from file.'
          call read_data(filename_wflux,'wflux',wflux(isc:iec,jsc:jec), Domain) 
       else
          call mpp_error(FATAL,'==> Error from ocean_freesurf_mod: Expected file INPUT/wflux.nc does not exist')
       endif
    endif

    id_eta = register_diag_field ('ocean_model', 'eta',(/Grid%id_xt, Grid%id_yt/), &
                                  Time%Model_time, 'free surface height', 'm',     &
                                  missing_value=-999.0, range=(/-999.0,999.0/) )

    !--- set newtonian cooling time scale
    newt_coeff = 1.0/(newtonian_timescale*secpday) 

    module_is_initialized = .true.

  end subroutine ocean_freesurf_init
  ! </SUBROUTINE>

  !###################################################################
! <SUBROUTINE NAME="update_ocean_freesurf">
!
! <DESCRIPTION>
! Time step the freesurf fields using a leap-frog scheme. 
! </DESCRIPTION>
!
  subroutine update_ocean_freesurf(Velocity, Freesurf)
    type(ocean_velocity_type),    intent(in) :: Velocity
    type(ocean_freesurf_type), intent(inout) :: Freesurf

    integer :: i, j, taum1, tau, taup1
    logical :: used
    real    :: twodt

    twodt = 2.0*Time%dt
    taum1 = Time%taum1
    tau   = Time%tau
    taup1 = Time%taup1

    !--- construct quantities:  Freesurf%ud, Freesurf%vd
    if(nonlinear) then
       do j = jsc, jec
          do i = isc, iec
             Freesurf%ud(i,j) = 0.5*(Grid%ht(i,j) + Freesurf%eta(i,j,tau)            &
                                + Grid%ht(i-1,j) + Freesurf%eta(i-1,j,tau))*Velocity%u(i,j,tau)
             Freesurf%vd(i,j) = 0.5*(Grid%ht(i,j) + Freesurf%eta(i,j,tau)            &
                                + Grid%ht(i,j-1) + Freesurf%eta(i,j-1,tau))*Velocity%v(i,j,tau)*Grid%cosc(j)
          enddo
       enddo
    else          
       do j = jsc, jec
          do i = isc, iec
             Freesurf%ud(i,j) = 0.5*(Grid%ht(i,j) + Grid%ht(i-1,j))*Velocity%u(i,j,tau)
             Freesurf%vd(i,j) = 0.5*(Grid%ht(i,j) + Grid%ht(i,j-1))*Velocity%v(i,j,tau)*Grid%cosc(j)
          enddo
       enddo
    endif

    call mpp_update_domains(Freesurf%ud, Domain)
    call mpp_update_domains(Freesurf%vd, Domain)

    !--- solve free surface (eta must remain zero on land points)

    if (Time%Model_time > Time_free_surf ) then
       do j = jsc, jec
          do i = isc, iec
             Freesurf%eta(i,j,taup1) = Freesurf%eta(i,j,taum1) + Grid%tmask(i,j)*twodt*( &
                                     - (Grid%rcosdxt(i,j)*(Freesurf%ud(i+1,j)-Freesurf%ud(i,j))&
                                     + Grid%rcosdyt(i,j)*(Freesurf%vd(i,j+1)-Freesurf%vd(i,j)))&
                                     - newt_coeff*Freesurf%eta(i,j,taum1))
          enddo
       enddo

       !--- adding freshwater
       if (freshwater_flux_on) then
          Freesurf%eta(isc:iec,jsc:jec,taup1) = Freesurf%eta(isc:iec,jsc:jec,taup1) + Grid%tmask(i,j)*twodt*&
                                                wflux(isc:iec,jsc:jec)
       endif

       !--- adding lateral diffusion
       if (horz_diffuse_on) then
          Freesurf%eta(isc:iec,jsc:jec,taup1) = Freesurf%eta(isc:iec,jsc:jec,taup1) + Grid%tmask(i,j)*twodt*&
                                                ocean_horz_diffuse(Freesurf%eta(:,:,taum1) )
       endif

       call mpp_update_domains(Freesurf%eta(:,:,taup1), Domain)
    endif

    !--- send out data
    if (id_eta > 0) used = send_data (id_eta, Freesurf%eta(isc:iec,jsc:jec,taup1), &
                                      Time%Model_time, rmask=Grid%tmask(isc:iec,jsc:jec))


  end subroutine update_ocean_freesurf
  ! </SUBROUTINE>

  !#####################################################################
  ! <SUBROUTINE NAME="ocean_freesurf_end">
  !
  ! <DESCRIPTION>
  !  Write ocean freesurf restart file.
  ! </DESCRIPTION>
  !
  subroutine ocean_freesurf_end(Freesurf)
    type(ocean_freesurf_type), intent(inout) :: Freesurf

    character(len=32) :: filename = 'RESTART/ocean_freesurf.res.nc'
    integer           :: tau, taup1

    tau   = Time%tau
    taup1 = Time%taup1

    call write_data(filename, 'eta', Freesurf%eta(isc:iec,jsc:jec,tau),   domain=Domain)
    call write_data(filename, 'eta', Freesurf%eta(isc:iec,jsc:jec,taup1), domain=Domain)

    write(stdout(),*) 'Ending free surface chksum (tau) ==>', &
               mpp_chksum(Freesurf%eta(isc:iec,jsc:jec,tau)*Grid%tmask(isc:iec,jsc:jec))
    write(stdout(),*) 'Ending free surface chksum (taup1) ==>', &
               mpp_chksum(Freesurf%eta(isc:iec,jsc:jec,taup1)*Grid%tmask(isc:iec,jsc:jec))

    deallocate(Freesurf%eta, Freesurf%ud, Freesurf%vd)
    if (freshwater_flux_on) then
      deallocate(wflux)
    endif
  end subroutine ocean_freesurf_end
  ! </SUBROUTINE>

  !#####################################################################

end module ocean_freesurf_mod
