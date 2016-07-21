module ocean_model_mod
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

  !<OVERVIEW>
  ! Time step the ocean model.

  !</OVERVIEW>
  !<DESCRIPTION>
  ! This module steps the shallow water equations forward in time using a leap-frog time scheme. 
  ! shallow water equations of R. Sadourny (JAS,1975 vol 32 pp 2103,2110) implemented 
  ! with geometry and topography on a sphere by Ron C. Pacanowski   12/92 
  !</DESCRIPTION>
  use mpp_mod,            only: mpp_npes, mpp_pe, mpp_error, FATAL, stdout, stdlog
  use mpp_mod,            only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_min
  use mpp_domains_mod,    only: mpp_define_domains, mpp_define_layout
  use mpp_domains_mod,    only: domain2d, CYCLIC_GLOBAL_DOMAIN, mpp_get_compute_domain
  use mpp_domains_util_mod, only: mpp_domains_set_stack_size
  use fms_mod,            only: open_namelist_file, check_nml_error, close_file
  use fms_mod,            only: write_version_number
  use constants_mod,      only: GRAV
  use time_manager_mod,   only: time_type, get_time, increment_time, operator( + )
  use ocean_types_mod,    only: ocean_grid_type, ocean_time_type
  use ocean_types_mod,    only: ocean_velocity_type, ocean_freesurf_type 
  use ocean_grid_mod,     only: ocean_grid_init, set_ocean_grid_size
  use ocean_grid_mod,     only: set_ocean_grid, ocean_grid_end
  use ocean_velocity_mod, only: ocean_velocity_init, update_ocean_velocity, ocean_velocity_end
  use ocean_freesurf_mod, only: ocean_freesurf_init, update_ocean_freesurf, ocean_freesurf_end 
  use ocean_horz_diffuse_mod, only: ocean_horz_diffuse_init

  implicit none
  private

  !--- public interface ----
  public :: ocean_model_init, update_ocean_model, ocean_model_end

  !--- namelist interface ----
  !<NAMELIST NAME="ocean_model_nml">
  ! <DATA NAME="nonlinear" TYPE="logical" DEFAULT="TRUE" >
  !  If true, the nonlinear term is included in the shallow water equation. 
  !  If false, nonlinear term is not included.
  ! </DATA> 
  ! <DATA NAME="robert" TYPE="real" DEFAULT="0.001">
  !   coefficient for time filter  
  ! </DATA>
  !</NAMELIST>

  logical           :: nonlinear = .true.
  real              :: robert    = 0.001
  namelist /ocean_model_nml/ nonlinear, robert
                             
  type(ocean_time_type), target, save :: Time
  type(domain2D),        target, save :: Domain  
  type(ocean_grid_type), target, save :: Grid
  type(ocean_freesurf_type),     save :: Freesurf
  type(ocean_velocity_type),     save :: Velocity

  real, parameter :: secpday = 86400.0 
  integer         :: id_ocean              ! id for the clock to time update_ocean_model
  logical         :: module_is_initialized = .FALSE.

  !--- version information variables ----
  character(len=128) :: version = '$Id: ocean_model.f90,v 1.1.2.3.2.4 2005/03/01 16:17:56 z1l Exp $'
  character(len=128) :: tagname = '$Name:  $'

contains

  !#####################################################################
  ! <SUBROUTINE NAME="ocean_model_init">
  !
  ! <DESCRIPTION>
  ! Initialize the ocean model 
  ! </DESCRIPTION>

  subroutine ocean_model_init(Time_init, Model_time, Time_step)
    type(time_type), intent(in) :: Time_init, Model_time, Time_step

    integer         :: unit, io_status, ierr
    integer         :: npes, layout(2) = (/1,0/)
    integer         :: seconds, days, icg, jcg, isc, iec, jsc, jec, i, j
    real            :: gridmin, max_dt_for_gravity_wave, dtcg,  max_dt_for_gravity_wave0

    if ( module_is_initialized ) then
       call mpp_error(FATAL, '==>Error from ocean_grid_mod (ocean_grid_init): module has been initialized')
    endif

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

    id_ocean = mpp_clock_id('ocean')

    !--- read the namelist variables ---
    unit = open_namelist_file('input.nml')
    read  (unit, ocean_model_nml,iostat=io_status)
    write (stdout(), ocean_model_nml)  
    write (stdlog(), ocean_model_nml)
    ierr = check_nml_error(io_status,'ocean_model_nml')
    call close_file (unit)

    !--- write out some namelist information
    if(nonlinear) then
       write(stdout(),*)'==> NOTE: The shallow water model nonlinear term includes '
    else
       write(stdout(),*)'==> NOTE: The shallow water model nonlinear term is neglect '
    endif

    Time%Time_init  = Time_init
    Time%Model_Time = Model_time
    Time%Time_step  = Time_step
    Time%taum1      = 1
    Time%tau        = 2
    Time%taup1      = 3
    Time%itt        = 0

    !--- define time step, if user didn't specify time step, 
    !--- it will be computed based on CFL criteria.
    call get_time(Time_step, seconds, days)
    Time%dt = days*secpday + seconds

    !--- call ocean_grid initialization routine
    call ocean_grid_init

    !--- get number of grids
    call set_ocean_grid_size(Grid)

    !--- domain decompsition -------------------------------------------
    npes = mpp_npes()
! Swathi.. To make it run on 6 procs or less
    if (npes .le. 6) call mpp_domains_set_stack_size(1000000)
!
    call mpp_define_layout((/1,Grid%ni,1,Grid%nj/),npes,layout)
    if(Grid%cyclic) then
       call mpp_define_domains((/1,Grid%ni,1,Grid%nj/),layout, Domain, &
                               xhalo=1, yhalo=1, xflags = CYCLIC_GLOBAL_DOMAIN)
    else
       call mpp_define_domains((/1,Grid%ni,1,Grid%nj/),layout, Domain, xhalo=1, yhalo=1)
    endif

    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)

    call set_ocean_grid(Domain, Grid)

    !--- check the CFL condition ----------
    icg = isc; jcg = jsc; gridmin = 1.0e20; max_dt_for_gravity_wave = 1.0e20
    do j=jsc,jec
       do i=isc,iec
          if (Grid%tmask(i,j) > 0.0) then
              gridmin = min(gridmin,Grid%dyt(i,j),Grid%dxt(i,j)*Grid%cost(j) )
              dtcg = 0.5*gridmin/sqrt(GRAV*Grid%ht(i,j))          ! 0.5 arises from use of leap-frog 
              if (dtcg < max_dt_for_gravity_wave) then
                  max_dt_for_gravity_wave = dtcg; icg  = i; jcg  = j
              endif
          endif
       enddo
    enddo

    max_dt_for_gravity_wave = nint(max_dt_for_gravity_wave)
    max_dt_for_gravity_wave = max_dt_for_gravity_wave + 0.001*mpp_pe() ! to separate redundancies
    max_dt_for_gravity_wave0 = max_dt_for_gravity_wave
    call mpp_min (max_dt_for_gravity_wave)

! show the most unstable location for internal gravity waves
  
    if (max_dt_for_gravity_wave == max_dt_for_gravity_wave0 ) then
        if (Time%dt > max_dt_for_gravity_wave) then
            write (stdout(),'(/a,i4,a,i4,a,f6.2,a,f6.2,a)')'=>Error: Stability is violated at eta-cell (i,j) = (',&
                 icg,',',jcg,'), (lon,lat) = (',Grid%xt(icg),',',Grid%yt(jcg),').'
            write(stdout(),'(a,e12.6)') ' The dxt grid distance (m) at this point is ',Grid%dxt(icg,jcg) 
            write(stdout(),'(a,e12.6)') ' The dyt grid distance (m) at this point is ',Grid%dyt(icg,jcg) 
            call mpp_error(FATAL,'==>Error: time step instability detected for gravity wave in ocean_model_mod')
         else
            write(stdout(),'(a,f5.2,a)') '=>Note: the timestep was specified = ',Time%dt,' sec.'
            write(stdout(),'(a,f5.2,a)') &
            '  Linear stability theory restricts the max dt to < ',max_dt_for_gravity_wave,' sec.'
        endif 
    endif

    !--- call ocean_velocity_init
    call ocean_velocity_init(Domain, Grid, Time, nonlinear, Velocity)

    !--- call ocean_freesurf_init
    call ocean_freesurf_init(Domain, Grid, Time, nonlinear, Freesurf)

    !--- call ocean_horz_diffuse_init
    call ocean_horz_diffuse_init(Grid, Domain)

    module_is_initialized = .true.

  end subroutine ocean_model_init
  !</SUBROUTINE>

  !#####################################################################
  ! <SUBROUTINE NAME="update_ocean_model">
  !
  ! <DESCRIPTION>
  ! Update in time the ocean model fields. 
  ! </DESCRIPTION>
  !
  subroutine update_ocean_model
    integer :: taum1, tau, taup1

    call mpp_clock_begin(id_ocean)

    !  increment ocean time
    Time%Model_time = Time%Model_Time + Time%Time_step
    Time%itt   = Time%itt+1
    taum1      = mod(Time%itt  ,3) + 1
    tau        = mod(Time%itt+1,3) + 1
    taup1      = mod(Time%itt+2,3) + 1   
    Time%taum1 = taum1
    Time%tau   = tau
    Time%taup1 = taup1

    call update_ocean_freesurf(Velocity, Freesurf)

    call update_ocean_velocity(Velocity, Freesurf)

    ! apply robert time filter to all predicted variables

    Velocity%u(:,:,tau)   = Velocity%u(:,:,tau) + robert*(Velocity%u(:,:,taup1)     &
                          -2.0*Velocity%u(:,:,tau)+Velocity%u(:,:,taum1))
    Velocity%v(:,:,tau)   = Velocity%v(:,:,tau) + robert*(Velocity%v(:,:,taup1)     &
                          -2.0*Velocity%v(:,:,tau)+Velocity%v(:,:,taum1))
    Freesurf%eta(:,:,tau) = Freesurf%eta(:,:,tau) + robert*(Freesurf%eta(:,:,taup1) &
                          -2.0*Freesurf%eta(:,:,tau)+Freesurf%eta(:,:,taum1))

    call mpp_clock_end(id_ocean)

  end subroutine update_ocean_model
  !</SUBROUTINE>

  !#####################################################################
  ! <SUBROUTINE NAME="ocean_model_end">
  !
  ! <DESCRIPTION>
  ! Close down the ocean model 
  ! </DESCRIPTION>
  !
  subroutine  ocean_model_end
    character(len=32) :: filename = 'RESTART/ocean_model.res'

    call ocean_velocity_end(Velocity)

    call ocean_freesurf_end(Freesurf)

    call ocean_grid_end(Grid)

    module_is_initialized = .false.


  end subroutine  ocean_model_end
  !</SUBROUTINE>

  !#####################################################################


end module ocean_model_mod
