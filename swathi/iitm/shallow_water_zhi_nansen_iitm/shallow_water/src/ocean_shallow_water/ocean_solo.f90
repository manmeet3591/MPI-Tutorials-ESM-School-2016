program ocean_solo 
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
  ! Driver for the ocean-only simulations.
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  ! Driver for the ocean-only simulations. Similar to the FMS coupler, but 
  ! allows one to run the ocean model without compiling  other models. 
  ! Much simpler than the full FMS coupler. 
  !</DESCRIPTION>
  !
  !<NAMELIST NAME="ocean_solo_nml">
  !
  ! <DATA NAME="date_init"  TYPE="integer, dimension(6)"  DEFAULT="0">
  !   The date that the current integration starts with. If the restart file
  !    ocean_solo.res is present, date_init will be taken from there.
  ! </DATA>
  ! <DATA NAME="calendar"  TYPE="character(maxlen=17)"  DEFAULT="''">
  !   The calendar type used by the current integration. Valid values are consistent 
  !   with the time_manager module: 'julian', 'noleap', or 'thirty_day'. The value 
  !   'no_calendar' can not be used because the time_manager's date  function are used. 
  !   All values must be lowercase.
  ! </DATA>
  ! <DATA NAME="months "  TYPE="integer"  DEFAULT="0">
  !   The number of months that the current integration will be run for. 
  ! </DATA>
  ! <DATA NAME="days "  TYPE="integer"  DEFAULT="0">
  !   The number of days that the current integration will be run for. 
  ! </DATA>
  ! <DATA NAME="hours"  TYPE="integer"  DEFAULT="0">
  !   The number of hours that the current integration will be run for. 
  ! </DATA>
  ! <DATA NAME="minutes "  TYPE="integer"  DEFAULT="0">
  !   The number of minutes that the current integration will be run for. 
  ! </DATA>
  ! <DATA NAME="seconds"  TYPE="integer"  DEFAULT="0">
  !   The number of seconds that the current integration will be run for. 
  ! </DATA>
  ! <DATA NAME="dt_ocean"  TYPE="integer"  DEFAULT="0">
  !   Ocean model time step in seconds. 
  ! </DATA>
  !
  !</NAMELIST>
  !
  !<NOTE>
  !  <PRE>
  !   1.The actual run length will be the sum of months, 
  !     days, hours, minutes, and seconds. A run length of zero
  !     is not a valid option. 
  !   2.The run length must be an integral multiple of the timestep dt_ocean. 
  !  </PRE>
  !</NOTE>

  use mpp_mod,          only: mpp_pe, mpp_npes, mpp_root_pe, mpp_error, FATAL, stdout, stdlog
  use mpp_io_mod,       only: mpp_open, mpp_close
  use mpp_io_mod,       only: MPP_ASCII, MPP_RDONLY, MPP_APPEND, MPP_SINGLE   
  use fms_mod,          only: fms_init, fms_end, file_exist
  use fms_mod,          only: open_namelist_file, check_nml_error, close_file
  use fms_io_mod,       only: fms_io_exit
  use time_manager_mod, only: time_type, set_date, set_time, get_date, increment_date
  use time_manager_mod, only: operator( > ), operator( * ), operator( /= )
  use time_manager_mod, only: operator( - ), operator( / ), operator( + ), month_name
  use time_manager_mod, only: set_calendar_type, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use diag_manager_mod, only: diag_manager_init, diag_manager_end
  use constants_mod,    only: constants_init
  use ocean_model_mod,  only: ocean_model_init, update_ocean_model, ocean_model_end

  implicit none

  !--- namelist interface

  integer           :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0
  integer           :: date_init(6) = 0
  character(len=64) :: calendar = 'julian' 
  integer           :: dt_ocean
  namelist /ocean_solo_nml/ date_init, calendar, years, months, days, hours, minutes, seconds, dt_ocean

  !--- some other variables
  integer :: unit, io_status, ierr, calendar_type, n
  integer :: num_ocean_calls, date(6)
  integer :: yr, mon, day, hr, min, sec
  type(time_type) :: Time, Time_init, Time_start, Time_end, Time_step, Run_len
  character(len=32) :: month

  ! initialize shared modules

  call fms_init()

  call constants_init

  !--- reading namelist -------
  unit = open_namelist_file('input.nml')
  read  (unit, ocean_solo_nml,iostat=io_status)
  write (stdout(), ocean_solo_nml)  
  write (stdlog(), ocean_solo_nml)
  ierr = check_nml_error(io_status,'ocean_solo_nml')
  call close_file (unit)

  ! set the calendar 
  if (calendar(1:6) == 'julian') then
     calendar_type = JULIAN
  else if (calendar(1:6) == 'noleap') then
     calendar_type = NOLEAP
  else if (calendar(1:10) == 'thirty_day') then
     calendar_type = THIRTY_DAY_MONTHS
  else if (calendar(1:11) == 'no_calendar') then
     calendar_type = NO_CALENDAR
  else if (calendar(1:1) /= ' ') then
     call mpp_error (FATAL,'==>Error from ocean_solo: invalid namelist value for calendar')
  else
     call mpp_error (FATAL,'==>Error from ocean_solo: no namelist value for calendar')
  endif


  !--- get ocean_solo restart: his can override settings from namelist
  if (file_exist('INPUT/ocean_solo.res')) then
     call mpp_open(unit,'INPUT/ocean_solo.res',form=MPP_ASCII,action=MPP_RDONLY)
      read(unit,*) calendar_type
      read(unit,*) date_init
      read(unit,*) date
     call mpp_close(unit)
  endif

  call set_calendar_type (calendar_type)
  call diag_manager_init()

  if (sum(date_init) <= 0) then
     call mpp_error(FATAL,'==>Error from ocean_solo: '// &
          'date_init must be set either in ocean_solo.res or in ocean_solo_nml')
  else
     Time_init  = set_date(date_init(1),date_init(2), date_init(3), &
           date_init(4),date_init(5),date_init(6))
  endif

  if (file_exist('INPUT/ocean_solo.res')) then
      Time_start =  set_date(date(1),date(2),date(3),date(4),date(5),date(6))
  else
      Time_start = Time_init
      date = date_init
  endif

  Time_end         = increment_date(Time_start, years, months, days, hours, minutes, seconds)
  Run_len          = Time_end - Time_start
  Time_step        = set_time(dt_ocean,0)
  num_ocean_calls  = Run_len / Time_step
  Time             = Time_start

!----- initial time cannot be greater than current time -------
  if ( Time_init > Time ) call mpp_error (FATAL, 'program ocean_solo:initial time is greater than current time')

!----- make sure run length is a multiple of ocean time step ------
  if ( num_ocean_calls * Time_step /= Run_len ) call mpp_error (FATAL, &
         'program ocean_solo: run length must be multiple of ocean time step')

  call mpp_open (unit, 'time_stamp.out', form=MPP_ASCII, action=MPP_APPEND,threading=MPP_SINGLE)

  month = month_name(date(2))
  if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

  call get_date (Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
  month = month_name(date(2))
  if ( mpp_pe() == mpp_root_pe() ) write (unit,'(6i4,2x,a3)') date, month(1:3)

  call close_file (unit)  

  call ocean_model_init(Time_init, Time, Time_step) ! initialize ocean model & set boundary conditions

  ! loop over the calls 
  do n = 1, num_ocean_calls
     call update_ocean_model() ! integrate ocean model for desired time
     Time = Time + Time_step
  enddo

  call ocean_model_end()
  call diag_manager_end(Time)

  ! write restart file
  call mpp_open( unit, 'RESTART/ocean_solo.res', nohdrs=.TRUE. )
  if ( mpp_pe().EQ.mpp_root_pe() )then
        write( unit, '(i6,8x,a)' )calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
        call get_date(Time_init,yr,mon,day,hr,min,sec)
        write( unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
             'Model start time:   year, month, day, hour, minute, second'
     call get_date(Time_end ,yr,mon,day,hr,min,sec)
     write( unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
          'Current model time: year, month, day, hour, minute, second'
  end if
  call mpp_close(unit)

  call fms_io_exit
  call fms_end

end program ocean_solo
