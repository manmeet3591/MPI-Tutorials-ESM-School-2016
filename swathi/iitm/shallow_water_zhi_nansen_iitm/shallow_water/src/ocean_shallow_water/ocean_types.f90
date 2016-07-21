module ocean_types_mod
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
  ! This module contains type declarations for ocean model.
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  ! This module contains type declarations for ocean model.
  !</DESCRIPTION>
  !

  use time_manager_mod, only : time_type

  implicit none
  private


  public :: ocean_grid_type, ocean_time_type, ocean_velocity_type, ocean_freesurf_type

!<PUBLICTYPE >
  type ocean_grid_type
     integer                       :: ni, nj            ! grid size
     integer                       :: id_xt, id_yt      ! axis ids for fms diag_manager
     integer                       :: id_xc, id_yc      ! axis ids for fms diag_manager
     logical                       :: cyclic            ! cyclic boundary condition
     real                          :: totarea           ! total area of working domain
     real, dimension(:),   pointer :: xt       =>NULL() ! longitude of the T-cell center in degrees.
     real, dimension(:),   pointer :: xc       =>NULL() ! longitude of the C-cell center  
                                                        ! (southwest corner of T-cell) in degrees.
     real, dimension(:),   pointer :: yt       =>NULL() ! latitude  of the T-cell center in degrees.
     real, dimension(:),   pointer :: yc       =>NULL() ! latitude  of the C-cell center 
                                                        ! (southwest corner of T-cell)in degrees.
     real, dimension(:),   pointer :: cost     =>NULL() ! cosine of T-cell latitude
     real, dimension(:),   pointer :: cosc     =>NULL() ! cosine of C-cell latitude
     real, dimension(:,:), pointer :: dxt      =>NULL() ! longitudinal width of T-cells at grid point (m) 
                                                        ! without multiplying cos(yt)
     real, dimension(:,:), pointer :: dxc      =>NULL() ! longitudinal width of C-cells at grid point (m) 
                                                        ! without multiplying cos(yc)
     real, dimension(:,:), pointer :: dyt      =>NULL() ! latitudinal width of T-cells at grid point (m) 
     real, dimension(:,:), pointer :: dyc      =>NULL() ! latitudinal width of C-cells at grid point (m) 
     real, dimension(:,:), pointer :: dxte     =>NULL() ! longitudinal width between grid points at i+1 and i 
                                                        ! in T-cells (m) without multiplying cos(yt)
     real, dimension(:,:), pointer :: dytn     =>NULL() ! latitudinal width between grid points at j+1 and j in T-cells (m) 
     real, dimension(:,:), pointer :: rdxt     =>NULL() ! 1/dxt
     real, dimension(:,:), pointer :: rdyt     =>NULL() ! 1/dyt
     real, dimension(:,:), pointer :: rdxc     =>NULL() ! 1/dxc
     real, dimension(:,:), pointer :: rdyc     =>NULL() ! 1/dyc
     real, dimension(:,:), pointer :: rcosdxt  =>NULL() ! 1/(dxt*cost)
     real, dimension(:,:), pointer :: rcosdyt  =>NULL() ! 1/(dyt*cost)
     real, dimension(:,:), pointer :: rcosdxte =>NULL() ! 1/(dxte*cost)
     real, dimension(:,:), pointer :: rdytn    =>NULL() ! 1/(dytn)
     real, dimension(:,:), pointer :: area     =>NULL() ! area of T-cell
     real, dimension(:),   pointer :: f        =>NULL() ! rotation rate  
     real, dimension(:),   pointer :: rc8   =>NULL()    ! inverse of ( 8*cosine latitude )
     real, dimension(:,:), pointer :: ht    =>NULL()    ! depth at T-cell
     real, dimension(:,:), pointer :: tmask =>NULL()    ! land/sea mask for Tracer
     real, dimension(:,:), pointer :: umask =>NULL() ! land/sea mask for U (west face center)
     real, dimension(:,:), pointer :: vmask =>NULL() ! land/sea mask for V (south face center)
     real, dimension(:,:), pointer :: cmask =>NULL() ! land/sea mask for Z (southwest corner)
  end type ocean_grid_type
!</PUBLICTYPE>

!<PUBLICTYPE >
  type ocean_time_type
     type(time_type) :: Time_init         ! initial time 
     type(time_type) :: Model_time        ! model time
     type(time_type) :: Time_step         ! time step
     real            :: dt                ! time step in seconds
     integer         :: itt               ! timestep counter
     integer         :: taum1, tau, taup1 ! time level indices
  end type ocean_time_type
!</PUBLICTYPE>

!<PUBLICTYPE >
  type ocean_velocity_type
     real, dimension(:,:,:), pointer :: u => NULL() ! zonal velocity
     real, dimension(:,:,:), pointer :: v => NULL() ! meridinal velocity
  end type ocean_velocity_type
!</PUBLICTYPE>

!<PUBLICTYPE >
  type ocean_freesurf_type
     real, dimension(:,:,:), pointer :: eta=>NULL() ! surface height on tracer cell center (m) 
     real, dimension(:,:),   pointer :: ud =>NULL() ! vertically integrated zonal velocity (m^2/s)
     real, dimension(:,:),   pointer :: vd =>NULL() ! vertically integrated meridinal velocity (m^2/s)
  end type ocean_freesurf_type
!</PUBLICTYPE>


end module ocean_types_mod
