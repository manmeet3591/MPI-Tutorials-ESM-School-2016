module ocean_grid_mod
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
  ! Set up the ocean model grid spacing 
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  ! This module sets up the ocean model grid based on information read in 
  ! from the grid_spec.nc file (when read_my_grid is true) or namelist 
  ! setup(when read_my_grid is false). 
  !</DESCRIPTION>

  use mpp_mod,         only: mpp_error, FATAL, stdout, stdlog
  use mpp_io_mod,      only: mpp_open, mpp_close, mpp_read, mpp_get_atts
  use mpp_io_mod,      only: mpp_get_axes, mpp_get_info, mpp_get_axis_data
  use mpp_io_mod,      only: mpp_get_fields, axistype, fieldtype, atttype 
  use mpp_io_mod,      only: MPP_RDONLY,MPP_NETCDF,MPP_MULTI,MPP_SINGLE
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only: mpp_update_domains, domain2D 
  use mpp_domains_mod, only: mpp_global_sum, BITWISE_EXACT_SUM 
  use fms_mod,         only: open_namelist_file, check_nml_error, close_file
  use fms_mod,         only: write_version_number, read_data
  use diag_manager_mod,only: diag_axis_init
  use axis_utils_mod,  only: nearest_index  
  use constants_mod,   only: PI, RADIUS, EPSLN, OMEGA
  use ocean_types_mod, only: ocean_grid_type

  implicit none
  private

  !--- public interface ----
  public :: ocean_grid_init, set_ocean_grid_size, set_ocean_grid, ocean_grid_end

  !--- namelist interface ----
  !<NAMELIST NAME="ocean_grid_nml">  
  ! <DATA NAME="read_my_grid" TYPE="logical" DEFAULT="FALSE">
  !   when true, read grid and topography from grid file INPUT/grid_spec.nc. 
  !   when false, the grid and topography are specified by nml nlon, nlat,
  !   slon, slat, dxdeg, dydeg, cyclic and mean_depth.  
  ! </DATA> 
  ! <DATA NAME="nlon" TYPE="integer" DEFAULT="100">
  !  number of longitude grid points. 
  !  This value will be used when read_my_grid is false. 
  ! </DATA>
  ! <DATA NAME="nlat" TYPE="integer" DEFAULT="80">
  !  number of latitude grid points. 
  !  This value will be used when read_my_grid is false.
  ! </DATA>
  ! <DATA NAME="slon" TYPE="real" DEFAULT="0.0">
  !   starting (westernmost) longitude of "eta" grid point (degrees)
  ! </DATA>
  ! <DATA NAME="slat" TYPE="real" DEFAULT="-20.25">
  !   starting (southernmost) latitude of "eta" grid point (degrees)
  ! </DATA>  
  ! <DATA NAME="dxdeg" TYPE="real" DEFAULT="1.0">
  !   longitudinal width of a grid box (degrees)
  ! </DATA>
  ! <DATA NAME="dydeg" TYPE="real" DEFAULT="0.5">
  !   latitudinal width of a grid box (degrees)
  ! </DATA>
  ! <DATA NAME="cyclic" TYPE="logical" DEFAULT="TRUE">
  !  If true, cyclic boundary condition is applied along zonal direction. 
  !  If false, the domain is bounded by a solid wall.
  ! </DATA> 
  ! <DATA NAME="mean_depth" TYPE="real" DEFAULT="0.2">
  !   mean ocean depth (meter) disregarding any topography. It can be set to
  !   the actual mean ocean depth or the equivalent depth of a single 
  !   baroclinic mode (default value is for the 2nd baroclinic mode). 
  ! </DATA> 
  !</NAMELIST>
  logical           :: read_my_grid     = .FALSE.
  integer           :: nlon             = 100
  integer           :: nlat             = 80
  real              :: slon             = 0.0
  real              :: slat             = -20.25
  real              :: dxdeg            = 1.0
  real              :: dydeg            = 0.5
  logical           :: cyclic           = .true.
  real              :: mean_depth             = 0.2

  namelist /ocean_grid_nml/ read_my_grid, nlon, nlat, slon, slat, dxdeg, dydeg, &
                            cyclic, mean_depth

  integer :: isc, iec, jsc, jec  ! compute domain decomposition
  integer :: isd, ied, jsd, jed  ! data domain decomposition
  integer :: ni, nj              ! grid size
  real    :: D2R, D2M  
  logical :: is_new_grid
  logical :: module_is_initialized = .false.
  character(len=32) :: grid_file = "INPUT/grid_spec.nc"

  !--- version information variables ----
  character(len=128) :: version = '$Id: ocean_grid.f90,v 1.1.2.2.2.4 2005/03/01 15:52:41 z1l Exp $'
  character(len=128) :: tagname = '$Name:  $'

contains


  !#######################################################################
  ! <SUBROUTINE NAME="ocean_grids_init">
  !
  ! <DESCRIPTION>
  ! Initialize the grids module. 
  ! </DESCRIPTION>
  !
  subroutine ocean_grid_init

    integer :: ioun, io_status, ierr

    if ( module_is_initialized ) then
       call mpp_error(FATAL, '==>Error from ocean_grid_mod (ocean_grid_init): module has been initialized')
    endif

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

    ! provide for namelist over-ride of defaults 
    ioun = open_namelist_file()
    read (ioun,ocean_grid_nml,IOSTAT=io_status)
    write (stdout(),'(/)')
    write (stdout(),ocean_grid_nml)  
    write (stdlog(),ocean_grid_nml)
    ierr = check_nml_error(io_status, 'ocean_grid_nml')
    call close_file(ioun)

    D2R  = PI/180.0
    D2M = PI/180.*RADIUS

  end subroutine ocean_grid_init
  ! </SUBROUTINE> NAME="ocean_grid_init"

  !#####################################################################
  ! <SUBROUTINE NAME="set_ocean_grid_size">
  !
  ! <DESCRIPTION>
  ! Set the ocean grid size.  when read_my_grid is true, will get grid size
  ! from grid file 'INPUT/grid_spec.nc'. otherwise grid size will be determined
  ! by the namelist variable nlon and nlat  
  ! </DESCRIPTION>
  subroutine set_ocean_grid_size(Grid)
    type(ocean_grid_type), intent(inout) :: Grid

    integer                              :: unit, ndim, nvar, natt, ntime, i, len
    type(axistype),          allocatable :: axes(:)
    type(atttype),           allocatable :: global_atts(:)
    character(len=64)                    :: name

    if(read_my_grid) then
       ni = 0; nj = 0
       call mpp_open(unit,trim(grid_file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
       call mpp_get_info(unit, ndim, nvar, natt, ntime)
       allocate(axes(ndim))
       call mpp_get_axes(unit,axes)
       !--- determine if it is new_grid or old grid
       do i=1, ndim
          call mpp_get_atts(axes(i),name=name)  
          select case(trim(name))
          case ('grid_x_T')
             is_new_grid = .true.
          case ('gridlon_t')
             is_new_grid = .false.
          end select
       enddo

       if(is_new_grid) then
          write(stdout(),*)'==>NOTE, running with is_new_grid=.true.'
       else
          write(stdout(),*)'==>NOTE, running with is_new_grid=.false.'
       endif

       do i=1, ndim
          call mpp_get_atts(axes(i),name=name,len=len)       
          if(is_new_grid) then
             select case(trim(name))
             case ('grid_y_T')
                nj = len
             case ('grid_x_T')
                ni = len
             end select
          else
             select case(trim(name))
             case ('gridlat_t')
                nj = len
             case ('gridlon_t')
                ni = len
             end select
          endif
       enddo

       if (ni == 0 .or. nj == 0 ) call mpp_error(FATAL,   &
            '==>Error reading grid information from grid file (read_my_grid=true).Are you sure file exists?')

       !--- get the boundary condition
       Grid%cyclic = .false.
       allocate(global_atts(natt))
       call mpp_get_atts(unit,global_atts)
       do i=1,natt
          select case (trim(global_atts(i)%name))
          case ('x_boundary_type')
             if (trim(global_atts(i)%catt) == 'cyclic') Grid%cyclic = .true.
          case ('y_boundary_type')
             if (trim(global_atts(i)%catt) == 'fold_north_edge') &
                  call mpp_error(FATAL,'==>Error from ocean_model_mod: The grid should be simple spheric grid.')
          end select
       end do
      call mpp_close(unit)
    else
       if(nlon .le. 0 .or. nlat .le. 0) call mpp_error(FATAL, '==>Error from ocean_model_mod: '&
            //'when read_my_grid is false, nml ni and nj should be set to positive')
       ni     = nlon
       nj     = nlat
       Grid%cyclic = cyclic
    endif

    if(Grid%cyclic) then
       write(stdout(),*)'==> NOTE: The shallow water model use cyclic conditions in zonal direction.'
    else
       write(stdout(),*)'==> NOTE: The shallow water model domain is bounded by solid walls.'
    endif

    Grid%ni = ni; Grid%nj = nj

  end subroutine set_ocean_grid_size
  ! </SUBROUTINE> NAME="set_ocean_grid_size"

  !#####################################################################
  ! <SUBROUTINE NAME="set_ocean_grid">
  !
  ! <DESCRIPTION>
  !  define ocean horizontal grid and topography. when read_my_grid is true,
  !  grid information will be read from file 'INPUT/grid_spec.nc'. otherwise,
  !  grid will be specified by namelist variable nlon, nlat, slon, slat, 
  !  dxdeg, dydeg, cyclic, and mean_depth
  ! </DESCRIPTION>
  subroutine set_ocean_grid(Domain, Grid)
    type(domain2D),        intent(inout) :: Domain
    type(ocean_grid_type), intent(inout) :: Grid

    call mpp_get_compute_domain(Domain,isc,iec,jsc,jec)        
    call mpp_get_data_domain   (Domain,isd,ied,jsd,jed)

    allocate(Grid%xt(0:ni+1), Grid%yt(0:nj+1), Grid%xc(0:ni+1), Grid%yc(0:nj+1) )
    allocate(Grid%f(0:nj+1), Grid%cosc(0:nj+1),  Grid%cost(0:nj+1), Grid%rc8(0:nj+1) )
    allocate(Grid%ht(isd:ied,jsd:jed),       Grid%tmask(isd:ied,jsd:jed) )
    allocate(Grid%umask(isd:ied,jsd:jed),    Grid%vmask(isd:ied,jsd:jed) )
    allocate(Grid%cmask(isd:ied,jsd:jed) )
    allocate(Grid%dxt(isd:ied,jsd:jed),      Grid%dyt(isd:ied,jsd:jed) )
    allocate(Grid%dxc(isd:ied,jsd:jed),      Grid%dyc(isd:ied,jsd:jed) )
    allocate(Grid%rdxt(isd:ied,jsd:jed),     Grid%rdyt(isd:ied,jsd:jed) )
    allocate(Grid%rdxc(isd:ied,jsd:jed),     Grid%rdyc(isd:ied,jsd:jed) )
    allocate(Grid%dxte(isd:ied,jsd:jed),     Grid%dytn(isd:ied,jsd:jed) )
    allocate(Grid%rcosdxt(isd:ied,jsd:jed),  Grid%rcosdyt(isd:ied,jsd:jed) )
    allocate(Grid%rcosdxte(isd:ied,jsd:jed), Grid%rdytn(isd:ied,jsd:jed) )
    allocate(Grid%area(isd:ied,jsd:jed))
    Grid%ht    = 0.0
    Grid%tmask = 0.0
    Grid%umask = 0.0
    Grid%vmask = 0.0
    Grid%cmask = 0.0
    !--- set up grid, topography 
    if(read_my_grid) then
       call read_grid(Domain,Grid)
    else
       call define_grid(Domain,Grid)
    endif

    !--- write out grid informaiton 
    write(stdout(), *)'Grid%xt(i),i=1,ni', Grid%xt(1:ni)
    write(stdout(), *)'Grid%yt(j),j=1,nj', Grid%yt(1:nj)
    write(stdout(), *)'Grid%xc(i),i=1,ni', Grid%xc(1:ni)
    write(stdout(), *)'Grid%yc(j),j=1,nj', Grid%yc(1:nj)

    !--- define grid other varialbe, such as land/sea mask, grid coefficient.
    call define_grid_other_var(Domain, Grid)

    ! --- set up diag_axis for diagnostic output
    Grid%id_xt = diag_axis_init ('xt',Grid%xt(1:Grid%ni),'degrees_E','x','T-cell longitude', Domain2=Domain)
    Grid%id_yt = diag_axis_init ('yt',Grid%yt(1:Grid%nj),'degrees_N','y','T-cell latitude', Domain2=Domain)
    Grid%id_xc = diag_axis_init ('xc',Grid%xc(1:Grid%ni),'degrees_E','x','C-cell longitude', Domain2=Domain)
    Grid%id_yc = diag_axis_init ('yc',Grid%yc(1:Grid%nj),'degrees_N','y','C-cell latitude', Domain2=Domain)


  end subroutine set_ocean_grid
  ! </SUBROUTINE> NAME="set_ocean_grid"

  !###################################################################
  ! <SUBROUTINE NAME="ocean_grid_end">
  !
  ! <DESCRIPTION>
  !  destructor routine: release memory.
  ! </DESCRIPTION>
  subroutine ocean_grid_end(Grid)
    type(ocean_grid_type), intent(inout) :: Grid

    deallocate(Grid%xt, Grid%yt, Grid%xc, Grid%yc, Grid%f, Grid%cosc, Grid%cost )
    deallocate(Grid%dxt, Grid%dyt, Grid%dxc, Grid%dyc, Grid%dxte, Grid%dytn )
    deallocate(Grid%rdxt, Grid%rdyt, Grid%rdxc, Grid%rdyc, Grid%rcosdxte, Grid%rdytn )
    deallocate(Grid%rcosdxt, Grid%rcosdyt, Grid%area)
    deallocate(Grid%rc8, Grid%ht )
    deallocate(Grid%tmask, Grid%umask, Grid%vmask, Grid%cmask )

    module_is_initialized = .false.

  end subroutine ocean_grid_end
  ! </SUBROUTINE> NAME="ocean_grid_end"

  !#####################################################################
  !--- define longitude and latitude, topography and grid mask.
  !--- we can read these information from a grid spec file if needed.
  subroutine define_grid(Domain, Grid)
    type(domain2D),        intent(inout) :: Domain
    type(ocean_grid_type), intent(inout) :: Grid

    integer                              :: i, j
    integer                              :: island, ieland, jsland, jeland
    integer                              :: istopo, ietopo, jstopo, jetopo
    real                                 :: xtopo, ytopo, xisle, yisle
    real                                 :: hheight, hwidth, hwidt, hheightt
    real                                 :: degecm, dx, dy

    !--- define grid
    Grid%xt(0) = slon - dxdeg
    Grid%yt(0) = slat - dydeg
    Grid%xc(0) = Grid%xt(0) - 0.5*dxdeg
    Grid%yc(0) = Grid%yt(0) - 0.5*dydeg

    do i=1, ni+1
       Grid%xt(i) = Grid%xt(i-1) + dxdeg
       Grid%xc(i) = Grid%xc(i-1) + dxdeg
    enddo

    do j=1, nj+1
       Grid%yt(j) = Grid%yt(j-1) + dydeg
       Grid%yc(j) = Grid%yc(j-1) + dydeg
    enddo

    do j = jsc, jec
       do i = isc, iec
          Grid%ht(i,j) = mean_depth
       enddo
    enddo

    call mpp_update_domains(Grid%ht, Domain)

  end subroutine define_grid

  !#####################################################################
  subroutine read_grid(Domain,Grid)
    type(domain2D),        intent(inout) :: Domain
    type(ocean_grid_type), intent(inout) :: Grid

    real,            allocatable         :: tmp(:,:,:)
    real                                 :: small = 1.0e-10
    real                                 :: dx, dy
    integer                              :: i,j

    if(is_new_grid) then
       allocate(tmp(ni,nj,4))
       call read_data(grid_file, 'grid_x_T', Grid%xt(1:ni) )
       call read_data(grid_file, 'grid_y_T', Grid%yt(1:nj) )
       call read_data(grid_file, 'x_vert_T', tmp)
       Grid%xc(1:ni) = tmp(1:ni,1,1)
       Grid%xc(ni+1) = tmp(ni,1,2)
       call read_data(grid_file, 'y_vert_T', tmp)
       Grid%yc(1:nj) = tmp(1,1:nj,1)
       Grid%yc(nj+1) = tmp(1,nj,4)
       deallocate(tmp)
    else
       call read_data(grid_file, 'gridlon_t', Grid%xt(1:ni) )
       call read_data(grid_file, 'gridlat_t', Grid%xt(1:nj) )
       call read_data(grid_file, 'gridlon_vert_t', Grid%xc(1:ni+1) )
       call read_data(grid_file, 'gridlat_vert_t', Grid%yc(1:nj+1) )
    endif

    Grid%xt(0)    = 2*Grid%xc(1)  - Grid%xt(1)
    Grid%yt(0)    = 2*Grid%yc(1)  - Grid%yt(1)
    Grid%xc(0)    = 2*Grid%xt(0)  - Grid%xc(1)
    Grid%yc(0)    = 2*Grid%yt(0)  - Grid%yc(1)
    Grid%xt(ni+1) = 2*Grid%xc(ni+1)  - Grid%xt(ni)
    Grid%yt(nj+1) = 2*Grid%yc(nj+1)  - Grid%yt(nj)

    !--- right now the model assume dx and dy are uniform. We may remove
    !--- this assumption in the future.
    dx = Grid%xt(1)-Grid%xt(0)
    do i = 1, ni
       if( abs(Grid%xt(i+1) - Grid%xt(i) - dx) > small ) then
           write(stdout(),*) "Grid%xt(0:ni+1) is ", Grid%xt
           call mpp_error(FATAL,"ocean_grid: the longitude is not uniform.")
       endif
       if( abs(Grid%xc(i+1) - Grid%xc(i) - dx) > small ) then
           write(stdout(),*) "Grid%xc(0:ni+1) is ", Grid%xc
           call mpp_error(FATAL,"ocean_grid: the longitude is not uniform.")
       endif
    enddo
   
    dy = Grid%yt(1)-Grid%yt(0)
    do j = 1, nj
       if( abs(Grid%yt(j+1) - Grid%yt(j) - dy) > small ) then
           write(stdout(),*) "Grid%yt(0:nj+1) is ", Grid%yt
           call mpp_error(FATAL,"ocean_grid: the latitude is not uniform.")
       endif
       if( abs(Grid%yc(j+1) - Grid%yc(j) - dy) > small ) then
           write(stdout(),*) "Grid%yc(0:nj+1) is ", Grid%yc
           call mpp_error(FATAL,"ocean_grid: the latitude is not uniform.")
       endif
    enddo

    if(is_new_grid) then
       call read_data(grid_file, 'depth_t', Grid%ht(isc:iec,jsc:jec), Domain)
    else
       call read_data(grid_file, 'ht', Grid%ht(isc:iec,jsc:jec), Domain)
    endif

    call mpp_update_domains(Grid%ht, Domain)

  end subroutine read_grid

  !#####################################################################
  subroutine define_grid_other_var(Domain, Grid)
    type(domain2D),        intent(inout) :: Domain
    type(ocean_grid_type), intent(inout) :: Grid

    integer                              :: i, j
    real, allocatable                    :: tmp(:,:)

    !--- define some coefficents related to grid
    do j = 0, nj+1
       Grid%cosc(j)    = cos(Grid%yc(j)*D2R)
       Grid%cost(j)    = cos(Grid%yt(j)*D2R)
       Grid%f(j)       = 2.0*omega*sin(Grid%yc(j)*D2R)
       Grid%rc8(j)     = 0.125/Grid%cost(j)
    enddo

    !--- define grid length ( in meters )
    Grid%dxt = 0.0; Grid%dyt = 0.0; Grid%dxc = 0.0; Grid%dyc = 0.0
    Grid%dxte = 0.0; Grid%dytn = 0.0

    do j = jsc, jec
       do i = isc, iec    
          Grid%dxt(i,j)  = ( Grid%xc(i+1)-Grid%xc(i) )*D2M
          Grid%dxc(i,j)  = ( Grid%xt(i)-Grid%xt(i-1) )*D2M
          Grid%dyt(i,j)  = ( Grid%yc(j+1)-Grid%yc(j) )*D2M
          Grid%dyc(i,j)  = ( Grid%yt(j)-Grid%yt(j-1) )*D2M
          Grid%dxte(i,j) = ( Grid%xt(i+1)-Grid%xt(i) )*D2M
          Grid%dytn(i,j) = ( Grid%yt(j+1)-Grid%yt(j) )*D2M
       enddo
    enddo

    call mpp_update_domains(Grid%dxt, Domain)
    call mpp_update_domains(Grid%dyt, Domain)
    call mpp_update_domains(Grid%dxc, Domain)
    call mpp_update_domains(Grid%dyc, Domain)
    call mpp_update_domains(Grid%dxte, Domain)
    call mpp_update_domains(Grid%dytn, Domain)
          

    do j = jsd, jed
       do i = isd, ied
          Grid%rdxt(i,j)     = 1./(Grid%dxt(i,j)+epsln)
          Grid%rdyt(i,j)     = 1./(Grid%dyt(i,j)+epsln)
          Grid%rdxc(i,j)     = 1./(Grid%dxc(i,j)+epsln)
          Grid%rdyc(i,j)     = 1./(Grid%dyc(i,j)+epsln)
          Grid%rcosdxt(i,j)  = 1./(Grid%dxt(i,j)*Grid%cost(j)+epsln)
          Grid%rcosdyt(i,j)  = 1./(Grid%dyt(i,j)*Grid%cost(j)+epsln)
          Grid%rcosdxte(i,j) = 1./(Grid%dxte(i,j)*Grid%cost(j)+epsln)
          Grid%rdytn(i,j)    = 1./(Grid%dytn(i,j)+epsln)
          Grid%area(i,j)     = Grid%dxt(i,j)*Grid%cost(j)*Grid%dyt(i,j) 
       enddo
    enddo

    !--- calculate mask( 1 means ocean ) 
    do j=jsc, jec
       do i=isc, iec
          if (Grid%ht(i,j) .gt. 0.0) Grid%tmask(i,j) = 1.0
       enddo
    enddo

    call mpp_update_domains(Grid%tmask, Domain)

    do j = jsc, jec
       do i = isc, iec
          Grid%umask(i,j) = min(Grid%tmask(i,j), Grid%tmask(i-1,j) )
          Grid%vmask(i,j) = min(Grid%tmask(i,j), Grid%tmask(i,j-1) )
          Grid%cmask(i,j) = min(Grid%tmask(i,j), Grid%tmask(i-1,j), &
                                Grid%tmask(i,j-1), Grid%tmask(i-1,j-1) )
       enddo
    enddo
    call mpp_update_domains(Grid%umask, Domain)
    call mpp_update_domains(Grid%vmask, Domain)

    !--- calculate total ocean area for integrals.
    allocate(tmp(isc:iec,jsc:jec) )
    tmp = 0.0
    do j = jsc, jec
       do i = isc, iec
          if(Grid%ht(i,j) > 0 ) tmp(i,j) = Grid%dxt(i,j)*Grid%dyt(i,j)*Grid%cost(j)
       enddo
    enddo
    Grid%totarea = mpp_global_sum(Domain,tmp(:,:), BITWISE_EXACT_SUM)
    write(stdout(),*)' total area of the region is ', Grid%totarea
    deallocate(tmp)

  end subroutine define_grid_other_var

  !#######################################################################

end module ocean_grid_mod
