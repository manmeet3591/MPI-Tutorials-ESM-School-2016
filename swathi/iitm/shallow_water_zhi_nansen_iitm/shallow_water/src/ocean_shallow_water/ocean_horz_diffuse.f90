module ocean_horz_diffuse_mod
!-----------------------------------------------------------------------
!                   GNU General Public License                        
! This file is a part of MOM.                                                                 
!                                                                      
! MOM is free software; you can redistribute it and/or modify it and  
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
! time tendency from horizontal laplacian diffusion
!</OVERVIEW>

!<NAMELIST NAME="ocean_horz_diffuse_nml">
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM diffusivity. 
!  </DATA> 
!</NAMELIST>
!

use fms_mod,         only: stdout, stdlog, check_nml_error, write_version_number
use fms_mod,         only: open_namelist_file, close_file
use mpp_domains_mod, only: domain2D, mpp_update_domains
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
use ocean_types_mod, only: ocean_grid_type

implicit none
private

public :: ocean_horz_diffuse_init, ocean_horz_diffuse

!--- namelist interface
real    :: vel_micom       = 1.0e-2    ! constant velocity scale (m/s) for setting micom diffusivity  

namelist /ocean_horz_diffuse_nml/ vel_micom


real, dimension(:,:),      allocatable :: diff_ce ! horizontal diffusivity for eastern face  cell (m^2/s)
real, dimension(:,:),      allocatable :: diff_cn ! horizontal diffusivity for northern face cell (m^2/s)
real, dimension(:,:),      allocatable :: fx, fy
type(domain2D),                pointer :: Domain => NULL()
type(ocean_grid_type),         pointer :: Grid   => NULL()

integer                                :: isc, iec, jsc, jec
integer                                :: isd, ied, jsd, jed

  !--- version information variables ----
  character(len=128) :: version = '$Id: ocean_horz_diffuse.f90,v 1.1.2.2 2005/02/23 19:33:20 z1l Exp $'
  character(len=128) :: tagname = '$Name:  $'

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_horz_diffuse_init">
!
! <DESCRIPTION>
! Initialize the horizontal laplacian diffusion module.
! </DESCRIPTION>
!
  subroutine ocean_horz_diffuse_init(Grid_in, Domain_in)

    type(ocean_grid_type), target, intent(in)  :: Grid_in
    type(domain2D),        target, intent(in)  :: Domain_in

    integer :: ioun, io_status, i, j, k, n, num, ierr
    real    :: tmp

    ! provide for namelist over-ride of defaults 
    ioun  = open_namelist_file()
    read (ioun,ocean_horz_diffuse_nml,IOSTAT=io_status)
    write (stdout(),'(/)')
    write (stdout(),ocean_horz_diffuse_nml)  
    write (stdlog(),ocean_horz_diffuse_nml)
    ierr = check_nml_error(io_status, 'ocean_horz_diffuse_nml')
    call close_file(ioun)

    Domain => Domain_in
    Grid   => Grid_in

    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain, isd, ied, jsd, jed)

    call write_version_number( version, tagname )

    if(vel_micom > 0 ) then 
          write (stdout(),'(/1x,a)') ' ==> Note: USING horizontal Laplacian diffusion'
    endif

    allocate (diff_ce(isd:ied,jsd:jed), diff_cn(isd:ied,jsd:jed))
    allocate (fx(isd:ied,jsd:jed))
    allocate (fy(isd:ied,jsd:jed))  

    fx = 0.0
    fy = 0.0

    ! Micom background diffusivity Space scale set by grid size
    ! Velocity scale input via namelist 
    diff_ce = 0.0;
    diff_cn = 0.0;

    do j = jsc, jec
       do i = isc, iec
          tmp = vel_micom*2.0*Grid%dxt(i,j)*Grid%dyt(i,j)/(Grid%dxt(i,j)+Grid%dyt(i,j))
          diff_ce(i,j) = tmp
          diff_cn(i,j) = tmp
       enddo
    enddo

    call mpp_update_domains(diff_ce, Domain)
    call mpp_update_domains(diff_cn, Domain)

  end subroutine ocean_horz_diffuse_init
! </SUBROUTINE>  NAME="ocean_horz_diffuse_init"

!#######################################################################
! <FUNCTION NAME="ocean_horz_diffuse">
!
! <DESCRIPTION>
! This function computes the thickness weighted time tendency for tracer 
! from horizontal laplacian diffusion. 
! </DESCRIPTION>
!
  function ocean_horz_diffuse(fld)
    real, dimension(isd:ied,jsd:jed), intent(in) :: fld

    real, dimension(isc:iec,jsc:jec) :: ocean_horz_diffuse 
    integer                          :: i, j
    real                             :: diff_x, diff_y

    ! fx = flux component through "eastern"  face of T-cells at level k
    ! fy = flux component through "northern" face of T-cells at level k

    !--- forward scheme
    do j = jsd, jec
       do i = isd, iec
          fx(i,j) = diff_ce(i,j)*(fld(i+1,j) - fld(i,j))*Grid%rcosdxte(i,j) * Grid%tmask(i+1,j)*Grid%tmask(i,j) 
          fy(i,j) = diff_cn(i,j)*(fld(i,j+1) - fld(i,j))*Grid%rdytn(i,j) * Grid%tmask(i,j+1)*Grid%tmask(i,j) 
       enddo
    enddo

!    call mpp_update_domains(fx, Domain)
!    call mpp_update_domains(fy, Domain)

    !--- backward scheme     
    do j = jsc, jec
       do i = isc, iec
          diff_x = (fx(i,j) - fx(i-1,j)) * Grid%rcosdxt(i,j)
          diff_y = (fy(i,j) - fy(i,j-1)) * Grid%rdyt(i,j)
          ocean_horz_diffuse(i,j) = Grid%tmask(i,j) * (diff_x + diff_y)
       enddo
    enddo

    return

  end function ocean_horz_diffuse
! </FUNCTION> NAME="ocean_horz_diffusion"

  !#####################################################################


end module ocean_horz_diffuse_mod
      
