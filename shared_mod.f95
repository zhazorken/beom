! Copyright 2016 Pierre St-Laurent.
! pierrestlau@globetrotter.net

! The parameters defining the calculation are to be set in this file.
! This module also defines public procedures and variables
!   such as the precision of floats and the length of strings.

! This file is part of beom.

! beom is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! beom is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with beom. If not, see <http://www.gnu.org/licenses/>.

module shared_mod
  implicit none
  private

! Parameter rw defines the `Working precision' for Real numbers.
  integer, parameter, public ::           &
    i4   = selected_int_kind ( R =   9 ), &
    r4   = selected_real_kind( P =   6 ), &
    r8   = selected_real_kind( P =  12 ), &
    rw   = r8,                            &
!   Length of Long  STRings.
    lstr = 999,                           &
!   Length of Short STRings.
    sstr =  99

!<=============BEGINNING OF USER-MODIFIABLE SECTION============================>
! i.e., parameters that must be set for a given calculation:

  integer,     parameter, public :: &
    lm         =  125,              & ! Grid size in x (zonal) direction.
    mm         =  501,              & ! Grid size in y (merid) direction.
    nlay       =    2,              & ! Number of layers.
    ndeg       =  63252                ! Grid size after vectorization.
  real( rw ),  parameter, public :: &
    dl         = 400,             & ! Mesh size (meters).
    cext       =  82.8,              & ! c_{external} = sqrt(g H_{max}) (m/s).
    f0         = 1.409e-4,             & ! Coriolis parameter (rad/s).
    rhon(nlay) = (/1027.47,1027.75/),         & ! Neutral density of layers (kg/m3).
    topl(nlay) = (/0.00,0.2143/),            & ! Position of layer interfaces (non-dim).
    dt_s       =  30.,             & ! Duration of Simulation           (days).
    dt_o       = .01,                & ! Time between model Outputs       (days).
    dt_r       = 0.,                & ! Time over which winds are Ramped (days).
    dt3d       = 0.,                & ! Time between updates of stress   (days).
    bvis       = 0.,                & ! Background viscosity (m2/s).
    dvis       = 0.9,              & ! Dyn. Visc. (0.< dvis < 1., non-dim).
    bdrg       = 0.,            & ! Bottom Drag (non-dim quadr; m/s linear).
    hmin       = 0.3,          & ! Minimum layer thickness allowed   (m).
    hsbl       =  5.,               & ! Thickn. surf bound lay if outcrop (m).
    hbbl       =  5.,               & ! Thickn. bott bound lay if outcrop (m).
    g_fb       = 1.,                & ! 1. = Use Generalized Forward Backward.
    uadv       = 1.,                & ! 1. = Use momentum advection.
    qdrg       = 1.,                & ! 1. = Use quadratic bottom drag.
    ocrp       = 1.,                & ! 1. = Allow isopycnal outcrops.
    rsta       = 0.,                & ! 1. = Restart from previous calc.
    xper       = 0.,                & ! 1. = Domain periodic along x.
    yper       = 0.,                & ! 1. = Domain periodic along y.
    diag       = 0.,                & ! 1. = Output diag fields.
    rgld       = 0.,                & ! 1. = Rigid lid 
    mcbc       = 1.                   ! 1. = Modified closed boundary conditions with 'eta' nudging.
  complex(rw), parameter, public :: &
    tauw       = (0.00, 0.0)         ! Wind stress taux,tauy (Pascals).
  character(len = sstr), public ::  &
    idir       = '/home/zhazorken/Desktop/BEOM/runs/h350day30cbc/',           & ! Path to input  files.
    odir       = '/home/zhazorken/Desktop/BEOM/runs/h350day30cbc/',           & ! Path to output files.
    desc       = 'Test-case: 3D sill exchange'

!<=============END OF USER-MODIFIABLE SECTION==================================>

! Other constants used in the program (not to be modified).

  real( rw ), parameter, public ::       &
    dt   = 0.5_rw * dl / cext,           & ! Model timestep (seconds).
    hsal = 10._rw * hmin,                & ! Salmon thickness if outcrop   (m).
    hdry = 1.e-3,                        & ! Depth threshold for dry areas (m).
    tole = 1.e-6,                        & ! Tolerance for h_0 if outcrop  (m).
    pi   = 3.1415927,                    &
    grav = 9.8,                          & ! (m/s2)
    rho0 = rhon( nlay ),                 & ! Density reference (kg/m3).
    beta = 0.281105,                     & ! Coefficients for Generalized FB,
    epsi = 0.013,                        & !   taken from Shchepetkin &
    gamm = 0.088,                        & !   McWilliams (2005),
    del1 = 0.5_rw + gamm + 2._rw * epsi, & !   Ocean Modelling,
    del2 = 1._rw  - del1 - gamm - epsi     !   page 374.

  real( r8 ), parameter, public ::       &
!   Successive-Over-Relaxation during calculation of layers' thickness.
    sor  = 1.9                             ! 1. <= sor <= 2.

  integer,    parameter, public ::       &
!   itmx is involved in the calculation of the layers' thickness.
!   More iterations are required when the number of layers (nlay) increases.
    itmx = 99999,                        & ! Max. number of iterations allowed.
    nsal = 4,                            & ! Salmon's exponent for outcrop.
    ix_n = 1,                            & ! Index to surface elevation.
    ix_u = 2,                            & !          zonal   velocity.
    ix_v = 3,                            & !          merid.  velocity.
    iosi = 5,                            & ! Unit for Standard input.
    ioso = 6,                            & ! Unit for Standard output.
    iose = 0                               ! Unit for Standard error.

  integer,           public :: errc        ! Error code.
  character( lstr ), public :: errm        ! Error message.

  public quit, get_un                      ! Public procedures.

contains

subroutine quit()   ! Ensure a graceful exit.
  implicit none
  integer           :: iunit
  logical           :: this_unit_is_opened
  character( lstr ) :: action_of_unit

  if ( errc /= 0 ) then ! An exception flag was raised.

    write( iose, * ) '  *** ERROR CODE = ', errc, ' ***  '
    write( iose, * ) trim( errm )

!   Make sure that all messages were sent to the user.
!   CAREFUL, `flush' is not available in strict Fortran 95.
!   flush( iose )
  end if

! Ensure that all data were written to disk
!   (for debugging purposes, in case an exception was raised).

  do iunit = 99, 7, -1
    inquire(   unit = iunit, opened = this_unit_is_opened )

    if ( this_unit_is_opened ) then
      inquire( unit = iunit, action = action_of_unit )

      if ( trim( action_of_unit ) == 'WRITE' .or. &
           trim( action_of_unit ) == 'write' .or. &
           trim( action_of_unit ) == 'Write' ) then
!       CAREFUL, `flush' is not available in strict Fortran 95.
!       flush( iunit )
      end if

      close( iunit )
    end if
  end do

  stop
end subroutine quit

function get_un() ! Get a unit number. Raise exception if none available.
  implicit none
! Units (0,5,6) are usually reserved on UNIX-like systems:
! 0 : stderr (standard error )
! 5 : stdin  (standard input )
! 6 : stdout (standard output)
  integer, parameter :: base_unit =  7, &
                        max_unit  = 99
  integer            :: unum, get_un, lerm
  logical            :: l1

  get_un = - 999 ! Initialization.

  if ( errc /= 0 ) then
    return
  else
    lerm = len_trim( errm )
    errm = trim(errm) // ' inside function get_un from module shared_mod.f95,'
  end if

  unum = base_unit
  do
    inquire( unit = unum, opened = l1 )
    if ( l1 ) then
      unum = unum + 1
    else
      get_un = unum
      errm   = errm(1 : lerm)
      exit
    end if
    if ( unum > max_unit ) then
      errc = unum
      errm = errm(1 : len_trim(errm)) &
           // ' no more free I/O unit numbers available.'
      call quit()
    end if
  end do
end function get_un

end module shared_mod
