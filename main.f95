! Copyright 2016 Pierre St-Laurent.
! pierrestlau@globetrotter.net

! Main program for `beom'. The coupling between the model components
! (ocean and e.g. SWAN or CICE) is to be defined in this file.
! If there is no coupling, the main program just calls/runs the ocean model
! and then comes to a stop.
! The ocean model itself is in 2 modules: shared_mod.f95 and private_mod.f95.
! The parameters defining the ocean model are to be set in shared_mod.f95.

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

program main
  use shared_mod,  only: errc, errm, quit
  use private_mod, only: run
  implicit none

  errc = 0              ! Initialize error code.
  errm = 'In main.f95,' ! Initialize error message.

  call run ()           ! Initialize and run the ocean model.
  call quit()           ! Ensure graceful exit if exception is raised.

end program main
