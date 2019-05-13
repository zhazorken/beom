! Copyright 2016 Pierre St-Laurent.
! pierrestlau@globetrotter.net

! Module containing the procedures and arrays involved in a calculation.
! This module is entirely private except for subroutine run().

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

module private_mod
  use shared_mod
  implicit none
  private

  integer ::              &
!   Connectivity table, providing the vector index of neighbors.
!   1- Right, 2- Upper right, 3- Up,    4- Upper left,
!   5- Left,  6- Lower left,  7- Lower, 8- Lower right.
    neig( 8, 0 : ndeg ),  &
!   (i,j) Sub-indices corresponding to vectorized grid Cells [0, ndeg].
    subc( 0 : ndeg, 2 )

  real( rw ) ::           &
    h_u (0 : ndeg, nlay), & ! h*u                 (m**2 s**(-1)).
    h_v (0 : ndeg, nlay), & ! h*v                 (m**2 s**(-1)).
    u   (0 : ndeg, nlay), & ! Zonal velocity      (m    s**(-1)).
    v   (0 : ndeg, nlay), & ! Meridional velocity (m    s**(-1)).
    UU4 (0 : ndeg, nlay), & ! Biharmonic viscosity term
    VV4 (0 : ndeg, nlay), & ! Biharmonic viscosity term 
    delu(0 : ndeg, nlay), & ! Laplacian of u      (m**(-1) s**(-1))
    delv(0 : ndeg, nlay), & ! Laplacian of v      (m**(-1) s**(-1)) 
    tt3d(0 : ndeg, 2, nlay),& ! Surface stress forcing (Pascals), cell-centered.
    tb3d(0 : ndeg, 2, nlay),& ! Bottom  stress         (Pascals), at u/v points.
    tu3d(0 : ndeg, 2, nlay),& ! Top  stress         (Pascals), at u/v points.
    taus(0 : ndeg, 2   ), & ! Surface stress forcing (Pascals), cell-centered.
    rvor(0 : ndeg      ), & ! Relative  vorticity   (        s**(-1)).
    pvor(0 : ndeg      ), & ! Potential vorticity   (m**(-1) s**(-1)).
    dive(0 : ndeg      ), & ! Horizontal divergence (        s**(-1)).
    v_cc(0 : ndeg, nlay), & ! Viscosity `nu' at cell-center (m**2 s**(-1)).
    v_ll(0 : ndeg, nlay), & ! Viscosity `nu' at lower-left  (m**2 s**(-1)).
    fcor(0 : ndeg      ), & ! Coriolis parameter (rad s**(-1)).
    mk_u(0 : ndeg      ), & ! Mask for `u'   points (0. = no zonal      flux).
    mk_v(0 : ndeg      ), & ! Mask for `v'   points (0. = no meridional flux).
    mk_n(0 : ndeg      ), & ! Mask for `eta' points (0. = dry cell          ).
    mkpe(0 : ndeg      ), & ! Mask for `psi' points *excluding coastline*.
    mkpi(0 : ndeg      ), & ! Mask for `psi' points *including coastline*.
    mont(0 : ndeg      ), & ! Montgomery potential / grav (m).
    d2hx(0 : ndeg      ), & ! Correction term for upstream-biased scheme.
    d2hy(0 : ndeg      ), & ! Correction term for upstream-biased scheme.
    h_bo(0 : ndeg      ), &    ! Bottom topography (meters, >0).
    h_to(0 : ndeg      ), &   ! Surface topography (meters, >0).
    h_th(0 : ndeg      ), &   ! Water column thickness (meters, >0).
    Ow(  0 : ndeg      ) = 0._rw,&
    Os(  0 : ndeg      ) = 0._rw,&
    Osum_(0: ndeg     ) = 0._rw
    
  real( rw ) ::                   &
    ctim,                         & ! Current simulation time (days).
    invf,                         & ! Absolute inverse Coriolis parameter (sec).
    ramp,                         & ! Wind forcing is ramped over `dt_r' days.
    gene,                         & ! Switch when GFB is ready to be used.
    bodf(             nlay,  2),  & ! Body force along x,y        (m s**(-2)).
!   External fields toward which model solution is relaxed at open boundaries.
!     (the same fields are used as the initial condition).
    fnud(   0 : ndeg, nlay,  3),  & ! 1=h, 2=u, 3=v.
!   Rate at which the model variables are relaxed at open boundaries.
    nudg(   0 : ndeg,        3),  & ! (non-dimensional).
    hdot(   0 : ndeg, nlay    ),  & ! Buoyancy forcing   dh/dt    (m s**(-1)).
    rs_h(2, 0 : ndeg, nlay    ),  & ! Right-hand Side of dh/dt    (m s**(-1)).
    dmdx(3, 0 : ndeg, nlay    ),  & ! Montgomery gradient along x (m s**(-2)).
    dmdy(3, 0 : ndeg, nlay    ),  & ! Montgomery gradient along y (m s**(-2)).
    tide(2, 1, 0 : ndeg,     3),  & ! ampl+phase, nconst, grid points, eta+u+v.
    w_ti(   1                 )     ! Frequency of tidal constituents (rad/day).

  real( r8 )           ::         &
    tres,                         &    ! Time at restart (days).
    hlay(0 : ndeg, nlay),      &               ! Layer thickness (meters).
        pi_rhs(0: ndeg),              &     ! Pressure correction
    pi_s(0: ndeg)                   !Pressure field
  integer, allocatable :: segm(:,:)    ! Segments of nudged open boundaries.
  logical              :: flag_nudging ! If active nudging at open boundaries.

  public run           ! Allow the procedure to be called from main program.

contains

subroutine run()
  implicit none
  call read_input_data()
  call integrate_time ()
end subroutine run

subroutine read_input_data()
  implicit none
  real   ( r4 ), allocatable :: ior4(:,:)
  real   ( rw )              :: h_2d( - 1 : lm + 2, - 1 : mm + 2 )
  real   ( r8 )              :: dmin, dmax, habv, hbel, &
                                h_0(  0 : ndeg, nlay )
  integer                    :: ipnt, lerm, ilay, i, j, lrec, unum

  lerm = len_trim( errm )
  errm = trim(     errm ) // &
         ' in subroutine read_input_data from module private_mod.f95,'

  call initialize_variables()

  h_0( :,:           ) = 0._r8
  h_2d(:,:           ) = 0._rw
  h_2d(1 : lm, 1 : mm) = cext**2._rw / grav

  call check_consistency_options()

! Get realistic bathymetry if available.

  call read_input_file( keyw = 'h_bo', h_2d = h_2d )

  call index_grid_points( h_2d ) ! Index velocity and scalar grid points.

! Set h_0, the equilibrium thickness of layers
!   (i.e. the thickness for an ocean at rest).

  dmin = minval( real(h_2d, r8), mask = h_2d > hdry   ) ! Min. depth.
  dmax = maxval( real(h_2d, r8)                       ) ! Max. depth.

  if   ( ocrp < 0.5_rw .and. nlay > 1 ) then
    if ( ocrp < 0.5_rw .and. &
         (real(topl(nlay), r8) * dmax + 10._r8 * real(hmin, r8)) >= dmin ) then
      errc = - 1
      errm = trim( errm ) // ' Please modify topl so that bathymetry is ' // &
                             'contained within lower layer.'
      call quit()
    end if
  elseif ( ocrp < 0.5_rw .and. nlay == 1 ) then
    if ( dmin <= 10._r8 * real(hmin, r8) ) then
      errc = - 1
      errm = trim(errm) // &
             ' Please adjust h_bo or hmin so that min(h_bo) > 10. * hmin.'
      call quit()
    end if
  end if

  if ( ocrp < 0.5_rw ) then ! First assume no outcrop.
    do ipnt = 1, ndeg
      if ( mk_n(ipnt) > 0.5_rw ) then
        i = subc(ipnt, 1)
        j = subc(ipnt, 2)
  
        do ilay = nlay, 1, - 1
          habv = 0._r8         ! Cumulated thickness of layers above ilay.
          hbel = 0._r8         ! Cumulated thickness of layers below ilay.
  
          if ( ilay > 1    ) then
            habv = dmax * real(topl(ilay), r8)
          end if
          if ( ilay < nlay ) then
            hbel = sum( h_0(ipnt, ilay + 1 : nlay) )
          end if
    
          h_0(ipnt, ilay) = real(h_2d(i, j), r8) - habv - hbel
        end do
      end if
    end do
  end if

  if ( ocrp > 0.5_rw ) then
    write( ioso, * ) 'Calculating equilibrium thickness h_0 of layers...'

    call get_equilibrium_thickness_h_0( h_2d = h_2d, h_0 = h_0 )

    write( ioso, * ) 'Completed the calculation of h_0.'
  end if

  unum = get_un()
  allocate( ior4( ndeg, nlay ) )
  inquire( iolength = lrec ) ior4(:,:)
  ior4(:,:) = real( h_0( 1 :, : ), r4 )
  open( unit = unum, form   = 'unformatted', access = 'direct', &
        recl = lrec, status = 'replace',     action = 'write',  &
        file = trim( odir ) // 'h_0.bin' )
    write( unum, rec = 1 ) ior4(:,:)
  close( unit = unum )
  deallocate( ior4 )

! Assume for now that the calculation starts from rest.

  do ilay = 1, nlay
    hlay( :, ilay ) = h_0( :, ilay ) * real( mk_n( : ), rw )
  end do

  ! Get nudging coefficients if avail.
  flag_nudging = .false.       ! Initialization.
  call read_input_file( keyw = 'nudg', h_2d = h_2d )

  call read_input_file( keyw = 'init' ) ! Initial/nudging fields.

  call read_input_file( keyw = 'bodf' ) ! Body force.

  call read_input_file( keyw = 'hdot' ) ! Buoyancy forcing (dh/dt).

  call read_input_file( keyw = 'taus' ) ! Surface stress (winds).

  call read_input_file( keyw = 'tide' ) ! Tidal forcing.

  call read_input_file( keyw = 'fcor' ) ! Coriolis parameter.

! Inverse Coriolis parameter: invf = 1/fcor.
!   Only used to prescribe Ekman transport at nudged open boundaries.
!   A domain-averaged value should be good enough.
!   If the value is within 5degrees from Equator, set invf to zero.

  invf = sum( fcor(:) ) / real(size(fcor(:)), rw)

  if ( abs( invf ) > 1.25e-5_rw ) then
    invf = 1._rw / invf
  else
    invf = 0._rw
  end if

  if ( rsta < 0.5_rw ) call save_metadata()

  write( ioso, * ) 'lm = ', lm
  write( ioso, * ) 'mm = ', mm

! If restart, read last output record and store in prognostic arrays.

  if ( rsta > 0.5_rw ) call read_restart_record

  if ( rsta < 0.5_rw ) then
!   Initialize output files, and write the initial condition on disk.
    call write_outputs()
  end if

  if ( errc == 0 ) then
    errm = errm(1 : lerm)
  else
    call quit()
  end if
end subroutine read_input_data

subroutine initialize_variables()
  implicit none

  tide(:,:,:,:) = 0._rw
  dmdx(:,:,:  ) = 0._rw
  dmdy(:,:,:  ) = 0._rw
  rs_h(:,:,:  ) = 0._rw
  fnud(:,:,:  ) = 0._rw
  nudg(:,:    ) = 0._rw
  bodf(:,:    ) = 0._rw
  h_u( :,:    ) = 0._rw
  h_v( :,:    ) = 0._rw
  u(   :,:    ) = 0._rw
  v(   :,:    ) = 0._rw
  delu(:,:    ) = 0._rw
  delv(:,:    ) = 0._rw
  UU4( :,:    ) = 0._rw
  VV4( :,:    ) = 0._rw
  hdot(:,:    ) = 0._rw
  w_ti(:      ) = 0._rw
  mkpe(:      ) = 0._rw
  rvor(:      ) = 0._rw
  pvor(:      ) = 0._rw
  dive(:      ) = 0._rw
  v_cc(:,:    ) = bvis
  v_ll(:,:    ) = bvis
  mk_n(:      ) = 0._rw
  mk_u(:      ) = 0._rw ! 0 everywhere except for offshore velocity pts.
  mk_v(:      ) = 0._rw ! 0 everywhere except for offshore velocity pts.
  mont(:      ) = 0._rw
  mkpi(:      ) = 0._rw
  d2hx(:      ) = 0._rw
  d2hy(:      ) = 0._rw
  h_bo(:      ) = 0._rw
  h_to(:      ) = 0._rw
  h_th(:      ) = 0._rw
  Ow(  :      ) = 0._rw
  Os(  :      ) = 0._rw
  Osum_(:     ) = 0._rw
  
  subc(:,:    ) = 0
  neig(:,:    ) = 0

  hlay(:,:    ) = 0._r8
  pi_rhs(:    ) = 0._r8
  pi_s(:        ) = 0._r8
  ctim          = 0._rw
  tres          = 0._r8

  fcor(:      ) = f0
  taus(:,1    ) = real(  tauw )
  taus(:,2    ) = aimag( tauw )
  tt3d(:,:,:  ) = 0._rw
  tb3d(:,:,:  ) = 0._rw
  tu3d(:,:,:  ) = 0._rw
end subroutine initialize_variables

subroutine get_equilibrium_thickness_h_0( h_2d, h_0 )
! The number of iterations needed for convergence is
!   approximately proportional to the number of layers (nlay).
! The convergence can be slow. Successive-over-relaxation (SOR) helps a bit.
  implicit none
  real     (rw  ), intent(in   ) :: h_2d( - 1 : lm + 2, - 1 : mm + 2 )
  real     (r8  ), intent(inout) :: h_0(  0 : ndeg, nlay )
  character(sstr)                :: strg
  logical                        :: conv(ndeg), flag
  integer                        :: ilay, ipnt, lerm, i, j, k, l, iter, imax, &
                                     c__1, c__3, c__5, c__7
  real     (r8  )                :: gues(nlay), rho8(nlay), dmax, hbot, habv, &
                                    cons(nlay), func(nlay), hbel, maxv, resu, &
                                    maug(nlay, nlay + 1), line(nlay + 1), thre, &
                                    Osum(0:ndeg), h_west(0:ndeg), h_sout(0:ndeg)       

  lerm = len_trim( errm )
  errm = trim(     errm ) // ' in subroutine get_equilibrium_thickness_h_0 ' //&
         'from module private_mod.f95,'

! Tolerance for calculation of equilibrium thickness h_0 (meters).

  thre = real( tole, r8 )
  write( unit = ioso, fmt = * ) 'Tolerance = ', tole, ' meters.'

  dmax = real( maxval( h_2d ), r8 ) ! Max. depth.
  rho8 = real( rhon(:), r8 )

! Calculate the thickness of each layer at the deepest point of domain.
!   Note that array `gues' will be overwritten afterward.
  do ilay = 1, nlay
    gues(ilay) = dmax * ( 1._r8 - real(topl(ilay), r8) )
    if ( ilay < nlay ) then
      gues(ilay) = gues(ilay)                        &
                 - dmax * ( 1._r8 - real(topl(ilay + 1),  r8) )
    end if
  end do

! Calculate the constant for each layer.

  do ilay = 1, nlay
    cons(ilay) = dmax * (- 1._r8) + sum( gues(:) )
    do k = 1, ilay - 1
      cons(ilay) = cons(ilay)                              &
                 - ( rho8(ilay) - rho8(k) ) * gues(k) / rho8(ilay)
    end do
  end do

  conv(:) = .true.
  flag    = .false.

!$OMP PARALLEL DO PRIVATE(ipnt,hbot,ilay,habv,hbel,gues,iter,i,j,k,l,imax,maxv,func,maug,line,resu) SHARED(flag)
  do ipnt = 1, ndeg
    if ( mk_n(ipnt) < 0.5_rw .or. flag ) then
      cycle
    end if

    hbot = real( h_2d( subc(ipnt, 1), subc(ipnt, 2) ), r8 )

!   Calculate an initial guess for the layer thickness.

    do ilay = nlay, 1, -1
      habv = 0._r8 ! Cumulated thickness of layers above ilay.
      hbel = 0._r8 ! Cumulated thickness of layers below ilay.
    
      habv = dmax * real(topl(ilay), r8)
      if ( ilay < nlay ) then
        hbel = sum( real(gues(ilay + 1 : nlay), r8) )
      end if
      gues(ilay) = hbot - habv - hbel
      gues(ilay) = max( gues(ilay), real(hsal, r8) )
    end do
      
    do iter = 1, itmx
      do i = 1, nlay
        func(i) = ( hbot - sum(gues(:)) )                                   &
                + 1._r8 / real(nsal - 1, r8)                                &
                * real(hsal, r8) * (real(hsal, r8) / gues(i))**(nsal - 1)   &
                + cons(i)
        func(i) = func(i) * (-1._r8)
        do j = 1, i - 1
          func(i) = func(i)                         &
                  - (rho8(i) - rho8(j)) * gues(j) / rho8(i)
        end do
      end do
    
      if ( iter == itmx ) then
!       Convergence was not achieved.
!       Throw an exception (flag) and end (cycle) the loop.
        conv(ipnt) = .false.
        flag       = .true.
        cycle
      end if
    
      if ( all( abs(func(:)) < thre ) ) then
!       Convergence was achieved. Move on to next grid point.
!       Diagnostics will be provided below in the serial part.
        h_0( ipnt, :) = gues(:)
        exit
      end if
      
      do i = 1, nlay
        do j = 1, nlay
!         Jacobian_i,j = df_i/dh_j
          maug(i, j) = min( rho8(i), rho8(j) ) / rho8(i)
          if ( i == j ) then
            maug(i, j) = maug(i, j) &
                       + (real(hsal, r8) / gues(j))**nsal
          end if
        end do
      end do
      
      maug(:, nlay + 1) = func(:) * (- 1._r8)

!     Solve system of linear equations Ax=b by Gaussian elimination.
!     maug is the augmented matrix [A|b].

      do k = 1, nlay ! Loop over columns.
!       Find pivot for column k.

        imax = 0
        maxv = 0._r8

        do ilay = k, nlay
          if ( abs( maug(ilay, k) ) > maxv ) then
            maxv = abs( maug(ilay, k) )
            imax = ilay
          end if
        end do

        if ( imax /= k ) then
!         Swap rows(k, imax).
          line(      :) = maug(   k, :)
          maug(   k, :) = maug(imax, :)
          maug(imax, :) = line(      :)
        end if

        do ilay = k + 1, nlay   ! Loop over all rows below pivot.
!         Do for all remaining elements in current row:
          do l = k, nlay + 1 ! Loop over remaining elements in current row.
            maug(ilay, l) = maug(ilay, l) &
                          - maug(   k, l) * (maug(ilay, k) / maug(k, k))
          end do
!         Fill lower triangular matrix with zeros.
          maug(ilay, k) = 0._r8
        end do
      end do

!     Back substitution.

      do ilay = nlay, 1, - 1
        resu = 0._r8
        do j = ilay + 1, nlay
          resu = resu + maug(ilay, j) * maug(j, nlay + 1)
        end do
        maug( ilay, nlay + 1 ) = (maug(ilay, nlay + 1) - resu) &
                               /  maug(ilay, ilay    )
      end do

      gues(:) = (1._r8 - sor) * gues(:) + sor * (maug(:, nlay + 1) + gues(:))

      if ( any( gues(:) <= thre ) ) then
        gues(:) = max( gues(:), thre )
      end if
    end do

  end do
!$OMP END PARALLEL DO

  if ( flag ) then ! If convergence was not achieved, provide diags.
    do ipnt = 1, ndeg
      if ( .not. conv( ipnt ) ) then ! Give diag on first point that failed.
        hbot = real( h_2d( subc( ipnt, 1 ), subc( ipnt, 2 ) ), r8 )
        errc = int( ipnt ) * (- 1)
        errm = trim( errm ) // ' calculation of h_layers did not converge,'
        write( strg, * ) thre
        errm = trim( errm ) // ' tolerance (meters) was ' // strg
        write( strg, * ) hbot
        errm = trim( errm ) // ', local depth (meters) is ' // strg
        errm = trim( errm ) // ', current value (meters) is'
        do i = 1, nlay
          write( strg, * ) h_0( ipnt, i )
          errm = trim( errm ) // strg
        end do
        call quit()
      end if
    end do
  end if

  if ( errc == 0 ) then
    errm = errm( 1 : lerm )
  else
    call quit()
  end if
  
    
  if (rgld > 0.5_rw) then
	   
! convert top layer equilibrium thickness to starting surface pressure	     

      pi_s(:)= (sum(h_0(:,:),dim=2)-h_th(:))*grav

! calculate operators Os, Ow, Osum, Osum_ based on h_0	  	  
	  Ow(:)=0._rw
	  Os(:)=0._rw
	  Osum(:)=0._rw
	  Osum_(:)=0._rw
  
	  do ipnt = 1, ndeg
		
		if (1<subc(ipnt,1) .and. subc(ipnt,1)<lm+1 .and. 1<subc(ipnt,2) .and. subc(ipnt,2)<mm+1) then
                        c__5=neig(5,ipnt)
                        c__7=neig(7,ipnt)
                        h_west(ipnt)= 0.5_rw*real( h_th( ipnt ) + h_th( c__5 ), rw )
                        h_sout(ipnt)= 0.5_rw*real( h_th( ipnt ) + h_th( c__7 ), rw )  
                        Ow(ipnt)=h_west(ipnt)/dl**2
                        Os(ipnt)=h_sout(ipnt)/dl**2

		elseif (subc(ipnt,1)==1 .and. 1<subc(ipnt,2) .and. subc(ipnt,2)<mm+1) then					
			c__7=neig(7,ipnt)
			h_sout(ipnt)= 0.5_rw*real( h_th( ipnt ) + h_th( c__7 ), rw )
			Ow(ipnt)=0._rw
			Os(ipnt)=h_sout(ipnt)/dl**2
		elseif (1<subc(ipnt,1) .and. subc(ipnt,1)<lm+1 .and. subc(ipnt,2)==1) then	
			c__5=neig(5,ipnt)
                        h_west(ipnt)=0.5_rw*real( h_th( ipnt ) + h_th( c__5 ), rw )
                        Ow(ipnt)=h_west(ipnt)/dl**2
                        Os(ipnt)=0._rw

                end if

          end do
  
          do ipnt = 1,ndeg
                if (subc(ipnt,1)<lm .and. subc(ipnt,2)<mm) then
                        c__1=neig(1,ipnt)
                        c__3=neig(3,ipnt)
                        Osum(ipnt) = Ow(ipnt) + Ow(c__1) + Os(ipnt) + Os(c__3);
                elseif (subc(ipnt,1)==lm .and. subc(ipnt,2)<mm) then
                        c__3=neig(3,ipnt)
                        Osum(ipnt)= Ow(ipnt)+ Os(ipnt)+Os(c__3)
                elseif (subc(ipnt,2)==mm .and. subc(ipnt,1)<lm) then
                        c__1=neig(1,ipnt)
                        Osum(ipnt)=Ow(ipnt)+Os(ipnt)+Ow(c__1)
                else
                        Osum(ipnt)=Ow(ipnt)+Os(ipnt)
                end if


                if (subc(ipnt,1)>0 .and. subc(ipnt,1)<lm+1 .and. subc(ipnt,2)>0 .and. subc(ipnt,2)<mm+1) then
                        Osum_(ipnt) = 1 / Osum(ipnt)
                end if
          end do

    end if
  
end subroutine get_equilibrium_thickness_h_0

subroutine index_grid_points( h_2d )
  implicit none
  real   (rw), intent(in)  :: h_2d( - 1 : lm + 2, - 1 : mm + 2 )
  character(sstr)          :: str1
  integer(i4), allocatable :: ioi4(:)
  integer                  :: lrec, unum, i_c, i, j, lerm,        &
! Position of each grid Cell, with posc = i + max_i * (j - 1).
!   Used for unpack under Matlab/Octave, e.g. u_ij(posc(:)) = vector_u(:).
                              indc( - 1 : lm + 2, - 1 : mm + 2 ), &
! Index of grid Cells within the vector.
                              posc( ndeg )

  lerm = len_trim( errm )
  errm = trim(     errm ) // &
         ' in subroutine index_grid_points from module private_mod.f95,'

! Note: By design, h_2d should be 0._rw outside the computational domain.

  indc(:,:) = 0
  posc(:  ) = 0

  i_c  = 0 ! Backhaus 2008 (Ocean Modelling) vectorized indexing.
           !   keep only wet ij cells,
           !   (extended to cover all eta,u,v,pvor points).

  do j = 0, mm + 1
    do i = 0, lm + 1
      if (     h_2d(i,         j        ) > hdry  .or. & ! eta  points.
           any(h_2d(i - 1 : i, j        ) > hdry) .or. & ! u    points.
           any(h_2d(i,         j - 1 : j) > hdry) .or. & ! v    points.
           any(h_2d(i - 1 : i, j - 1 : j) > hdry) ) then ! pvor points.
        i_c       = i_c + 1
        indc(i,j) = i_c
      end if
    end do
  end do

  if ( i_c /= ndeg ) then
    write(str1, '(1i8)') i_c
    errc = min( - 1, int( i_c ) * (- 1) )
    errm = trim(errm) // ' wrong input parameter! Please set ndeg = ' &
                      // trim(str1) // ' inside file shared_mod.f95.'
    call quit()
  end if

! Set periodic open boundaries.

  if ( xper > 0.5_rw ) then
    i = 1
    do j = 1, mm
      if ( h_2d(i, j) > hdry .and. &
           h_2d(lm,j) > hdry ) then
        indc(0,      j) = indc(lm, j)
        indc(lm + 1, j) = indc(1,  j)
        mk_u(indc(i,j)) = 1._rw
      end if

      if ( j > 1 .and. j <= mm ) then
        if ( all( h_2d(i,  j - 1 : j) > hdry ) .and. &
             all( h_2d(lm, j - 1 : j) > hdry ) ) then
          mkpe(indc(i,j)) = 1._rw
        end if
      end if

      if ( j == mm ) then
        if ( h_2d(i,  j) > hdry .and. &
             h_2d(lm, j) > hdry ) then
!         For psi point located along northern boundary.
          indc(0,      mm + 1) = indc(lm, mm + 1)
          indc(lm + 1, mm + 1) = indc(1,  mm + 1)
        end if
      end if
    end do
  end if

  if ( yper > 0.5_rw ) then
    j = 1
    do i = 1, lm
      if ( h_2d(i, j) > hdry .and. &
           h_2d(i,mm) > hdry ) then
        indc(i, 0     ) = indc(i, mm)
        indc(i, mm + 1) = indc(i,  1)
        mk_v(indc(i,j)) = 1._rw
      end if

      if ( i > 1 .and. i <= lm ) then
        if ( all( h_2d(i - 1 : i,  j) > hdry ) .and. &
             all( h_2d(i - 1 : i, mm) > hdry ) ) then
          mkpe(indc(i,j)) = 1._rw
        end if
      end if

      if ( i == lm ) then
        if ( h_2d(i, mm) > hdry .and. &
             h_2d(i,  j) > hdry ) then
!         For psi point located along eastern boundary.
          indc(lm + 1,      0) = indc(lm + 1, mm)
          indc(lm + 1, mm + 1) = indc(lm + 1,  1)
        end if
      end if
    end do
  end if

! Special treatment for corners in double-periodic domain.

  if ( xper > 0.5_rw .and. yper > 0.5_rw ) then
    if ( h_2d( 1, 1 ) > hdry .and. h_2d(lm, 1 ) > hdry &
                             .and. h_2d( 1, mm) > hdry ) then
      indc(   0, 0     ) = indc(lm, mm) ! Lower-left  corner.
      mkpe( indc(1, 1) ) = 1._rw
      indc(   0, mm + 1) = indc(lm, 1 ) ! Upper-left  corner.
    end if

    if ( h_2d(lm, mm) > hdry .and. h_2d( 1, mm) > hdry &
                             .and. h_2d(lm, 1 ) > hdry ) then
      indc(lm + 1, 0     ) = indc( 1, mm) ! Lower-right corner.
      indc(lm + 1, mm + 1) = indc( 1, 1 ) ! Upper-right corner.
    end if
  end if

! Here we explicitly assume that all vertical layers share the same horizontal
!   indexing.

  i_c = 0 ! Initialization.

  do j = 0, mm + 1
    do i = 0, lm + 1

      if (     h_2d(i,         j        ) > hdry  .or. & ! eta  points.
           any(h_2d(i - 1 : i, j        ) > hdry) .or. & ! u    points.
           any(h_2d(i,         j - 1 : j) > hdry) .or. & ! v    points.
           any(h_2d(i - 1 : i, j - 1 : j) > hdry) ) then ! pvor points.
        i_c = i_c + 1

        if ( h_2d(i,j) > hdry ) mk_n(i_c) = 1._rw

        if ( all(h_2d(i - 1 : i, j) > hdry) ) then
          mk_u(i_c) = 1._rw
        end if
        if ( all(h_2d(i, j - 1 : j) > hdry) ) then
          mk_v(i_c) = 1._rw
        end if
        if ( all( h_2d(i - 1 : i, j - 1 : j ) > hdry ) ) then
          mkpe(i_c) = 1._rw
        end if
        if ( any( h_2d(i - 1 : i, j - 1 : j ) > hdry ) ) then
          mkpi(i_c) = 1._rw
        end if
!       Matlab/Octave does not support negative/zero indices so they start at 1.
        posc(i_c   ) = i + 1 + j * (lm + 2)
        subc(i_c, 1) = i
        subc(i_c, 2) = j
        neig(1, i_c) = indc(i + 1, j    ) ! Counter-clockwise sweep,
        neig(2, i_c) = indc(i + 1, j + 1) !   starting from the right.
        neig(3, i_c) = indc(i,     j + 1)
        neig(4, i_c) = indc(i - 1, j + 1)
        neig(5, i_c) = indc(i - 1, j    )
        neig(6, i_c) = indc(i - 1, j - 1)
        neig(7, i_c) = indc(i,     j - 1)
        neig(8, i_c) = indc(i + 1, j - 1)
      end if

    end do
  end do

  unum = get_un()
  allocate( ioi4( size(posc(:)) ) )
  inquire( iolength = lrec ) ioi4(:)
  ioi4(:) = int( posc(:), i4 )
  open( unit = unum, form   = 'unformatted', access = 'direct', &
        recl = lrec, status = 'replace',     action = 'write',  &
        file = trim( odir ) // 'grid.bin' )
    write( unum, rec = 1 ) ioi4(:)
    ioi4(:) = nint( mk_n(1 :), i4 )
    write( unum, rec = 2 ) ioi4(:)
    ioi4(:) = nint( mk_u(1 :), i4 )
    write( unum, rec = 3 ) ioi4(:)
    ioi4(:) = nint( mk_v(1 :), i4 )
    write( unum, rec = 4 ) ioi4(:)
    ioi4(:) = nint( mkpi(1 :), i4 )
    write( unum, rec = 5 ) ioi4(:)
  close( unit = unum )
  deallocate( ioi4  )

! Store the bathymetric and top topography grid in indexed format.

  do i_c = 0, ndeg
    i           = subc( i_c, 1 )
    j           = subc( i_c, 2 )
    h_th( i_c ) = h_2d( i,   j )
  end do  

  if ( errc == 0 ) then
    errm = errm( 1 : lerm )
  else
    call quit()
  end if
end subroutine index_grid_points

subroutine read_input_file( keyw, h_2d )
  implicit none
  character(4), intent(in   )           :: keyw
  real   (rw),  intent(inout), optional :: h_2d( - 1 : lm + 2, - 1 : mm + 2 )
  logical                               :: is_e
  integer                               :: lrec, unum, i, j, ilay, ipnt
  real   (r4),  allocatable             :: ior4(:,:), ior4_(:,:), ior4_3d(:,:,:), &
                                           ior4_4d(:,:,:,:), ior4_5d(:,:,:,:,:)

  inquire( iostat = errc, exist = is_e, file = trim(idir) // keyw // '.bin' )
  if ( .not. is_e ) return
  if     ( keyw == 'h_bo' .or. keyw == 'fcor' ) then
    allocate( ior4   ( 0 : lm + 1, 0 : mm + 1 ) )
    inquire( iolength = lrec ) ior4(:,:)
    if (topt > 0.5) then
      allocate ( ior4_ ( 0: lm + 1, 0 : mm + 1 ) )
      inquire(iolength = lrec ) ior4_(:,:)
    end if
  elseif ( keyw == 'nudg' ) then
    allocate( ior4_3d( 0 : lm + 1, 0 : mm + 1, 3 ) )
    inquire( iolength = lrec ) ior4_3d(:,:,:)
  elseif ( keyw == 'init' ) then
    allocate( ior4_4d( 0 : lm + 1, 0 : mm + 1, nlay, 3 ) )
    inquire( iolength = lrec ) ior4_4d(:,:,:,:)
  elseif ( keyw == 'bodf'                     ) then
    allocate( ior4   ( nlay, 2 ) )
    inquire( iolength = lrec ) ior4(:,:)
  elseif ( keyw == 'hdot' ) then
    allocate( ior4_3d( 0 : lm + 1, 0 : mm + 1, nlay ) )
    inquire( iolength = lrec ) ior4_3d(:,:,:)
  elseif ( keyw == 'tide' ) then
    allocate( ior4_5d( 2, size( w_ti ), 0 : lm + 1, 0 : mm + 1, 3 ) )
    inquire( iolength = lrec ) ior4_5d(:,:,:,:,:)
  elseif ( keyw == 'taus' ) then
    allocate( ior4_3d( 0 : lm + 1, 0 : mm + 1, 2 ) )
    inquire( iolength = lrec ) ior4_3d(:,:,:)
  end if
  unum = get_un()
  open( unit = unum, status = 'old',    iostat = errc, action = 'read', &
        recl = lrec, access = 'direct', form   = 'unformatted',         &
        file = trim( idir ) // keyw // '.bin' )
  if     ( keyw == 'h_bo' .or. keyw == 'fcor' .or. keyw == 'bodf' ) then
    read( unit = unum, iostat = errc, rec = 1 ) ior4(:,:)
    if (topt>0.5) then
      open( unit = unum, status = 'old',    iostat = errc, action = 'read', &
        recl = lrec, access = 'direct', form   = 'unformatted',         &
        file = trim( idir ) // 'h_to.bin' )
      read( unit = unum, iostat = errc, rec = 1 ) ior4_(:,:)
    end if
  elseif ( keyw == 'init'                                         ) then
    read( unit = unum, iostat = errc, rec = 1 ) ior4_4d(:,:,:,:)
  elseif ( keyw == 'hdot' .or. keyw == 'taus' .or. keyw == 'nudg' ) then
    read( unit = unum, iostat = errc, rec = 1 ) ior4_3d(:,:,:)
  elseif ( keyw == 'tide'                                         ) then
    read( unit = unum, iostat = errc, rec = 1 ) ior4_5d(:,:,:,:,:)
  end if
  if ( errc /= 0 ) then ! Fatal error.
    errm = trim(errm) // ' could not open/read file ' // keyw &
        // '.bin from directory ' // trim(idir)
    call quit()
  end if
  if     ( keyw == 'h_bo' ) then
    h_2d(  :,          :       ) = 0._rw
    h_2d(0 : lm + 1, 0 : mm + 1) = real( ior4, rw )
    if (topt > 0.5) then
      h_2d(0 : lm + 1, 0 : mm + 1) = real( ior4 - ior4_, rw )
      deallocate( ior4_)
    end if
    deallocate( ior4 )
    where( h_2d(:,:) < hdry ) h_2d(:,:) = 0._rw
    h_2d(     0,      :) = 0._rw ! Enforce margins of dry cells.
    h_2d(     :,      0) = 0._rw
    h_2d(lm + 1,      :) = 0._rw
    h_2d(     :, mm + 1) = 0._rw
  elseif ( keyw == 'bodf' ) then
    bodf(:,:) = real( ior4(:,:), rw )
    deallocate( ior4 )
  elseif ( keyw == 'nudg' ) then
!   Nudging coefficients inside nudg.bin are cell-centered.
!   This is fine for `eta', but interpolation is necessary for u,v points.
!$OMP PARALLEL DO PRIVATE(ipnt,i,j)
    do ipnt = 1, ndeg
      i = subc( ipnt, 1 )
      j = subc( ipnt, 2 )

      nudg( ipnt, ix_n ) = real( ior4_3d( i,     j,     ix_n ), rw )

      if ( all( ior4_3d( i - 1 : i, j, ix_u ) > 1.e-9_r4 ) ) then
        nudg( ipnt, ix_u ) = real( ior4_3d( i,     j,     ix_u ), rw ) * 0.5_rw &
                           + real( ior4_3d( i - 1, j,     ix_u ), rw ) * 0.5_rw
      else
        nudg( ipnt, ix_u ) = 0._rw
      end if

      if ( all( ior4_3d( i, j - 1 : j, ix_v ) > 1.e-9_r4 ) ) then
        nudg( ipnt, ix_v ) = real( ior4_3d( i,     j,     ix_v ), rw ) * 0.5_rw &
                           + real( ior4_3d( i,     j - 1, ix_v ), rw ) * 0.5_rw
      else
        nudg( ipnt, ix_v ) = 0._rw
      end if
    end do
!$OMP END PARALLEL DO
    if ( any( nudg(:,:) > 1.e-9_rw ) ) then
      flag_nudging = .true.
      call index_boundary_points( nudg = ior4_3d, h_2d = h_2d )
    end if
    deallocate( ior4_3d )

    do ilay = 1, nlay
      do ipnt = 1, ndeg
        i = subc( ipnt, 1 )
        j = subc( ipnt, 2 )

        fnud( ipnt, ilay, ix_n ) = real( hlay( ipnt, ilay ), rw )
      end do
    end do
  elseif ( keyw == 'init' ) then
    do ilay = 1, nlay
      do ipnt = 1, ndeg
        i = subc( ipnt, 1 )
        j = subc( ipnt, 2 )

        if ( ilay < nlay ) then
          fnud( ipnt, ilay, ix_n )                              &
            = real( hlay( ipnt, ilay )                          &
                  + real( ior4_4d( i, j, ilay,     ix_n ), r8 ) &
                  - real( ior4_4d( i, j, ilay + 1, ix_n ), r8 ), rw )
        else
          fnud( ipnt, ilay, ix_n )                              &
            = real( hlay( ipnt, ilay )                          &
                  + real( ior4_4d( i, j, ilay,     ix_n ), r8 ), rw )
        end if

        fnud( ipnt, ilay, ix_n ) = fnud( ipnt, ilay, ix_n ) * mk_n( ipnt )
        fnud( ipnt, ilay, ix_u ) = real( ior4_4d(i, j, ilay, ix_u), rw )
        fnud( ipnt, ilay, ix_v ) = real( ior4_4d(i, j, ilay, ix_v), rw )

        if ( rsta < 0.5_rw ) then
          hlay( ipnt, ilay ) = real( fnud( ipnt, ilay, ix_n ) * mk_n( ipnt ), r8 )
          u(    ipnt, ilay ) = fnud( ipnt, ilay, ix_u )
          v(    ipnt, ilay ) = fnud( ipnt, ilay, ix_v )
        end if
      end do
    end do
    deallocate( ior4_4d )
  elseif ( keyw == 'hdot' ) then
    do ilay = 1, nlay
      do ipnt = 1, ndeg
        i = subc( ipnt, 1 )
        j = subc( ipnt, 2 )
        hdot( ipnt, ilay ) = real( ior4_3d( i, j, ilay ), rw )
      end do
    end do
    deallocate( ior4_3d )
  elseif ( keyw == 'taus' ) then
    taus(:,:) = 0._rw
!   Index cell-centered wind stress forcing.
!$OMP PARALLEL DO PRIVATE(ipnt,i,j)
    do ipnt = 1, ndeg
      i = subc(ipnt, 1)
      j = subc(ipnt, 2)
      taus(ipnt,1) = real( ior4_3d(i,j,1), rw )
      taus(ipnt,2) = real( ior4_3d(i,j,2), rw )
    end do
!$OMP END PARALLEL DO
    deallocate( ior4_3d )
  elseif ( keyw == 'fcor' ) then
    fcor(0) = real( sum( ior4 ) / real( size( ior4 ), r4 ), rw )
!   fcor(0) = real( ior4(1,1), rw )
!$OMP PARALLEL DO PRIVATE(ipnt,i,j)
    do ipnt = 1, ndeg
      i = subc( ipnt, 1 )
      j = subc( ipnt, 2 )
!     Interpolate the 2-D array at the lower-left corner (psi point).
      if ( i > 0 .and. j > 0 ) then
        fcor( ipnt ) = real( ior4( i,     j     ) * 0.25_r4 &
                           + ior4( i - 1, j     ) * 0.25_r4 &
                           + ior4( i,     j - 1 ) * 0.25_r4 &
                           + ior4( i - 1, j - 1 ) * 0.25_r4, rw )
      else
        fcor( ipnt ) = real( ior4( i, j ), rw )
      end if
    end do
!$OMP END PARALLEL DO
    deallocate( ior4 )
  elseif ( keyw == 'tide' ) then
    do ipnt = 1, size( w_ti )
      w_ti( ipnt ) = real( ior4_5d( 1, ipnt, 0, 0, 1 ), rw )
      write( unit = ioso, fmt = * ) &
        'Tidal constituent (rad/day) = ', w_ti( ipnt )
    end do
!$OMP PARALLEL DO PRIVATE(ipnt,i,j)
    do ipnt = 1, ndeg
      i = subc( ipnt, 1 )
      j = subc( ipnt, 2 )
      tide( :, :, ipnt, : ) = real( ior4_5d( :, :, i, j, : ), rw )
    end do
!$OMP END PARALLEL DO
    deallocate( ior4_5d )
  end if
  close( unit = unum, iostat = errc, status = 'keep' )
end subroutine read_input_file

subroutine check_consistency_options()
  implicit none
  integer :: lerm
  logical :: is_e(9)

  lerm = len_trim( errm )
  errm = trim(     errm ) // ' in subroutine check_consistency_options ' // &
         'from module private_mod.f95,'

! Ensure that paths are left-justified.

  idir = adjustl( idir )
  odir = adjustl( odir )

! Ensure that idir and odir exist.

  inquire( iostat = errc, exist = is_e(1), file = trim(idir) )
  inquire( iostat = errc, exist = is_e(2), file = trim(odir) )

  if ( .not. is_e(1) ) then
    errc = - 1
    errm = trim(errm) // ' idir is set to ' // trim(idir) &
         // ' but this directory does not exist. Program stopped.'
    call quit()
  end if

  if ( .not. is_e(2) ) then
    errc = - 1
    errm = trim(errm) // ' odir is set to ' // trim(odir) &
         // ' but this directory does not exist. Program stopped.'
    call quit()
  end if

! Ensure that paths end with the '/' character.

  if ( idir(len_trim(idir) : len_trim(idir)) /= '/' ) then
    idir(len_trim(idir) + 1 : len_trim(idir) + 1) = '/'
  end if
  if ( odir(len_trim(odir) : len_trim(odir)) /= '/' ) then
    odir(len_trim(odir) + 1 : len_trim(odir) + 1) = '/'
  end if

  if ( lm < 1 .or. mm < 1 ) then
    errc = - 1
    errm = trim(errm) // ' grid dimensions (lm,mm) should be >= 1.'
  end if

  if ( dl < 1.e1_rw ) then
    errc = - 1
    errm = trim(errm) // ' mesh size (dl) should be >= 10 meters.'
  end if

  if ( abs(f0) > 2.e-4_rw ) then
    errc = - 1
    errm = trim(errm) // ' Coriolis parameter (f0, in s**(-1)) should ' // &
           'be within: -2x10**(-4) < f0 < 2x10**(-4).'
  end if

  if ( dvis < 0._rw .or. dvis > 5._rw ) then
    errc = - 1
    errm = trim(errm) // ' Viscosity coefficient should be within: ' // &
           '0 <= dvis < 5.0.'
  end if

  if ( (bdrg < 0._rw .or. bdrg > 15.e-3_rw) .and. qdrg > 0.5_rw ) then
    errc = - 1
    errm = trim(errm) // ' quadratic bottom drag coefficient bdrg ' // &
           'should be within: 0 <= bdrg < 5x10**(-3).'
  elseif ( bdrg < 0._rw .or. bdrg > 5.e-2_rw ) then
    errc = - 1
    errm = trim(errm) // ' linear bottom drag coefficient bdrg has units ' // &
           'of m s**(-1) and should be within: 0 <= bdrg < 5x10**(-3) x u_max.'
  end if
  
  if ( (tdrg < 0._rw .or. tdrg > 15.e-3_rw) .and. qdrg > 0.5_rw ) then
    errc = - 1
    errm = trim(errm) // ' quadratic top drag coefficient tdrg ' // &
           'should be within: 0 <= tdrg < 5x10**(-3).'
  elseif ( tdrg < 0._rw .or. tdrg > 5.e-2_rw ) then
    errc = - 1
    errm = trim(errm) // ' linear top drag coefficient tdrg has units ' // &
           'of m s**(-1) and should be within: 0 <= tdrg < 5x10**(-3) x u_max.'
  end if

  if ( errc == 0 ) then
    errm = errm(1 : lerm)
  else
    call quit()
  end if
end subroutine check_consistency_options

subroutine index_boundary_points( nudg, h_2d )
  implicit none
  real( r4 ), intent( in ) :: nudg(   0 : lm + 1,   0 : mm + 1, 3 )
  real( rw ), intent( in ) :: h_2d( - 1 : lm + 2, - 1 : mm + 2    )
  integer                  :: nseg, i, j, irep, lerm, i_c, nnud
  integer,    allocatable  :: indc(:,:)

  lerm = len_trim( errm )
  errm = trim(     errm ) // ' in subroutine index_boundary_points ' // &
         'from module private_mod.f95,'

! Note: by default, h_2d is 0. outside the computational domain.
! Note: values inside of array `nudg' are cell-centered;
!       i.e. they must be interpolated for velocity points.

  allocate( indc(-1 : lm + 2, -1 : mm + 2) )
  indc(:,:) = 0
  i_c       = 0

  do irep = 1, 2

    if ( irep == 2 ) then
      allocate( segm( nseg, 18 ) )
      segm(:,:) = 0
    end if
    nseg = 0
    nnud = 0

    do j = 0, mm + 1
      do i = 0, lm + 1
        if ( irep == 1 ) then
          if (     h_2d(i,         j        ) > hdry  .or. & ! rho points.
               any(h_2d(i - 1 : i, j        ) > hdry) .or. & ! u   points.
               any(h_2d(i,         j - 1 : j) > hdry) .or. & ! v   points.
               any(h_2d(i - 1 : i, j - 1 : j) > hdry) ) then ! psi points.
            i_c       = i_c + 1
            indc(i,j) = i_c
          end if
          if ( any(nudg(i,j,:) > 1.e-9_r4) .and. indc(i,j) > 0 ) then
            nnud = nnud + 1
          end if
        else
          if ( any(nudg(i,j,:) > 1.e-9_r4) .and. indc(i,j) > 0 ) then
            nnud       = nnud + 1
          end if
        end if
        if   (        h_2d(i,     j) > hdry .and. & ! Western boundary.
               .not. (h_2d(i - 1, j) > hdry) ) then
          if ( nudg( i,     j, ix_u ) > tiny( 0._r4 ) .and. &
               nudg( i - 1, j, ix_u ) > tiny( 0._r4 ) .and. &
               xper < 0.5_rw ) then            ! Western nudged open
            nseg = nseg + 1                    !   boundary.
            if ( irep == 2 ) then
              segm( nseg, 1 ) = indc(i,j)      ! Index of the velocity point
                                               !   normal to the boundary.
              segm( nseg, 2 ) = i              ! i     of the ...
              segm( nseg, 3 ) = j
              segm( nseg, 4 ) = 1              ! This is a zonal velocity point.
              segm( nseg, 6 ) = 1
              segm( nseg, 7 ) = indc(i - 1, j) ! The dry cell.
              segm( nseg, 8 ) = i - 1          ! The dry cell.
              segm( nseg, 9 ) = j              ! The dry cell.
              segm( nseg,10 ) = indc(i,     j) ! The wet cell.
              segm( nseg,11 ) = i              ! The wet cell.
              segm( nseg,12 ) = j              ! The wet cell.
              segm( nseg,13 ) = indc(i + 1, j) ! The interior, normal velo. pnt.
              segm( nseg,14 ) = i + 1          ! The interior, normal velo. pnt.
              segm( nseg,15 ) = j              ! The interior, normal velo. pnt.
              segm( nseg,16 ) = indc(i + 1, j) ! The interior cell.
              segm( nseg,17 ) = i + 1          ! The interior cell.
              segm( nseg,18 ) = j              ! The interior cell.
            end if
          end if
        end if
        if   ( .not. (h_2d(i,     j) > hdry) .and. & ! Eastern boundary.
                      h_2d(i - 1, j) > hdry ) then
          if ( nudg( i - 1, j, ix_u ) > tiny(0._r4) .and. &
               nudg( i,     j, ix_u ) > tiny(0._r4) .and. &
               xper < 0.5_rw ) then            ! Eastern nudged open
                                               !   boundary.
            nseg = nseg + 1
            if ( irep == 2 ) then
              segm( nseg, 1 ) = indc(i,j)      ! Index of the velocity point
                                               !   normal to the boundary.
              segm( nseg, 2 ) = i              ! i     of the ...
              segm( nseg, 3 ) = j
              segm( nseg, 4 ) = 1              ! This is a zonal velocity point.
              segm( nseg, 6 ) = - 1
              segm( nseg, 7 ) = indc(i,     j) ! The dry cell.
              segm( nseg, 8 ) = i              ! The dry cell.
              segm( nseg, 9 ) = j              ! The dry cell.
              segm( nseg,10 ) = indc(i - 1, j) ! The wet cell.
              segm( nseg,11 ) = i - 1          ! The wet cell.
              segm( nseg,12 ) = j              ! The wet cell.
              segm( nseg,13 ) = indc(i - 1, j) ! The interior, normal velo. pnt.
              segm( nseg,14 ) = i - 1          ! The interior, normal velo. pnt.
              segm( nseg,15 ) = j              ! The interior, normal velo. pnt.
              segm( nseg,16 ) = indc(i - 2, j) ! The interior cell.
              segm( nseg,17 ) = i - 2          ! The interior cell.
              segm( nseg,18 ) = j              ! The interior cell.
            end if
          end if
        end if
        if   (        h_2d(i, j    ) > hdry .and. & ! Southern boundary.
               .not. (h_2d(i, j - 1) > hdry) ) then
          if ( nudg( i, j,     ix_v ) > tiny(0._r4) .and. &
               nudg( i, j - 1, ix_v ) > tiny(0._r4) .and. &
               yper < 0.5_rw ) then            ! Southern nudged open
                                               !   boundary.
            nseg = nseg + 1
            if ( irep == 2 ) then
              segm( nseg, 1 ) = indc(i,j)      ! Index of the velocity point
                                               !   normal to the boundary.
              segm( nseg, 2 ) = i              ! i     of the ...
              segm( nseg, 3 ) = j
              segm( nseg, 5 ) = 1              ! This is a meridional velocity
                                               !   point.
              segm( nseg, 6 ) = 1
              segm( nseg, 7 ) = indc(i, j - 1) ! The dry cell.
              segm( nseg, 8 ) = i              ! The dry cell.
              segm( nseg, 9 ) = j - 1          ! The dry cell.
              segm( nseg,10 ) = indc(i, j    ) ! The wet cell.
              segm( nseg,11 ) = i              ! The wet cell.
              segm( nseg,12 ) = j              ! The wet cell.
              segm( nseg,13 ) = indc(i, j + 1) ! The interior, normal velo. pnt.
              segm( nseg,14 ) = i              ! The interior, normal velo. pnt.
              segm( nseg,15 ) = j + 1          ! The interior, normal velo. pnt.
              segm( nseg,16 ) = indc(i, j + 1) ! The interior cell.
              segm( nseg,17 ) = i              ! The interior cell.
              segm( nseg,18 ) = j + 1          ! The interior cell.
            end if
          end if
        end if
        if   ( .not. (h_2d(i, j    ) > hdry) .and. & ! Northern boundary.
                      h_2d(i, j - 1) > hdry ) then
          if ( nudg( i, j - 1, ix_v ) > tiny(0._r4) .and. &
               nudg( i, j,     ix_v ) > tiny(0._r4) .and. &
               yper < 0.5_rw ) then            ! Northern nudged open
                                               !   boundary.
            nseg = nseg + 1
            if ( irep == 2 ) then
              segm( nseg, 1 ) = indc(i,j)      ! Index of the velocity point
                                               !   normal to the boundary.
              segm( nseg, 2 ) = i              ! i     of the ...
              segm( nseg, 3 ) = j
              segm( nseg, 5 ) = 1              ! This is a meridional velocity
                                               !   point.
              segm( nseg, 6 ) = - 1
              segm( nseg, 7 ) = indc(i,     j) ! The dry cell.
              segm( nseg, 8 ) = i              ! The dry cell.
              segm( nseg, 9 ) = j              ! The dry cell.
              segm( nseg,10 ) = indc(i, j - 1) ! The wet cell.
              segm( nseg,11 ) = i              ! The wet cell.
              segm( nseg,12 ) = j - 1          ! The wet cell.
              segm( nseg,13 ) = indc(i, j - 1) ! The interior, normal velo. pnt.
              segm( nseg,14 ) = i              ! The interior, normal velo. pnt.
              segm( nseg,15 ) = j - 1          ! The interior, normal velo. pnt.
              segm( nseg,16 ) = indc(i, j - 2) ! The interior cell.
              segm( nseg,17 ) = i              ! The interior cell.
              segm( nseg,18 ) = j - 2          ! The interior cell.
            end if
          end if
        end if
      end do
    end do

    if ( nseg == 0 ) then
      errc = - 1
      errm = trim( errm ) // ' the nudged open boundary segments ' // &
             'could not be identified.'
      call quit()
    end if
  end do ! Repetition loop (1st time: get nseg, allocate,
         !   then 2nd time: store information).

  if ( errc == 0 ) then
    errm = errm(1 : lerm)
  else
    call quit()
  end if
end subroutine index_boundary_points

subroutine save_metadata()
  implicit none
  integer :: unum

! Write a copy of the parameters in directory `odir'.

  if ( rsta < 0.5_rw ) then
    unum = get_un()
    open(   unit = unum, iostat = errc, action = 'write', status = 'replace', &
            file = trim( odir ) // 'param_basin.txt' )
      write(unit = unum, fmt = *) 'lm             = ',  lm,              ';'
      write(unit = unum, fmt = *) 'mm             = ',  mm,              ';'
      write(unit = unum, fmt = *) 'nlay           = ',  nlay,            ';'
      write(unit = unum, fmt = *) 'ndeg           = ',  ndeg,            ';'
      write(unit = unum, fmt = *) 'dl             = ',  dl,              ';'
      write(unit = unum, fmt = *) 'cext           = ',  cext,            ';'
      write(unit = unum, fmt = *) 'f0             = ',  f0,              ';'
      write(unit = unum, fmt = *) 'rhon           = [', rhon(1 : nlay), '];'
      write(unit = unum, fmt = *) 'topl           = [', topl(1 : nlay), '];'
      write(unit = unum, fmt = *) 'dt_s           = ',  dt_s,            ';'
      write(unit = unum, fmt = *) 'dt_o           = ',  dt_o,            ';'
      write(unit = unum, fmt = *) 'dt_r           = ',  dt_r,            ';'
      write(unit = unum, fmt = *) 'dt3d           = ',  dt3d,            ';'
      write(unit = unum, fmt = *) 'bvis           = ',  bvis,            ';'
      write(unit = unum, fmt = *) 'dvis           = ',  dvis,            ';'
      write(unit = unum, fmt = *) 'svis           = ',  svis,            ';' 
      write(unit = unum, fmt = *) 'bdrg           = ',  bdrg,            ';'
      write(unit = unum, fmt = *) 'tdrg           = ',  tdrg,            ';'
      write(unit = unum, fmt = *) 'tole           = ',  tole,            ';'
      write(unit = unum, fmt = *) 'nsal           = ',  nsal,            ';'
      write(unit = unum, fmt = *) 'hsal           = ',  hsal,            ';'
      write(unit = unum, fmt = *) 'hmin           = ',  hmin,            ';'
      write(unit = unum, fmt = *) 'hdry           = ',  hdry,            ';'
      write(unit = unum, fmt = *) 'hsbl           = ',  hsbl,            ';'
      write(unit = unum, fmt = *) 'hbbl           = ',  hsbl,            ';'
      write(unit = unum, fmt = *) 'g_fb           = ',  g_fb,            ';'
      write(unit = unum, fmt = *) 'uadv           = ',  uadv,            ';'
      write(unit = unum, fmt = *) 'qdrg           = ',  qdrg,            ';'
      write(unit = unum, fmt = *) 'ocrp           = ',  ocrp,            ';'
      write(unit = unum, fmt = *) 'tauwx          = ',  real ( tauw ),   ';'
      write(unit = unum, fmt = *) 'tauwy          = ',  aimag( tauw ),   ';'
      write(unit = unum, fmt = *) 'rsta           = ',  rsta,            ';'
      write(unit = unum, fmt = *) 'xper           = ',  xper,            ';'
      write(unit = unum, fmt = *) 'yper           = ',  yper,            ';'
      write(unit = unum, fmt = *) 'diag           = ',  diag,            ';'
      write(unit = unum, fmt = *) 'rgld           = ',  rgld,            ';'
      write(unit = unum, fmt = *) 'mcbc           = ',  mcbc,            ';'
      write(unit = unum, fmt = *) 'topt           = ',  topt,            ';'
      write(unit = unum, fmt = *) 'idir           = ', '''',              &
                                                  trim( idir ), '''',    ';'
      write(unit = unum, fmt = *) 'desc           = ', '''',              &
                                                  trim( desc ), '''',    ';'
      write(unit = unum, fmt = *) 'dt             = ',  dt,              ';'
    close( unit = unum )
  end if
end subroutine save_metadata

subroutine read_restart_record
  implicit none
  real   (r8) :: date
  integer     :: unum, lerm, irec

  lerm = len_trim( errm )
  errm = trim(     errm ) // ' in subroutine read_restart_record ' // &
                             'from module private_mod.f95,'

! Get the number of records (from file time.txt).

  unum = get_un()
  irec = 0
  open( unit = unum, form = 'formatted',   iostat = errc, action = 'read', &
        file = trim( odir ) // 'time.txt', status = 'old' )
    do
      read( unit = unum, iostat = errc, fmt = * ) date
      if ( errc == 0 ) then
        irec = irec + 1
        tres = date
      else
        errc = 0
        exit
      end if
    end do
  close( unit = unum, iostat = errc, status = 'keep' )

  write( unit = ioso, fmt = * ) '*** Restarting from record number ', irec, &
                                ' at time = ', real( tres, rw )

! Read u,v,eta.

  call read_array( path = odir, var = 'u___', irec = irec )
  call read_array( path = odir, var = 'v___', irec = irec )
  call read_array( path = odir, var = 'eta_', irec = irec )

! Get ready for next outputs (increment irec within subroutine write_outputs).

  call write_outputs( restart_record = irec )

  if ( errc == 0 ) then
    errm = errm(1 : lerm)
  else
    call quit()
  end if
end subroutine read_restart_record

subroutine read_array( path, var, irec )
! This subroutine is used to restart the model from a previous calculation.
  implicit none
  character( 4),   intent(in)  :: var
  character(sstr), intent(in)  :: path
  integer,         intent(in)  :: irec
  real     (r4),   allocatable :: ior4(:,:), h_0(:,:)
  integer :: lrec(2), unum(2), lerm, ipnt, ilay

  lerm = len_trim( errm )
  errm = trim(     errm ) &
       // ' in subroutine read_array from module private_mod.f95,'

  allocate( ior4( ndeg, nlay ) )

  inquire( iolength = lrec(1) ) ior4(:,:)
  unum(1) = get_un()
  open( unit = unum(1), form = 'unformatted', iostat = errc,     &
        file = trim( path ) // var // '.bin', access = 'direct', &
        recl = lrec(1), status = 'old',       action = 'read' )
  read( unit = unum(1), iostat = errc, rec = irec ) ior4(:,:)
  if ( errc /= 0 ) then ! Fatal error.
    errm = trim( errm ) // ' could not open/read file ' // var // &
                           '.bin from directory ' // trim( path )
    call quit()
  else
    if     ( var == 'u___' ) then
      u(1 : ndeg, :) = real( ior4(:,:), kind = rw )
    elseif ( var == 'v___' ) then
      v(1 : ndeg, :) = real( ior4(:,:), kind = rw ) 
    elseif ( var == 'eta_' ) then

!     You need to read the isopycnal elevation `eta' and convert it into a
!       layer thickness. Such inversion requires the undisturbed thickness
!       `h_0' from a file.

      unum(2) = get_un()
      allocate( h_0( ndeg, nlay ) )
      inquire( iolength = lrec(2) ) h_0(:,:)
      open( unit = unum(2), form   = 'unformatted', access = 'direct', &
            recl = lrec(2), status = 'old',         action = 'read',   &
            file = trim( odir ) // 'h_0.bin' )
      read( unum(2), rec = 1 ) h_0(:,:)
      close( unit = unum(2) )

      do ilay = 1, nlay
        do ipnt = 1, ndeg

          if ( ilay < nlay ) then
            hlay( ipnt, ilay )                     &
              = real( h_0(  ipnt, ilay     ), r8 ) &
              + real( ior4( ipnt, ilay     ), r8 ) &
              - real( ior4( ipnt, ilay + 1 ), r8 )
          else
            hlay( ipnt, ilay )                     &
              = real( h_0(  ipnt, ilay     ), r8 ) &
              + real( ior4( ipnt, ilay     ), r8 )
          end if
          hlay( ipnt, ilay ) = hlay( ipnt, ilay ) * real( mk_n( ipnt ), r8 )

        end do
      end do
      deallocate( h_0 )

    end if
    deallocate( ior4 )
    close( unit = unum(1), iostat = errc, status = 'keep' )
  end if

  if ( errc == 0 ) then
    errm = errm( 1 : lerm )
  else
    call quit()
  end if
end subroutine read_array

subroutine update_u( ilay )
  implicit none
  integer, intent(in) :: ilay
  integer             :: ipnt, c__3, c__4, c__5
  real ( rw )         :: vcor, tauw, ufor, mask, hcen, &
                         i__h, rhsi, uold, i_dl, i_r0, i_r1, dmd4
  i_dl = 1._rw / dl
  i_r0 = 1._rw / rho0
  i_r1 = 1._rw / rhon(1)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ipnt,c__3,c__4,c__5,tauw,vcor,ufor,mask,hcen,i__h,rhsi,uold,dmd4)
  do ipnt = 1, ndeg
    c__3 = neig( 3,    ipnt )
    c__4 = neig( 4,    ipnt )
    c__5 = neig( 5,    ipnt )
    mask = mk_u(       ipnt )
    hcen = real( hlay( c__5, ilay ) + hlay( ipnt, ilay ), rw ) / (1._rw + mask)
    i__h = 1._rw / (hcen + 1._rw - mask)
    uold = u( ipnt, ilay )

    dmd4 = ( mont( c__5 ) - mont( ipnt ) ) * i_dl * grav * mask

    tauw = 0.5_rw * ( tt3d( c__5, 1, ilay )                           &
                    + tt3d( ipnt, 1, ilay ) ) * ramp
    vcor = 0.25_rw * ( h_v( ipnt, ilay )                              &
                     + h_v( c__3, ilay )                              &
                     + h_v( c__4, ilay )                              &
                     + h_v( c__5, ilay ) )
    ufor = fnud( ipnt, ilay, ix_u )                                   &
         + 0.5_rw * ( tt3d( ipnt, 2, ilay ) + tt3d( c__5, 2, ilay ) ) &
         * i_r1   * invf * i__h * ramp                                &
         + ramp * tide( 1, 1, ipnt, ix_u )                            &
         * cos(   tide( 2, 1, ipnt, ix_u ) - w_ti(1) * ctim )

    rhsi = dmd4 * (1._rw - gene)                                      &
         + 0.25_rw * pvor( ipnt ) * ( h_v( ipnt, ilay )               &
                                    + h_v( c__5, ilay ) )             &
         + 0.25_rw * pvor( c__3 ) * ( h_v( c__3, ilay )               &
                                    + h_v( c__4, ilay ) )             &
         + tauw            * i_r0 * i__h                              &
         - tb3d( ipnt, 1, ilay ) * i_r0 * i__h                        &
         - tu3d( ipnt, 1, ilay ) * i_r0 * i__h                        &
         + bodf( ilay, 1 )                                            &
         + ( del1 * dmd4                                              &
           + del2 * dmdx( 3, ipnt, ilay )                             &
           + gamm * dmdx( 2, ipnt, ilay )                             &
           + epsi * dmdx( 1, ipnt, ilay )                             &
           ) * gene

    if (svis>0._rw) then    
      rhsi = rhsi -  svis * i_dl* real( UU4(ipnt,ilay)                 &
          - UU4(c__5,ilay) + VV4(c__3,ilay) - VV4(ipnt,ilay))*i__h

    else   
      rhsi = rhsi +  ( v_cc( ipnt, ilay ) * dive( ipnt )             &
           - v_cc( c__5, ilay ) * dive( c__5 ) ) * i_dl               &
         - ( v_ll( c__3, ilay ) * rvor( c__3 )                        &
           - v_ll( ipnt, ilay ) * rvor( ipnt ) ) * i_dl               
    end if
    uold = uold + rhsi * mask * dt
     
    uold = ufor *          nudg( ipnt, ix_u )                         &
           + uold * (1._rw - nudg( ipnt, ix_u ))
   
    u( ipnt, ilay ) = uold

!   3rd-order upstream-biased advection of layer thickness.
!   See Shchepetkin & McWilliams, 2005, p.394.

    if (rgld < 0.5_rw) then
            h_u( ipnt, ilay ) = 0.5_rw * ( uold + abs( uold ) )      & ! max( u, 0. )
                                * ( hcen - 0.16667_rw * d2hx( c__5 ) ) &
                                + 0.5_rw * ( uold - abs( uold ) )      & ! min( u, 0. )
                                * ( hcen - 0.16667_rw * d2hx( ipnt ) )
    end if

    dmdx( 1, ipnt, ilay ) = dmdx( 2, ipnt, ilay )
    dmdx( 2, ipnt, ilay ) = dmdx( 3, ipnt, ilay )
    dmdx( 3, ipnt, ilay ) = dmd4
  end do
!$OMP END PARALLEL DO
end subroutine update_u

subroutine update_v( ilay )
  implicit none
  integer, intent(in) :: ilay
  integer             :: ipnt, c__1, c__7, c__8
  real ( rw )         :: ucor, tauw, vfor, mask, hcen, &
                         i_dl, i_r0, i_r1, i__h, rhsi, vold, dmd4
  i_dl = 1._rw / dl
  i_r0 = 1._rw / rho0
  i_r1 = 1._rw / rhon(1)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ipnt,c__1,c__7,c__8,tauw,ucor,vfor,mask,hcen,i__h,rhsi,vold,dmd4)
  do ipnt = 1, ndeg
    c__1 = neig( 1,    ipnt )
    c__7 = neig( 7,    ipnt )
    c__8 = neig( 8,    ipnt )
    mask = mk_v(       ipnt )
    hcen = real( hlay( ipnt, ilay ) + hlay( c__7, ilay ), rw ) / (1._rw + mask)
    
    
    i__h = 1._rw / ( hcen + 1._rw - mask )
    vold = v( ipnt, ilay )

    dmd4 = ( mont( c__7 ) - mont( ipnt ) ) * i_dl * grav * mask

    tauw = 0.5_rw * ( tt3d( c__7, 2, ilay )                           &
                    + tt3d( ipnt, 2, ilay ) ) * ramp
    ucor = 0.25_rw * ( h_u( ipnt, ilay )                              &
                     + h_u( c__1, ilay )                              &
                     + h_u( c__7, ilay )                              &
                     + h_u( c__8, ilay ) )
    vfor = fnud( ipnt, ilay, ix_v )                                   &
         - 0.5_rw * ( tt3d( ipnt, 1, ilay ) + tt3d( c__7, 1, ilay ) ) &
         * i_r1   * invf * i__h * ramp                                &
         + ramp * tide( 1, 1, ipnt, ix_v )                            &
         * cos(   tide( 2, 1, ipnt, ix_v ) - w_ti(1) * ctim )

    rhsi = dmd4 * (1._rw - gene)                                      &
         - 0.25_rw * pvor( ipnt ) * ( h_u( ipnt, ilay )               &
                                    + h_u( c__7, ilay ) )             &
         - 0.25_rw * pvor( c__1 ) * ( h_u( c__1, ilay )               &
                                    + h_u( c__8, ilay ) )             &
         + tauw                  * i_r0 * i__h                        &
         - tb3d( ipnt, 2, ilay ) * i_r0 * i__h                        &
         - tu3d( ipnt, 2, ilay ) * i_r0 * i__h                        &
         + bodf( ilay, 2 )                                            &
         + ( del1 * dmd4                                              &
           + del2 * dmdy( 3, ipnt, ilay )                             &
           + gamm * dmdy( 2, ipnt, ilay )                             &
           + epsi * dmdy( 1, ipnt, ilay )                             &
           ) * gene
    if (svis>0._rw) then         
      rhsi = rhsi -  svis * i_dl* ( VV4(c__1,ilay)                    &
       - VV4(ipnt,ilay) - UU4(ipnt,ilay) + UU4(c__7,ilay))*i__h      


    else   
       rhsi = rhsi+ ( v_cc( ipnt, ilay ) * dive( ipnt )               &
           - v_cc( c__7, ilay ) * dive( c__7 ) ) * i_dl               &
         + ( v_ll( c__1, ilay ) * rvor( c__1 )                        &
           - v_ll( ipnt, ilay ) * rvor( ipnt ) ) * i_dl               
    end if
   
    vold = vold + rhsi * mask * dt
    
    vold = vfor *          nudg( ipnt, ix_v )                         &
            + vold * (1._rw - nudg( ipnt, ix_v ))

    v( ipnt, ilay ) = vold

!   3rd-order upstream-biased advection of layer thickness.
!   See Shchepetkin & McWilliams, 2005, p.394.

    if (rgld < 0.5_rw) then
       h_v( ipnt, ilay ) = 0.5_rw * ( vold + abs( vold ) )        & ! max( v, 0. )
                           * ( hcen - 0.16667_rw * d2hy( c__7 ) ) &
                           + 0.5_rw * ( vold - abs( vold ) )      & ! min( v, 0. )
                           * ( hcen - 0.16667_rw * d2hy( ipnt ) )
    end if

    dmdy( 1, ipnt, ilay ) = dmdy( 2, ipnt, ilay )
    dmdy( 2, ipnt, ilay ) = dmdy( 3, ipnt, ilay )
    dmdy( 3, ipnt, ilay ) = dmd4
  end do  

!$OMP END PARALLEL DO

end subroutine update_v

subroutine update_h()
  implicit none
  integer    :: ipnt, ilay, c__1, c__3
  real( rw ) :: i_dl, hfor, rs_3, vecl( nlay ), rhsi, hfor1, hfor2
  real( r8 ) :: hold

  i_dl    = 1._rw / dl
  vecl(:) = 0._rw
  vecl(1) = 1._rw

          
  do ilay = nlay, 1, - 1
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ipnt,c__1,c__3,hold,hfor,rs_3,rhsi)
    do ipnt = 1, ndeg
      c__1 = neig( 1, ipnt )
      c__3 = neig( 3, ipnt )

      hold = hlay( ipnt, ilay )

      rs_3 = (                                                       &
               h_u( ipnt, ilay )                                     &
             - h_u( c__1, ilay )                                     &
             ) * i_dl                                                &
           + (                                                       &
               h_v( ipnt, ilay )                                     &
             - h_v( c__3, ilay )                                     &
             ) * i_dl                                                &
           + hdot(  ipnt, ilay )

      rs_3 = rs_3 * mk_n( ipnt )

      rhsi = ( (1.5_rw +         beta) * rs_3                        &
             - (0.5_rw + 2._rw * beta) * rs_h( 2, ipnt, ilay )       &
             +                   beta  * rs_h( 1, ipnt, ilay )       &
             ) * dt * gene                                           &
           + rs_3 * dt * (1._rw - gene)

      hold = hold + real( rhsi, r8 )

      hfor = fnud( ipnt, ilay, ix_n )                                &
           + ramp * tide( 1, 1, ipnt, ix_n ) * vecl( ilay )          &
           * cos(   tide( 2, 1, ipnt, ix_n ) - w_ti(1) * ctim )
            
      hlay ( ipnt, ilay ) = hold
! add water mass transformation to layer 1 and 2
      hfor1= 800._rw  ! fixed nudging top (not IC)
      hfor2= 0._rw  ! fixed nudging bot (not IC)
      if (hlay( ipnt, 2) > 20*hsal .and. subc(ipnt,1)>lm/2) then
         if (ilay == 1) then
           ! if (hfor1 > hlay (ipnt, ilay) ) then 
                hlay( ipnt, ilay) = hlay(ipnt, ilay) + 0*real(nudg(ipnt, ix_n), r8) &
                    +  max( real( hfor1 * nudg(ipnt, ix_n ), r8 )  + real( - nudg(ipnt, ix_n ), r8) * hlay(ipnt, ilay) , 0._rw )   
           !end if
        elseif (ilay ==2) then   
           ! if (hfor2 < hlay (ipnt, ilay) ) then    
                hlay( ipnt, ilay) = hlay(ipnt, ilay) - 0*real(nudg(ipnt, ix_n), r8) & 
                    +  min( real( hfor2 * nudg(ipnt, ix_n ), r8 )  + real( - nudg(ipnt, ix_n ), r8) * hlay(ipnt, ilay) , 0._rw )  
           ! end if
        end if
      end if
 ! both bot and top layers must be nudged
    if (subc(ipnt,1)<lm/2) then
         hlay( ipnt, ilay) = real( hfor * nudg( ipnt, ix_n ), r8 )   &
                       + real(1._rw - nudg( ipnt, ix_n ), r8 ) * hold
    !  elseif (ilay == 1 .and. subc(ipnt,1)> lm/2 .and. hfor1 > hlay(ipnt,ilay) ) then
    !        hlay( ipnt, ilay ) = real( hfor1  * nudg( ipnt, ix_n ), r8 )    &
    !                     + real(1._rw -  nudg( ipnt, ix_n ), r8 ) * hlay(ipnt, ilay) 
    !  elseif (ilay == 2 .and. subc(ipnt,1)> lm/2 .and. hfor2 < hlay(ipnt, ilay) ) then
    !        hlay( ipnt, ilay ) = real( hfor2 * nudg( ipnt, ix_n ), r8 )    &
    !                     + real(1._rw - nudg( ipnt, ix_n ), r8 ) * hlay(ipnt, ilay)
      end if
    
      rs_h( 1, ipnt, ilay ) = rs_h( 2, ipnt, ilay )
      rs_h( 2, ipnt, ilay ) = rs_3
    end do
!$OMP END PARALLEL DO
  end do
  
  if (rgld>0.5_rw) then
          do ipnt=1, ndeg

          !tiny correction to top layer to ensure sum of layers is htop-hbot to machine precision

          hlay(ipnt, 1) = hlay(ipnt,1)-0.5*real(sum(hlay(ipnt,:)) - real((h_th(ipnt)),r8))
          hlay(ipnt, 2) = hlay(ipnt,2)-0.5*real(sum(hlay(ipnt,:)) - real((h_th(ipnt)),r8))

!if (ipnt == 50) then
!print *, hlay(50,1)
!print *, hlay(50,2)
!print *, sum(hlay(ipnt,:))
!end if


!     hlay(ipnt, 1) = hlay(ipnt,1)-real(hlay(ipnt,1)/sum(hlay(ipnt,:)))*(real(sum(hlay(ipnt,:)) - real(h_bo(ipnt),r8)))

 !    hlay(ipnt, 2) = hlay(ipnt,2)-real(hlay(ipnt,2)/sum(hlay(ipnt,:)))*(real(sum(hlay(ipnt,:)) - real(h_bo(ipnt),r8)))


!     hlay(ipnt, 1) = hlay(ipnt,1)-.9*(real(sum(hlay(ipnt,:)) - real(h_bo(ipnt),r8)))

 !  hlay(ipnt, 2) = hlay(ipnt,2)-.1*(real(sum(hlay(ipnt,:)) - real(h_bo(ipnt),r8)))
  ! if (hlay(ipnt,1)>700 .or. hlay(ipnt,1)<0) then
  !         print *, 'eee', hlay(ipnt,1)
  ! end if
   
   !if (hlay(ipnt,2)>700 .or. hlay(ipnt,2)<0) then
   !        print *, 'dddd', hlay(ipnt,2)
   !end if
   
 !if (sum(hlay(ipnt,:)) > 700.00001) then
  !       print *, 'aaaa', sum(hlay(ipnt,:))
  !       print *, 'bbbb', hlay(ipnt, 1)
  !       print *, 'cccc', hlay(ipnt,2)
  !       !read(*,*) ilay
 !end if    



               !do ilay= nlay,1, -1                 
               !     hlay(ipnt, ilay) = hlay(ipnt,ilay) - hlay(ipnt,ilay)/real(sum(hlay(ipnt,:)))   &
               !			*real(sum(hlay(ipnt,:)) - real(h_bo(ipnt),r8))
               !end do
               if (sum(real(hlay(ipnt, :),r8)) > h_th(ipnt)) then
                       errm = trim(errm) //' The water column depth is greater than htop-hbot\n '//  &
                               'from module private_mod.f95,'
               elseif (sum(real(hlay(ipnt, :),r8)) < h_th(ipnt)) then
                       errm = trim(errm) //' The water column depth is less than htop-hbot\n '//  &
                               'from module private_mod.f95,'        
               end if
               end do
  end if
  
end subroutine update_h


subroutine surf_pressure()
  implicit none
  real(rw)     :: rp, maxdiff, temp, diff, pi_tol, pi_rhs(1:ndeg), pi_prev(1:ndeg)             
  integer     :: iters, i, j,k,im1,ip1,jp1,jm1,maxiters, c__1, c__3, c__5, c__7, ilay, ipnt
  logical     :: hasConverged
  
  pi_rhs(: ) = 0._rw
  rp=1.000
  iters=0
  diff=0
  hasConverged=.true.
  pi_tol=1.e-5_rw
  maxiters=1000
  

!  Set right-hand side of Poisson equation

  do ilay = nlay, 1, -1
  
    ! Add contribution due to x-volume fluxes
    do ipnt=1, ndeg
        if (subc(ipnt,1)>1) then
        c__5= neig(5,ipnt)
         
       
        pi_rhs(ipnt) = pi_rhs(ipnt)- h_u(ipnt,ilay) / (dl*dt)
        pi_rhs(c__5) = pi_rhs(c__5)+ h_u(ipnt,ilay) / (dl*dt)

        end if
    end do
    

    ! Add contribution due to y-volume fluxes
    do ipnt=1, ndeg


      if (subc(ipnt,2)>1) then
        c__7=neig(7,ipnt)

        ! already calculated h_v = v(ipnt,ilay)*h_south(ipnt,ilay)
        pi_rhs(ipnt) =pi_rhs(ipnt) - h_v(ipnt,ilay) / (dl*dt)
        pi_rhs(c__7) =pi_rhs(c__7)+ h_v(ipnt,ilay) / (dl*dt)
        
        
      end if
    end do
    
  end do

  
  !  Perform SOR iteration
  maxdiff = pi_tol + 1
  iters = 0
  do while (maxdiff > pi_tol .and. iters < maxiters)
  
    !print *, 'The maxdiff is', maxdiff
    !print *, 'The iteration is', iters
    
    maxdiff = 0
    
   ! PSOR default
      do ipnt=1, ndeg

          ! Store current grid value of pi_s
          pi_prev(ipnt) = pi_s(ipnt)
          pi_s(ipnt) = (1-rp)*pi_s(ipnt) - rp * Osum_(ipnt) * pi_rhs(ipnt)
 
      
          if (subc(ipnt,1)<lm) then
                  c__1= neig(1, ipnt)
                  pi_s(ipnt) = pi_s(ipnt)+ rp * Osum_(ipnt) * Ow(c__1)*pi_s(c__1)
          end if
          if (subc(ipnt,2)<mm) then
                  c__3= neig(3, ipnt)
                  pi_s(ipnt) = pi_s(ipnt)+ rp * Osum_(ipnt) * Os(c__3)*pi_s(c__3)
          end if
          if (subc(ipnt,1)>1) then
                  c__5= neig(5, ipnt)
                  pi_s(ipnt) = pi_s(ipnt)+ rp * Osum_(ipnt) * Ow(ipnt)*pi_s(c__5)
          end if
          if (subc(ipnt,2)>1) then
                  c__7= neig(7, ipnt)
                  pi_s(ipnt) = pi_s(ipnt)+ rp * Osum_(ipnt) * Os(ipnt)*pi_s(c__7)
          end if

          !  operator Os, Ow are defined to be 0 at walls and H/dl^2 elsewhere

          end do
      
      !  Calculate the absolute difference between iterations (pointwise convergence)
      do ipnt=1, ndeg
        diff = abs(pi_s(ipnt)-pi_prev(ipnt))
        if (diff > maxdiff) then
          maxdiff = diff
        end if
      end do

    iters=iters+1
  end do

  
  !  Correct u-velocity
  do ilay=1, nlay
  
    
    do ipnt = 1, ndeg
        
        if (subc(ipnt,1)>1 .and. subc(ipnt,1)<lm+1) then
                c__5= neig( 5, ipnt)
                u(ipnt,ilay) =u(ipnt,ilay) - dt/dl*pi_s(ipnt)
                u(ipnt,ilay) =real(u(ipnt,ilay) +dt/dl*pi_s(c__5), r8)
        end if
    end do
    
  end do
  
  
  !   Correct v-velocity
  do ilay=1, nlay
    
    do ipnt = 1, ndeg
    if (subc(ipnt,2)>1 .and. subc(ipnt,2)<mm+1) then
        c__7 = neig(7, ipnt)
        v(ipnt,ilay) = v(ipnt,ilay)- dt/dl*pi_s(ipnt)
        v(ipnt, ilay) = real(v(ipnt, ilay) + dt/dl* pi_s(c__7), r8)
    end if
    end do
    
  end do
  
 ! insert statement to state if haven't converged within max iterations?
        
! print *, maxval(pi_s)  
end subroutine surf_pressure

subroutine integrate_time()
  implicit none
  real   ( r8 ) :: dtd8
  integer       :: lerm, nstp, notp, n_3d, tstp
  logical       :: upst

  lerm = len_trim( errm )
  errm = trim(     errm ) &
      // ' in subroutine integrate_time from module private_mod.f95,'

  write( ioso, * ) 'dl = ', dl, ' meters.'
  write( ioso, * ) 'dt = ', dt, ' seconds.'

  dtd8 = real( dt, r8 ) / 24._r8 / 3600._r8
  nstp = nint(       real( dt_s, r8 ) / dtd8 )
  notp = max(  nint( real( dt_o, r8 ) / dtd8 ), 1 )
  n_3d = max(  nint( real( dt3d, r8 ) / dtd8 ), 1 ) ! Every `dt3d'.

  ramp = 1._rw ! Wind forcing is ramped over a period of `dt_r' days.
  gene = 0._rw ! Generalized Forward-Backward cannot be activated until
               !   at least 3 time-steps are completed.
  tstp = 1
  ctim = real( tres + dtd8 * real(tstp, r8), rw )
  call distribute_stress()
  if ( rsta < 0.5_rw .and. ctim < dt_r ) then
    ramp = ctim / dt_r
  end if
  call first_three_timesteps( tstp )

  tstp = 2
  ctim = real( tres + dtd8 * real(tstp, r8), rw )
  call first_three_timesteps( tstp )

  tstp = 3
  ctim = real( tres + dtd8 * real(tstp, r8), rw )
  call first_three_timesteps( tstp )

  gene = g_fb ! Activate Generalized Forward-Backward
              !   (if it has been selected in shared_mod.f95).
              
  if (gene > 0.5_rw .and. rgld > 0.5_rw) then
        gene = 0._rw
        errm = trim(errm) //' Cant use multistep method with rigid lid\n '//  &
                   'from module private_mod.f95,'
  end if
              
  do tstp = 4, nstp
    ctim = real( tres + dtd8 * real(tstp, r8), rw )

    upst = .false.
    if ( mod( tstp, n_3d ) == 0 ) then
      upst = .true.
    end if

    if ( upst ) then
      call distribute_stress()
    end if

    ramp = 1._rw
    if ( rsta < 0.5_rw .and. ctim < dt_r ) then
      ramp = ctim / dt_r
    end if

!   Generalized Forward-Backward time-stepping.
!   Shchepetkin & McWilliams (2005, Ocean Modelling, Eq. 2.49).

    call gener_forward_backward( tstp, upst )

    if ( mod( tstp, notp ) == 0 ) then
      call write_outputs()
    end if

  end do

  if ( errc == 0 ) then
    errm = errm(1 : lerm)
  else
    call quit()
  end if
end subroutine integrate_time

subroutine distribute_stress()
  implicit none
  integer    :: ilay, ipnt, klay, c__1, c__3, c__4, c__5, c__7, c__8
  real( rw ) :: hcum, sofar, uatv, vatu, &
                layt( 0 : ndeg, nlay ),  &
                layb( 0 : ndeg, nlay ),  &
                layu( 0 : ndeg, nlay ),  &
                taub( 0 : ndeg, 2    ),  &
                taum( 0 : ndeg, 2    )
  real( r8 ) :: hs_8

  hs_8 = real( hsal, r8 )

! Corner cases to watch for:
! - Grid points shallower than hsbl/hbbl;
! - Grid points that are dry.
! Following Salmon 2002, we allow wind stress to be transmitted to
!   bedrock wherever sum(h(1 : nlay)) < hsbl.
! More precisely, wind stress is only applied over the top hsbl meters.
!   Layers that are thinner than `hsal' are considered `outcropped' and
!   no wind stress is applied to them.
! At the bottom, no bottom stress is applied to the ocean wherever
!   sum(h(1 : nlay)) -> 0.

  if ( any( abs(taus(:,:)) > 1.e-7_rw ) .and. ocrp > 0.5_rw ) then
    do ilay = 1, nlay ! Top -> Down
!$OMP PARALLEL DO PRIVATE(ipnt,hcum,sofar,klay)
      do ipnt = 0, ndeg
        layt(ipnt, ilay) = 0._rw ! Initialization.
        hcum             = 0._rw
        sofar            = sum( layt(ipnt, 1 : ilay) ) ! 0. <= sofar <= 1.
        do klay = 1, ilay
          hcum = hcum + max(0._rw, real( hlay(ipnt, klay), rw ) - 1.5_rw * hsal)
        end do
        layt(ipnt, ilay) = min( hcum, hsbl ) / hsbl - sofar
        layt(ipnt, ilay) = max( layt(ipnt, ilay), 0._rw ) ! Round-off errors.
      end do
!$OMP END PARALLEL DO
    end do
  elseif ( any( abs(taus(:,:)) > 1.e-7_rw )) then
!$OMP PARALLEL DO PRIVATE(ipnt)
    do ipnt = 0, ndeg
      layt( ipnt, : ) = 0._rw
      layt( ipnt, 1 ) = 1._rw
    end do
!$OMP END PARALLEL DO
  end if

  if ( bdrg > 1.e-7_rw .and. ocrp > 0.5_rw ) then ! Save some cycles.
    do ilay = nlay, 1, - 1 ! Bottom -> Up
!$OMP PARALLEL DO PRIVATE(ipnt,hcum,sofar)
      do ipnt = 0, ndeg
        layb(ipnt, ilay) = 0._rw ! Initialization.
        sofar            = sum( layb(ipnt, ilay : nlay) ) ! 0. <= sofar <= 1.
        hcum             = sum( real( hlay(ipnt, ilay : nlay), rw ) ) ! >= 0.
        layb(ipnt, ilay) = min( hcum, hbbl ) / hbbl - sofar
        layb(ipnt, ilay) = max( layb(ipnt, ilay), 0._rw ) ! Round-off errors.
      end do
!$OMP END PARALLEL DO
    end do
  elseif ( bdrg > 1.e-7_rw ) then
!$OMP PARALLEL DO PRIVATE(ipnt)
    do ipnt = 0, ndeg
      layb( ipnt, :    ) = 0._rw
      layb( ipnt, nlay ) = 1._rw
    end do
!$OMP END PARALLEL DO
  end if


  if ( tdrg > 1.e-7_rw  .and. ocrp > 0.5_rw ) then ! Save some cycles
    do ilay = 1, nlay ! Top -> Down
!$OMP PARALLEL DO PRIVATE(ipnt,hcum,sofar,klay)
      do ipnt = 0, ndeg
        layu(ipnt, ilay) = 0._rw ! Initialization.
        hcum             = 0._rw
        sofar            = sum( layu(ipnt, 1 : ilay) ) ! 0. <= sofar <= 1.
        do klay = 1, ilay
          hcum = hcum + max(0._rw, real( hlay(ipnt, klay), rw ) - 1.5_rw * hsal)
        end do
        layu(ipnt, ilay) = min( hcum, hsbl ) / hsbl - sofar
        layu(ipnt, ilay) = max( layu(ipnt, ilay), 0._rw ) ! Round-off errors.
      end do
!$OMP END PARALLEL DO
    end do
  elseif ( tdrg > 1.e-7_rw ) then
!$OMP PARALLEL DO PRIVATE(ipnt)
    do ipnt = 0, ndeg
      layu( ipnt, : ) = 0._rw
      layu( ipnt, 1 ) = 1._rw
    end do
!$OMP END PARALLEL DO
  end if

  if ( bdrg > 1.e-7_rw ) then ! Need to calculate/update bottom stress.
!$OMP PARALLEL DO PRIVATE(ipnt,ilay,klay,c__1,c__3,c__4,c__5,c__7,c__8,vatu,uatv)
    do ipnt = 0, ndeg
      ilay = nlay ! Default layer for bottom stress.
      if ( ocrp > 0.5_rw ) then
!       Bottom stress is calculated from the layer that is closest to the
!         bottom while having a thickness > (2 * hsal).
        do klay = nlay, 1, - 1 ! Bottom -> Up.
          if ( hlay( ipnt, klay ) > (0._r8 * hs_8) ) then
            ilay = klay
            exit
          end if
        end do
      end if
      c__1 = neig( 1, ipnt )
      c__3 = neig( 3, ipnt )
      c__4 = neig( 4, ipnt )
      c__5 = neig( 5, ipnt )
      c__7 = neig( 7, ipnt )
      c__8 = neig( 8, ipnt )
      vatu = 0.25_rw * v( ipnt, ilay ) &
           + 0.25_rw * v( c__3, ilay ) &
           + 0.25_rw * v( c__4, ilay ) &
           + 0.25_rw * v( c__5, ilay )
      uatv = 0.25_rw * u( ipnt, ilay ) &
           + 0.25_rw * u( c__1, ilay ) &
           + 0.25_rw * u( c__7, ilay ) &
           + 0.25_rw * u( c__8, ilay )
      taub( ipnt, 1 ) = u( ipnt, ilay ) * bdrg * rhon( ilay )         &
                      * ( qdrg * sqrt( u( ipnt, ilay )**2 + vatu**2 ) &
                        + 1._rw - qdrg )
      taub( ipnt, 2 ) = v( ipnt, ilay ) * bdrg * rhon( ilay )         &
                      * ( qdrg * sqrt( v( ipnt, ilay )**2 + uatv**2 ) &
                        + 1._rw - qdrg )
    end do
  
!$OMP END PARALLEL DO
!print *, 'min h ' , minval(hlay(: ,2)), minloc(hlay(:,2)) , klay       
!print *, 'minbottaub=', minval(taub(:,2)), 'minbotv',minval(v(:,2)), 'bdrg*rhon', bdrg*rhon(2)
!print *, 'first five terms of taub v', taub(1:5,2)

    do ilay = 1, nlay
!$OMP PARALLEL DO PRIVATE(ipnt,c__5,c__7)
      do ipnt = 1, ndeg
        c__5 = neig( 5, ipnt )
        c__7 = neig( 7, ipnt )
!       tb3d, and taub are defined at u/v points..
!       layb is cell-centered (needs interpolation).
        tb3d( ipnt, 1, ilay ) =   taub( ipnt, 1    ) * 0.5_rw &
                              * ( layb( ipnt, ilay )          &
                                + layb( c__5, ilay ) )
        tb3d( ipnt, 2, ilay ) =   taub( ipnt, 2    ) * 0.5_rw &
                              * ( layb( ipnt, ilay )          &
                                + layb( c__7, ilay ) )
      end do
!$OMP END PARALLEL DO
    end do
  end if

 ! print *, '1111', maxval(taum(:,1)), maxval(taum(:,2)), maxval(taub(:,1)), maxval(taub(:,2))
  if ( tdrg > 1.e-7_rw ) then ! Need to calculate/update top stress.
!$OMP PARALLEL DO PRIVATE(ipnt,ilay,klay,c__1,c__3,c__4,c__5,c__7,c__8,vatu,uatv)
    do ipnt = 0, ndeg
      ilay = 1 ! Default layer for top stress.
      if ( ocrp > 0.5_rw ) then
!       Top stress is calculated from the layer that is closest to the
!         surface while having a thickness > (2 * hsal).
        do klay = 1, nlay ! Top -> Down.
          if ( hlay( ipnt, klay ) > (2._r8 * hs_8) ) then
            ilay = klay
            exit
          end if
        end do
      end if
      c__1 = neig( 1, ipnt )
      c__3 = neig( 3, ipnt )
      c__4 = neig( 4, ipnt )
      c__5 = neig( 5, ipnt )
      c__7 = neig( 7, ipnt )
      c__8 = neig( 8, ipnt )
      vatu = 0.25_rw * v( ipnt, ilay ) &
           + 0.25_rw * v( c__3, ilay ) &
           + 0.25_rw * v( c__4, ilay ) &
           + 0.25_rw * v( c__5, ilay )
      uatv = 0.25_rw * u( ipnt, ilay ) &
           + 0.25_rw * u( c__1, ilay ) &
           + 0.25_rw * u( c__7, ilay ) &
           + 0.25_rw * u( c__8, ilay )
      taum( ipnt, 1 ) = u( ipnt, ilay ) * tdrg * rhon( ilay )         &
                      * ( qdrg * sqrt( u( ipnt, ilay )**2 + vatu**2 ) &
                        + 1._rw - qdrg )
      taum( ipnt, 2 ) = v( ipnt, ilay ) * tdrg * rhon( ilay )         &
                      * ( qdrg * sqrt( v( ipnt, ilay )**2 + uatv**2 ) &
                        + 1._rw - qdrg )
    end do
!$OMP END PARALLEL DO
!print *, '2222', maxval(taum(:,1)) , maxval(taum(:,2))
!print *, 'maxtoptaum=', maxval(taum(:,1)), 'maxtopv',maxval(v(:,1)), 'tdrg*rhon(top)', tdrg*rhon(1)
! print *, 'frist five terms of taum v', taum(1:5, 2)

!
    do ilay = 1, nlay
!$OMP PARALLEL DO PRIVATE(ipnt,c__5,c__7)
      do ipnt = 1, ndeg
        c__5 = neig( 5, ipnt )
        c__7 = neig( 7, ipnt )
!       tu3d, taum are defined at u/v points..
!       layu is cell-centered (needs interpolation).

        tu3d( ipnt, 1, ilay ) =   taum( ipnt, 1    ) * 0.5_rw &
                              * ( layu( ipnt, ilay )          &
                                + layu( c__5, ilay ) )
        tu3d( ipnt, 2, ilay ) =   taum( ipnt, 2    ) * 0.5_rw &
                              * ( layu( ipnt, ilay )          &
                                + layu( c__7, ilay ) )

      end do
!$OMP END PARALLEL DO
    end do
  end if

  if ( any( abs(taus(:,:)) > 1.e-7_rw ) ) then
    do ilay = 1, nlay
!$OMP PARALLEL DO PRIVATE(ipnt)
      do ipnt = 1, ndeg
!       tt3d, layt and taus are all cell-centered.
        tt3d( ipnt, 1, ilay ) = taus( ipnt, 1 ) * layt( ipnt, ilay )
        tt3d( ipnt, 2, ilay ) = taus( ipnt, 2 ) * layt( ipnt, ilay )
      end do
!$OMP END PARALLEL DO
    end do
  end if

 ! print *, 'maxtoptu3d and minbottb3d',  maxval(tu3d(:,1,1)), minval(tb3d(:,1,2))
end subroutine distribute_stress

subroutine first_three_timesteps( tstp )
  implicit none
  integer, intent( in ) :: tstp
  integer               :: ilay, ipnt, c__5, c__7

! Standard Forward-Backward scheme
!   Shchepetkin & McWilliams (2005, Ocean Modelling), Eq. 2.13.
!   This scheme is used for the first three time-steps.
!   It is not parallelized.
!   After the first three time-steps, the code keeps using the standard
!     scheme (g_fb=0.) or switches to the "Generalized" FB scheme if
!     g_fb is set to 1. in file shared_mod.f95.

! Make sure that h*u,h*v is up-to-date in prevision of forward step.

  do ilay = 1, nlay
    do ipnt = 1, ndeg
      c__5 = neig( 5, ipnt )
      c__7 = neig( 7, ipnt )
      h_u( ipnt, ilay ) = u( ipnt, ilay ) * real( hlay( ipnt, ilay )       &
                                                + hlay( c__5, ilay ), rw ) &
                        / ( 1._rw + mk_u( ipnt ) )
      h_v( ipnt, ilay ) = v( ipnt, ilay ) * real( hlay( ipnt, ilay )       &
                                                + hlay( c__7, ilay ), rw ) &
                        / ( 1._rw + mk_v( ipnt ) )
    end do
  end do

! Forward step: eta is updated to time-level `n+1', using velocities from `n'.

  call update_h()

  do ilay = 1, nlay
!   Update Montgomery potential to use the updated eta(t=n+1).
!   Also, calculate auxiliary fields required for update of u,v.

    call update_mont_rvor_pvor_dive_kine( ilay )
    call update_viscosity( ilay )

!   Backward step: u,v are updated to `n+1' using eta(t=n+1).
!   Alternate the order in which velocity components are updated.
!     (Bleck and Smith JGR 1990 vol.95 no.C3, p.3276).
        if ( mod( tstp, 2 ) == 0 ) then ! If n is even, ...
                call update_u( ilay )
                call update_v( ilay )
        else                            ! If n is odd,  ...
                call update_v( ilay )
                call update_u( ilay )
        end if
    
    if ( flag_nudging .and. mcbc < 0.5_rw ) then
!     Apply vanishing normal derivative at non-periodic open boundaries.
      call no_gradient_obc( ilay )
    end if
  end do
  
  if (rgld > 0.5_rw) then
  do ilay = 1, nlay
    do ipnt = 1, ndeg
      c__5 = neig( 5, ipnt )
      c__7 = neig( 7, ipnt )
      h_u( ipnt, ilay ) = u( ipnt, ilay ) * real( hlay( ipnt, ilay )       &
                                                + hlay( c__5, ilay ), rw ) &
                        / ( 1._rw + mk_u( ipnt ) )
      h_v( ipnt, ilay ) = v( ipnt, ilay ) * real( hlay( ipnt, ilay )       &
                                                + hlay( c__7, ilay ), rw ) &
                        / ( 1._rw + mk_v( ipnt ) )
    end do
  end do
    call surf_pressure()
  end if
  
end subroutine first_three_timesteps

subroutine gener_forward_backward( tstp, upst )
  implicit none
  logical, intent( in ) :: upst
  integer, intent( in ) :: tstp
  integer               :: ilay, ipnt,  c__5, c__7
  real ( rw )         :: mask, hcen
! `Generalized' Forward-Backward
! Shchepetkin & McWilliams (2005, Ocean Modelling, Eq. 2.49)
! The `standard' scheme is used if g_fb=0. inside file shared_mod.f95.


!This is the proper place to update hu, hv for a rigid lid
  if (rgld > 0.5_rw) then
      do ilay = 1, nlay
          do ipnt = 1, ndeg
          c__5 = neig( 5, ipnt )
          c__7 = neig( 7, ipnt )
          mask = mk_u(       ipnt )
          hcen = real( hlay( c__5, ilay ) + hlay( ipnt, ilay ), rw ) / (1._rw + mask)
          h_u( ipnt, ilay ) = 0.5_rw * ( u(ipnt,ilay) + abs( u(ipnt,ilay) ) )      & ! max( u, 0. )
                             * ( hcen - 0.16667_rw * d2hx( c__5 ) ) &
                             + 0.5_rw * ( u(ipnt,ilay) - abs( u(ipnt,ilay) ) )      & ! min( u, 0. )
                             * ( hcen - 0.16667_rw * d2hx( ipnt ) )
          mask = mk_v(       ipnt )
          hcen = real( hlay( c__7, ilay ) + hlay( ipnt, ilay ), rw ) / (1._rw + mask)
          h_v( ipnt, ilay ) = 0.5_rw * ( v(ipnt,ilay) + abs( v(ipnt,ilay) ) )      & ! max( v, 0. )
                             * ( hcen - 0.16667_rw * d2hy( c__7 ) ) &
                             + 0.5_rw * ( v(ipnt,ilay) - abs( v(ipnt,ilay) ) )      & ! min( v, 0. )
                             * ( hcen - 0.16667_rw * d2hy( ipnt ) )

          end do
      end do
  end if
  
  call update_h() ! 14%

  do ilay = 1, nlay

!   Update Montgomery potential to use the updated eta(t=n+1).
!   Also, calculate auxiliary fields required for the update of u,v.

    call update_mont_rvor_pvor_dive_kine( ilay ) ! 19%

    if ( (dvis > 1.e-3_rw .and. upst) .or. svis > 0 ) then
      call update_viscosity( ilay )              ! 10%
    end if

!   Backward step: u(n) becomes u(n+1), using eta(n+1).
!   Alternate the order in which velocity components are updated.
!     (Bleck and Smith JGR 1990 vol.95 no.C3, p.3276).

    if ( mod( tstp, 2 ) == 0 ) then ! If n is even, ...
            call update_u( ilay ) ! 28%
            call update_v( ilay ) ! 28%
    else                            ! If n is odd,  ...
            call update_v( ilay )
            call update_u( ilay )
    end if


    if ( flag_nudging .and. mcbc < 0.5_rw ) then
!     Apply vanishing normal derivative at non-periodic open boundaries.
      call no_gradient_obc( ilay )
    end if

  end do ! loop on layers
  
  if (rgld > 0.5_rw) then
       do ilay = 1, nlay
           do ipnt = 1, ndeg
                c__5 = neig( 5, ipnt )
                c__7 = neig( 7, ipnt )
                mask = mk_u(       ipnt )
                hcen = real( hlay( c__5, ilay ) + hlay( ipnt, ilay ), rw ) / (1._rw + mask)
                h_u( ipnt, ilay ) = 0.5_rw * ( u(ipnt,ilay) + abs( u(ipnt,ilay) ) )      & ! max( u, 0. )
                                   * ( hcen - 0.16667_rw * d2hx( c__5 ) ) &
                                   + 0.5_rw * ( u(ipnt,ilay) - abs( u(ipnt,ilay) ) )      & ! min( u, 0. )
                                   * ( hcen - 0.16667_rw * d2hx( ipnt ) )
                mask = mk_v(       ipnt )
                hcen = real( hlay( c__7, ilay ) + hlay( ipnt, ilay ), rw ) / (1._rw + mask)
                h_v( ipnt, ilay ) = 0.5_rw * ( v(ipnt,ilay) + abs( v(ipnt,ilay) ) )      & ! max( v, 0. )
                                   * ( hcen - 0.16667_rw * d2hy( c__7 ) ) &
                                   + 0.5_rw * ( v(ipnt,ilay) - abs( v(ipnt,ilay) ) )      & ! min( v, 0. )
                                   * ( hcen - 0.16667_rw * d2hy( ipnt ) )
            end do
        end do

        call surf_pressure()

  end if
  
end subroutine gener_forward_backward

subroutine update_mont_rvor_pvor_dive_kine( ilay )
  implicit none
  integer, intent( in ) :: ilay
  real( rw ), parameter :: i_dl = 1._rw / dl, &
                           i_gr = 1._rw / grav
  real( rw ) :: i_ns, u_le,                   &
                u_ri, v_bo, v_to, have, hloc
  real( r8 ) :: hcol, mpot, i_rn( nlay ), hs_8
  integer    :: ipnt, c__1, c__3, c__5, c__6, c__7, i

  i_ns = 1._rw / real( nsal - 1, rw )
  i_rn = 1._r8 / real( rhon(:),  r8 )
  hs_8 = real( hsal, r8 )

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ipnt,c__1,c__3,c__5,c__6,c__7,i,u_le,u_ri,v_bo,v_to,have,hcol,hloc,mpot)
  do ipnt = 1, ndeg
    c__1 = neig( 1, ipnt ) ! Right
    c__3 = neig( 3, ipnt ) ! Top
    c__5 = neig( 5, ipnt )
    c__6 = neig( 6, ipnt )
    c__7 = neig( 7, ipnt )

    u_le = u( ipnt, ilay )
    u_ri = u( c__1, ilay )
    v_bo = v( ipnt, ilay )
    v_to = v( c__3, ilay )

    hloc = real( hlay( ipnt, ilay ), rw ) ! `Local' layer thickness.
 
!   Update `mont' (the Montgomery potential divided by gravity) to time n+1.

!   1- Salmon's contribution (if `outcrop' is activated).

    mpot = hlay( ipnt, ilay ) + real( hmin * ( 1._rw - mk_n( ipnt ) ), r8 )
    mpot = ( real( hsal, r8 ) / mpot )**(nsal - 1)
    mpot = mpot * real( - ocrp * i_ns * hsal * mk_n( ipnt ), r8 )

!   2- Baroclinic contribution.
 mpot=mpot-h_to(ipnt)
    do i = 1, ilay - 1
      mpot = mpot                                                &
           - real( rhon( ilay ) - rhon( i ), r8 ) * i_rn( ilay ) &
           * hlay( ipnt, i )
    end do

!   3- Barotropic contribution.

    if (rgld<0.5_rw) then

      hcol = 0._r8

      do i = 1, nlay
         hcol = hcol + hlay( ipnt, i )
      end do

      mpot = hcol - real( h_th(ipnt), r8 ) + mpot

    end if
!   4- Kinetic energy (Montgomery potential becomes Bernoulli potential).
!      kine = 0.5 * ave_x(u**2) + 0.5 * ave_y(v**2), defined at cell-center.
!        e.g., Arakawa & Lamb MWR 1981, Sadourny MWR 1975.

    mont( ipnt ) = real( mpot, rw )              &
                 + 0.25_rw * uadv * i_gr         &
                 * ( u_ri * u_ri + u_le * u_le   & ! Arithmetic mean.
                   + v_to * v_to + v_bo * v_bo )   ! Meters.

!   Relative vorticity along land boundaries is set to zero (Free Slip, FSM)
!     (see Ketefian & Jacobson 2009, and Dupont et al. 2003).

    rvor( ipnt ) = ( v_bo - v( c__5, ilay )   &
                   - u_le + u( c__7, ilay ) ) * i_dl * mkpe( ipnt ) ! s**(-1).

!   Second spatial derivative of layer thickness. Zero unless all cells wet.
!     (from Shchepetkin and McWilliams, 2005, p.394).

    d2hx( ipnt ) = real( hlay( c__1, ilay )               & ! Meters.
                       + hlay( c__5, ilay )               &
                       - hlay( ipnt, ilay ) * 2._r8, rw ) &
                 * mk_n( c__1 ) * mk_n( c__5 ) * mk_n( ipnt )

    d2hy( ipnt ) = real( hlay( c__3, ilay )               & ! Meters.
                       + hlay( c__7, ilay )               &
                       - hlay( ipnt, ilay ) * 2._r8, rw ) &
                 * mk_n( c__3 ) * mk_n( c__7 ) * mk_n( ipnt )

    if ( ocrp > 0.5_rw ) then ! Do not interfere with Salmon's method.
!     Consider x,y dimensions separately (for cases with 2-D domain).
      if ( hlay( c__1, ilay ) < 2._r8 * hs_8 .or.         &
           hlay( c__5, ilay ) < 2._r8 * hs_8 .or.         &
           hlay( ipnt, ilay ) < 2._r8 * hs_8 ) then
        d2hx( ipnt ) = 0._rw
      end if
      if ( hlay( c__3, ilay ) < 2._r8 * hs_8 .or.         &
           hlay( c__7, ilay ) < 2._r8 * hs_8 .or.         &
           hlay( ipnt, ilay ) < 2._r8 * hs_8 ) then
        d2hy( ipnt ) = 0._rw
      end if
    end if

!   Thickness along land boundaries is defined as in
!     Ketefian & Jacobson 2009 J.Comput.Phys. (their Appendix A).

    have = real( hlay( ipnt, ilay ) &
               + hlay( c__5, ilay ) &
               + hlay( c__6, ilay ) &
               + hlay( c__7, ilay ), rw )
  
    pvor( ipnt ) = ( fcor( ipnt )                         &
                   + rvor( ipnt ) * uadv )                &
                 * mkpi(   ipnt )                         &
                 * ( mk_n( ipnt )                         &
                   + mk_n( c__5 )                         &
                   + mk_n( c__6 )                         &
                   + mk_n( c__7 ) )                       &
                 / have             ! m**(-1) s**(-1).

    dive( ipnt ) = ( u_ri - u_le &
                   + v_to - v_bo ) * i_dl
  end do
!$OMP END PARALLEL DO
end subroutine update_mont_rvor_pvor_dive_kine

subroutine update_viscosity( ilay )
  implicit none
  integer, intent( in ) :: ilay
  integer    :: ipnt, c__1, c__2, c__3, c__5, c__6, c__7
  real( rw ) :: r_bl, r_br, r_tr, r_tl, rbll, rbbl, &
                d_cc, d_ri, d_to, d_le, d_bl, d_bo, hh_q

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ipnt,c__1,c__2,c__3,c__5,c__6), &
!$OMP& PRIVATE(c__7,r_bl,r_br,r_tr,r_tl,rbll,rbbl,d_cc,d_ri,d_to,d_le,d_bl,d_bo)
  do ipnt = 1, ndeg
    c__1 = neig( 1, ipnt )
    c__2 = neig( 2, ipnt )
    c__3 = neig( 3, ipnt )
    c__5 = neig( 5, ipnt )
    c__6 = neig( 6, ipnt )
    c__7 = neig( 7, ipnt )

    r_bl = rvor( ipnt )  ! rvor at bottom-left  corner of cell.
    r_br = rvor( c__1 )  !         bottom-right
    r_tr = rvor( c__2 )  !         top   -right
    r_tl = rvor( c__3 )  !         top   -left
    rbll = rvor( c__5 )  !         bottom-left-left
    rbbl = rvor( c__7 )  !         bottom-bottom-left

    d_cc = dive( ipnt )  ! dive at center of cell.
    d_ri = dive( c__1 )  !         cell on right.
    d_to = dive( c__3 )  !         cell on top.
    d_le = dive( c__5 )  !         cell on left.
    d_bl = dive( c__6 )  !         cell at bottom-left.
    d_bo = dive( c__7 )  !         cell at bottom.
 
!   Leith's viscosity `nu' (1996, Physica D).
!     `v_ll' is defined at Lower-Left corner.
!     `v_cc' is defined at Center of Cell.
!     Units are m**2 s**(-1).

    v_ll( ipnt, ilay ) = ( r_br - r_bl ) * ( r_br - r_bl )           &
                       + ( r_bl - rbll ) * ( r_bl - rbll )           &
                       + ( r_tl - r_bl ) * ( r_tl - r_bl )           &
                       + ( r_bl - rbbl ) * ( r_bl - rbbl )           &
!   `Modified Leith viscosity'.
!   See Fox-Kemper and Menemenlis, 2008, Ocean Modeling in an Eddying Regime,
!   their Eq.39.
                       + ( d_cc - d_le ) * ( d_cc - d_le )           &
                       + ( d_bo - d_bl ) * ( d_bo - d_bl )           &
                       + ( d_cc - d_bo ) * ( d_cc - d_bo )           &
                       + ( d_le - d_bl ) * ( d_le - d_bl )
    v_ll( ipnt, ilay ) = sqrt( v_ll( ipnt, ilay ) ) * dvis * dl * dl &
                       + bvis 

!   Leith's viscosity `nu' (1996, Physica D).
    v_cc( ipnt, ilay ) = ( r_br - r_bl ) * ( r_br - r_bl )           &
                       + ( r_tr - r_tl ) * ( r_tr - r_tl )           &
                       + ( r_tl - r_bl ) * ( r_tl - r_bl )           &
                       + ( r_tr - r_br ) * ( r_tr - r_br )           &
!   `Modified Leith viscosity'.
                       + ( d_ri - d_cc ) * ( d_ri - d_cc )           &
                       + ( d_cc - d_le ) * ( d_cc - d_le )           &
                       + ( d_to - d_cc ) * ( d_to - d_cc )           &
                       + ( d_cc - d_bo ) * ( d_cc - d_bo )
    v_cc( ipnt, ilay ) = sqrt( v_cc( ipnt, ilay ) ) * dvis * dl * dl &
                       + bvis
! Biharmonic viscosity.               
! Calculated using thickness weighted grad^4(u,v).
! See Schepetkin and O'Brien '96.
! First calculate Laplacians of velocity (adding component by component).
! Laplacian component is 0 in direction of boundary.
    if ( svis > 0._rw ) then
      delu(ipnt,ilay) = 0._rw
      delv(ipnt,ilay) = 0._rw
      
   !   if (mk_u(ipnt) > 0.5) then !not west or east bdy
   !     delu(ipnt,ilay) = delu(ipnt,ilay) + 1._rw/dl**2 * (u(c__1,ilay) - 2*u(ipnt,ilay) + u(c__5,ilay))
   !     delv(ipnt,ilay) = delv(ipnt,ilay) + 1._rw/dl**2 * (v(c__1,ilay) - 2*v(ipnt,ilay) + v(c__5,ilay))
   !     if (mk_v(ipnt) < 0.5 .and. mk_n(ipnt) > 0.5) then ! if south
   !       delu(ipnt,ilay) = delu(ipnt,ilay) + 1._rw/dl**2 * (u(c__3,ilay) - u(ipnt,ilay))
   !     elseif (mk_v(ipnt) < 0.5 .and. mk_n(ipnt) < 0.5) then ! if north
   !       delu(ipnt,ilay) = delu(ipnt,ilay) + 1._rw/dl**2 * (-u(ipnt,ilay) + u(c__7,ilay))
   !     end if
   !   end if

   !   if (mk_v(ipnt) > 0.5) then !not north or south bdy
   !     delu(ipnt,ilay) = delu(ipnt,ilay) + 1._rw/dl**2 * (u(c__3,ilay) - 2*u(ipnt,ilay) + u(c__7,ilay))
   !     delv(ipnt,ilay) = delv(ipnt,ilay) + 1._rw/dl**2 * (v(c__3,ilay) - 2*v(ipnt,ilay) + v(c__7,ilay))
   !     if (mk_u(ipnt) < 0.5 .and. mk_n(ipnt) > 0.5) then ! if west
   !       delv(ipnt,ilay) = delv(ipnt,ilay) + 1._rw/dl**2 * (v(c__1,ilay) - v(ipnt,ilay))
   !     elseif (mk_u(ipnt) < 0.5 .and. mk_n(ipnt) < 0.5) then ! if east
   !       delv(ipnt,ilay) = delv(ipnt,ilay) + 1._rw/dl**2 * (-v(ipnt,ilay) + v(c__5,ilay))
   !     end if
   !   end if
       if (mk_u(ipnt) > 0.5) then
         delu(ipnt,ilay) = delu(ipnt,ilay) + 1._rw/dl**2 * &
           (mk_u(c__1)*u(c__1,ilay)+ mk_u(c__3)*u(c__3,ilay)+ mk_u(c__5)*u(c__5,ilay)+ &
           mk_u(c__7)*u(c__7,ilay) )
       
        delu(ipnt,ilay) = delu(ipnt,ilay) - 1._rw/dl**2 * &
          (mk_u(c__1)+mk_u(c__3)+mk_u(c__5)+mk_u(c__7))*u(ipnt,ilay)  
       
       end if

       if (mk_v(ipnt) > 0.5) then 
         delv(ipnt,ilay) = delv(ipnt,ilay) + 1._rw/dl**2 * &
           (mk_v(c__1)*v(c__1,ilay)+ mk_v(c__3)*v(c__3,ilay)+ mk_v(c__5)*v(c__5,ilay)+ &
           mk_v(c__7)*v(c__7,ilay) )
         delv(ipnt,ilay) = delv(ipnt,ilay) - 1._rw/dl**2 * &
           (mk_v(c__1)+mk_v(c__3)+mk_v(c__5)+mk_v(c__7))*v(ipnt,ilay)

       end if
   
     end if
  end do

!$OMP END PARALLEL DO

! Calculate thickness weighted grad^4(u,v) using delu and delv.
  if (svis > 0._rw ) then
    do ipnt = 1, ndeg
        c__1 = neig( 1, ipnt )
        c__3 = neig( 3, ipnt )
        c__5 = neig( 5, ipnt )
        c__6 = neig( 6, ipnt )
        c__7 = neig( 7, ipnt )
        UU4(ipnt,ilay) = 0._rw
        VV4(ipnt,ilay) = 0._rw
        hh_q= real(hlay(ipnt,ilay) +mk_n(c__5)*hlay(c__5,ilay)             &
              +mk_n(c__6)*hlay(c__6,ilay) + mk_n(c__7)*hlay(c__7,ilay))    &
              / (1._rw+mk_n(c__5)+mk_n(c__6)+mk_n(c__7)) 

        UU4(ipnt,ilay)=UU4(ipnt,ilay)-1._rw/dl*hlay(ipnt,ilay)             &
                  *delu(ipnt,ilay)+1._rw/dl*hlay(ipnt,ilay)*delv(ipnt,ilay)
        VV4(ipnt,ilay)=VV4(ipnt,ilay)+1._rw/dl*hh_q*delu(ipnt,ilay)        &
                  +1._rw/dl*hh_q*delv(ipnt,ilay)

        if (subc(ipnt,1)<=lm-1) then!not east
        !if (mk_u(c__1)> 0.5) then
          UU4(ipnt,ilay) = UU4(ipnt,ilay) +  1._rw/dl                      &
                           *hlay(ipnt,ilay)*delu(c__1,ilay)
        end if
        if (subc(ipnt,2)<=mm-1) then !not north
        !if (mk_v(c__3)> 0.5) then
          UU4(ipnt,ilay) = UU4(ipnt,ilay) - 1._rw/dl                       &
                           *hlay(ipnt,ilay)*delv(c__3,ilay)
        end if
        if (subc(ipnt,1)>1) then !not west
        !if (mk_v(c__5)> 0.5) then
          VV4(ipnt,ilay) = VV4(ipnt,ilay) - 1._rw/dl*hh_q*delv(c__5,ilay)
        end if
        if (subc(ipnt,2)>1) then !not south
        !if (mk_u(c__7)> 0.5) then
          VV4(ipnt,ilay) = VV4(ipnt,ilay) -  1._rw/dl*hh_q*delu(c__7,ilay)
        end if
    
       ! if (subc(ipnt,1)>=lm+1 .or. subc(ipnt,1)<=2 .or.                   &
       !                subc(ipnt,2)==mm+2 .or. subc(ipnt,2)==2) then
       if (mk_u(ipnt)*mk_v(ipnt)<0.5) then !this is effectively mk_vort
          VV4(ipnt,ilay) = 0._rw
        end if 
    end do
  end if

!print *, 'u',u(1:ndeg,2)
!print *, 'v',v(1:ndeg,2)

!print *, 'delu',  delu(1:ndeg,2)
!print *, 'delv', delv(1:ndeg,2)
!print *, 'UU4', UU4(1:ndeg,2)
!print *, 'VV4', VV4(1:ndeg,2)
!print *, 'hlay', hlay(1:ndeg,2)
!read(*,*) 

end subroutine update_viscosity

subroutine no_gradient_obc( ilay )
  implicit none
  integer, intent(in) :: ilay
  integer             :: iseg, ipnt

! The component of velocity that is parallel to nudged open boundaries is
!   typically not nudged. Its value at the boundary is obtained by assuming
!   vanishing normal derivative (Lavelle & Thacker 2008, Ocean Modelling).
! Note that `at the boundary' means the wet cell closest to the boundary;
!   see Fig.2 of Herzfeld (2009, Ocean Modelling).

  do iseg = 1, size( segm, dim = 1 )
    ipnt = segm( iseg, 10 )
    if     ( segm( iseg, 5 ) == 1 ) then ! Northern or Southern open boundary.
      if ( mk_u( ipnt ) > 0.5 ) then
!       ipnt,i,j=16,17,18: the interior cell (i.e. the cell just inside of the
!                                             wet cell).
!       ipnt,i,j=10,11,12: the wet      cell (i.e. the cell on wet side of
!                                             boundary).
        u(   ipnt, ilay ) = u(    segm(iseg, 16), ilay       )      &
                          - fnud( segm(iseg, 16), ilay, ix_u )      &
                          + fnud( ipnt,           ilay, ix_u )
        h_u( ipnt, ilay ) = u( ipnt, ilay )                         &
                          * real( hlay( ipnt,          ilay )       &
                                + hlay( neig(5, ipnt), ilay ), rw ) &
                          / (1._rw + mk_u( ipnt ))
      end if
    elseif ( segm( iseg, 4 ) == 1 ) then ! Eastern  or Western  open boundary.
      if ( mk_v( ipnt ) > 0.5 ) then
        v(   ipnt, ilay ) = v(    segm(iseg, 16), ilay       )      &
                          - fnud( segm(iseg, 16), ilay, ix_v )      &
                          + fnud( ipnt,           ilay, ix_v )
        h_v( ipnt, ilay ) = v( ipnt, ilay )                         &
                          * real( hlay( ipnt,          ilay )       &
                                + hlay( neig(7, ipnt), ilay ), rw ) &
                          / (1._rw + mk_v( ipnt ))
      end if
    end if
  end do

! Apply a similar `no-gradient' condition on velocity component perpendicular
!   to nudged open boundaries (Lavelle & Thacker 2008 Ocean Modelling).
! `eta' at the boundary is directly calculated from these velocities.

  do iseg = 1, size( segm, dim = 1 )
    ipnt = segm(iseg, 1)
    if     ( segm(iseg, 5) == 1 ) then ! Northern or Southern open boundary.
!     ipnt,i,j=13,14,15: the interior, normal velocity point.
!     ipnt,i,j= 1, 2, 3: the normal velocity point at the boundary.
      v  ( ipnt, ilay ) = v   ( segm(iseg, 13), ilay       )      &
                        - fnud( segm(iseg, 13), ilay, ix_v )      &
                        + fnud( ipnt,           ilay, ix_v )
      h_v( ipnt, ilay ) = v( ipnt, ilay )                         &
                        * real( hlay( ipnt,          ilay )       &
                              + hlay( neig(7, ipnt), ilay ), rw ) &
                        / (1._rw + mk_v( ipnt ))
    elseif ( segm(iseg, 4) == 1 ) then ! Eastern  or Western  open boundary.
      u  ( ipnt, ilay ) = u   ( segm(iseg, 13), ilay       )      &
                        - fnud( segm(iseg, 13), ilay, ix_u )      &
                        + fnud( ipnt,           ilay, ix_u )
      h_u( ipnt, ilay ) = u( ipnt, ilay )                         &
                        * real( hlay( ipnt,          ilay )       &
                              + hlay( neig(5, ipnt), ilay ), rw ) &
                        / (1._rw + mk_u( ipnt ))
    end if
  end do
end subroutine no_gradient_obc

subroutine write_outputs( restart_record )
  implicit none
  integer, intent(in), optional :: restart_record
  logical, save                 :: is_init = .false.
  integer, save                 :: irec
  integer                       :: lerm, ilay
  integer                       :: unum
  real( r8 )                    :: date
  character( 9 )                :: str1

  lerm = len_trim( errm )
  errm = trim( errm ) // &
         ' in subroutine write_outputs from module private_mod.f95,'

! How to read output files:
!   1- Time-axis is in its own file (`time.txt'),
!        readable with `load' under Octave/Matlab/Octave.
!      Each entry in time.txt corresponds to an output record.
!   2- The size of the grid (lm x mm) is obtained from file param_basin.txt.
!   3- Each variable is saved as its own file, e.g., u___.bin for u velocity.
!   4- Only wet grid cells are saved in output files.
!      The (i,j) location of the wet grid cells is in file grid.bin.

  if ( present(restart_record) ) then
    is_init = .true.
    irec    = restart_record + 1 ! Ready for next writing.
    errm    = errm(1 : lerm)
    return  ! Do not write anything (initial record is last record of
            !   previous run).
  end if

  if ( is_init ) then

    call write_array( var = 'eta_', irec = irec, is_init = is_init )
    call write_array( var = 'u___', irec = irec, is_init = is_init )
    call write_array( var = 'v___', irec = irec, is_init = is_init )

    if ( diag > 0.5_rw ) then
      call write_array( var = 'pvor', irec = irec, is_init = is_init )
      call write_array( var = 'mont', irec = irec, is_init = is_init )
      call write_array( var = 'v_cc', irec = irec, is_init = is_init )
    end if

!   The time.txt file is only updated if (i.e. after) ALL arrays have been
!     successfully written on disk. In case of a restart, the model looks
!     for the record that corresponds to the last entry in file `time.txt'.
!   For instance, if for some reason the code stops while writing outputs
!     (e.g. the disk becomes full), the run can still be restarted since
!     the code will seek the last complete record and ignore the partially
!     written records.

    date = real(ctim, r8) !+ real(itim, r8)
    unum = get_un()
    open( unit   = unum,  form   = 'formatted', action   = 'write',  &
          status = 'old', iostat = errc,        position = 'append', &
          file   = trim( odir ) // 'time.txt' )
      write( unum, fmt = * ) date
    close( unit = unum )

    irec = irec + 1 ! Get ready for next writing.

  else ! If the outputs are not initialized yet...

    irec = 1

    call write_array( var = 'eta_', irec = irec, is_init = is_init )
    call write_array( var = 'u___', irec = irec, is_init = is_init )
    call write_array( var = 'v___', irec = irec, is_init = is_init )

    if ( diag > 0.5_rw ) then
      call write_array( var = 'pvor', irec = irec, is_init = is_init )
      call write_array( var = 'mont', irec = irec, is_init = is_init )
      call write_array( var = 'v_cc', irec = irec, is_init = is_init )
    end if

    date = real(ctim, r8) !+ real(itim, r8)
    unum = get_un()
    open( unit   = unum,      form = 'formatted', action = 'write', &
          status = 'replace', iostat = errc,                        &
          file   = trim( odir ) // 'time.txt' )
      write( unum, fmt = * ) date
    close( unit = unum )

    irec    = irec + 1 ! Get ready for next writing.

    is_init = .true.
  end if

  write( ioso, * ) 'ctim = ', ctim, ' days; dt_s = ', dt_s, &
                   ' days; record = ', irec - 1

  do ilay = 1, nlay
    write( ioso, * ) 'min/max h', ilay, '= ',                     &
      minval( real( hlay(:, ilay), rw ), mask = mk_n(:) > 0.5 ),  &
      maxval( real( hlay(:, ilay), rw ), mask = mk_n(:) > 0.5 )
    if ( any( mk_u(:) > 0.5 ) ) then
      write(ioso,*) 'min/max u', ilay, '= ',                      &
                     minval(u   (:, ilay), mask = mk_u(:) > 0.5), &
                     maxval(u   (:, ilay), mask = mk_u(:) > 0.5)
    else
      write(ioso,*) 'min/max u', ilay, '= ',                      &
                     minval(u   (:, ilay)),                       &
                     maxval(u   (:, ilay))
    end if
    if ( any( mk_v(:) > 0.5 ) ) then
      write(ioso,*) 'min/max v', ilay, '= ',                      &
                     minval(v   (:, ilay), mask = mk_v(:) > 0.5), &
                     maxval(v   (:, ilay), mask = mk_v(:) > 0.5)
    else
      write(ioso,*) 'min/max v', ilay, '= ',                      &
                     minval(v   (:, ilay)),                       &
                     maxval(v   (:, ilay))
    end if
  end do

! Check for layer thickness < hmin.

  do ilay = 1, nlay
    if ( any( mk_n( 1 :       ) >       0.5_rw                &
        .and. hlay( 1 :, ilay ) < real( 0.5_rw * hmin, r8 ) ) ) then
      write( str1, '(1i8)' ) ilay
      errc = int(  ilay ) * (- 1)
      errm = trim( errm ) // ' layer number ' // trim( str1 ) &
                          // ' has its thickness < hmin;'     &
                          // ' Calculation halted.'
      call quit()
    end if
  end do

  if ( errc == 0 ) then
    errm = errm( 1 : lerm )
  else
    call quit()
  end if
end subroutine write_outputs

subroutine write_array( var, irec, is_init )
  implicit none
  character(4), intent(in)  :: var
  integer,      intent(in)  :: irec
  logical,      intent(in)  :: is_init
  real   (r4),  allocatable :: ior4(:,:), h_0( :,:)
  real   (rw),  allocatable :: wrk1(:,:), wrk2(:,:)
  real   (rw) :: zeta
  integer     :: lrec, unum, lerm, ilay, ipnt, i, &
                 c__1, c__2, c__3, c__5, c__6, c__7

  lerm = len_trim( errm )
  errm = trim(     errm ) &
       // ' in subroutine write_array from module private_mod.f95,'

  allocate( ior4( ndeg, nlay ) )

  if     ( var == 'eta_' ) then ! Meters.

!   You must invert the layer thickness `hlay' to obtain the isopycnal
!     elevation `eta'. This operation requires the array `h_0'. 

    unum = get_un()
    allocate( h_0( ndeg, nlay ) )
    inquire( iolength = lrec ) h_0(:,:)
    open( unit = unum, form   = 'unformatted', access = 'direct', &
          recl = lrec, status = 'old',         action = 'read',   &
          file = trim( odir ) // 'h_0.bin' )
    read( unum, rec = 1 ) h_0(:,:)
    close( unit = unum )

    do ilay = nlay, 2, - 1
      do ipnt = 1, ndeg
        if ( ilay == nlay ) then
          ior4( ipnt, ilay )                           &
            = real(       hlay( ipnt, ilay     )       &
                  - real( h_0(  ipnt, ilay     ), r8 ), r4 )
        else
          ior4( ipnt, ilay )                           &
            = real(       hlay( ipnt, ilay     )       &
                  - real( h_0(  ipnt, ilay     ), r8 ) &
                  + real( ior4( ipnt, ilay + 1 ), r8 ), r4 )
        end if
      end do
    end do
    
    ilay = 1
    if (rgld < 0.5_rw) then
       do ipnt=1, ndeg
          ior4( ipnt, ilay )                           &
            = real(       hlay( ipnt, ilay     )       &
                  - real( h_0(  ipnt, ilay     ), r8 ) &
                  + real( ior4( ipnt, ilay + 1 ), r8 ), r4 )
       end do
    else
       do ipnt=1, ndeg
       ior4(ipnt, ilay) = real(pi_s(ipnt) ) !+ real (ior4(ipnt, ilay+1), r8), r4)
       end do
    end if 
    
    
    deallocate( h_0 )

  elseif ( var == 'u___' ) then ! m s**(-1).
    ior4(:,:) = real( u(1 :, :), r4 )
  elseif ( var == 'v___' ) then ! m s**(-1).
    ior4(:,:) = real( v(1 :, :), r4 )
  elseif ( var == 'v_cc' ) then ! m**2 s**(-1).
    allocate( wrk1(0 : ndeg, nlay) )
    allocate( wrk2(0 : ndeg, nlay) )
    ior4(:,:) = 0._r4
    wrk1(:,:) = 0._rw
    wrk2(:,:) = 0._rw

    do ilay = 1, nlay
!$OMP PARALLEL DO PRIVATE(ipnt,c__1,c__3,c__5,c__7)
      do ipnt = 1, ndeg
        c__1 = neig(1, ipnt)
        c__3 = neig(3, ipnt)
        c__5 = neig(5, ipnt)
        c__7 = neig(7, ipnt)
        wrk1( ipnt, ilay ) = ( ( v(ipnt, ilay) - v(c__5, ilay) ) / dl   &
                             - ( u(ipnt, ilay) - u(c__7, ilay) ) / dl ) &
                           * mkpe( ipnt ) ! s**(-1).
        wrk2( ipnt, ilay ) = ( u(c__1, ilay) - u(ipnt, ilay) ) / dl     &
                           + ( v(c__3, ilay) - v(ipnt, ilay) ) / dl
      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(ipnt,c__1,c__2,c__3,c__5,c__7)
      do ipnt = 1, ndeg
        c__1 = neig(1, ipnt)
        c__2 = neig(2, ipnt)
        c__3 = neig(3, ipnt)
        c__5 = neig(5, ipnt)
        c__7 = neig(7, ipnt)

        ior4( ipnt, ilay )                                               &
          = real( bvis                                                   &
                + dvis * dl**2                                           &
                * sqrt( ( wrk1( c__1, ilay ) - wrk1( ipnt, ilay ) )**2   &           
                      + ( wrk1( c__2, ilay ) - wrk1( c__3, ilay ) )**2   &
                      + ( wrk1( c__3, ilay ) - wrk1( ipnt, ilay ) )**2   &
                      + ( wrk1( c__2, ilay ) - wrk1( c__1, ilay ) )**2   &
                      + ( wrk2( c__1, ilay ) - wrk2( ipnt, ilay ) )**2   &
                      + ( wrk2( ipnt, ilay ) - wrk2( c__5, ilay ) )**2   &
                      + ( wrk2( c__3, ilay ) - wrk2( ipnt, ilay ) )**2   &
                      + ( wrk2( ipnt, ilay ) - wrk2( c__7, ilay ) )**2 ), r4 )
      end do
!$OMP END PARALLEL DO
    end do
    deallocate( wrk1 )
    deallocate( wrk2 )
  elseif ( var == 'mont' ) then ! Meters.
    ior4(:,:) = 0._r4
    do ilay = 1, nlay
!$OMP PARALLEL DO PRIVATE(ipnt,i)
      do ipnt = 1, ndeg
        ior4( ipnt, ilay )                                               &
          = real( - ocrp / real( nsal - 1, rw ) * hsal * mk_n( ipnt )    &
                  * ( hsal / ( hmin * ( 1._rw - mk_n( ipnt ) )           &
                    + real( hlay( ipnt, ilay ), rw ) ) )**(nsal - 1), r4 )
        do i = 1, ilay - 1
          ior4( ipnt, ilay ) = ior4( ipnt, ilay )                        &
                             - real( ( rhon( ilay ) - rhon( i ) )        &
                                     * real( hlay( ipnt, i ), rw )       &
                                     / rhon( ilay ), r4 )
        end do
        ior4(ipnt, ilay) = ior4(ipnt, ilay)                              &
                         + real( sum(       hlay( ipnt, : )       )      &
                                    - real( h_th(ipnt), r8 ), r4 )
      end do
!$OMP END PARALLEL DO
    end do
  elseif ( var == 'pvor' ) then ! m**(-1) s**(-1).
    ior4(:,:) = 0._r4
    do ilay = 1, nlay
!$OMP PARALLEL DO PRIVATE(ipnt,c__5,c__6,c__7,zeta)
      do ipnt = 1, ndeg
        c__5 = neig(5, ipnt)
        c__6 = neig(6, ipnt)
        c__7 = neig(7, ipnt)
        zeta = ( ( v(ipnt, ilay) - v(c__5, ilay) ) / dl       &
               - ( u(ipnt, ilay) - u(c__7, ilay) ) / dl )     &
             * mkpe(ipnt)
        ior4(ipnt, ilay) = real( ( fcor(ipnt) + zeta * uadv ) &
                                 * mkpi(ipnt)                 &
                                 * ( mk_n(ipnt)               &
                                   + mk_n(c__5)               &
                                   + mk_n(c__7)               &
                                   + mk_n(c__6) )             &
                                 / real( hlay(ipnt, ilay)     &
                                       + hlay(c__5, ilay)     &
                                       + hlay(c__6, ilay)     &
                                       + hlay(c__7, ilay), rw ), r4 )
      end do
!$OMP END PARALLEL DO
    end do
  end if

  inquire( iolength = lrec ) ior4(:,:)
  unum = get_un()

  if ( is_init ) then
    open( unit = unum, status = 'old',     access = 'direct', iostat = errc, &
          recl = lrec, action = 'write',   form   = 'unformatted',           &
          file = trim( odir ) // var // '.bin' )
  else
    open( unit = unum, status = 'replace', access = 'direct', iostat = errc, &
          recl = lrec, action = 'write',   form   = 'unformatted',           &
          file = trim( odir ) // var // '.bin' )
  end if

  if ( errc == 0 ) then
    write( unit = unum, rec = irec, iostat = errc ) ior4(:,:)
    deallocate( ior4 )
  end if
  close( unit = unum )

  if ( errc == 0 ) then
    errm = errm(1 : lerm)
  else
    call quit()
  end if
end subroutine write_array

end module private_mod
