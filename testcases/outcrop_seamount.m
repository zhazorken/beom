
% Test-case for isopycnal outcropping over a seamount and sloppy coast.
%   Verifies that the initial state corresponds to the state of rest (du/dt=0).
%   See Salmon, R. (2002), J. Mar. Res., vol.60, p.605-638.

  clear;
  close all;
  more  off;

  outd   = '/tmp/';        % Path to input/output files.
  lx     = 600.e3;         % Length of domain in zonal direction (meters).
  dl     =  5.e3;          % Mesh size (meters).
  nlay   = 5;              % Number of layers.
  fcor   = + 1.e-4;        % Coriolis parameter (rad/s).
  if ( nlay == 1 )
    rhon = 1000.;          % Density of layers (kg/m3).
    topl = 0.;             % Top     of layers (fraction of Hmax).
  else
    rhon = 1000. + 30. * (0. : 1. / (nlay - 1) : 1.);
    topl = (0 : nlay - 1) ./ nlay;
  end
  grav   = 9.8;            % m/s2
  hmax   = 300.;           % Max basin depth before adding seamount (meters).
  lm     = round(lx / dl); % Grid size in zonal (x) direction.
  if ( mod(lm, 2) == 0 )
    lm = lm + 1;           % Make sure it is odd.
  end
  mm     = 1;              % For 2-D calculation (x-z plane).
% mm     = lm;             % For 3-D calculation.

  [yy, xx] = meshgrid( 1 : mm + 2, 1 : lm + 2 );
  xx = xx - mean(xx(:));
  yy = yy - mean(yy(:));
  xx = xx * dl;
  yy = yy * dl;

% Create a closed basin with depth increasing from zero at the edges to hmax
%   at the center. Then, add a Gaussian seamount at the center.
%   The steepest bottom slopes are found on the sides of the seamount.

  ix_0 = ceil(0.5 * (mm + 2)); % Center of grid.

  h_bo = sqrt(xx.^2 + yy.^2) / dl;
  h_bo = h_bo ./ max(h_bo(:, ix_0)) * hmax;
  h_bo = max(h_bo(:, ix_0)) - h_bo;
  smnt = hmax - exp(- (xx.^2 + yy.^2) ./ (50.e3)^2) * 0.75 * hmax;
  h_bo = min([h_bo(:)'; smnt(:)'])'; clear smnt;
  h_bo = reshape(h_bo, lm + 2, mm + 2);
  h_bo( find( h_bo < min(h_bo(:, ix_0)) ) ) = 0.;

% Bottom slope (non-dimensional).

  dhdx = zeros( lm + 2, mm + 2 );
  dhdy = zeros( lm + 2, mm + 2 );
  dhdx( 2 : end - 1, : ) = (h_bo( 3 : end, : ) - h_bo( 1 : end - 2, : )) / (2. * dl);
  dhdy( :, 2 : end - 1 ) = (h_bo( :, 3 : end ) - h_bo( :, 1 : end - 2 )) / (2. * dl);
  slop = sqrt( dhdx.^2 + dhdy.^2 );
  clear dhdx dhdy;

% hsal is the `Salmon's thickness' (in meters).
% `Baroclinic' version of Salmon's criterion (Eq.6.5, Salmon 2002 JMR):
%    hsal > (slope * dl) * (delta rho / rho_0)

  hsal = 1.0 * max( slop( find( h_bo > 1.e-3 ) ) ) * dl ...
       * (rhon( 2 ) - rhon( 1 )) / rhon( 2 );

% `hmin' is the minimum layer thickness allowed during a calculation.
%   By default the code assumes `hsal = 10 * hmin'.

  hmin = hsal / 10.;                     % Meters.

% There should be no `wetting-drying' in this test-case.
% Therefore, you must ensure that hlay(1) is always >> hsal by drying-up all
%   cells of insufficient depth. Note that the top layer has to accommodate all
%   the deeper layers that are outcropped. Layers that are outcropped have
%   a thickness approximately equal to `hsal'.

  dryd = 10. * hsal + (nlay - 1) * hsal; % Meters.

  h_bo( find( h_bo < dryd ) ) = 0.;      % Everything under `dryd' is dryied-up.

  h_bo(  1, :) = 0.;                     % Enforce dry margins.
  h_bo(end, :) = 0.;
  h_bo(:, end) = 0.;
  h_bo(:,   1) = 0.;

  cext = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).

  ndeg = get_nbr_deg_freedom( h_bo );

% Write input files over disk.

  [fid, msg] = fopen([outd 'h_bo.bin'], 'w', 'ieee-le');
  cnt = fwrite(fid, h_bo, 'real*4', 0,       'ieee-le');
  fclose(fid); clear fid msg cnt;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, fcor, rhon, topl, 15.,  ...
                0.2,  0.,   0.,   0.,   0.,   0.,   hmin, 10.,  10.,  1.,   ...
                1.,   0.,   1.,   0.,   0.,   0.,   0.,   [0,0],outd, outd, ...
                'Test-case for state of rest allowing isopycnal outcrop' );

% See beom documentation on how to compile and execute the program.

  disp(['Now compile and execute the program in a different terminal. ' ...
        'Press a key when it is over.']);
  pause

  disp('Reading the output data, please wait...');

% Read metadata.

  [taxi, h_0, f0, dl, rhon, desc] = get_metadata(strtrim(outd));
  nrec = length(taxi);         % Nbr of records in output files.
  lm   = size(h_0, 1) - 2;     % Size of grid.
  mm   = size(h_0, 2) - 2;
  nlay = size(h_0, 3);

  kine = nan(nrec, nlay);

  for irec = 1 : 1 : nrec
    eta = get_field('eta_', irec, outd);
    u   = get_field('u___', irec, outd);
    v   = get_field('v___', irec, outd);

    if ( any( h_bo(2 : lm + 1, 2 : mm + 1) > 0.5 ...
            & isnan( eta(2 : lm + 1, 2 : mm + 1, 1) ) ) )
      if ( irec > 1 )
        eta = get_field('eta_', irec - 1, outd);
        u   = get_field('u___', irec - 1, outd);
        v   = get_field('v___', irec - 1, outd);
      end
      stop
    end

    zlay = nan(lm + 2, mm + 2, nlay); % Distance from geoid z=0 (m).
    hlay = nan(lm + 2, mm + 2, nlay); % Thickness (m).

    for ilay = 1 : nlay
      zlay(:,:,ilay) = h_bo - squeeze(sum(h_0(:,:,ilay : nlay), 3)) - eta(:,:,ilay);
    end; clear ilay;

    hlay(:,:,nlay) = h_0(:,:,nlay) + eta(:,:,nlay);
    for ilay = nlay - 1 : -1 : 1
      hlay(:,:,ilay) = h_0(:,:,ilay) + eta(:,:,ilay) - eta(:,:,ilay + 1);
    end; clear ilay;

    if ( mm == 1 )
      kine(irec, :) = nansum( ...
                        reshape( ...
                          0.5 * hlay .* (u.^2       ), [(lm + 2) * (mm + 2), nlay] ), 1 );
    else
      kine(irec, :) = nansum( ...
                        reshape( ...
                          0.5 * hlay .* (u.^2 + v.^2), [(lm + 2) * (mm + 2), nlay] ), 1 );
    end

    if ( irec == 1 || irec == nrec )
      figure;
      subplot(1,2,1);
      hand = plot( xx(:,ix_0) / 1.e3, [squeeze(zlay(:,ix_0,:)), h_bo(:,ix_0)] );
      legend( num2str( (1 : nlay)' ) );
      set( hand(end), 'color', [0 0 0] );
      xlabel('Distance from seamount (km)');
      ylabel('Depth (m)')
      title(['Position of layer interface after ' num2str(taxi(irec)) ' days']);
      axis ij;
      axis([min(xx(:)) / 1.e3, max(xx(:)) / 1.e3, - 1., max(h_bo(:))]);
  
      subplot(1,2,2);
      plot( xx(:,ix_0) / 1.e3, squeeze(hlay(:,ix_0,:)) );
      legend( num2str( (1 : nlay)' ) );
      xlabel('Distance from seamount (km)');
      ylabel('Thickness (m)');
      title(['Layer thickness after ' num2str(taxi(irec)) ' days']);
      axis([min(xx(:)) / 1.e3, max(xx(:)) / 1.e3, 0, max(h_0(:)) + 1.]);
      sleep(0.1);
    end
  end; clear irec;

figure;
  hand = plot( taxi, kine );
  legend( num2str( (1 : nlay)' ) );
  xlabel('Time (days)');
  ylabel('Total layer kinetic energy (m^3 s^{-2})');
  title('Timeseries of kinetic energy');
  axis([0, ceil(taxi(end))]);
