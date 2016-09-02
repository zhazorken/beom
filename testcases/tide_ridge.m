
% Test-case of a semidiurnal tidal flow over a ridge (2-D).

  clear;
  close all;
  more  off;

  ocrp = 1;                                    % Allow isopycnal outcrops.
  movi = 0;                                    % Create frames for a movie.
  outd = '/tmp/';                              % Path to input/output files.
  hmax = 30.;                                  % Maximum depth of basin (m).
  dt_s = 3.;                                   % Duration of calc.     (days).
  dt_r = 0.75 * dt_s;                          % Ramp on tidal forcing (days).
  bdrg = 0.;                                   % Bottom drag coefficient.

  if ( ocrp == 1 )
    nlay = 7;                                  % Number of layers.
  else
    nlay = 3;                                  % Number of layers.
  end

  fcor = 0.;                                   % Coriolis parameter (rad/s).
  dl   = 100.;                                 % Mesh size (m).
  mm   = 1;                                    % Grid size in across-stream direction.
  lm   = 500;                                  % Grid size in along -stream direction.
  if ( mod(lm, 2) == 0 ); lm = lm + 1; end;    % Make it odd.
  npts = 15;                                   % Number of grid points in the sponge.
  grav = 9.8;                                  % m/s2

% `rhop' is top-to-bottom neutral density profile (kg/m3).

  rhop = 1028. - 6. * exp( - (0. : 1. : hmax)' ./ (0.2 * hmax) );
  drho = ( max( rhop ) - min( rhop ) ) / nlay;
  rhon = min( rhop ) + (0.5 : nlay - 0.5) * drho;

  [xx, yy] = meshgrid( (1 : lm + 2)' - 1.5, (1 : mm + 2) - 1.5 );
  xx   = xx' * dl;
  yy   = yy' * dl;
  xx   = xx - mean(xx(:)); % Centered on ridge.
  yy   = yy - mean(yy(:));

  h_bo = 0.75 * hmax * exp( - (xx ./ (75. * dl)).^2 );
  h_bo = hmax - h_bo;

  h_bo(:,   1) = 0.;
  h_bo(:, end) = 0.;
  h_bo(  1, :) = 0.;
  h_bo(end, :) = 0.;

  ndeg = get_nbr_deg_freedom( h_bo );
  cext = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).

  dhdx = zeros( lm + 2, mm + 2 );
  dhdx( 3 : end - 2, : ) = h_bo( 4 : end - 1, : ) ...
                         - h_bo( 2 : end - 3, : );
  dhdx = dhdx / (2. * dl);

% Salmon's thickness hsal (m).

  hsal = 20. * max( dhdx(:) ) * dl * (rhon(end) - rhon(1)) / rhon(1);
  hmin = hsal / 10.;            % Min. layer thickness allowed (m).

  topl = interp1( rhop, (0. : 1. : hmax), min( rhop ) + (0 : nlay - 1) * drho );

  if ( topl(2) < 10. * hsal )
    topl( 2 : end ) = topl( 2 : end ) + (10. * hsal - topl(2));
  end

  topl = topl / hmax;

% Array `tide' has dimensions [2 x ncon x (lm + 2) x (mm + 2) x 3]
%   and stores the tidal harmonics of the entire grid:
%   1st dim. is for amplitude (m or m/s) and phase (rad),
%   2nd dim. is for the tidal constituents (here 1 constituent),
%   3rd dim. is for zonal  grid points,
%   4th dim. is for merid. grid points,
%   5th dim. is for variables `eta' (interface displacement), `u', and 'v'.

  tide = nan( 2, 1, lm + 2, mm + 2, 3 );

% Assume an oscillating current starting from rest
%   (i.e. the phase of `u' is pi/2 which leads to cosine()=0 at t=0).

  tide( 1, 1, :, :, 1 ) = 0.;      % Ampl. eta (meter).
  tide( 1, 1, :, :, 2 ) = 0.1;     % Ampl. u   (meter per second).
  tide( 1, 1, :, :, 3 ) = 0.;      % Ampl. v   (meter per second).
  tide( 2, 1, :, :, 1 ) = 0.;      % Phase eta (rad).
  tide( 2, 1, :, :, 2 ) = pi / 2.; % Phase u   (rad).
  tide( 2, 1, :, :, 3 ) = pi / 2.; % Phase v   (rad).

% Insert the constituent frequency in the first element of the array.

  tide( 1, 1, 1, 1, 1 ) = 2. * pi / (12.4206012 / 24.); % M2 (rad/days).

% Add wave sponges.

  nudg = zeros( lm + 2, mm + 2, 3    ); % eta, u, v.

  east = zeros( lm + 2, mm + 2, 3 );
  west = zeros( lm + 2, mm + 2, 3 );
  dt   = 0.5  * dl / cext; % Model timestep (seconds).
  widt = npts * dl;        % Width of sponge (meters).

  for i = 1 : lm + 2 % Eastern boundary.
%  `xpos' is the position within the sponge zone (in units of grid points).
    xpos             = i - 1.5 + npts - lm;
    xpos             = max([xpos;        0. ]);
    xpos             = min([xpos; npts - 0.5]);
%   Lavelle and Thacker 2008, Ocean Modelling, vol.20 p.270-292:
%     for eastern boundary, relax `eta' (1) and `u' (2), but not `v' (3).
%   The nudging coefficients (non-dimenzionalized by timestep dt)
%     vary as Eq.29 of Modave et al. 2010 Ocean Dynamics vol.60 p.65-79.
    east(i,:,[1, 2]) = dt * cext / widt * xpos / (npts - xpos);
  end; clear i xpos;

  for i = 1 : lm + 2 % Western boundary.
    xpos             = npts - (i - 1.5);
    xpos             = max([xpos;        0. ]);
    xpos             = min([xpos; npts - 0.5]);
    west(i,:,[1, 2]) = dt * cext / widt * xpos / (npts - xpos);
  end; clear i xpos;

% For each variable (eta,u,v),
%   take the maximum of eastern/western sides
%   (e.g., Fig.2, Lavelle & Thacker 2008).

  for ivar = 1 : size( nudg, 3 )
    temp = max( [reshape(nudg(:,:,ivar), [1, (lm + 2) * (mm + 2)]); ...
                 reshape(east(:,:,ivar), [1, (lm + 2) * (mm + 2)]); ...
                 reshape(west(:,:,ivar), [1, (lm + 2) * (mm + 2)])] );
    nudg(:,:,ivar) = reshape( temp', lm + 2, mm + 2 );
    clear temp;
  end; clear ivar;

% Write input files over disk.

  [fid, msg] = fopen([outd 'h_bo.bin'], 'w', 'ieee-le');
  cnt = fwrite(fid, h_bo,            'real*4', 0, 'ieee-le');
  fclose(fid); clear fid msg cnt;

  [fid, msg] = fopen([outd 'nudg.bin'], 'w', 'ieee-le');
  cnt = fwrite( fid, nudg, 'real*4',     0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

  [fid, msg] = fopen([outd 'tide.bin'], 'w', 'ieee-le');
  cnt = fwrite( fid, tide, 'real*4',     0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, fcor, rhon, topl, dt_s, ...
                0.01, dt_r, 0.,   0.,   0.5,  bdrg, hmin, 5.,   5.,   1.,   ...
                1.,   1.,   ocrp, 0.,   0.,   0.,   0.,   [0,0],outd, outd, ...
                'Test-case for tidal flow over a ridge' );

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
  h_bo( find( h_bo < 1.e-3 ) ) = nan; % For plotting purposes.

  ix_0 = ceil(0.5 * (mm + 2));

  for irec = 1 : nrec
    n    = get_field('eta_', irec, outd);

    zlay = nan(lm + 2, mm + 2, nlay);

    for ilay = 1 : nlay
      zlay(:,:,ilay) = h_bo - squeeze(sum(h_0(:,:,ilay : nlay), 3)) - n(:,:,ilay);
    end; clear ilay;

    hand = plot(        xx( 1 + npts : end - npts, ix_0 ) / 1.e3, ...
                      h_bo( 1 + npts : end - npts, ix_0 ), ...
                        xx( 1 + npts : end - npts, ix_0 ) / 1.e3, ...
             squeeze( zlay( 1 + npts : end - npts, ix_0, :) ) );
    if ( irec == 1 )
      for ilin = 1 : length( hand )
        colo( ilin, : ) = get( hand( ilin ), 'color' );
      end; clear ilin;
      colo = [ 0.5, 0.5, 0.5; colo(1 : end - 1, :) ];
    end
    for ilin = 1 : length( hand )
      set( hand( ilin ), 'color', colo( ilin, : ), 'linewidth', 1.5 );
    end; clear ilin;
    xlabel('Distance from ridge (km)');
    ylabel('Depth (m)');
    if ( irec == 1 )
      lege = [ 'Topography' ];
      for ilay = 1 : nlay
        lege = [ lege; 'Interface ' num2str( ilay ) ];
      end; clear ilay;
    end
    legend( lege, 'location', 'south');
    axis ij;
    axis([xx( 1 + npts, ix_0 ) / 1.e3, xx( end - npts, ix_0 ) / 1.e3, - 1., hmax]);
    axis square;
    sleep( 1.e-6 );
    if ( movi )
      myft = '/usr/share/fonts/type1/gsfonts/n019003l.pfb';
      print('-dpng', [outd 'frame' num2str(irec, '%04u') '.png'], '-S1024,768', ['-F' myft]);
    end
  end; clear irec;
