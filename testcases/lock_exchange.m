
% Test-case of 'lock-exchange'.
%   See Shchepetkin et al. 2015, Ocean Modelling, vol.91,    p.38-69,
%         doi:10.1016/j.ocemod.2015.03.006,
%       Ilicak et al.      2012, Ocean Modelling, vol.45-46, p.37-58,
%         doi:10.1016/j.ocemod.2011.10.003.

  clear;
  close all;
  more  off;

  outd = '/tmp/';                       % Path to input/output files.
  hmax =  20.;                          % Depth of basin (m).
  nlay = 2;                             % Number of layers.
  rhon = [1025., 1030.];                % Density of layers (kg/m3).
  topl = [0., 0.5];                     % Top     of layers (fraction of Hmax).
  grav = 9.8;                           % m/s2
  dl   = 400.;                          % Mesh size (m).
  lx   = 64.e3;                         % Length of domain (m).
  mm   = 1;                             % Grid size in across-stream direction.
  hmin = 0.005;                         % Minimum layer thickness allowed   (m).
  dt_o = 1. / 24.;                      % Time between model Outputs (days).

  hsal = 10. * hmin;                    % Salmon's thickness                (m).
  lm   = round( lx / dl );              % Grid size in along -stream direction.

  [xx, yy] = meshgrid( (1 : lm + 2)' - 1.5, (1 : mm + 2) - 1.5 );
  xx   = xx' * dl;
  yy   = yy' * dl;
  xx   = xx - mean(xx(:));
  yy   = yy - mean(yy(:));

  h_bo = hmax * ones( lm + 2, mm + 2 );

  h_bo(:,   1) = 0.;
  h_bo(:, end) = 0.;
  h_bo(  1, :) = 0.;
  h_bo(end, :) = 0.;

  ndeg = get_nbr_deg_freedom( h_bo );
  cext = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).
  cint = 0.5 * sqrt( grav * hmax * (rhon(2) - rhon(1)) / rhon(2) );

  n    = zeros( lm + 2, mm + 2, nlay );
  u    = zeros( lm + 2, mm + 2, nlay );
  v    = zeros( lm + 2, mm + 2, nlay );

  n( 1 : round( 0.5 * (lm + 2) ),           :, 2 ) =   0.5 * hmax - 4. * hsal;
  n(     round( 0.5 * (lm + 2) ) + 1 : end, :, 2 ) = - 0.5 * hmax + 4. * hsal;

% Write input files over disk.

  [fid, msg] = fopen([outd 'init.bin'],       'w', 'ieee-le');
  cnt = fwrite(fid, cat(4, n, u, v), 'real*4', 0,  'ieee-le');
  fclose(fid); clear fid msg cnt n u v;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, 0.,   rhon, topl, 5.,   ...
                dt_o, 0.,   0.,   0.,   0.03, 0.,   hmin, 10.,  10.,  1.,   ...
                1.,   0.,   0.,   0,    0.,   0.,   0.,   [0,0],outd, outd, ...
                'Test-case for lock-exchange' );

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

  xcor = [          xx( 2 : end - 1, ix_0 ); ...
           flipdim( xx( 2 : end - 1, ix_0 ), 1 ); ...
                    xx( 2,           ix_0 ) ] / 1.e3;

  for irec = 1 : nrec
    n    = get_field( 'eta_', irec, outd );

%   Depth (vertical position, >0) of the top of the layers.
    dlay = nan( lm + 2, mm + 2, nlay );

    for ilay = 1 : nlay
      dlay(:,:,ilay) = h_bo ...
                     - squeeze(sum(h_0(:,:,ilay : nlay), 3)) - n(:,:,ilay);
    end; clear ilay;

    hand = fill( xcor, ...
                 [         dlay( 2 : end - 1, ix_0, 1 );      ...
                  flipdim( dlay( 2 : end - 1, ix_0, 2 ), 1 ); ...
                           dlay( 2,           ix_0, 1 )], [0.5 0 0], ...
                 xcor, ...
                 [         dlay( 2 : end - 1, ix_0, 2 ); ...
                           h_bo( 2 : end - 1, ix_0    ); ...
                           dlay( 2,           ix_0, 1 )], [0 0 0.5] );
    title( ['Lock-exchange experiment after ' ...
             num2str( taxi( irec ) ) ' days'] );
    set( hand, 'edgecolor', 'none' );
    ylabel('Depth (m)');
    xlabel('Distance (km)');
    axis ij;
    axis([xx(1, 2) / 1.e3, xx(end, 2) / 1.e3, - 0.5, max(h_bo(:))]);

    if ( taxi( irec ) <= 18. / 24. )
      hold on;
        hand = plot( cint * taxi( irec ) * 24. * 3600. / 1.e3 * [- 1.; 1.], ...
                     [0.25; 0.75] * hmax, '+' );
        set( hand, 'color', [1 1 1] );
      hold off;
    end

    sleep(0.01);
  end; clear irec;
