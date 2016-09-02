
% Test-case for barotropic instability, periodic boundary conditions,
%   momentum advection, and non-linear viscosity.

  clear;
  close all;
  more  off;

  movi = 0;       % Create a movie (optional).
  dl   = 15.e3;   % Mesh size (meters).
  hshf =   5;     % Depth of basin (meters).
  lx   = 3000.e3; % Zonal      length of basin (meters).
  ly   = 4000.e3; % Meridional length of basin (meters)
  rhon = 1030.;
  lm   = round(lx / dl); if ( mod(lm, 2) == 0 ); lm = lm + 1; end
  mm   = round(ly / dl); if ( mod(mm, 2) == 0 ); mm = mm + 1; end
  nlay = 1;
  grav = 9.8;
  fcor = 0.5e-4;  % Coriolis parameter (rad/s).
  outd = '/tmp/'; % Path to input/output files.

  [yy, xx] = meshgrid( (1 : mm + 2) - 1.5, (1 : lm + 2) - 1.5 );
  xx       = xx - mean(xx(:));
  yy       = yy - mean(yy(:));

  h_bo         = hshf + 0.1 * hshf * cos(4. * pi * xx / lm);
  h_bo(find( h_bo < 1. )) = 0.;
  h_bo(  1, :) = 0.;
  h_bo(end, :) = 0.;
  h_bo(:, 1  ) = 0.;
  h_bo(:, end) = 0.;

  cext = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).

  ndeg = get_nbr_deg_freedom( h_bo );

  disp(['Rossby radius of deformation = ' ...
         num2str(sqrt(grav * max(h_bo(:))) / abs(fcor) / 1.e3) ' km.']);

  n        = zeros(lm + 2, mm + 2, nlay);
  u        = zeros(lm + 2, mm + 2, nlay);
  v        = zeros(lm + 2, mm + 2, nlay);
  n(:,:,1) = 1.0 * exp( - yy.^2 / (0.1 * mm).^2 );

  for iy = 2 : mm + 1
    u(:, iy, :) = (n(:, iy + 1, :) - n(:, iy - 1, :)) / (2. * dl) * grav / abs(fcor) * (-1.);
  end; clear iy;

% Write input files on disk.

  [fid, msg] = fopen([outd 'h_bo.bin'], 'w', 'ieee-le');
  cnt = fwrite( fid, h_bo, 'real*4',     0,  'ieee-le');
  fclose(fid); clear fid msg cnt tmp;

  [fid, msg] = fopen([outd 'init.bin'],        'w', 'ieee-le');
  cnt = fwrite( fid, cat(4, n, u, v), 'real*4', 0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, fcor, rhon, 0.,   50.,  ...
                2.,   0.,   0.,   0.,   0.2,  0.,   0.1,  10.,  10.,  1.,   ...
                1.,   0.,   0.,   0.,   1.,   1.,   1.,   [0,0],outd, outd, ...
                'Test-case for barotropic instability' );

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

  for irec = 1 : nrec
    pvor = get_field( 'pvor', irec, outd );
%   u    = get_field( 'u___', irec, outd );
%   v    = get_field( 'v___', irec, outd );

%   u(   1, :) = u( end - 1, : ); % Periodic along `x'.
%   v(   1, :) = v( end - 1, : ); % Periodic along `x'.
%   u(   :, 1) = u( :, end - 1 ); % Periodic along `y'.
%   v(   :, 1) = v( :, end - 1 ); % Periodic along `y'.

%   rvor = nan(lm + 2, mm + 2);
%   rvor(2 : end, 2 : end, :) = v(2 : end,     2 : end,     :) ...
%                             - v(1 : end - 1, 2 : end,     :) ...
%                             + u(2 : end,     1 : end - 1, :) ...
%                             - u(2 : end,     2 : end,     :);
%   rvor = rvor / dl;

    t_pv( irec, 1 ) = mean( reshape( pvor( 2 : end - 1, ...
                                           2 : end - 1 ), ...
                                     lm * mm, 1 ) );
    t_pv( irec, 2 ) = min(  reshape( pvor( 2 : end - 1, ...
                                           2 : end - 1 ), ...
                                     lm * mm, 1 ) );
    t_pv( irec, 3 ) = max(  reshape( pvor( 2 : end - 1, ...
                                           2 : end - 1 ), ...
                                     lm * mm, 1 ) );

    figure(1);
    hand = imagesc( [ xx( 1, 1 ), xx( end, 1 ) ] * dl / 1.e3, ...
                    [ yy( 1, 1 ), yy( 1, end ) ] * dl / 1.e3, ...
                    pvor' );
%                   rvor' ./ abs( fcor ) );
    set( gca, 'xticklabel', [], ...
              'yticklabel', [] );
    axis xy; axis equal;
    title([ 'pvor '    'after ' num2str( round( taxi( irec ) ) ) ' days' ]);
%   title([ 'rvor / |f| after ' num2str( round( taxi( irec ) ) ) ' days' ]);
    hand = colorbar;
    set( get( hand, 'title' ), 'string', '(m s)^{-1}' );
    sleep(0.1);
    if ( movi )
      myft = '/usr/share/fonts/type1/gsfonts/n019003l.pfb';
      print('-dpng', [outd 'frame' num2str(irec, '%04u') '.png'], '-S1024,768', ['-F' myft]);
    end
  end; clear irec;

  figure;
    hand = plot( taxi, t_pv * 1.e6 );
    legend('Average', 'Min', 'Max');
    title('Evolution of potential vorticity');
    xlabel('Days of simulation');
    ylabel('10^{-6} (m s)^{-1}');
