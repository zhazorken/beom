
% Test-case for `wetting & drying'.
% Taken from:
%   G.F. Carrier and H.P. Greenspan,
%   Water waves of finite amplitude on a sloping beach,
%   JFM 1958, vol.4, p.97-109.

  clear;
  close all;
  more  off;

  movi = 0;                       % Create a movie from model results.
  alph = 1.e-3;                   % alph = - h* / x*  (p.97).
  l_0  = 3.e3;                    % Lengthscale       (meters).
  epsi = 0.1;                     % Non-dim. height of mound.
  dl   = l_0 / 50.;               % Mesh size         (meters).
  hhti = 3. * epsi * alph * l_0;  % Highest High Tide (meters).
% hhti = 2. * epsi * alph * l_0;  % Highest High Tide (meters).
  l_x  = 10. * l_0 ...            % Length of domain  (meters)
       + hhti / alph;             %   including dry area.
  lm   = round(l_x / dl);
  mm   = 1;
  hsal = 0.2 * alph * dl;         % Salmon's thickness (m).
  hmin = hsal / 10.;              % Min. layer thickness allowed (m).
  dt_s = 0.08;                    % Duration of calculation (days).
  ix_0 = ceil(0.5 * (mm + 2));    % Center of domain.
  nlay = 1;                       % Number of layers.
  grav = 9.8;                     % m/s2
  v_0  = sqrt(grav * l_0 * alph); % Velocity scale (m/s).
  T    = v_0 / (alph * grav);     % Timescale (s) (p.97).
  p    = 1. / 8. / (1. + epsi);   % Page 104.
  npts = 15;                      % Number of grid points in the sponge.
  outd = '/tmp/';                 % Path to input/output files.

% Bathymetry is a sloping beach.

  h_bo = repmat( (lm + 1 : - 1 : 0)' * dl * alph, [1, mm + 2] );
  h_bo(1 : npts, ix_0) = h_bo(npts, ix_0); % Enforce a flat bottom within sponge.
  h_bo( :,  1 ) = 0.;
  h_bo( :, end) = 0.;
  h_bo(end, : ) = 0.;
  h_bo( 1,  : ) = 0.;

  ndeg = get_nbr_deg_freedom( h_bo );

  disp(['Mesh size is '                 num2str(    dl      ) ' meters.']);
  disp(['Maximum depth over domain is ' num2str(max(h_bo(:))) ' meters.']);

  i_sl = find( abs(h_bo(:, ix_0)' - hhti) == min( abs(h_bo(:, ix_0)' - hhti) ) );
  i_sl = i_sl(1);
  hhti = h_bo( i_sl, ix_0 );
  topl = hhti / max( h_bo(:) );       % Position of mean sea level.
  cext = sqrt( grav * max(h_bo(:)) ); % Celerity of surface gravity waves (m/s).
  xref = (1 : lm + 2)' * dl;
  xref = xref - xref( i_sl );         % Cross-shore pos. referenced to shoreline (m).

% Initial condition for surface elevation `eta'.
% Eqs. 3.14, 3.15.

  sigm = (10. : - 1.e-1 : 0.)';

  x    = 0.25 * epsi * exp(2) * p^2 * sigm.^4 .* exp(- sigm.^2 * p) ...
       - sigm.^2 / 16.;
  eta  = 0.25 * epsi * p^2 * exp(2) * sigm.^4 .* exp(- sigm.^2 * p);
  n    = interp1( x * l_0, eta * alph * l_0, xref );
  n    = repmat( n, [1, mm + 2] );
  temp = nan(lm + 2, mm + 2, nlay);
  temp(:,:,1) = n; clear temp;
  n( find( isnan(n) ) ) = 0.;

  figure;
  hand = fill( [ xref(1 + npts),          ...
                 xref(end - 1),           ...
                 xref(end - 1),           ...
                 xref(1 + npts) ] / 1.e3, ...
               [ h_bo(1 + npts, ix_0),    ...
                 h_bo(end - 1,  ix_0),    ...
                 h_bo(1 + npts, ix_0),    ...
                 h_bo(1 + npts, ix_0) ] * (-1.) + hhti, [0.8 0.8 0] );
  set( hand, 'edgecolor', 'none' );
  hold on;
    hand = plot( [xref(1), 0.], [0., 0.], ...
                 x * l_0 / 1.e3, eta * alph * l_0 );
    set( hand(1), 'color', [0 0 0] );
    set( hand(2), 'color', [0 0 1] );
    set( hand(2), 'linewidth', 3 );
    legend('Sloping beach', 'Mean sea level', ...
           'Sea surface', 'location', 'southwest');
    title('Wave at t = 0 hours');
    xlabel('Distance from shoreline (km)');
    ylabel('Vertical displacement (meters)');
    axis([- 20., xref(end) / 1.e3, - hhti, hhti]);
  hold off;
  clear sigm x eta;

% Time-evolution of the shoreline.
% Eqs. 3.23 - 3.29.

  dlam = 1.e-3;
  lamb = (0. : dlam : 20.)';
  lmid = 0.5 * (lamb(1 : end - 1) + lamb(2 : end));
  E_la = [0.; cumsum( exp(0.25 * lmid.^2) ) * dlam]; % Eq.3.29.
  f_la = - lamb.^2 + 0.5 * lamb.^4     ...           % Eq.3.27.
       + exp(- 0.25 * lamb.^2) .* E_la ...
      .* (lamb + lamb.^3 - 0.25 * lamb.^5);
  f_la = f_la / 16.;
  dfdl = - lamb + 3. * lamb.^3 - 0.25 * lamb.^5 ...  % Eq.3.28.
       + exp(- 0.25 * lamb.^2) .* E_la          ...
      .* (1. + 2.5 * lamb.^2 - 7. * lamb.^4 / 4. + 0.125 * lamb.^6);
  dfdl = dfdl / 16.;

  v_sl = sqrt(pi * p) * epsi * exp(2) * dfdl;    % Eq.3.23.
  t_sl = 0.25 * lamb / sqrt(p) ...               % Eq.3.24.
       - sqrt(pi * p) * epsi * exp(2) * dfdl;
  x_sl = - v_sl.^2 / 16.       ...               % Eq.3.26.
       + epsi * exp(2) * sqrt(pi) * 0.25 * f_la;

  u    = zeros(lm + 2, mm + 2, nlay );
  v    = zeros(lm + 2, mm + 2, nlay );
  nudg = zeros(lm + 2, mm + 2, 3    ); % eta, u, v.

% Add wave sponge on western open boundary.

  west = zeros( lm + 2, mm + 2, 3 );
  dt   = 0.5  * dl / cext; % Model timestep (seconds).
  widt = npts * dl;        % Width of sponge (meters).

  for i = 1 : lm + 2 % Western boundary.
%  `xpos' is the position within the sponge zone (in units of grid points).
    xpos             = npts - (i - 1.5);
    xpos             = max([xpos;        0. ]);
    xpos             = min([xpos; npts - 0.5]);
%   Lavelle and Thacker 2008, Ocean Modelling, vol.20 p.270-292:
%     for eastern boundary, relax `eta' (1) and `u' (2), but not `v' (3).
%   The nudging coefficients (non-dimenzionalized by timestep dt)
%     vary as Eq.29 of Modave et al. 2010 Ocean Dynamics vol.60 p.65-79.
    west(i,:,[1, 2]) = dt * cext / widt * xpos / (npts - xpos);
  end; clear i xpos;

% For each variable (eta,u,v),
%   take the maximum of all sides
%   (e.g., Fig.2, Lavelle & Thacker 2008).

  for ivar = 1 : size( nudg, 3 )
    temp = max( [reshape(nudg(:,:,ivar), [1, (lm + 2) * (mm + 2)]); ...
                 reshape(west(:,:,ivar), [1, (lm + 2) * (mm + 2)])] );
    nudg(:,:,ivar) = reshape( temp', lm + 2, mm + 2 );
    clear temp;
  end; clear ivar;

% Write input files over disk.

  [fid, msg] = fopen([outd 'h_bo.bin'], 'w', 'ieee-le');
  cnt = fwrite(fid, h_bo, 'real*4',      0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

  [fid, msg] = fopen([outd 'init.bin'],       'w', 'ieee-le');
  cnt = fwrite(fid, cat(4, n, u, v), 'real*4', 0,  'ieee-le');
  fclose(fid); clear fid msg cnt n u v;

  [fid, msg] = fopen([outd 'nudg.bin'], 'w', 'ieee-le');
  cnt = fwrite( fid, nudg, 'real*4',     0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, 0.,   1030, topl, dt_s, ...
                5.e-4,0.,   0.,   0.,   0.,   0.,   hmin, 10.,  10.,  1.,   ...
                1.,   0.,   1.,   0.,   0.,   0.,   0.,   [0,0],outd, outd, ...
                'Test-case for wave on sloping beach' );

% See beom documentation on how to compile and execute the program.

  disp(['Now compile and execute the program in a different terminal. ' ...
        'Press a key when it is over.']);
  pause

  disp('Reading the output data, please wait...');

% Read metadata.

  [taxi, h_0, f0, dl, rhon, desc] = get_metadata( strtrim(outd) );
  nrec = length(taxi);         % Nbr of records in output files.
  lm   = size(h_0, 1) - 2;     % Size of grid.
  mm   = size(h_0, 2) - 2;
  nlay = size(h_0, 3);
  xmod = nan(nrec, 1);
  umod = nan(nrec, 1);

  if ( movi )
    figure;
  end
  for irec = 1 : nrec
    n    = get_field('eta_', irec, outd);
    u    = get_field('u___', irec, outd);
    hlay = h_0 + n;

    i_sl = find( hlay(:, ix_0) < hsal );
    xmod(irec) = interp1( hlay(i_sl(1) - 1 : i_sl(1), ix_0), ...
                          xref(i_sl(1) - 1 : i_sl(1)      ), ...
                          hsal );
    umod(irec) = interp1( hlay(i_sl(1) - 1 : i_sl(1), ix_0), ...
                          u(   i_sl(1) - 1 : i_sl(1), ix_0), ...
                          hsal );
    xitp = interp1( t_sl(:) * T, x_sl(:) * l_0, taxi(irec) * 24. * 3600. );
    zitp = interp1( xref,      - h_bo(:, ix_0) + hhti, xitp );
    clear i_sl;

    if ( movi )
      hand = fill(  [ xref(1 + npts),          ...
                      xref(end - 1),           ...
                      xref(end - 1),           ...
                      xref(1 + npts) ] / 1.e3, ...
                    [ h_bo(1 + npts, ix_0),    ...
                      h_bo(end - 1,  ix_0),    ...
                      h_bo(1 + npts, ix_0),    ...
                      h_bo(1 + npts, ix_0) ] * (-1.) + hhti, [0.8 0.8 0] );
      set( hand, 'edgecolor', 'none' );
      title(['Wave at t = ' num2str(taxi(irec) * 24., '%2.2f') ' hours']);
      hold on;
        hand = plot( [xref(1), 0.], [0., 0.], ...
                     xitp / 1.e3, zitp, 'x',  ...
                     xref / 1.e3, - h_bo(:, ix_0) + hhti + hlay(:, ix_0) );
        set( hand(1),       'color', [0 0 0] );
        set( hand(2),       'color', [1 0 0] );
        set( hand(3),       'color', [0 0 1] );
        set( hand(2 : end), 'linewidth', 3 );
        legend('Sloping beach', 'Mean sea level', 'Shoreline (theory)', ...
               'Sea surface (model)', 'location', 'southwest');
        xlabel('Distance from mean shoreline (km)');
        ylabel('Vertical displacement (m)');
        axis([- 20., xref(end) / 1.e3, - hhti, hhti]);
      hold off;
      sleep( 1.e-3 );
      myft = '/usr/share/fonts/type1/gsfonts/n019003l.pfb';
      print('-dpng', [outd 'frame' num2str(irec, '%04u') '.png'], '-S1024,768', ['-F' myft]);
    end
  end; clear irec;

  figure;
  hand = plot( [taxi(1), taxi(end)] * 24., [0., 0.], ...
               t_sl * T / 3600., x_sl * l_0, ...
               taxi(:) * 24., xmod(:) );
  legend('', 'Shoreline (theory)', 'Shoreline (model)');
  xlabel('Time (hours)');
  ylabel('Horizontal displacement (meters)');
  set( hand(2 : end), 'linewidth', 3 );
  set( hand(1), 'color', [0 0 0] );
  set( hand(2), 'color', [1 0 0] );
  set( hand(3), 'color', [0 0 1] );
  axis([0., round(dt_s * 24.), - 300., 500.]);
% grid on;
  axis square;
