
% Test-case for upwelling driven by longshore wind, taken from:
%   Morel, Y.G., D.S. Darr, and C. Talandier,
%   Possible Sources Driving the Potential Vorticity Structure and
%     Long-Wave Instability of Coastal Upwelling and Downwelling Currents,
%   J. Phys. Oceanogr., 2006, vol.36, p.875-896.

  clear;
  close all;
  more  off;

% Domain is periodic along `x', and coastline is at y=y_max.

  dl     = 1.e3;           % Mesh size (m)
  nlay   = 2;              % Number of layers.
  topl   = [   0.,   0.5]; % Top     of layers (fraction of Hmax).
  rhon   = [1015., 1030.]; % Density of layers (kg/m3).
  hfla   =  50.;           % Depth of flat basin (meters).
  hsal   = 0.5;            % Salmon's thickness  (meters).
  f0     = 1.e-4;          % Coriolis parameter (rad/s).
  grav   = 9.8;            % Earth's gravity (m/s2).
  rext   = sqrt(grav * hfla) ...
         / abs(f0);        % Barotropic Rossby radius (m).
  ly     = 1. * rext;      % Length of domain in across-shore direction (y, meters  ).
  mm     = round(ly / dl); % Grid size        in across-shore direction (y          ).
  lm     = 1;              % Grid size        in along -shore direction (x, periodic).
  outd   = '/tmp/';        % Path to input/output files.

  h_bo   = zeros(lm + 2, mm + 2);     % Depth of grid (meters).
  h_bo(2 : end - 1, 2 : end - 1) = hfla;

  ndeg   = get_nbr_deg_freedom( h_bo );
  cext   = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).

  hmin   = hsal / 10.;                % Minimum layer thickness allowed (m).

  H_1    =       topl(2)  * hfla;     % Undisturbed thickness of top   layer (m).
  H_2    = (1. - topl(2)) * hfla;     % Undisturbed thickness of lower layer (m).
  R_d    = sqrt(grav * (rhon(2) - rhon(1)) / rhon(2) * H_1 * H_2 / hfla) ...
         / abs(f0);                   % Baroclinic Rossby radius of deformation (m).
  delt   = H_1 / H_2;                 % Ratio of layer thicknesses.
  T_w    = 0.05 / rhon(1) / H_1;      % Wind force applied to top layer (m/s2).
  t_o    = abs(f0) * R_d * (1. + delt) ...
         / T_w;                       % Time when interface reaches surface (s).
  dt_s   = round( 3. * t_o / 3600. / 24. ); % Duration of simulation (days).

% Add body force F_x that mimics the westerly wind stress.

  bodf       = zeros(nlay, 2);        % F_x, F_y (m/s2).
  bodf(1, 1) = T_w;                   % Apply only over top layer.

% Write input files over disk.

  [fid, msg] = fopen([outd 'bodf.bin'], 'w', 'ieee-le');
  cnt = fwrite( fid, bodf, 'real*4',     0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
         0.05 * dt_s, 0.,   0.,   0.,   0.,   0.,   hmin, 10.,  10.,  1.,   ...
                0.,   0.,   1.,   0.,   1.,   0.,   1.,   [0,0],outd, outd, ...
                'Test-case: Upwelling in presence of outcrop' );

% See beom documentation on how to compile and execute the program.

  disp(['Now compile and execute the program in a different terminal. ' ...
        'Press a key when it is over.']);
  pause

  disp('Reading the output data, please wait...');

% Get results from calculation.

% Read metadata.

  [taxi, h_0, f0, dl, rhon, desc] = get_metadata(strtrim(outd));
  nrec = length(taxi);           % Nbr of records in output files.
  lm   = size(h_0, 1) - 2;       % Size of grid.
  mm   = size(h_0, 2) - 2;
  nlay = size(h_0, 3);
  dt_o = taxi(2) - taxi(1);      % Days.

  ix_0 = round(0.5 * (lm + 2));  % Center of grid.
  [xx_r, yy_r] = meshgrid( (     1 : lm + 2) - 1.5, ...
                           (mm + 2 : -1 : 1) - 1.5 );
  xx_r = xx_r' * dl;             % m
  yy_r = yy_r' * dl * (-1.);     % m

  h_bo( find( h_bo < 1.e-3 ) ) = nan;

  nday = round( 0.2 * dt_s ); % Plot results every `nday' days.

  figure;
  subplot(1,2,1);

  for irec = 1 : nrec
    n = get_field('eta_', irec, outd);

    % Layer thickness hlay (meters).

    hlay           = nan(lm + 2, mm + 2, nlay);
    hlay(:,:,nlay) = h_0(:,:,nlay) + n(:,:,nlay);
    for ilay = 1 : nlay - 1
      hlay(:,:,ilay) = h_0(:,:,ilay) + n(:,:,ilay) - n(:,:,ilay + 1);
    end; clear ilay;

    if ( taxi(irec) < (t_o / 24. / 3600.) ) % Before theoretical outcrop.
      U_c  = T_w * taxi(irec) * 24. * 3600. / (1. + delt);        % m/s
      U_b  = T_w * taxi(irec) * 24. * 3600. / (1. + delt) * delt; % m/s
      Y    = 0.;                                                  % Pos. of outcrop (m).
      eta2 = H_1 * U_c * exp(yy_r(ix_0, :) ./ R_d) / f0 / R_d;    % m
      U_1  =       U_c * exp(yy_r(ix_0, :) ./ R_d)        + U_b;  % m/s
      U_2  =     - U_c * exp(yy_r(ix_0, :) ./ R_d) * delt + U_b;  % m/s
    else
      Y    = (-1.) * T_w * (taxi(irec) * 24. * 3600. - t_o) / f0 / (1. + delt);   % m
      U_c  = f0 * R_d * exp(- Y / R_d);                                           % m/s
      U_b  = delt * f0 * (R_d - Y);                                               % m/s
      eta2 = H_1 - H_1 * (1. - U_c * exp(yy_r(ix_0, :) ./ R_d) / f0 / R_d);
      eta2( find( yy_r(ix_0, :) >= Y ) ) = H_1;
      U_1  = U_c * exp(yy_r(ix_0, :) ./ R_d) + U_b;                               % m/s
      U_1(  find( yy_r(ix_0, :) >= Y ) ) = nan;
      U_2  = - delt * U_c * exp(yy_r(ix_0, :) ./ R_d) + U_b;                      % m/s
      tmp  = - f0 * delt * (yy_r(ix_0, :) - Y) - delt * U_c * exp(Y / R_d) + U_b; % m/s
      U_2(  find( yy_r(ix_0, :) >= Y ) ) = tmp( find( yy_r(ix_0, :) >= Y ) );
      clear tmp;
    end

    if ( irec == 1 || mod(irec, round(nday / dt_o)) == 0 )
%     Plot the theoretical and modeled pycnoclines.
%       The theoretical value assumes a rigid lid at the surface.
%       The modeled     value is referenced to the free surface
%         (hence the "+ n(ix_0, :, 1)").
      hand = plot( yy_r(ix_0, :) / 1.e3,                          ...
                   n(ix_0, :, 1) * (- 1.),                        ...
                   yy_r(ix_0, :) / 1.e3,                          ...
                   h_bo(ix_0, :) - h_0(ix_0,:,2) - eta2,          ...
                   yy_r(ix_0, :) / 1.e3,                          ...
                   h_bo(ix_0, :) - h_0(ix_0,:,2) - n(ix_0, :, 2) );
      axis([R_d * (-3.) / 1.e3, 0., 0., 0.6 * hfla]);
      axis ij;
      axis square;
      set(hand(2 : end), 'linewidth', 2);
      if ( irec == 1 )
        legend('Surface (model)', 'Interface (theory)', 'Interface (model)', 'location', 'southwest');
        xlabel('Distance y from coast (km)');
        ylabel('Depth (m)');
        title(['Position of layer interface every ' num2str(nday) ' days.']);
        hold on;
      end
      sleep(0.1);
    end
    if ( irec == nrec )
      hold off;
    end
  end; clear irec;

  u = get_field('u___', nrec, outd);
  v = get_field('v___', nrec, outd);
  m = get_field('mont', nrec, outd);

  subplot(1,2,2);
  hand = plot( [yy_r(ix_0, 2), yy_r(ix_0, end)] / 1.e3, [0., 0.], ...
               yy_r(ix_0, 2 : end) / 1.e3, ...
               hlay(2, 2 : end, 1) * f0 .* u(2, 2 : end, 1) * 1.e3, ...
               yy_r(ix_0, 2 : end) / 1.e3, ...
               hlay(2, 2 : end, 1) * grav / dl ...
                .* ( n(2, 2 : end, 1) - n(2, 1 : end - 1, 1) ) * 1.e3, ...
               yy_r(ix_0, 2 : end) / 1.e3, ...
               hlay(2, 2 : end, 1) * grav / dl ...
                .* ( m(2, 2 : end,     1) - n(2, 2 : end,     1) ...
                   - m(2, 1 : end - 1, 1) + n(2, 1 : end - 1, 1) ) * 1.e3 );
  title(['Dynamical balance within outcrop (t = ' num2str(nrec * dt_o) ' days )']);
  for icur = length( hand(:) ) : - 1 : 2
    set( hand(icur), 'color', get( hand(icur - 1), 'color' ) );
  end; clear icur;
  set( hand(1), 'color', [0 0 0] );
  set( hand(2 : end), 'linewidth', 2 );
  xlabel('Distance y from coast (km)');
  ylabel('10^{-3} m^2 s^{-2}');
  axis([R_d * (-3.) / 1.e3, 0., - 5., 5.]);
  axis square;
  legend('', 'h f u', 'h g d(eta)/dy', 'h g d(M - eta)/dy');
