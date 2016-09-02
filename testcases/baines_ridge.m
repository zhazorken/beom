
% Test-case of a rotating flow over a ridge
%   Taken from Baines & Leonard, 1989, Q.J.R. Meteorol. Soc., vol.115, p.293-308.

  clear;
  close all;
  more  off;

  outd = '/tmp/';                              % Path to input/output files.
  hmax = 110.;                                 % Maximum depth of basin (m).
  U_0  = 1.2;                                  % Barotropic speed upstream (m/s).
  nlay = 2;                                    % Number of layers.
  rhon = [1025., 1030.];                       % Density of layers (kg/m3).
  topl = [0., 1. / 1.1];                       % Top     of layers (fraction of Hmax).
  fcor = 1.e-4;                                % Coriolis parameter (rad/s).
  grav = 9.8;                                  % m/s2
  gp   = grav * (rhon(2) - rhon(1)) / rhon(2); % Reduced gravity (m/s2).
  d_0  = hmax * (1. - topl(2));                % Thickness of lower layer (m).
  Lros = sqrt(gp * d_0) / abs(fcor);           % Baroclinic Rossby radius (m).
  lx   = 150. * Lros;                          % Length of domain (m).
  dl   = Lros / 5.;                            % Mesh size (m).
  lm   = round(lx / dl);                       % Grid size in along -stream direction.
  if ( mod(lm, 2) == 0 ); lm = lm + 1; end;    % Make it odd.
  npts = 15;                                   % Number of grid points in the sponge.
  mm   = 1;                                    % Grid size in across-stream direction.

  disp(['Rossby radius = ' num2str(Lros / 1.e3) ' km.']);

  H_m  = 0.1;                                  % Non-dimensional height of topogr.
  F_0  = U_0 / sqrt(gp * d_0);                 % Froude number.
  dksi = dl / Lros;                            % Non-dimensional mesh size.

  disp(['Froude number = ' num2str(F_0)]);

  [xx, yy] = meshgrid( (1 : lm + 2)' - 1.5, (1 : mm + 2) - 1.5 );
  xx   = xx' * dl;
  yy   = yy' * dl;
  xx   = xx - mean(xx(:));                     % Centered on ridge.
  yy   = yy - mean(yy(:));
  X    = xx(:, 2) / Lros;                      % Non-dim. along-stream position.

  h_bo = H_m * d_0 * cos(pi * xx / (10. * Lros)); % Cosine ridge.
  h_bo( find( xx < (- 5. * Lros) | xx > (5. * Lros) ) ) = 0.;
  H    = h_bo(:, 2) / d_0;                     % Non-dimensional topography.
  h_bo = hmax - h_bo;
  h_bo(:,   1) = 0.;
  h_bo(:, end) = 0.;
  h_bo(  1, :) = 0.;
  h_bo(end, :) = 0.;

  ndeg = get_nbr_deg_freedom( h_bo );
  cext = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).

  n    = zeros(lm + 2, mm + 2, nlay   );
  u    = ones( lm + 2, mm + 2, nlay   ) * U_0;
  v    = zeros(lm + 2, mm + 2, nlay   );
  bodf = zeros(                nlay, 2); % F_x, F_y (m/s2).
  nudg = zeros(lm + 2, mm + 2,       3); % eta, u, v.

% Add body force F_y that balances fcor*u.

  bodf(:,2) = fcor * U_0;

% Add wave sponges.

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

% Obtain analytical results (Eq.5.1 to 5.5 in Baines&Leonard 1989).

  Dtil = zeros(lm + 2, 1);
  if     ( F_0 < 1. )
    for ix = 1 : lm + 2
      xval     = X(ix);
      Dtil(ix) = - H(ix) / (1. - F_0^2) + 0.5 * (1. - F_0^2)^(- 1.5) ...
               * ( dksi * sum( exp(   (xval - X(ix : end)) / sqrt(1. - F_0^2) ) .* H(ix : end) ) ...
                 + dksi * sum( exp( - (xval - X( 1 : ix )) / sqrt(1. - F_0^2) ) .* H( 1 : ix ) ) );
      clear xval;
    end; clear ix;
  elseif ( F_0 > 1. )
    for ix = 1 : lm + 2
      xval     = X(ix);
      Dtil(ix) = H(ix) ./ (F_0^2 - 1.) ...
               + (F_0^2 - 1.)^(- 1.5) * dksi ...
               * sum( H(1 : ix) .* sin( (X(1 : ix) - xval) / sqrt(F_0^2 - 1.) ) );
      clear xval;
    end; clear ix;
  end

  disp(['min(D / F_0^(2./3.)) = ' num2str(min((1. + Dtil) / F_0^(2. / 3.)))]);
  disp(['max(D / F_0^(2./3.)) = ' num2str(max((1. + Dtil) / F_0^(2. / 3.)))]);

% Convert analytical results to dimensional form (d_an, u_an, v_an).

  d_an = (1. + Dtil(:)) * d_0; % meters.
  u_an = U_0 * d_0 ./ d_an(:); % m/s.
  V    = zeros(lm + 2, 1);
  V(2 : end - 1) = F_0 * u_an(2 : end - 1) ./ U_0 ...
                .* (u_an(3 : end) - u_an(1 : end - 2)) ./ (2. * dksi) ./ U_0 ...
                 + F_0^(-1.) * ( d_an(3 : end    ) ./ d_0 + H(3 : end    )   ...
                               - d_an(1 : end - 2) ./ d_0 - H(1 : end - 2) ) / (2. * dksi);
  v_an = V(:) ./ U_0;
  clear Dtil H V X dksi d_0 gp;

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

  [fid, msg] = fopen([outd 'bodf.bin'], 'w', 'ieee-le');
  cnt = fwrite( fid, bodf, 'real*4',     0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, fcor, rhon, topl, 10.,  ...
                0.2,  0.,   0.,   0.,   0.,   0.,   0.1,  10.,  10.,  1.,   ...
                1.,   0.,   0.,   0,    0.,   1.,   0.,   [0,0],outd, outd, ...
                'Test-case for flow over a ridge' );

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

    hand = plot( xx(:, 2) / 1.e3, h_bo(:, ix_0), ...
                 xx(:, 2) / 1.e3, h_bo(:, ix_0) - d_an, ...
                 xx(:, 2) / 1.e3, zlay(:, ix_0, 2) );
    xlabel('Distance from ridge (km)');
    ylabel('Depth (m)');
    set(hand, 'linewidth', 1.5);
    set(hand(1), 'color', [0 0 0]);
    legend('Sea floor topography', 'Position of interface (analytical)', 'Position of interface (model)', ...
           'location', 'east');
    axis ij;
    axis([xx(1, 2) / 1.e3, xx(end, 2) / 1.e3, hmax * topl(2) - 2., max(h_bo(:))]);
    axis square;
    sleep(0.01);
  end; clear irec;
