
% Test-case for beta effect, taken from:
% H. Stommel, 1948, The westward intensification of wind-driven ocean currents,
% Transactions of the American Geophysical Union, Volume 29, Number 2.

  clear;
  close all;
  more  off;

% R   is in m   s^(-1) <=> qdrg * u = 0.02 cm s^(-1) = 2.e-4 m s^(-1)
% psi is in m^2 s^(-1)
% See Stommel 1948 for definition of symbols.

  outd    = '/tmp/';         % Path to input/output files.
  F_PLANE = 0;               % Binary flag (0/1); 1 for uniform Earth rotation.

  dl      = 100.e3;          % Mesh size (meters).
  lambda  = 1.e7;            % Zonal      extent of domain (meters).
  b       = 2. * pi * 1.e6;  % Meridional extent of domain (meters).
  lm      = round(lambda / dl);
  mm      = round(b      / dl);
  rhon(1) = 1027.;           % kg m^(-3)
  hfla    = 200;             % Depth of basin (flat bottom, meters).
  F       = 0.1 / rhon(1);   % m^2 s^(-2)
  R       = 2.e-4;           % Linear bottom drag coefficient (m s^(-1)).
  if ( F_PLANE )
    beta  = 0.;              % m^(-1) s^(-1)
    fmin  = 0.25e-4;         % s^(-1)
%   fmin  = 0.;              % s^(-1)
  else
    beta  = 1.e-11;          % m^(-1) s^(-1)
    fmin  = 0.;              % s^(-1)
  end
  grav    = 9.8;             % m s^(-2)
  alpha   = hfla * beta / R; % m^(-1)
  gamma   = F * pi / R / b;  % s^(-1)
  A       = - 0.5 * alpha + sqrt(0.25 * alpha^2 + (pi / b)^2);
  B       = - 0.5 * alpha - sqrt(0.25 * alpha^2 + (pi / b)^2);
  p       = (1. - exp(B * lambda)) ./ (exp(A * lambda) - exp(B * lambda));
  q       = 1. - p;
  [xx,yy] = meshgrid( (0.5 : lm - 0.5) * dl, (0.5 : mm - 0.5) * dl );
  xx      = xx'; yy = yy';         % x as first dimension.
  h_bo    = zeros(lm + 2, mm + 2); % Depth of grid (meters).
  h_bo(2 : end - 1, 2 : end - 1) = hfla;

  ndeg = get_nbr_deg_freedom( h_bo );

  fcor = nan( lm, mm );
  for j = 1 : mm
    fcor(:, j) = fmin + (j - 0.5) * dl * beta;
  end; clear j;

  tausx = - rhon(1) * F * cos(pi * yy ./ b);
  tausy = zeros( size(tausx,1), size(tausx,2) );

  tmp   = nan(lm + 2, mm + 2);
  tmp(2 : end - 1, 2 : end - 1) = fcor;
  tmp(1,  :) = tmp(2,:);
  tmp(:,  1) = tmp(:,2);
  tmp(end,:) = tmp(end - 1,:);
  tmp(:,end) = tmp(:,end - 1);

  [fid, msg] = fopen('/tmp/fcor.bin', 'w', 'ieee-le');
  cnt = fwrite( fid, tmp, 'real*4',    0,  'ieee-le');
  fclose(fid); clear fid msg cnt tmp;

  tmp = nan(lm + 2, mm + 2, 2);
  tmp(2 : end - 1, 2 : end - 1, 1) = tausx(:,:);
  tmp(2 : end - 1, 2 : end - 1, 2) = tausy(:,:);
  tmp(1,  :,:) = tmp(2,      :,:);
  tmp(:,  1,:) = tmp(:,      2,:);
  tmp(end,:,:) = tmp(end - 1,:,:);
  tmp(:,end,:) = tmp(:,end - 1,:);

  [fid, msg] = fopen('/tmp/taus.bin', 'w', 'ieee-le');
  cnt = fwrite( fid, tmp, 'real*4',    0,  'ieee-le');
  fclose(fid); clear fid msg cnt tmp;

  cext = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   1,    ndeg, dl,   cext, 0.,   rhon, 0.,   40.,  ...
                1.,   0.,   0.,   0.,   0.,   R,    1.,   10.,  10.,  1.,   ...
                0.,   0.,   0.,   0.,   0.,   0.,   0.,   [0,0],outd, outd, ...
                'Test-case for Stommel 1948' );

% See beom documentation on how to compile and execute the program.

  disp(['Now compile and execute the program in a different terminal. ' ...
        'Press a key when it is over.']);
  pause

  disp('Reading the output data, please wait...');

  psi = gamma * (b / pi)^2 * sin(pi * yy / b) ...
     .* (p * exp(A * xx) + q * exp(B * xx) - 1.); % m^2 s^(-1)

  eta = - F / grav / hfla * (exp(A * xx) * p / A + exp(B * xx) * q / B) ...
      - (b / pi)^2 * F / grav / hfla * (p * A * exp(A * xx) + q * B * exp(B * xx)) .* (cos(pi * yy / b) - 1.) ...
      - (fcor * gamma / grav * (b / pi)^2 .* sin(pi * yy / b) + 1*beta * gamma / grav * (b / pi)^3 * (cos(pi * yy / b) - 1.)) ...
     .* (p * exp(A * xx) + q * exp(B * xx) - 1.);

figure(1);
  [ccon, hcon] = contourf(xx' / dl, yy' / dl, psi' / 1.e4, 20. );
  title('Analytical streamfunction (10^4 m^2 s^{-1})');
  xlabel('Zonal grid points');
  ylabel('Meridional grid points');
  axis xy;
  axis equal;
  colorbar

figure(2);
  eta = eta - eta(1,1);
  [ccon, hcon] = contourf(xx' / dl, yy' / dl, eta' * 1.e2, 20 );
  title('Analytical sea surface height (cm)');
  xlabel('Zonal grid points');
  ylabel('Meridional grid points');
  axis([0, size(xx,1) + 1, 0, size(xx,2) + 1]);
  axis xy;
  axis equal;
  colorbar

% Read metadata.

  [taxi, h_0, f0, dl, rhon, desc] = get_metadata(strtrim(outd));
  nrec = length(taxi);         % Nbr of records in output files.
  lm   = size(h_0, 1) - 2;     % Size of grid.
  mm   = size(h_0, 2) - 2;
  nlay = size(h_0, 3);

% Animate the sea surface height.

  for irec = 1 : nrec
    eta  = get_field('eta_', irec, outd);
    figure(3);
    hand = imagesc( squeeze(eta(:,:,1))' * 1.e2 );
    axis xy; axis equal;
    title(['Sea surface elevation (cm) after ' num2str(taxi(irec) - taxi(1)) ' days']);
    xlabel('Zonal grid points');
    ylabel('Meridional grid points');
    colorbar
    sleep(0.01);
  end

figure;
  eta(:,:,1) = eta(:,:,1) - eta(2,2,1); % Reference to value at southwest corner.
  [xx, yy] = meshgrid( (1 : size(eta,1)) - 1.5, (1 : size(eta,2)) - 1.5 );
  [ccon, hcon] = contourf( xx(2 : end - 1, 2 : end - 1), ...
                           yy(2 : end - 1, 2 : end - 1), ...
                  squeeze(eta(2 : end - 1, 2 : end - 1, 1))' * 1.e2, 20 );
  axis([0, size(eta,1) - 1, 0, size(eta,2) - 1]);
  axis xy; axis equal;
  title(['Sea surface elevation (cm) after ' num2str(taxi(irec) - taxi(1)) ' days']);
  xlabel('Zonal grid points');
  ylabel('Meridional grid points');
  colorbar
