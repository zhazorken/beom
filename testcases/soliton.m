
% Test-case of an equatorial soliton, taken from
% J.W. Lavelle & W.C. Thacker, 2008, Ocean Modelling, vol.20, p.270-292.

  clear;
  close all;
  more  off;

  outd  = '/tmp/';                 % Path to input/output files.
  rhon  = 1029.;                   % Density of layers (kg/m3).
  nlay  = 1;                       % Number of layers.
  lx    = 0.5 * 12.24e6;           % Zonal      extent of domain (meters).
  ly    = 1.5 * 2.04e6;            % Meridional extent of domain (meters).
  dl    = 20.e3;                   % Mesh size (meters).
  H     = 1.;                      % Equivalent barotropic depth (meters).
  grav  = 9.8;                     % m s^(-2)
  a_rd  = 6371.e3;                 % Earth's radius (meters).
  omeg  = 2. * pi / (24. * 3600.); % Earth's rotation rate (rad s^(-1)).
  x_0   = 0.;                      % Initial position of soliton (meters).
  parB  = 0.394;
  parA  = 0.772 * parB^2;
  Elam  = 4. * omeg^2 * a_rd^2 / grav / H;
  L_ls  = a_rd / Elam^(0.25);      % Horizontal length scale  (meters).

  lm    = round(lx / dl);
  mm    = round(ly / dl);
  if ( mod(lm, 2) == 0 )
    lm = lm + 1; % Make it odd.
  end
  if ( mod(mm, 2) == 0 )
    mm = mm + 1; % Make it odd.
  end

  h_bo    = zeros(lm + 2, mm + 2); % Depth of grid (meters).
  h_bo(2 : end - 1, 2 : end - 1) = H;

  ndeg = get_nbr_deg_freedom( h_bo );

  n        = zeros(lm + 2, mm + 2, nlay);
  u        = zeros(lm + 2, mm + 2, nlay);
  v        = zeros(lm + 2, mm + 2, nlay);

  [xx, yy] = meshgrid( ((1 : lm + 2) - 1.5) * dl, ...
                       ((1 : mm + 2) - 1.5) * dl );
  xx       = xx' - mean(xx(:));
  yy       = yy' - mean(yy(:));

  n(:,:,1) = parA             * H  * sech(parB * (xx - x_0) ./ L_ls).^2 ...
          .* (6. * yy.^2 + 3. * L_ls^2) ./ (4. * L_ls^2) ...
          .* exp(- yy.^2 ./ (2. * L_ls^2));
  u(:,:,1) = parA * sqrt(grav * H) * sech(parB * (xx - x_0) ./ L_ls).^2 ...
          .* (6. * yy.^2 - 9. * L_ls^2) ./ (4. * L_ls^2) ...
          .* exp(- yy.^2 ./ (2. * L_ls^2));
  v(:,:,1) = - 2. * parA * parB * sqrt(grav * H) * tanh(parB * (xx - x_0) ./ L_ls) ...
          .* sech(parB * (xx - x_0) ./ L_ls).^2 ...
          .* 2. .* yy ./ L_ls .* exp(- yy.^2 ./ (2. * L_ls^2));

  hlay = h_bo + squeeze( n );
  cext = sqrt( grav * max(hlay(:)) ); % Celerity of surface gravity waves (m/s).
  clear hlay;

  deld = 0.25 * ly / 40.e6;
  beta = 2. * omeg * sind(  deld) ...
       - 2. * omeg * sind(- deld);
  beta = beta / (2. * deld * 40.e6 / 360.);
  fcor = beta * yy; clear deld beta;

% Write input files on disk.

  [fid, msg] = fopen([outd 'h_bo.bin'], 'w', 'ieee-le');
  cnt = fwrite(fid, h_bo, 'real*4',      0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

  [fid, msg] = fopen([outd 'fcor.bin'], 'w', 'ieee-le');
  cnt = fwrite( fid, fcor, 'real*4',     0,  'ieee-le');
  fclose(fid); clear fid msg cnt tmp;

  [fid, msg] = fopen([outd 'init.bin'], 'w',       'ieee-le');
  cnt = fwrite( fid, cat(4, n, u, v), 'real*4', 0, 'ieee-le');
  fclose(fid); clear fid msg cnt;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, 0.,   rhon, 0.,   60.,  ...
                2.,   0.,   0.,   0.,   0.,   0.,   0.05, 10.,  10.,  1.,   ...
                1.,   0.,   0.,   0.,   1.,   0.,   0.,   [0,0],outd, outd, ...
                'Test-case for equatorial soliton' );

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

% Animate the sea surface height.

  for irec = 1 : nrec
    eta  = get_field('eta_', irec, outd);
    figure(1);
    hand = imagesc( eta' * 1.e2 );
    axis xy; axis equal;
    title(['Sea surface elevation (cm) after ' num2str(taxi(irec) - taxi(1)) ' days']);
    xlabel('Zonal grid points');
    ylabel('Meridional grid points');
    colorbar
    sleep(0.01);
  end

figure;
  [xx, yy] = meshgrid( (1 : size(eta,1)) - 1.5, (1 : size(eta,2)) - 1.5 );
subplot(2,1,1);
  eta0 = get_field('eta_', 1, outd);
  cint = (- 2. : 1. : nanmax(eta0(:) * 1.e2));
  [ccon, hcon] = contourf(  xx(2 : end - 1, 2 : end - 1), ...
                            yy(2 : end - 1, 2 : end - 1), ...
                  squeeze(eta0(2 : end - 1, 2 : end - 1, 1))' * 1.e2, cint );
  axis([0, size(eta,1) - 1, 0, size(eta,2) - 1]);
  axis xy; axis equal;
  title(['Sea surface elevation (cm) after ' num2str(taxi(1)) ' days']);
  ylabel('Meridional grid points');
subplot(2,1,2);
  [ccon, hcon] = contourf( xx(2 : end - 1, 2 : end - 1), ...
                           yy(2 : end - 1, 2 : end - 1), ...
                  squeeze(eta(2 : end - 1, 2 : end - 1, 1))' * 1.e2, cint );
  axis([0, size(eta,1) - 1, 0, size(eta,2) - 1]);
  axis xy; axis equal;
  title(['Sea surface elevation (cm) after ' num2str(taxi(irec) - taxi(1)) ' days']);
  xlabel('Zonal grid points');
  ylabel('Meridional grid points');
