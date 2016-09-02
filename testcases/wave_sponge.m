
% Test-case for wave sponge at open boundaries (Flow Relaxation Scheme).
%   See Lavelle, J.W. & W.C. Thacker (2008), Ocean Modelling, vol.20, p.270-292.
%       Martinsen, E.A. & H. Engedahl (1987), Coastal Engineering, vol.11, p.603-627.
%       Modave, A., E. Deleersnijder, E.J. Delhez, Ocean Dynamics, vol.60, p.65-79.
%       Davies, H.C. (1983), Monthly Weather Rev., vol.111, p.1002-1012.

  clear;
  close all;
  more  off;

  outd   = '/tmp/';
  lx     = 600.e3; % meters
  ly     = 600.e3; % meters
  dl     = 10.e3;
  npts   = 15;     % Number of grid points in the sponge.
  nlay   = 2;
  fcor   = 1.e-4;
  rhon   = [1000.0, 1030.0];
  topl   = [0.,       0.5];
  hfla   =  200.;  % meters
  grav   = 9.8;
  lm     = round(lx / dl);
  mm     = round(ly / dl);
  if ( mod(lm, 2) == 0 ); lm = lm + 1; end; % Make it odd.
  if ( mod(mm, 2) == 0 ); mm = mm + 1; end; % Make it odd.
  lm     = lm + 2 * npts;
  mm     = mm + 2 * npts;
  n      = zeros(lm + 2, mm + 2, nlay);
  u      = zeros(lm + 2, mm + 2, nlay);
  v      = zeros(lm + 2, mm + 2, nlay);
  nudg   = zeros(lm + 2, mm + 2, 3   ); % eta, u, v
  h_bo   = zeros(lm + 2, mm + 2      );
  h_bo(2 : end - 1, 2 : end - 1) = hfla;

  disp(['Rossby radius = ' ...
        num2str(sqrt(grav * (rhon(2) - rhon(1)) / rhon(2) ...
        * topl(2) * (1. - topl(2)) * hfla) / abs(fcor) / 1.e3) ' km.']);
  cext     = sqrt(grav * hfla); % (m s^(-1)).

  ndeg     = get_nbr_deg_freedom( h_bo );

  [xx, yy] = meshgrid( (1 : lm + 2) - 1.5, (1 : mm + 2) - 1.5 );
  xx       = xx' - median(xx(1,:) );
  yy       = yy' - median(yy(:,1)');
  xx       = xx * dl;
  yy       = yy * dl;

  n(:,:,1) = 1. * exp(- (xx.^2 + yy.^2) ./ (50.e3)^2);

% Add wave sponges.

  east = zeros( lm + 2, mm + 2, 3 );
  nort = zeros( lm + 2, mm + 2, 3 );
  west = zeros( lm + 2, mm + 2, 3 );
  sout = zeros( lm + 2, mm + 2, 3 );
  dt   = 0.5 * dl / cext; % Model timestep (seconds).
  widt = npts * dl;       % Width of sponge (meters).

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

  for j = 1 : mm + 2 % Northern boundary.
    xpos             = j - 1.5 + npts - mm;
    xpos             = max([xpos;        0. ]);
    xpos             = min([xpos; npts - 0.5]);
    nort(:,j,[1, 3]) = dt * cext / widt * xpos / (npts - xpos);
  end; clear j xpos;

  for j = 1 : mm + 2 % Southern boundary.
    xpos             = npts - (j - 1.5);
    xpos             = max([xpos;        0. ]);
    xpos             = min([xpos; npts - 0.5]);
    sout(:,j,[1, 3]) = dt * cext / widt * xpos / (npts - xpos);
  end; clear j xpos;

% For each variable (eta,u,v),
%   take the maximum of eastern/western/northern/southern sides
%   (e.g., Fig.2, Lavelle & Thacker 2008).

  for ivar = 1 : size( nudg, 3 )
    temp = max( [reshape(nudg(:,:,ivar), [1, (lm + 2) * (mm + 2)]); ...
                 reshape(east(:,:,ivar), [1, (lm + 2) * (mm + 2)]); ...
                 reshape(west(:,:,ivar), [1, (lm + 2) * (mm + 2)]); ...
                 reshape(nort(:,:,ivar), [1, (lm + 2) * (mm + 2)]); ...
                 reshape(sout(:,:,ivar), [1, (lm + 2) * (mm + 2)])] );
    nudg(:,:,ivar) = reshape( temp', lm + 2, mm + 2 );
    clear temp;
  end; clear ivar;

% Write input files over disk.

  [fid, msg] = fopen([outd 'init.bin'],        'w', 'ieee-le');
  cnt = fwrite( fid, cat(4, n, u, v), 'real*4', 0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

  [fid, msg] = fopen([outd 'nudg.bin'], 'w', 'ieee-le');
  cnt = fwrite( fid, nudg, 'real*4',     0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, fcor, rhon, topl, 0.25, ...
             4.17e-3, 0.,   0.,   0.,   0.,   0.,   1.,   10.,  10.,  1.,   ...
                1.,   0.,   0.,   0.,   0.,   0.,   0.,   [0,0],outd, outd, ...
                'Test-case for wave sponge' );

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

  pote = nan(nrec, nlay);
  kine = nan(nrec, nlay);

  for irec = 1 : 1 : nrec
    eta  = get_field('eta_', irec, outd);
    u    = get_field('u___', irec, outd);
    v    = get_field('v___', irec, outd);
    hlay(  :,:,2) = h_0(:,:,2) + eta(:,:,2);              % Thickness (meters).
    hlay(  :,:,1) = h_0(:,:,1) + eta(:,:,1) - eta(:,:,2); % Thickness (meters).
    pote(irec, 1) = 0.5 * rhon(1) * grav ...
                  * nansum( reshape(eta(:,:,1), (lm + 2) * (mm + 2), 1).^2 );
    kine(irec, 1) = nansum( reshape(u(:,:,1).^2 .* hlay(:,:,1), (lm + 2) * (mm + 2), 1) ) ...
                  + nansum( reshape(v(:,:,1).^2 .* hlay(:,:,1), (lm + 2) * (mm + 2), 1) );
    kine(irec, 2) = nansum( reshape(u(:,:,2).^2 .* hlay(:,:,2), (lm + 2) * (mm + 2), 1) ) ...
                  + nansum( reshape(v(:,:,2).^2 .* hlay(:,:,2), (lm + 2) * (mm + 2), 1) );
    kine(irec, 1) = kine(irec, 1) * rhon(1) * 0.5;
    kine(irec, 2) = kine(irec, 2) * rhon(2) * 0.5;
    clear hlay;

% Animate the sea surface height.
    if ( mod(irec, 5) == 0 )
      figure(1);
      subplot(1,2,1);
      hand = imagesc( squeeze(eta(:,:,1))' * 1.e2 );
      axis xy; axis equal;
      title(['Sea surface elevation (cm) after ' num2str(taxi(irec) - taxi(1)) ' days']);
      xlabel('Zonal grid points');
      ylabel('Meridional grid points');
      caxis([-10 10]);
      colorbar
      hold on;
        hand = plot( [npts, npts, (lm + 2) - npts, (lm + 2) - npts, npts], ...
                     [npts, (lm + 2) - npts, (lm + 2) - npts, npts, npts] );
        set(hand, 'color', [1 1 1]);
      hold off;

      subplot(1,2,2);
      hand = imagesc( squeeze(eta(:,:,2))' * 1.e2 );
      axis xy; axis equal;
      title(['Interface elevation (cm) after ' num2str(taxi(irec) - taxi(1)) ' days']);
      xlabel('Zonal grid points');
      ylabel('Meridional grid points');
      caxis([-50 50]);
      colorbar
      hold on;
        hand = plot( [npts, npts, (lm + 2) - npts, (lm + 2) - npts, npts], ...
                     [npts, (lm + 2) - npts, (lm + 2) - npts, npts, npts] );
        set(hand, 'color', [1 1 1]);
      hold off;
      sleep(0.01);
    end
  end

  disp('Plotting the results:');

figure;
  hand = plot( taxi * 24.,  pote(:,1)                          / 1.e3, ...
               taxi * 24., (kine(:,1) + kine(:,2)            ) / 1.e3, ...
               taxi * 24., (pote(:,1) + kine(:,1) + kine(:,2)) / 1.e3 );
  set(hand, 'linewidth', 6);
  legend('Potential energy', 'Kinetic energy', 'Sum');
  xlabel('Time (hours)');
  ylabel('10^3 J m^{-2}');
  grid on;
