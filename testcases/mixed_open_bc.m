
% Test-case for mixed open boundary conditions in case of an upwelling driven
%   by seaward wind. Analytical solution taken from
%   C. Millot and M. Crepon,
%   Inertial oscillations on the continental shelf of the Gulf of Lions
%   --Observations and Theory,
%   J. Phys. Oceanogr., 1981, vol.11, p.639-657.

  clear;
  close all;
  more  off;

  lm     = 200;       % Grid points in across-shore direction.
  mm     = 100;       % Grid points in along -shore direction.
  nlay   = 2;
  dl     = 1.e3;
  f0     = 1.e-4;
  dt_o   = 0.16667;   % Time between model outputs (days).
  hfla   = 40.;
  topl   = [0., 0.5];
  rhon   = [1028.95, 1030.];
  tauw   = [0.1, 0.]; % Blowing toward the East. Coast is along western edge.
  grav   = 9.8;
  outd   = '/tmp/';
  npts   = 15;        % Number of grid points in the sponge.

% `nudg' is rate at which the model variables are relaxed at open boundaries.
%  nudg  is non-dimensional and defined at the center of the cells

%  The dimensional relaxation time-scale is obtained as:
%    tau = dt / nudg, where dt is the model time-step (in seconds).

  nudg   = zeros(lm + 2, mm + 2, 3   ); % eta, u, v
  h_bo   = zeros(lm + 2, mm + 2);       % Depth of grid (meters).

  h_bo(2 : end - 1, 2 : end - 1) = hfla;
  ndeg   = get_nbr_deg_freedom( h_bo );

% Celerity of surface  gravity waves (m/s).
  cext   = sqrt(grav * max(h_bo(:)));

% Add sponge for long surface gravity waves at offshore open boundary (East).
% Add weak relaxation at the two other open boundaries (North and South).
% Coast is along western edge. Zero nudging is interpreted as a closed boundary.
% By default, model fields are relaxed toward the state of rest.

  east = zeros( lm + 2, mm + 2, 3 );
  nort = zeros( lm + 2, mm + 2, 3 );
  west = zeros( lm + 2, mm + 2, 3 );
  sout = zeros( lm + 2, mm + 2, 3 );

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

%   Add a weak relaxation on tangential velocities (for long integrations).
    east(i,:, 3    ) = dt / (31. * 24. * 3600.) * xpos / npts; % One month.
  end; clear i xpos;

  for j = 1 : mm + 2 % Northern boundary.
    xpos        = j - 1.5 + npts - mm;
    xpos        = max([xpos;        0. ]);
    xpos        = min([xpos; npts - 0.5]);
    nort(:,j,:) = dt / (31. * 24. * 3600.) * xpos / npts; % One month.
  end; clear j xpos;

  for j = 1 : mm + 2 % Southern boundary.
    xpos        = npts - (j - 1.5);
    xpos        = max([xpos;        0. ]);
    xpos        = min([xpos; npts - 0.5]);
    sout(:,j,:) = dt / (31. * 24. * 3600.) * xpos / npts; % One month.
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

  nudg(1,:,:) = 0.; % The western edge is dry (coast; no nudging).

% Write input files over disk.

  [fid, msg] = fopen([outd 'nudg.bin'], 'w', 'ieee-le');
  cnt = fwrite( fid, nudg, 'real*4',     0,  'ieee-le');
  fclose(fid); clear fid msg cnt;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, 6.,   ...
                dt_o, 4.,   0.,   0.,   0.,   0.,   0.001,1.,   1.,   1.,   ...
                0.,   0.,   0.,   0,    0.,   0.,   0.,   tauw, outd, outd, ...
                'Test-case: Upwelling seaward wind with mixed open boundary conditions' );

% See beom documentation on how to compile and execute the program.

  disp(['Now compile and execute the program in a different terminal. ' ...
        'Press a key when it is over.']);
  pause

  disp('Reading the output data, please wait...');

% Get results from calculation.

% Read metadata.

  [taxi, h_0, f0, dl, rhon, desc] = get_metadata(strtrim(outd));
  nrec = length(taxi);         % Nbr of records in output files.
  lm   = size(h_0, 1) - 2;     % Size of grid.
  mm   = size(h_0, 2) - 2;
  nlay = size(h_0, 3);
  ix_0 = round( 0.5 * (mm + 2) );

% Fields for analytical solution.

  gpri = grav * (rhon(2) - rhon(1)) / rhon(2);
  r_1  = sqrt( grav * hfla ) / f0;
  r_2  = sqrt( gpri * h_0(:,:,1) .* h_0(:,:,2) ./ hfla ) / abs(f0);
  r_2  = max(r_2(:));
  [xx_r, yy_r] = meshgrid( (1 : lm + 2) - 1.5, ...
                           (1 : mm + 2) - 1.5 );
  xx_r = xx_r';
  yy_r = yy_r';

% Calculate theoretical upwelling.

  t_eta = nan(lm + 2, nlay);
  t_v   = nan(lm + 2, nlay);
  t_eta(:, 2) = exp(- xx_r(:,ix_0) * dl / r_2) ...
              * (-1.) * tauw(1) / rhon(1) / (r_1 * f0^2) .* h_0(:,ix_0,2) / hfla ...
              * (-1.) * r_1 / r_2;
  t_eta(:, 1) = (-1.) .* tauw(1) / rhon(1) / (r_1 * f0^2) ...
              * ( exp(- xx_r(:,ix_0) * dl / r_1) ...
                + h_0(:,ix_0,2) ./ h_0(:,ix_0,1) * r_2 / r_1 ...
               .* exp(- xx_r(:,ix_0) * dl / r_2) ...
                );
  t_v(:, 1) = tauw(1) / rhon(1) / f0 ./ hfla .* h_0(:,ix_0,2) ./ h_0(:,ix_0,1) ...
           .* ( exp(- xx_r(:,ix_0) * dl / r_2) - 1. );
  t_v(:, 2) = t_v(:, 1) .* (h_0(:,ix_0,1) ./ h_0(:,ix_0,2)) * (-1.);

  m_eta = nan(lm + 2, nrec, nlay);
  m_u   = nan(lm + 2, nrec, nlay);
  m_v   = nan(lm + 2, nrec, nlay);

  for irec = 1 : nrec
    eta = get_field('eta_', irec, outd);
    u   = get_field('u___', irec, outd);
    v   = get_field('v___', irec, outd);
    for ilay = 1 : nlay
      m_eta(:, irec, ilay) = eta(:, ix_0, ilay);
      m_u(  :, irec, ilay) =   u(:, ix_0, ilay);
      m_v(  :, irec, ilay) =   v(:, ix_0, ilay);
    end; clear ilay;
  end; clear irec;

% Theory-model comparison.

figure(1);
  hand = plot( xx_r(2 : end, ix_0) * dl / 1.e3, t_eta(2 : end, 2), ...
               xx_r(  :,     ix_0) * dl / 1.e3, m_eta(  :,1 : round(0.5 / dt_o) : end, 2) );
  legend('Steady analytical solution', '12-hourly model');
  axis([0. 30. 0. 1.4]);
  set(hand(1), 'linewidth', 20, 'color', [0.8 0.8 0.8]);
  set(hand(2 : end), 'color', [0 0 0]);
  xlabel('Distance from shore (km)');
  ylabel('Interface displacement (m)');

figure;
  hand = plot( xx_r(2 : end, ix_0) * dl / 1.e3, t_v(2 : end, 2) * 1.e2, ...
               xx_r(  :,     ix_0) * dl / 1.e3, m_v(:,1 : round(0.5 / dt_o) : end, 2) * 1.e2 );
  legend('Steady analytical solution','12-hourly model');
  axis([0. 30. 0. 3.]);
  set(hand(1), 'linewidth', 20, 'color', [0.8 0.8 0.8]);
  set(hand(2 : end), 'color', [0 0 0]);
  xlabel('Distance from shore (km)');
  ylabel('v_2 (cm s^{-1})');
