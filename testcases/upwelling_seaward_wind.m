
% Test-case for upwelling driven by seaward wind, taken from
% C. Millot and M. Crepon,
% Inertial oscillations on the continental shelf of the Gulf of Lions--Observations and Theory,
% J. Phys. Oceanogr., 1981, vol.11, p.639-657.

  clear;
  close all;
  more  off;

  lm     = 200;
  mm     = 1;
  nlay   = 2;
  dl     = 1.e3;
  f0     = 1.e-4;
  dt_o   = 0.16667; % days
  hfla   = 40.;
  topl   = [0., 0.5];
  rhon   = [1028.95, 1030.];
  tauw   = [0.1, 0.];
  grav   = 9.8;
  outd   = '/tmp/';

  h_bo   = zeros(lm + 2, mm + 2); % Depth of grid (meters).
  h_bo(2 : end - 1, 2 : end - 1) = hfla;

  ndeg   = get_nbr_deg_freedom( h_bo );
  cext   = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, 6.,   ...
                dt_o, 4.,   0.,   0.,   0.,   0.,   1.,   10.,  10.,  1.,   ...
                0.,   0.,   0.,   0.,   0.,   1.,   0.,   tauw, outd, outd, ...
                'Test-case for upwelling seaward wind' );

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
  t_eta(:, 2) = exp(- xx_r(:,2) * dl / r_2) ...
              * (-1.) * tauw(1) / rhon(1) / (r_1 * f0^2) .* h_0(:,2,2) / hfla ...
              * (-1.) * r_1 / r_2;
  t_eta(:, 1) = (-1.) .* tauw(1) / rhon(1) / (r_1 * f0^2) ...
              * ( exp(- xx_r(:,2) * dl / r_1) ...
                + h_0(:,2,2) ./ h_0(:,2,1) * r_2 / r_1 ...
               .* exp(- xx_r(:,2) * dl / r_2) ...
                );
  t_v(:, 1) = tauw(1) / rhon(1) / f0 ./ hfla .* h_0(:,2,2) ./ h_0(:,2,1) ...
           .* ( exp(- xx_r(:,2) * dl / r_2) - 1. );
  t_v(:, 2) = t_v(:, 1) .* (h_0(:,2,1) ./ h_0(:,2,2)) * (-1.);

  m_eta = nan(lm + 2, nrec, nlay);
  m_u   = nan(lm + 2, nrec, nlay);
  m_v   = nan(lm + 2, nrec, nlay);

  for irec = 1 : nrec
    eta = get_field('eta_', irec, outd);
    u   = get_field('u___', irec, outd);
    v   = get_field('v___', irec, outd);
    for ilay = 1 : nlay
      m_eta(:, irec, ilay) = eta(:, 2, ilay);
      m_u(  :, irec, ilay) =   u(:, 2, ilay);
      m_v(  :, irec, ilay) =   v(:, 2, ilay);
    end; clear ilay;
  end; clear irec;

% Theory-model comparison.

figure(1);
  hand = plot( xx_r(2 : end, 2) * dl / 1.e3, t_eta(2 : end, 2), ...
               xx_r(  :,     2) * dl / 1.e3, m_eta(  :,1 : round(0.5 / dt_o) : end, 2) );
  legend('Steady analytical solution', '12-hourly model');
  axis([0. 30. 0. 1.4]);
  set(hand(1), 'linewidth', 20, 'color', [0.8 0.8 0.8]);
  set(hand(2 : end), 'color', [0 0 0]);
  xlabel('Distance from shore (km)');
  ylabel('Interface displacement (m)');

figure;
  hand = plot( xx_r(2 : end, 2) * dl / 1.e3, t_v(2 : end, 2) * 1.e2, ...
               xx_r(  :,     2) * dl / 1.e3, m_v(:,1 : round(0.5 / dt_o) : end, 2) * 1.e2 );
  legend('Steady analytical solution','12-hourly model');
  axis([0. 30. 0. 3.]);
  set(hand(1), 'linewidth', 20, 'color', [0.8 0.8 0.8]);
  set(hand(2 : end), 'color', [0 0 0]);
  xlabel('Distance from shore (km)');
  ylabel('v_2 (cm s^{-1})');
