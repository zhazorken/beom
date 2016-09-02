
% Test-case for integral conservation of volume, energy, potential enstrophy
%   and vorticity during the collapse of a Gaussian-shaped mound of water
%   in a rotating, stratified, square basin with bathymetry
%   contained inside the deepest layer.
% See Sadourny            (1975), J.Atmos.Sci., vol. 32, p.680-689,
%     Ketefian & Jacobson (2009), J.Comp.Phys., vol.228, p.  1- 32.

  clear;
  close all;
  more  off;

  outc = 0.;                      % Set to `1.' to include isopycnal outcropping.
  topo = 1;                       % Set to `1'  to include bottom topography.
  xper = 1;                       % `1' for periodic domain along `x' axis.
  yper = 1;                       % `1' for periodic domain along `y' axis.
  outd = '/tmp/';                 % Path to input/output files.
  lx   = 600.e3;                  % Length of domain (meters).
  dl   = 10.e3;                   % Mesh size        (meters).
  nlay = 2;                       % Number of layers.
  fcor = 1.e-4;                   % Coriolis parameter (rad/s).
  rhon = 1000. ...
       + 30. * (1 : nlay) / nlay; % Density of layers (kg/m3).
  topl = ((1 : nlay) - 1) / nlay; % Top     of layers (fraction of Hmax).
  hfla = 200.;                    % Depth of (flat) basin (meters).
  grav = 9.8;                     % Gravity (m/s2).
  lm   = round(lx / dl);
  if ( mod(lm, 2) == 0 )
    lm = lm + 1;                  % Make it odd.
  end
  mm   = lm;
% mm   =  1;

  n    = zeros(lm + 2, mm + 2, nlay);
  u    = zeros(lm + 2, mm + 2, nlay);
  v    = zeros(lm + 2, mm + 2, nlay);

  [xx, yy] = meshgrid( (1 : lm + 2) - 1.5, (1 : mm + 2) - 1.5 );
  xx       = xx' - mean(xx(1,:) );
  yy       = yy' - mean(yy(:,1)');
  xx       = xx * dl;
  yy       = yy * dl;

% Add topography (a seamount centered over the domain) to verify
%   that topography does not influence the conservation of properties.

  h_bo = hfla * ones(lm + 2, mm + 2);

  if ( topo )
    if ( outc )
%     Seamount rises half-way through layer (nlay - 1).
      h_bo = hfla - 0.5 * (2. - topl(end) - topl(end - 1)) * hfla ...
           * exp(- (xx.^2 + yy.^2) ./ (0.25 * lx)^2); % Seamount (m).
    else
%     Seamount rises half-way through the lowest layer.
      h_bo = hfla - 0.5  * (1. - topl( end )) * hfla ...
           * exp(- (xx.^2 + yy.^2) ./ (0.25 * lx)^2); % Seamount (m).
    end
  end
  h_bo(  1,  :) = 0.;
  h_bo(end,  :) = 0.;
  h_bo(  :,  1) = 0.;
  h_bo(  :,end) = 0.;

  ndeg = get_nbr_deg_freedom( h_bo );

  if ( nlay == 2 )
    disp(['Rossby radius = ' ...
          num2str(sqrt(grav * (rhon(2) - rhon(1)) / rhon(2) ...
          * topl(2) * (1. - topl(2)) * hfla) / abs(fcor) / 1.e3) ' km.']);
  end

% Add the mound of water at the top of the surface layer.

  n(:,:,1) = 1. * exp(- (xx.^2 + yy.^2) ./ (lx / 12.)^2); % (meters).

  cext     = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).

% Write input files over disk.

  if ( topo )
    [fid, msg] = fopen([outd 'h_bo.bin'], 'w', 'ieee-le');
    cnt = fwrite(fid, h_bo,            'real*4', 0, 'ieee-le');
    fclose(fid); clear fid msg cnt;
  end

  [fid, msg] = fopen([outd 'init.bin'], 'w', 'ieee-le');
  cnt = fwrite(fid, cat(4, n, u, v), 'real*4', 0, 'ieee-le');
  fclose(fid); clear fid msg cnt n u v;

% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, fcor, rhon, topl, 50.,  ...
                0.5,  0.,   0.,   0.,   0.,   0.,   1.,   10.,  10.,  0.,   ...
                1.,   0.,   outc, 0.,   xper, yper, 1.,   [0,0],outd, outd, ...
                'Test-case for integral conservation of properties' );

% See beom documentation on how to compile and execute the program.

  disp(['Now compile and execute the program in a different terminal. ' ...
        'Press a key when it is over.']);
  pause

  disp('Reading the output data, please wait...');

% Read metadata.

  [taxi, h_0, f0, dl, rhon, desc] = get_metadata(strtrim(outd));
  nrec = length(taxi);     % Nbr of records in output files.
  lm   = size(h_0, 1) - 2; % Size of grid.
  mm   = size(h_0, 2) - 2;
  nlay = size(h_0, 3);

  volu = nan(nrec, nlay); % Domain-averaged thickness of layers               (m).
  vstd = nan(nrec, nlay); % Standard-deviation from domain-averaged thickness (m).
  enst = nan(nrec, nlay); % Domain-averaged potential enstrophy (m^{-1} s^{-2}).
  estd = nan(nrec, nlay); % Standard-deviation from domain-averaged pot. enstr.
  rvor = nan(nrec, nlay); % Domain-averaged relative vorticity (s^{-1}).
  rstd = nan(nrec, nlay); % Standard-deviation from domain-averaged rvor.
  pote = nan(nrec, nlay); % Domain-integrated barotropic potential energy (J/m2).
  kine = nan(nrec, nlay); % Domain-integrated kinetic energy              (J/m2).

  for irec = 1 : nrec
    n    = get_field('eta_', irec, outd);
    u    = get_field('u___', irec, outd);
    v    = get_field('v___', irec, outd);
    pvor = get_field('pvor', irec, outd);

    % Layer thickness hlay (meters).

    hlay           = nan(lm + 2, mm + 2, nlay);
    hlay(:,:,nlay) = h_0(:,:,nlay) + n(:,:,nlay);
    for ilay = 1 : nlay - 1
      hlay(:,:,ilay) = h_0(:,:,ilay) + n(:,:,ilay) - n(:,:,ilay + 1);
    end; clear ilay;

    if ( irec == 1 )
      hl_0 = hlay;
    end

    volu(irec, :) = nanmean( reshape(     hlay,         (lm + 2) * (mm + 2), nlay) );
    vstd(irec, :) = nanmean( reshape( abs(hlay - hl_0), (lm + 2) * (mm + 2), nlay) );
%   vstd(irec, :) = nanstd(  reshape(hlay, (lm + 2) * (mm + 2), nlay) );

%   CAREFUL, we only consider the barotropic potential energy
%     (as in a shallow-water model).
    pote(irec, 1) = 0.5 * rhon(1) * grav * dl^2 ... % Joules
                  * nansum( reshape(n(:,:,1), (lm + 2) * (mm + 2), 1).^2 );

    if ( xper == 1 )
      hlay( 1, :,:) = hlay(lm + 1, :, :);
    end
    if ( yper == 1 )
      hlay(:, 1, :) = hlay(:, mm + 1, :);
    end

    hatu = zeros(lm + 2, mm + 2, nlay); % Thickness at `u'   grid points (left  edge).
    hatv = zeros(lm + 2, mm + 2, nlay); % Thickness at `v'   grid points (lower edge).
    hatp = zeros(lm + 2, mm + 2, nlay); % Thickness at `psi' grid points (lower-left).
    hzer = hlay; hzer( find( isnan(hzer) ) ) = 0.;
    uzer = u;    uzer( find( isnan(uzer) ) ) = 0.;
    vzer = v;    vzer( find( isnan(vzer) ) ) = 0.;

    mask = hlay; mask( find(~ isnan(mask)) ) = 1;
    mask( find( isnan(mask) ) ) = 0;
    mask(2 : end, 2 : end, :) = mask(1 : end - 1, 1 : end - 1, :) ...
                              + mask(1 : end - 1, 2 : end,     :) ...
                              + mask(2 : end,     1 : end - 1, :) ...
                              + mask(2 : end,     2 : end,     :);

    hatu(2 : end, :, :) = 0.5 * ( hzer(1 : end - 1, :, :) + hzer(2 : end, :, :) );
    hatv(:, 2 : end, :) = 0.5 * ( hzer(:, 1 : end - 1, :) + hzer(:, 2 : end, :) );

    utmp = uzer.^2 .* hatu;
    vtmp = vzer.^2 .* hatv;
    utmp = 0.5 * ( utmp(1 : end - 1, :, :) + utmp(2 : end, :, :) );
    vtmp = 0.5 * ( vtmp(:, 1 : end - 1, :) + vtmp(:, 2 : end, :) );
    utmp = reshape( utmp, (lm + 1) * (mm + 2), nlay );
    vtmp = reshape( vtmp, (lm + 2) * (mm + 1), nlay );
    kine(irec, :) = 0.5 * sum( utmp, 1 ) + 0.5 * sum( vtmp, 1 );
    kine(irec, :) = squeeze(kine(irec, :)) .* rhon(:)' * dl^2; % Joules.
    clear uzer vzer hatu hatv utmp vtmp;

    mask( find( mask == 0 ) ) = 1;
    hatp(2 : end, 2 : end, :) = ( hzer(1 : end - 1, 1 : end - 1, :) ...
                                + hzer(1 : end - 1, 2 : end,     :) ...
                                + hzer(2 : end,     1 : end - 1, :) ...
                                + hzer(2 : end,     2 : end,     :) ...
                                ) ...
                             ./ mask(2 : end, 2 : end, :);
    clear hzer mask;

    if ( xper == 1 )
      pvor(end,:,:) = nan; % Do not account twice for points along boundary.
    end
    if ( yper == 1 )
      pvor(:,end,:) = nan; % Do not account twice for points along boundary.
    end

    if ( irec == 1 )
      en_0 = pvor.^2 .* hatp * 0.5;
    end

    enst(irec, :) = nanmean( reshape(     pvor.^2 .* hatp * 0.5,         (lm + 2) * (mm + 2), nlay) );
    estd(irec, :) = nanmean( reshape( abs(pvor.^2 .* hatp * 0.5 - en_0), (lm + 2) * (mm + 2), nlay) );
%   estd(irec, :) = nanstd(  reshape(pvor.^2 .* hatp * 0.5,  (lm + 2) * (mm + 2), nlay) );
    rvor(irec, :) = nanmean( reshape(pvor    .* hatp - fcor, (lm + 2) * (mm + 2), nlay) );
    rstd(irec, :) = nanstd(  reshape(pvor    .* hatp - fcor, (lm + 2) * (mm + 2), nlay) );
  end

  disp('Plotting the results:');

figure;
subplot(3,1,1);
  hand = plot( [taxi(1); taxi(end)], [0.; 0.], ...
               taxi, (volu - repmat(volu(1,:), [nrec, 1])) );
  rgb  = get(hand, 'color');
  for i = length(hand) : - 1 : 2
    set(hand(i), 'color', cell2mat(rgb(i - 1, :)));
  end; clear i rgb;
  set(hand(1), 'color', [0 0 0], 'linewidth', 0.5);
  set(hand(2 : end), 'linewidth', 2);
  title('Volume normalized by area');
  legend( [' '; num2str( (1 : nlay)' )], 'orientation', 'horizontal' );
  ylabel('meters');
  set(gca, 'xticklabel', '');
  grid on;
subplot(3,1,2);
  hand = plot( taxi, (vstd - repmat(vstd(1,:), [nrec, 1])) );
  title('Volume anomaly normalized by area');
  set(hand, 'linewidth', 2);
  legend( num2str( (1 : nlay)' ), 'orientation', 'horizontal' );
  ylabel('meters');
  set(gca, 'xticklabel', '');
  grid on;
subplot(3,1,3);
  hand = plot( taxi, pote( :,  1),  ...
               taxi, sum(kine, 2), ...
               taxi, (pote(:,1) + sum(kine, 2)) );
  set(hand, 'linewidth', 2);
  title('Mechanical energy');
  legend('P', 'K', 'P + K', 'orientation', 'horizontal');
  xlabel('Time (days)');
  ylabel('Joules');
  grid on;

figure;
subplot(4,1,1);
  hand = plot( [taxi(1); taxi(end)], [0.; 0.], ...
               taxi, (enst - repmat(enst(1,:), [nrec, 1])) );
  rgb  = get(hand, 'color');
  for i = length(hand) : - 1 : 2
    set(hand(i), 'color', cell2mat(rgb(i - 1, :)));
  end; clear i rgb;
  set(hand(1), 'color', [0 0 0], 'linewidth', 0.5);
  set(hand(2 : end), 'linewidth', 2);
  title('Potential enstrophy normalized by area');
  legend( [' '; num2str( (1 : nlay)' )], 'orientation', 'horizontal' );
  ylabel('m^{-1} s^{-2}');
  set(gca, 'xticklabel', '');
  grid on;
subplot(4,1,2);
  hand = plot( taxi, (estd - repmat(estd(1,:), [nrec, 1])) );
  title('Potential enstrophy anomaly normalized by area');
  legend( num2str( (1 : nlay)' ), 'orientation', 'horizontal' );
  set(hand, 'linewidth', 2);
  ylabel('m^{-1} s^{-2}');
  set(gca, 'xticklabel', '');
  grid on;
subplot(4,1,3);
  hand = plot( [taxi(1); taxi(end)], [0.; 0.], ...
               taxi, (rvor - repmat(rvor(1,:), [nrec, 1])) );
  rgb  = get(hand, 'color');
  for i = length(hand) : - 1 : 2
    set(hand(i), 'color', cell2mat(rgb(i - 1, :)));
  end; clear i rgb;
  set(hand(1), 'color', [0 0 0], 'linewidth', 0.5);
  set(hand(2 : end), 'linewidth', 2);
  title('Vorticity normalized by area');
  legend( [' '; num2str( (1 : nlay)' )], 'orientation', 'horizontal' );
  ylabel('s^{-1}');
  set(gca, 'xticklabel', '');
  grid on;
subplot(4,1,4);
  hand = plot( taxi, (rstd - repmat(rstd(1,:), [nrec, 1])) );
  title('Vorticity anomaly normalized by area');
  set(hand, 'linewidth', 2);
  legend( num2str( (1 : nlay)' ), 'orientation', 'horizontal' );
  xlabel('Time (days)');
  ylabel('s^{-1}');
  grid on;
