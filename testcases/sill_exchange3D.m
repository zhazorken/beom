% Test-case of a sill exchange over a ridge (3-D).

  clear;
  close all;
  more  off;
  ocrp = 1;                                    % Allow isopycnal outcrops.
  movi = 1;                                    % Create frames for a movie.
  name='hsal5day30';                        % Name for movie file
  outd = '/tmp/';                              % Path to input/output files.
  hmax = 700.;                                 % Maximum depth of basin (m).
  dt_s = 30;                                    % Duration of calc.     (days).
  dt_r = 0.;                                   % Ramp on wind forcing (days).
  bdrg = 0.;                                   % Bottom drag coefficient.
  nlay = 2;                                    % Number of layers.
  fcor = 0.00014087;                           % Coriolis parameter (rad/s).
  lx   = 50.e3;                                % Length of domain (m).
  ly   = 200.e3;                               % Width of domain (m).
  dl   = 400.;                                 % Mesh size (m).
  lm   = round( lx / dl );                     % Grid size in along -stream direction.
  mm   = round( ly / dl );                     % Grid size in across-stream direction.
  if ( mod(lm, 2) == 0 ); lm = lm + 1; end;    % Make it odd.
  if ( mod(mm, 2) == 0 ); mm = mm + 1; end;    % Make it odd.
  npts = 15;                                   % Number of grid points in the sponge.
  grav = 9.8;                                  % m/s2
  rhon = [1027.47, 1027.75];                   % Density of layers (kg/m3).
  hsill= 400;                                  % Height of the sill maximum.   
  topl = [0.0, 0.1428];                        % Top of layers (fraction of Hmax).

  [xx, yy] = meshgrid( (1 : lm + 2)' - 1.5, (1 : mm + 2) - 1.5 );
  xx   = xx' * dl;
  yy   = yy' * dl;
  xx   = xx - mean(xx(:)); % Centered on ridge.
  yy   = yy - mean(yy(:));

  h_bo = hsill* exp( - (yy ./ (50. * dl)).^2 );
  h_bo = hmax - h_bo;

  ndeg = get_nbr_deg_freedom( h_bo );
  cext = sqrt(grav * max(h_bo(:))); % Celerity of surface gravity waves (m/s).
  cint = 0.5 * sqrt( grav * hmax * (rhon(2) - rhon(1)) / rhon(2) ); %Celerity of internal gravity waves (m/s).
    
  dhdy = zeros( lm + 2, mm + 2 ); %interface slope used to calculated salmon layer condition
  dhdy( :, 3 : end - 2 ) = h_bo( :, 4 : end - 1 ) ...
                         - h_bo( :, 2 : end - 3 );
  dhdy = dhdy / (2. * dl);

  U_0=1;
  n    = zeros(lm + 2, mm + 2, nlay   );
  u    = zeros( lm + 2, mm + 2, nlay   );
  v    = zeros(lm + 2, mm + 2, nlay   )*U_0;
  % Salmon's thickness hsal (m).

  hsal = 5; %20. * max( dhdy(:) ) * dl * (rhon(end) - rhon(1)) / rhon(1);
  hmin = hsal / 10.;            % Min. layer thickness allowed (m).


%    %initial layer height
%    n( 1 : round( 0.5 * (lm + 2) ),           :, 2 ) =   0.5 * hmax - 4. * hsal-100;
%                         %  - repmat( hmax - h_bo( 1 : round( 0.5 * (lm + 2) ) , 2 ), 1, 3 );
%    n(     round( 0.5 * (lm + 2) ) + 1 : end, :, 2 ) =   -0.5 * hmax + 4. * hsal+450 ...
%        + repmat( hmax - h_bo( round( 0.5 * (lm + 2) ) + 1 : end, 2 ), 1, 3 );
%     n(     round( 0.5 * (lm + 2) ) + 1 : end, :, 2 )= ...
%         min( n(     round( 0.5 * (lm + 2) ) + 1 : end, :, 2 ),0.5 * hmax - 4. * hsal-100);
n( :, 1 : round( 0.6 * (mm + 2) ), 2 ) =0;
n( :,     round( 0.6 * (mm + 2) ) + 1 : end, 2 ) = ...
    min(-h_bo(:, round( 0.6 * (mm + 2) ) + 1 : end)+hmax-100+4.*hsal,0);


  %h_bo(:,   1) = 0.;
%  h_bo(:, end) = 0.;
% h_bo(  1, :) = 0.;
 % h_bo(end, :) = 0.;
  
  %bodf = zeros(                nlay, 2); % F_x, F_y (m/s2).
  nudg = zeros(lm + 2, mm + 2,       3); % eta, u, v.
  

% Add wave sponges.

  nudg = zeros( lm + 2, mm + 2, 3    ); % eta, u, v.

  east = zeros( lm + 2, mm + 2, 3 );
  nort = zeros( lm + 2, mm + 2, 3 );
  west = zeros( lm + 2, mm + 2, 3 );
  sout = zeros( lm + 2, mm + 2, 3 );
  
  dt   = 0.5  * dl / cext; % Model timestep (seconds).
  widt = npts * dl;        % Width of sponge (meters).

 for j = 1 : mm + 2 % Northern boundary.
%  `xpos' is the position within the sponge zone (in units of grid points).
    xpos        = j - 1.5 + npts - mm;
    xpos        = max([xpos;        0. ]);
    xpos        = min([xpos; npts - 0.5]);
%   Lavelle and Thacker 2008, Ocean Modelling, vol.20 p.270-292:
%     for eastern boundary, relax `eta' (1) and `u' (2), but not `v' (3).
%   The nudging coefficients (non-dimenzionalized by timestep dt)
%     vary as Eq.29 of Modave et al. 2010 Ocean Dynamics vol.60 p.65-79.
%Here we relax u,v but not eta.
    nort(:,j,[2, 3]) = dt * cext / widt * xpos / (npts - xpos);

%   Add a weak relaxation on tangential velocities (for long integrations).
    nort(:, j, 1) = dt / (31. * 24. * 3600.) * xpos / npts; % One month.
  end; clear j xpos;


  for j = 1 : mm + 2 % Southern boundary.
    xpos        = npts - (j - 1.5);
    xpos        = max([xpos;        0. ]);
    xpos        = min([xpos; npts - 0.5]);
    sout(:,j,[2, 3]) = dt * cext / widt * xpos / (npts - xpos);
    sout(:,j,1) = dt / (31. * 24. * 3600.) * xpos / npts; % One month.
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
  nudg(end,:,:) = 0.; % The eastern edge is dry (coast; no nudging).


% Write input files over disk.

  [fid, msg] = fopen([outd 'init.bin'],       'w', 'ieee-le');
  cnt = fwrite(fid, cat(4, n, u, v), 'real*4', 0,  'ieee-le');
  fclose(fid); clear fid msg cnt n u v;
  
  [fid, msg] = fopen([outd 'h_bo.bin'], 'w', 'ieee-le');
  cnt = fwrite(fid, h_bo,            'real*4', 0, 'ieee-le');
  fclose(fid); clear fid msg cnt;

   [fid, msg] = fopen([outd 'nudg.bin'], 'w', 'ieee-le');
   cnt = fwrite( fid, nudg, 'real*4',     0,  'ieee-le');
   fclose(fid); clear fid msg cnt;


% print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, dt_s, ...
%               dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, hbbl, g_fb, ...
%               uadv, qdrg, ocrp, rsta, xper, yper, diag, tauw, idir, odir, ...
%               desc );
  print_params( lm,   mm,   nlay, ndeg, dl,   cext, fcor, rhon, topl, dt_s, ...
                0.01, dt_r, 0.,   0.,   0.9,  bdrg, hmin, 5.,   5.,   1.,   ...
                1.,   1.,   ocrp, 0.,   0.,   0.,   0.,   [0,0],outd, outd, ...
                'Test-case: 3D sill exchange' );

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

  iy_0 = ceil(0.5 * (lm + 2));
  
  if (movi)
  vidObj = VideoWriter(strcat(name,'.avi'));
  vidObj2 = VideoWriter(strcat(name,'2.avi'));
  open(vidObj);
  open(vidObj2);
  end 
  
  for irec = 1 :3:1500%nrec
    n    = get_field('eta_', irec, outd);

    u    = get_field( 'u___', irec, outd );
    v    = get_field( 'v___', irec, outd );
    %plot(u(:,2,1))
    %plot(v(:,2,1))
    %calculate total mass
    mass(irec)=sum(n(2:end-1,2,2))+100000;


    zlay = nan(lm + 2, mm + 2, nlay);

    for ilay = 1 : nlay
      zlay(:,:,ilay) = h_bo - squeeze(sum(h_0(:,:,ilay : nlay), 3)) - n(:,:,ilay);
    end; clear ilay;

    figure(1)

        hand= plot(        yy(iy_0, 1 + npts : end - npts ) / 1.e3, ...
                          h_bo(iy_0, 1 + npts : end - npts ), ...
                            yy(3, 1 + npts : end - npts ) / 1.e3, ...
                      squeeze( zlay(3, 1 + npts : end - npts, :) ), ... 
                            yy(iy_0, 1 + npts : end - npts ) / 1.e3, ...
                      squeeze( zlay(iy_0, 1 + npts : end - npts, 2) ), ...
                            yy(lm, 1 + npts : end - npts ) / 1.e3, ...
                      squeeze( zlay(lm, 1 + npts : end - npts, 2) ) );
             
             


        if ( irec == 1 )
          for ilin = 1 : length( hand )
            colo( ilin, : ) = get( hand( ilin ), 'color' );
          end; clear ilin;
          colo = [ 0.5, 0.5, 0.5; colo(1 : end - 1, :) ];
        end
        for ilin = 1 : length( hand )
          set( hand( ilin ), 'color', colo( ilin, : ), 'linewidth', 1.5 );
        end; clear ilin;

    

    if ( irec == 1 )
      lege = [ 'Topography         ' ];
      lege = [ lege; 'Interface ' num2str( 1     ),'        '];
      lege = [ lege; 'Interface ', num2str( 2 ), ': x=west' ];
      lege = [ lege; 'Interface ', num2str( 2 ), ': x=mid ' ];
      lege = [ lege; 'Interface ', num2str( 2 ), ': x=east' ];
    end


    days=sprintf('%.2f',taxi(irec));
    xlabel('Distance from ridge (km)');
    ylabel('Depth (m)');
    title(strcat('Time: ' ,num2str(days),' days'))
    legend(hand, lege, 'location', 'south');
    axis ij;
    axis([yy(iy_0, 1 + npts ) / 1.e3, yy(iy_0, end - npts ) / 1.e3, - 1., hmax]);
    axis square;
    
    figure(2)
    umin=-.4; umax=.4;
    subplot(1,3,1)
    contourf(xx(:,1),yy(1,:),n(:,:,2)'-topl(2)*hmax)
    shading interp, colorbar
    set(gca,'fontsize',14)
    xlabel('x (km)')
    ylabel('y (km)')
    title( strcat('T=' ,num2str(days), 'days, interface depth z(m) '))
    subplot(1,3,2)
    contourf(xx(:,1),yy(1,:),v(:,:,2)')
    caxis([umin umax])
    shading interp, colorbar
    set(gca,'fontsize',14)
    xlabel('x (km)')
    ylabel('y (km)')
    title('bot along-chan u_2(m/s)')
    subplot(1,3,3)
    contourf(xx(:,1),yy(1,:),v(:,:,1)')   
    caxis([umin umax])
    shading interp, colorbar
    set(gca,'fontsize',14)
    xlabel('x (km)')
    ylabel('y (km)')
    title('top along-chan u_1(m/s)')
    set(gcf,'units','points','position',[10,10,1400,700])
    
    pause( 1.e-6 );
    if ( movi )
      %myft = '/usr/share/fonts/type1/gsfonts/n019003l.pfb';
      %print('-dpng', [outd 'frame' num2str(irec, '%04u') '.png'], '-S1024,768', ['-F' myft]);
      figure(1)
      currFrame=getframe(gcf);
      writeVideo(vidObj,currFrame);
      figure(2)
      currFrame2=getframe(gcf);
      writeVideo(vidObj2,currFrame2);
      
    end
  end; clear irec currFrame currFrame2; 
  
  if (movi)
  close(vidObj); close(vidObj2);
  end 
  
  
