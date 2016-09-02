% -- External Function: fiel = get_field( vnam, irec, outd );
%
% Returns a model variable in gridded format (i,j,nlay).
%
% Function arguments:
% `vnam' (input ): 4 characters long string identifying the variable desired:
%                    'eta_' Interface/surface deviation from state of rest (meters),
%                    'u___' Zonal      velocity component  (m/s    ),
%                    'v___' Meridional velocity component  (m/s    ),
%                    'pvor' Potential  vorticity           (1/(m s)),
%                    'mont' Montgomery potential / gravity (m      ),
%                    'v_cc' Epipycnal  eddy viscosity      (m2/s   ).
% `irec' (input ): integer, Must be >=0 and <= the number of entries in `time.txt'.
%                           1 = Initial condition,
%                           0 = Last record available in output files.
% `outd' (input ): string, the path to the output directory.
% `fiel' (output): float array, the variable in gridded format (i,j,nlay).
function [fiel] = get_field(vnam, irec, outdir)

  vnam = strtrim(vnam);
  nvar = ['eta_'; ... %  1
          'u___'; ... %  2
          'v___'; ... %  3
          'pvor'; ... %  4
          'smcc'; ... %  5
          'bern'; ... %  6
          'bave'; ... %  7
          'smll'; ... %  8
          'layt'; ... %  9
          'layb'; ... % 10
          'rvor'; ... % 11
          'mont'; ... % 12
          'v_cc'];    % 13
  ivar = 0;
  for i = 1 : size(nvar, 1)
    if ( strcmp(vnam, nvar(i,:)) )
      ivar = i;
    end
  end; clear i;
  if ( ivar == 0 )
    error(['Unknown variable ' vnam '. See help.']);
  end

% Read metadata.

  [fid, msg] = fopen( [deblank(outdir) 'param_basin.txt'], 'r' );
  if ( fid <= 0 )
    error(msg);
  end
  not_eof = 1;
  while ( not_eof )
    line = strtrim( fgetl(fid) );
    if ( strcmp( line(1 : 4), 'desc' ) )
      not_eof = 0;
    end
    eval( line );
  end; clear not_eof;
  fclose(fid); clear fid msg;

  taxi = load([deblank(outdir) 'time.txt']);
  nrec = length(taxi); clear taxi;
  if ( irec == 0 )
    irec = nrec;
  end

  [fid, msg] = fopen([deblank(outdir) 'grid.bin' ], 'r', 'ieee-le');
% [fid, msg] = fopen([deblank(outdir) 'pos_c.bin'], 'r', 'ieee-le');
  [val, cnt] = fread(fid, Inf, 'int32', 0, 'ieee-le');
  pos_b = val(1 : round(length(val) / 5.));
  fclose(fid); clear fid msg cnt;

  if     ( ivar == 2 )
    mask = val(1 + 2 * round(length(val) / 5.) : 3 * round(length(val) / 5.));
%   [fid,  msg] = fopen([deblank(outdir) 'mask_u.bin'], 'r', 'ieee-le');
%   [mask, cnt] = fread(fid, Inf, 'int32', 0, 'ieee-le');
%   fclose(fid); clear fid msg cnt;
  elseif ( ivar == 3 )
    mask = val(1 + 3 * round(length(val) / 5.) : 4 * round(length(val) / 5.));
%   [fid,  msg] = fopen([deblank(outdir) 'mask_v.bin'], 'r', 'ieee-le');
%   [mask, cnt] = fread(fid, Inf, 'int32', 0, 'ieee-le');
%   fclose(fid); clear fid msg cnt;
  elseif ( (ivar == 4 || ivar == 8 || ivar == 11) )
%       && exist([deblank(outdir) 'mask_p.bin'], 'file') )
    mask = val(1 + 4 * round(length(val) / 5.) : 5 * round(length(val) / 5.));
%   [fid,  msg] = fopen([deblank(outdir) 'mask_p.bin'], 'r', 'ieee-le');
%   [mask, cnt] = fread(fid, Inf, 'int32', 0, 'ieee-le');
%   fclose(fid); clear fid msg cnt;
  else
    mask = val(1 + round(length(val) / 5.) : 2 * round(length(val) / 5.));
%   [fid,  msg] = fopen([deblank(outdir) 'mask_n.bin'], 'r', 'ieee-le');
%   [mask, cnt] = fread(fid, Inf, 'int32', 0, 'ieee-le');
%   fclose(fid); clear fid msg cnt;
  end
  clear val;

  mask( find( mask == 0 ) ) = nan;

  fiel = nan(lm + 2, mm + 2, nlay);

  [fid, msg] = fopen([deblank(outdir) nvar(ivar,:) '.bin'], 'r', 'ieee-le');
  if ( fid <= 0 )
%   error(msg);
    return;
  end

  if ( irec > 1 )
    fseek(fid, 4 * (irec(1) - 1) * length(pos_b(:)) * nlay, 'bof');
  end
  [vect, cnt] = fread(fid, [length(pos_b(:)), nlay], 'real*4', 0, 'ieee-le');
  for ilay = nlay : -1 : 1
    vari           = nan( lm + 2, mm + 2 );
    vect(:,ilay)   = squeeze(vect(:,ilay)) .* mask(:);
    vari(pos_b)    = vect(:,ilay);
    fiel(:,:,ilay) = vari(:,:);
  end; clear ilay;
  fiel = squeeze(fiel);
  fclose(fid);
