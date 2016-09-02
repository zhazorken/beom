% -- External Function: [taxi, h_0, f0, dl, rhon, desc] = get_metadata( outd );
%
% Returns metadata (information about the calculation and its output files).
%
% Function arguments:
% `outd' (input ): string, the path to the output directory.
% `taxi' (output): float vector, the time-axis corresponding to the records
%                    inside the output files (days).
% `h_0'  (output): float array of size (lm + 2, mm + 2, nlay),
%                    containing the thickness of the layers at rest (meters).
% `f0'   (output): float, Coriolis parameter entered in `shared_mod.f95' (rad/s ).
% `dl'   (output): float, mesh size          entered in `shared_mod.f95' (meters).
% `rhon' (output): float vector, the density entered in `shared_mod.f95' (kg/m3 ).
% `desc' (output): string, the description   entered in `shared_mod.f95'.
function [taxi, h_0, f0, dl, rhon, desc] = get_metadata(outdir)

% Read metadata.

  [fid, msg] = fopen( [deblank(outdir) 'param_basin.txt'], 'r' );
  if ( fid <= 0 )
    msg
    return
  end
  not_eof = 1;
  while ( not_eof )
    line = fgetl(fid);
    if ( length(line) <= 1 )
      not_eof = 0;
    else
      line = strtrim( line );
      eval( line );
    end
  end; clear not_eof;
  fclose(fid); clear fid msg;

  taxi = load([deblank(outdir) 'time.txt']);

  [fid,   msg] = fopen([deblank(outdir) 'grid.bin' ], 'r', 'ieee-le');
% [fid,   msg] = fopen([deblank(outdir) 'pos_c.bin'], 'r', 'ieee-le');
  [pos_b, cnt] = fread(fid, Inf, 'int32', 0, 'ieee-le');
  pos_b = pos_b(1 : round(length(pos_b) / 5.));
  fclose(fid); clear fid msg cnt;

  [fid,   msg] = fopen([deblank(outdir) 'h_0.bin'], 'r', 'ieee-le');
  [vh_0,  cnt] = fread(fid, [length(pos_b(:)), nlay], 'real*4', 0, 'ieee-le');
  fclose(fid); clear fid msg cnt;

  for ilay = 1 : nlay;
    hi_0 = zeros(lm + 2, mm + 2);
    hi_0(pos_b)   = vh_0(:,ilay);
    h_0(:,:,ilay) = hi_0(:,:);
    clear hi_0;
  end; clear ilay;
end
