% -- External Function: print_params( ...params... );
%
% Prints on the screen the parameters to be copy-pasted in file shared_mod.f95.
%   The function automatically formats the parameters in a smart way;
%     e.g. "dl = 100000." becomes "dl = 100.e3".
function print_params( lm,   mm,   nlay, ndeg, dl,   cext, f0,   rhon, topl, ...
                       dt_s, dt_o, dt_r, dt3d, bvis, dvis, bdrg, hmin, hsbl, ...
                       hbbl, g_fb, uadv, qdrg, ocrp, rsta, xper, yper, diag, ...
                       tauw, idir, odir, desc )
  disp(' ');
  disp('Please copy-paste these parameters inside shared_mod.f95:');

  disp(  ['lm         = ' sprintf( '%li', lm   )]);
  disp(  ['mm         = ' sprintf( '%li', mm   )]);
  disp(  ['nlay       = ' sprintf( '%i',  nlay )]);
  disp(  ['ndeg       = ' sprintf( '%li', ndeg )]);

  if ( dl < 1.e3 )
    disp(['dl         = ' sprintf( '%#0.0f', dl        )     ]);
  elseif ( (dl - floor( dl / 1.e3 ) * 1.e3) > 1. )
    disp(['dl         = ' sprintf( '%g',     dl / 1.e3 ) 'e3']);
  else % If an exact multiple of 1 kilometer,
    disp(['dl         = ' sprintf( '%#0.0f', dl / 1.e3 ) 'e3']);
  end

  disp(  ['cext       = ' sprintf( '%0.1f',  cext      )   ]);

  if     ( abs( f0 ) < 1.e-4    )
    disp(['f0         = ' sprintf( '%e',     f0         )    ]);
  elseif ( abs( f0 ) > 1.001e-4 )
    disp(['f0         = ' sprintf( '%0.3f',  f0 / 1.e-4 ) 'e-4']);
  else
    disp(['f0         = ' sprintf( '%#0.0f', f0 / 1.e-4 ) 'e-4']);
  end

  if ( all( (rhon(:) * 1.e3 - floor( rhon(:) ) * 1.e3 ) < 1. ) )
    strg = sprintf( '%#0.0f,', rhon );
  else
    strg = sprintf( '%0.3f,',  rhon );
  end
  disp(  ['rhon(nlay) = (/' strg( 1 : end - 1 ) '/)']);

  strg = sprintf( '%f,', topl );
  disp(  ['topl(nlay) = (/' strg( 1 : end - 1 ) '/)']);

  disp(  ['dt_s       = ' sprintf( '%#f',    dt_s )]);
  disp(  ['dt_o       = ' sprintf( '%#f',    dt_o )]);
  if ( dt_r > 0. )
    disp(['dt_r       = ' sprintf( '%#f',    dt_r )]);
  else
    disp(['dt_r       = ' sprintf( '%#0.0f', dt_r )]);
  end
  if ( dt3d > 0. )
    disp(['dt3d       = ' sprintf( '%#f',    dt3d )]);
  else
    disp(['dt3d       = ' sprintf( '%#0.0f', dt3d )]);
  end
  if ( bvis > 0. )
    disp(['bvis       = ' sprintf( '%#f',    bvis )]);
  else
    disp(['bvis       = ' sprintf( '%#0.0f', bvis )]);
  end
  if ( dvis > 0. )
    disp(['dvis       = ' sprintf( '%#0.3f', dvis )]);
  else
    disp(['dvis       = ' sprintf( '%#0.0f', dvis )]);
  end

  if ( bdrg > 0. )
    disp(['bdrg       = ' sprintf( '%e',     bdrg )]);
  else
    disp(['bdrg       = ' sprintf( '%#0.0f', bdrg )]);
  end

  disp(  ['hmin       = ' sprintf( '%#f',    hmin )]);
  disp(  ['hsbl       = ' sprintf( '%#0.0f', hsbl )]);
  disp(  ['hbbl       = ' sprintf( '%#0.0f', hbbl )]);
  disp(  ['g_fb       = ' sprintf( '%#0.0f', g_fb )]);
  disp(  ['uadv       = ' sprintf( '%#0.0f', uadv )]);
  disp(  ['qdrg       = ' sprintf( '%#0.0f', qdrg )]);
  disp(  ['ocrp       = ' sprintf( '%#0.0f', ocrp )]);
  disp(  ['rsta       = ' sprintf( '%#0.0f', rsta )]);
  disp(  ['xper       = ' sprintf( '%#0.0f', xper )]);
  disp(  ['yper       = ' sprintf( '%#0.0f', yper )]);
  disp(  ['diag       = ' sprintf( '%#0.0f', diag )]);

  strg = sprintf( '%#0.2f,', tauw );
  disp(  ['tauw       = (' strg( 1 : end - 1 ) ')']);

  disp(  ['idir       = ''' idir '''']);
  disp(  ['odir       = ''' odir '''']);
  disp(  ['desc       = ''' desc '''']);
  disp(' ');
end
