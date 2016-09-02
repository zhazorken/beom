% -- External Function: ndeg = get_nbr_deg_freedom( h_bo );
%
% Returns the size of the grid after discarding the grid points
%   located on land.
% Each model variable has size (0 : ndeg, nlay) in the code.
%
% Function arguments:
% `h_bo' (input ): float array of size (0 : lm + 1, 0 : mm + 1) (meters).
% `ndeg' (output): integer, must be entered in file `shared_mod.f95'.
function ndeg = get_nbr_deg_freedom( h_bo )

% h_bo( 0 : lm + 1,  0 : mm + 1)
% hext(-1 : lm + 2, -1 : mm + 2)

% Local depth must be >hdry for the cell to be considered wet.
% CAREFUL, this value must be consistent with the one in file shared_mod.f95.

  hdry = 1.e-3; % Threshold (meters).

% To be on the safe side with respect to roundoff error,
%   everything shallower than 2*hdry is declared `dry'.

  h_bo( find( h_bo < (2. * hdry) ) ) = 0.;

% Margins must be set to zero
% (even if there are open boundaries).

  h_bo(1,  :) = 0.;
  h_bo(end,:) = 0.;
  h_bo(:,  1) = 0.;
  h_bo(:,end) = 0.;

  lm   = size( h_bo, 1 ) - 2;
  mm   = size( h_bo, 2 ) - 2;

  mask = zeros( lm + 4, mm + 4 );
  hext = zeros( lm + 4, mm + 4 );
  hext(2 : end - 1, 2 : end - 1) = h_bo;
  mask( find( hext > hdry ) ) = 1;
  clear h_bo hext;

  ndeg = 0;

  neig = mask(2 : end - 1, 2 : end - 1) ... % rho point.
       + mask(1 : end - 2, 2 : end - 1) ... % u   point.
       + mask(2 : end - 1, 1 : end - 2) ... % v   point.
       + mask(1 : end - 2, 1 : end - 2);    % psi point.
  ndeg = length( find( neig > 0 ) );

% for j = 2 : mm + 3
%   for i = 2 : lm + 3
%     if (             hext(i,         j        )      > hdry  || ... % rho point.
%          any(        hext(i - 1 : i, j        )      > hdry) || ... % u   point.
%          any(        hext(i,         j - 1 : j)      > hdry) || ... % v   point.
%          any(reshape(hext(i - 1 : i, j - 1 : j),4,1) > hdry) )      % psi point.
%       ndeg = ndeg + 1;
%     end
%   end
% end

end
