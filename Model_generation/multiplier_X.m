function mX = multiplier_X( X, k )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% X = Xgen( 500, 1000, 0.7, 'ar1' );
% mX = multiplier_X( X, 500 );
% [ ~, ldscores_adjusted_X ] = ldscore_calc( X );
% [ ~, ldscores_adjusted_mX ] = ldscore_calc( mX );
% plot(ldscores_adjusted_X, ldscores_adjusted_mX, '*')
%--------------------------------------------------------------------------
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
n = size(X,1);
G = randn([k,n]);
mX = G*X/sqrt(n);

end

