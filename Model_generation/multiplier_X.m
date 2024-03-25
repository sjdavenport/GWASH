function mX = multiplier_X( X, k, use_matrix )
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
%
% X = Xgen( 500, 10000, 0.7, 'ar1' );
% tic; mX = multiplier_X( X, 500, 0 ); toc
% tic; mX = multiplier_X( X, 500, 1 ); toc
%--------------------------------------------------------------------------
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'use_matrix', 'var' )
    % Default value
    use_matrix = 1;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
n = size(X,1);

if use_matrix == 1
    G = randn([k,n]);
    mX = G*X/sqrt(n);
else
    % This way is the same but slower!
    m = size(X,2);
    mX = zeros([k,m]);
    for I = 1:k
        g = randn([1,n]);
        mX(I,:) = g*X/sqrt(n);
    end
end

end

