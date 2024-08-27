function [ ldscores, ldscores_adjusted ] = ldscore_calc( X, do_loader )
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
% [ chi2, X ] = gengenmodel( 1000, 1000, 0.2, 0.5, 'ar1', 1);
% [ ldscores, ldscores_adjusted ] = ldscore_calc( X, 0 );
%--------------------------------------------------------------------------
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'do_loader', 'var' )
   % Default value
   do_loader = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
s_X = size(X);
n = s_X(1);
m = s_X(2);

XXT = X*X';

ldscores = zeros(1,m);
for j = 1:m
    if do_loader
        loader(j,m, 'Computing ldscores:')
    end
    ldscores(j) = (1/n^2)*X(:,j)'*XXT*X(:,j);
end

ldscores_adjusted = ldscores - (m-ldscores)/(n-2);

end

