function [ chi2, X, phi ] = gengenmodel( n, m, h2, rho, method, do_standardize, distbn )
% GENGENMODEL Generate genetic model for LD score regression.
%
%   [chi2, X, phi] = GENGENMODEL(n, m, h2, rho, method) generates genetic 
%   model for LD score regression based on the given parameters and returns 
%   the chi-squared statistics (chi2), the genetic predictors (X), and the 
%   phenotype (phi).
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory:
%   n       - Number of individuals.
%   m       - Number of SNPs.
%   h2      - Heritability value.
%   rho     - Genetic correlation between SNPs.
%   method  - Method for generating genetic predictors.
%--------------------------------------------------------------------------
% OUTPUT
%   chi2    - Chi-squared statistics.
%   X       - Genetic predictors.
%   phi     - Phenotype.
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'rho', 'var' )
   % Default value
   rho = 0;
end

if ~exist( 'distbn', 'var' )
   % Default value
   distbn = 'norm';
end
%%  Main Function Loop
%--------------------------------------------------------------------------
X = Xgen( n, m, rho, method );
% X = X - mean(X);
% X = X./std(X,0,1);

if strcmp(distbn, 'norm')
% iid beta
    beta = ((h2/m)^(1/2))*randn(m,1); % similar to gwash sims, once you normalize
elseif strcmp(distbn(1), 't')
    df = str2double(distbn(2:end));
    beta = ((h2/m)^(1/2))*trnd(df, m,1)*sqrt((df-2)/df);
end

% Actually seems to improve if you makes the betas smooth!
% smoothing_parameter = 30;
% smooth_prevarscaled_betas = convfield(wfield(m,1), smoothing_parameter);
% global PIloc
% load([PIloc,'Variance/storevars'], 'allvars')
% beta = (smooth_prevarscaled_betas.field/sqrt(allvars(smoothing_parameter)))*((h2/m)^(1/2));

% beta = randn(m,1);
e = ((1-h2)^(1/2))*randn(n,1);

phi = X*beta + e;
if do_standardize
    phi = phi - mean(phi);
    phi = phi./std(phi,0,1);
    X = X - mean(X);
    X = X./std(X,0,1);
end

betahat = zeros(m,1);
for j = 1:m
    betahat(j) = X(:,j)'*phi/n; %/norm(X(:,j))^2 consider normalizing the Xs
end

chi2 = n*betahat.^2;

% Calculate the ld scores
% [ ldscores, ldscores_adjusted ] = ldscores_calc( X, do_loader );

end

