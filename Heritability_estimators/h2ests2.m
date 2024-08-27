function [ ldsc_full, ldsc_intercept1, ldsc_conditional, gwash] = ...
                          h2ests2( n, m, ldscores_adjusted, chi2, ldscores )
% H2ESTS Estimate heritability using LD score regression and related
% estimators.
%
%   [ldsc_full, ldsc_intercept1, ldsc_conditional, gwash, gwashmn] = 
%                       H2ESTS(n, m, ldscores_adjusted, chi2, ldscores) 
%   estimates heritability using LD score regression with different settings 
%   and returns the results.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory:
%   n                   - Number of individuals.
%   m                   - Number of SNPs.
%   ldscores_adjusted   - Adjusted LD scores.
%   chi2                - Chi-squared statistics.
%   ldscores            - LD scores.
%--------------------------------------------------------------------------
% OUTPUT
%   ldsc_full           - LD score regression.
%   ldsc_intercept1     - LD score regression with intercept set to 1.
%   ldsc_conditional    - Conditional LD score regression.
%   gwash               - GWAS heritability estimate.
%   gwashmn             - Conditional GWAS heritability estimate.
%--------------------------------------------------------------------------
% EXAMPLES
% n = 1000; m = 1000;
% [ chi2, X ] = gengenmodel( n, m, 0.2, [0.2,0.9], 'ar1mix', 1);
% [ ldscores, ldscores_adjusted ] = ldscore_calc( X, 0 );
% [ ldsc_full, ldsc_intercept1, ldsc_conditional, gwash] = h2ests2( n, m, ldscores_adjusted, chi2, ldscores )
% --------------------------------------------------------------------------
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
% Unconstrained ld score regression
% design = [(ldscores_adjusted)'*(n/m)];
% design = design - mean(design);
% design = [design, ones(m,1)];
design = [(ldscores_adjusted)'*(n/m), ones(m,1)];
% design = [(ldscores)'*(n/m), ones(m,1)];
% 
ldsc_full = (design'*design)^(-1)*design'*chi2;
% ldsc = ldsc(1);

% Ld score regression with the intercept set to 1
design = [(ldscores_adjusted)'*(n/m)];
ldsc_intercept1.orig = (design'*design)^(-1)*design'*(chi2-1);

% Conditional
if exist('ldscores', 'var')
    design = [(ldscores*(n/m)-1)'];
    ldsc_conditional = (design'*design)^(-1)*design'*(chi2-1);
else
    ldsc_conditional = NaN;
end

% GWASH (up to O(1/n)) see GWASH supplementary!
gwash.orig = (mean(chi2) - 1)/mean((ldscores_adjusted*(n/m))');

% Weighting
ld_adj_max1 = max(ldscores_adjusted,1);

% Weighted LDSC

% Calculate weights
W_ldsc_intercept1 = 1./(1+ldsc_intercept1.orig*n/m*ld_adj_max1).^2;
W1_ldsc_intercept1 = (1./(1+ldsc_intercept1.orig*n/m*ld_adj_max1).^2)./(ld_adj_max1);

% Implement weighted ldsc estimates
design_W = [((W_ldsc_intercept1.^(1/2)).*ldscores_adjusted)'*(n/m)];
ldsc_intercept1.W = (design_W'*design_W)^(-1)*design_W'*((chi2-1).*(W_ldsc_intercept1.^(1/2))');

design_W1 = [((W1_ldsc_intercept1.^(1/2)).*ldscores_adjusted)'*(n/m)];
ldsc_intercept1.W1 = (design_W1'*design_W1)^(-1)*design_W1'*((chi2-1).*(W1_ldsc_intercept1.^(1/2))');

% W_ldsc_full = 1/(1+ldsc_full(1)*n/m*ell_adj)^2;
% W1_ldsc_full = 1/(1+ldsc_full(1)*n/m*ell_adj)^2*1/(ell_adj);

% Weighted gwash

% Calculate weights
W_gwash = 1./(1+gwash.orig*n/m*ld_adj_max1).^2;
W1_gwash = (1./(1+gwash.orig*n/m*ld_adj_max1).^2)./(ld_adj_max1);

% Implement weighted gwash estimates
gwash.W = (mean((chi2-1).*W_gwash'))/mean(W_gwash'.*(ldscores_adjusted*(n/m))');
gwash.W1 = (mean((chi2-1).*W1_gwash'))/mean(W1_gwash'.*(ldscores_adjusted*(n/m))');

end

% condtional GWASH
% (mean(chi2) - 1)/mean((ldscores*(n/m)-1)')

% ldratiomn = (mean(chi2) - 1)*(m/n);
