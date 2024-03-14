function [ ldsc_full, ldsc_intercept1, ldsc_conditional, ldratio, ldratiomn] = ...
                          h2ests( n, m, ldscores_adjusted, chi2, ldscores )
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
%
%--------------------------------------------------------------------------
% OUTPUT
%   ldsc_full           - LD score regression.
%   ldsc_intercept1     - LD score regression with intercept set to 1.
%   ldsc_conditional    - Conditional LD score regression.
%   gwash               - GWAS heritability estimate.
%   gwashmn             - Conditional GWAS heritability estimate.
%
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
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Unconstrained ld score regression
design = [(ldscores_adjusted)'*(n/m), ones(m,1)];
% 
ldsc_full = (design'*design)^(-1)*design'*chi2;
% ldsc = ldsc(1);

% Ld score regression with the intercept set to 1
design = [(ldscores_adjusted)'*(n/m)];

ldsc_intercept1 = (design'*design)^(-1)*design'*(chi2-1);

% Conditional
if exist('ldscores', 'var')
    design = [(ldscores*(n/m)-1)'];
    ldsc_conditional = (design'*design)^(-1)*design'*(chi2-1);
else
    ldsc_conditional = NaN;
end

% GWASH (up to O(1/n)) see GWASH supplementary!
ldratio = (mean(chi2) - 1)/mean((ldscores_adjusted*(n/m))');

% condtional GWASH
% (mean(chi2) - 1)/mean((ldscores*(n/m)-1)')

ldratiomn = (mean(chi2) - 1)*(m/n);

end

