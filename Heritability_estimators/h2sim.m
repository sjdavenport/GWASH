function [ OG_h2, ref_h2, OG_h2_summary, ref_h2_summary ] = ...
            h2sim( n, m, nref, h2, rho, method, nsims, do_standardize )
% h2sim( n, m, nref, h2, rho, method, do_standardize )
%--------------------------------------------------------------------------
% H2SIM Estimate heritability using LD score regression in a simulation setup.
%
%   [OG_h2, ref_h2] = H2SIM(n, m, nref, h2, rho, method, do_standardize) 
%   runs heritability estimation simulations. It generates LD scores,
%   estimates heritability for the original sample (OG_h2) and a reference 
%   sample (ref_h2). It considers several different methods and returns 
%   the results.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory:
%   n                   - Number of individuals in the original sample.
%   m                   - Number of SNPs.
%   nref                - Number of individuals in the reference sample.
%   h2                  - Heritability value.
%   rho                 - Genetic correlation between the original and reference samples.
%   method              - Method for generating LD scores.
%   do_standardize      - Flag indicating whether to standardize the data. 
%                         (default: 0)
%--------------------------------------------------------------------------
% OUTPUT
%   OG_h2               - Structure containing heritability estimates for 
%                         the original sample.
%   ref_h2              - Structure containing heritability estimates for 
%                         the reference sample.
%--------------------------------------------------------------------------
% EXAMPLES
% [ OG_h2, ref_h2, OG_h2_summary, ref_h2_summary ] = h2sim( 1000, 1000, 1000, 0.2, 0.5, 'ar1', 100, 1 )
%--------------------------------------------------------------------------
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'do_standardize', 'var' )
   % Default value
   do_standardize = 0;
end

%% Initialize the results
OG_h2.ldsc_free = zeros(1,nsims);
OG_h2.ldsc_free_intercept = zeros(1,nsims);
OG_h2.ldsc_fixed_intercept = zeros(1,nsims);
OG_h2.ldsc_conditional = zeros(1,nsims);
OG_h2.ld_ratio = zeros(1,nsims);

ref_h2.ldsc_free = zeros(1,nsims);
ref_h2.ldsc_free_intercept = zeros(1,nsims);
ref_h2.ldsc_fixed_intercept = zeros(1,nsims);
ref_h2.ldsc_conditional = zeros(1,nsims);
ref_h2.ld_ratio = zeros(1,nsims);

%%  Main Function Loop
%-------------------------------------------------------------------------
for I = 1:nsims
    [ chi2, X ] = gengenmodel( n, m, h2, rho, method);

    if do_standardize
        X = X - mean(X);
        X = X./std(X,0,1);
    end

    [ ldscores, ldscores_adjusted ] = ldscore_calc( X, 1 );
    [ ldsc_free, ldsc_fixed_intercept, ldsc_conditional, ld_ratio] = ...
        h2ests( n, m, ldscores_adjusted, chi2, ldscores );

    OG_h2.ldsc_free(I) = ldsc_free(1);
    OG_h2.ldsc_free_intercept(I) = ldsc_free(2);
    OG_h2.ldsc_fixed_intercept(I)= ldsc_fixed_intercept;
    OG_h2.ldsc_conditional(I) = ldsc_conditional;
    OG_h2.ld_ratio(I) = ld_ratio;

    % [ chi2_othersample, X_othersample ] = gengenmodel( nref, m, h2, rho, method);
    % X_othersample = Xgen( nref, m, rho, method );
    % if do_standardize
    %     X_othersample = X_othersample - mean(X_othersample);
    %     X_othersample = X_othersample./std(X_othersample,0,1);
    % end
    % [ chi2_othersample, X_othersample ] = gengenmodel( nref, m, h2, rho, method);
    X_othersample = Xgen( nref, m, rho, method );
    if do_standardize
        X_othersample = X_othersample - mean(X_othersample);
        X_othersample = X_othersample./std(X_othersample,0,1);
    end
    [ ldscores_othersample, ldscores_othersample_adjusted ] = ldscore_calc( X_othersample, 1 );
    % histogram(ldscores_othersample)

    % Using original chi2
    [ ldsc_free, ldsc_fixed_intercept, ldsc_conditional, ld_ratio] = ...
        h2ests( n, m, ldscores_othersample_adjusted, chi2, ldscores_othersample );

    ref_h2.ldsc_free(I) = ldsc_free(1);
    ref_h2.ldsc_free_intercept(I) = ldsc_free(2);
    ref_h2.ldsc_fixed_intercept(I) = ldsc_fixed_intercept;
    ref_h2.ldsc_conditional(I) = ldsc_conditional;
    ref_h2.ld_ratio(I) = ld_ratio;
end

OG_h2_summary.ldsc_free.mean = mean(OG_h2.ldsc_free);
OG_h2_summary.ldsc_free.ss = sum(OG_h2.ldsc_free.^2);
OG_h2_summary.ldsc_free_intercept.mean = mean(OG_h2.ldsc_free_intercept);
OG_h2_summary.ldsc_free_intercept.ss = sum(OG_h2.ldsc_free_intercept.^2);
OG_h2_summary.ldsc_fixed_intercept.mean = mean(OG_h2.ldsc_fixed_intercept);
OG_h2_summary.ldsc_fixed_intercept.ss = sum(OG_h2.ldsc_fixed_intercept.^2);
OG_h2_summary.ldsc_conditional.mean = mean(OG_h2.ldsc_conditional);
OG_h2_summary.ldsc_conditional.ss = sum(OG_h2.ldsc_conditional.^2);
OG_h2_summary.ld_ratio.mean = mean(OG_h2.ld_ratio);
OG_h2_summary.ld_ratio.ss = sum(OG_h2.ld_ratio.^2);

ref_h2_summary.ldsc_free.mean = mean(ref_h2.ldsc_free);
ref_h2_summary.ldsc_free.ss = sum(ref_h2.ldsc_free.^2);
ref_h2_summary.ldsc_free_intercept.mean = mean(ref_h2.ldsc_free_intercept);
ref_h2_summary.ldsc_free_intercept.ss = sum(ref_h2.ldsc_free_intercept.^2);
ref_h2_summary.ldsc_fixed_intercept.mean = mean(ref_h2.ldsc_fixed_intercept);
ref_h2_summary.ldsc_fixed_intercept.ss = sum(ref_h2.ldsc_fixed_intercept.^2);
ref_h2_summary.ldsc_conditional.mean = mean(ref_h2.ldsc_conditional);
ref_h2_summary.ldsc_conditional.ss = sum(ref_h2.ldsc_conditional.^2);
ref_h2_summary.ld_ratio.mean = mean(ref_h2.ld_ratio);
ref_h2_summary.ld_ratio.ss = sum(ref_h2.ld_ratio.^2);

end

