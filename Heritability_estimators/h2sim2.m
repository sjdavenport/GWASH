function [ OG_h2, ref_h2, OG_h2_summary, ref_h2_summary ] = ...
            h2sim2( n, m, nref, h2, rho, method, nsims, do_standardize, distbn )
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
% [ OG_h2, ref_h2, OG_h2_summary, ref_h2_summary ] = h2sim2( 1000, 1000, 1000, 0.2, 0.5, 'ar1', 100, 1 )
% [ OG_h2, ref_h2, OG_h2_summary, ref_h2_summary ] = h2sim2( 1000, 1000, 1000, 0.2, [0.2,0.9], 'ar1mix', 100, 1 )
%--------------------------------------------------------------------------
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'do_standardize', 'var' )
   % Default value
   do_standardize = 0;
end

if ~exist( 'distbn', 'var' )
   % Default value
   distbn = 'norm';
end
%% Initialize the results
OG_h2.ldsc_free = zeros(1,nsims);
OG_h2.ldsc_free_intercept = zeros(1,nsims);
OG_h2.ldsc_fixed_intercept = zeros(1,nsims);
OG_h2.ldsc_conditional = zeros(1,nsims);
OG_h2.gwash = zeros(1,nsims);

ref_h2.ldsc_free = zeros(1,nsims);
ref_h2.ldsc_free_intercept = zeros(1,nsims);
ref_h2.ldsc_fixed_intercept = zeros(1,nsims);
ref_h2.ldsc_conditional = zeros(1,nsims);
ref_h2.gwash = zeros(1,nsims);

%%  Main Function Loop
%-------------------------------------------------------------------------
for I = 1:nsims
    I
    [ chi2, X ] = gengenmodel( n, m, h2, rho, method, do_standardize, distbn);

    [ ldscores, ldscores_adjusted ] = ldscore_calc( X, 0 );
    [ ldsc_free, ldsc_fixed_intercept, ldsc_conditional, gwash] = ...
        h2ests2( n, m, ldscores_adjusted, chi2, ldscores );

    OG_h2.ldsc_free(I) = ldsc_free(1);
    OG_h2.ldsc_free_intercept(I) = ldsc_free(2);
    OG_h2.ldsc_fixed_intercept(I)= ldsc_fixed_intercept.orig;
    OG_h2.ldsc_fixed_interceptW(I)= ldsc_fixed_intercept.W;
    OG_h2.ldsc_fixed_interceptW1(I)= ldsc_fixed_intercept.W1;
    OG_h2.ldsc_conditional(I) = ldsc_conditional;
    OG_h2.gwash(I) = gwash.orig;
    OG_h2.gwashW(I) = gwash.W;
    OG_h2.gwashW1(I) = gwash.W1;

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
    [ ldscores_othersample, ldscores_othersample_adjusted ] = ldscore_calc( X_othersample, 0 );
    % histogram(ldscores_othersample)

    % Using original chi2
    
    [ ldsc_free, ldsc_fixed_intercept, ldsc_conditional, gwash] = ...
        h2ests2( n, m, ldscores_othersample_adjusted, chi2, ldscores_othersample );

    ref_h2.ldsc_free(I) = ldsc_free(1);
    ref_h2.ldsc_free_intercept(I) = ldsc_free(2);
    ref_h2.ldsc_fixed_intercept(I)= ldsc_fixed_intercept.orig;
    ref_h2.ldsc_fixed_interceptW(I)= ldsc_fixed_intercept.W;
    ref_h2.ldsc_fixed_interceptW1(I)= ldsc_fixed_intercept.W1;
    ref_h2.ldsc_conditional(I) = ldsc_conditional;
    ref_h2.gwash(I) = gwash.orig;
    ref_h2.gwashW(I) = gwash.W;
    ref_h2.gwashW1(I) = gwash.W1;
end

reg_types = {'ldsc_free', 'ldsc_free_intercept', 'ldsc_fixed_intercept', 'ldsc_fixed_interceptW', 'ldsc_fixed_interceptW1', 'ldsc_conditional', 'gwash', 'gwashW', 'gwashW1'};

for K = 1:length(reg_types)
    type2use = reg_types{K};
    OG_h2_summary.(type2use).mean = mean(OG_h2.(type2use));
    OG_h2_summary.(type2use).std = std(OG_h2.(type2use));
    OG_h2_summary.(type2use).ss = sum(OG_h2.(type2use).^2);

    ref_h2_summary.(type2use).mean = mean(ref_h2.(type2use));
    ref_h2_summary.(type2use).std = std(ref_h2.(type2use));
    ref_h2_summary.(type2use).ss = sum(ref_h2.(type2use).^2);
end

end

