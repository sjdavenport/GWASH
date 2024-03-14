function X = Xgen( n, m, rho, method )
%XGEN Generate a matrix of size n-by-m according to specified method.
%   X = XGEN(n, m, rho, method) generates a matrix X of size n-by-m according
%   to the specified method. 
%--------------------------------------------------------------------------
%   ARGUMENTS:
%   - n: Number of rows of the matrix.
%   - m: Number of columns of the matrix.
%   - rho: Autoregressive parameter for 'ar1' method. For 'equi' methods, 
%          it specifies the correlation coefficient and for Gaussian the
%          FWHM of the smoothing kernel.
%   - method: Method of generation. Available options are:
%             - 'ar1': Generates an AR(1) process matrix.
%             - 'equi': Generates a matrix with equicorrelated noise.
%             - 'Gaussian': Generates a matrix with Gaussian noise.
%--------------------------------------------------------------------------
%   OUTPUT:
%   - X: n by m generated matrix according to the specified method.
%--------------------------------------------------------------------------
% EXAMPLES
% X = Xgen( 1000, 1000, 0.7, 'ar1' );
% mean(std(X,0,1))
%
% X = Xgen( 1000, 1000, 0.7, 'equi' );
% mean(std(X,0,1))
%
% X = Xgen( 1000, 1000, 0.7, 'Gaussian' );
% mean(std(X,0,1))
%--------------------------------------------------------------------------
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
if rho == 0
    X = randn(n,m);
    return
end

if strcmp(method, 'Gaussian')
    noise = randn([n,m]);
    [X, ss] = fconv(noise, rho, 1);
    X = X/sqrt(ss); % Standardize X
elseif strcmp(method, 'ar1')
    X = zeros(n,m);
    X(:,1) = randn(n,1);
    for I = 2:m
        %Generate noise to add
        epsilon = randn(n,1);
        %Calculate the updated column of the X matrix
        X(:,I) = rho*X(:,I-1) + epsilon;
    end
    %Standardize X
    X = X/sqrt(1/(1-rho^2));
elseif strcmp(method, 'equi')
    X = zeros(n,m);
    w = randn(1,n);
    for I = 1:n
        X(I,:) = sqrt(rho)*w(I) + randn(1, m)*sqrt(1-rho);
    end 
end

end
        % 
        % Generate noise to add
        % epsilon = randn(n,1);
        % Calculate the updated column of the X matrix
        % X(:,I) = rho*X(:,I-1) + epsilon;
        % Standardize X
        % X = X/sqrt((1-rho^2));

    %     for I = 1:n
    %         % Generate white noise
    %         epsilon = randn(1, m);
    %         % Generate AR(1) process
    %         X(I,:) = filter(1, [1, -rho], epsilon);
    %         % Standardize X
    %         X = X/sqrt((1-rho^2));
    % end