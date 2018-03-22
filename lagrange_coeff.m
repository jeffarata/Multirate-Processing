% Jeff Arata
% 2/20/18

function [ h ] = lagrange_coeff( Q, R )
% This function evaluates the lagrange polynomial coefficients used in
% interpolating an upsampled signal. The output coefficients are the FIR
% filter coefficients to apply to an upsampled signal to complete
% interpolation. Uses information from the paper "A DSP Approach to
% Interpolation" by Ronald Schafer and Lawrence Rabiner, 1973.
%
% Inputs:
%
% Q:        the number of samples being used to evaluate the coefficients
%           Q must be even!
% R:        the conversion rate of the upsampler/downsampler that we need
%           to filter to make an interpolator/ddecimator. Also the the
%           upper bound of the range over which to evaluate each block of
%           lagrange coefficients where we evaluate:
%               h(t) at t = n/R with 0<=n<R;
%
% Outputs:
%
% h:        the lagrange coefficients/the filter coefficients

k = -Q/2+1:Q/2;
n = 0:R-1;
t = n/R;

h = zeros(1, Q*R);      % Initialize filter coefficient vector

for ii = 1:length(k);
    product = 1;        % Initialize product as 1
    for jj = 1:Q
        product = product.*(t+Q/2-jj);
    end
    % Numerator for filter coefficient calculation, without 'product'
    numerator = ((-1)^(k(ii) + Q/2));
    % Denominator for filter coefficient calculation
    denominator = (factorial((Q-2)/2 + k(ii)) * factorial(Q/2 - k(ii)) * (t-k(ii)));
    
    h( (length(k)-ii)*R+1:(length(k)-ii+1)*R ) = (numerator ./ denominator) .* product;
    
end

h = h(2:end);
h(floor(length(h)/2)+1) = 1;    % Ensures center coefficient is 1

end
