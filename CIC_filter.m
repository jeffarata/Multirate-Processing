% Jeff Arata
% 1/22/18

function [ y ] = CIC_filter( x, R, N )
% This function returns the filtered output of applying a CIC filter in N
% stages at a rate conversion of R
%
% Input:
%
% x -       the input signal
% R -       the rate conversion of the interpolator/decimator
% N -       the number of stages/number of times to filter the signal,
%           defaults to 2 (if it were 1, we would get zero order hold
%           interpolation i.e. a step function)
%
% Output:
%
% y -       CIC filtering output

if nargin < 3           % Defaults to 2 stages
   N = 2; 
end

b = zeros(1, 2*R);      % Builds CIC filter coefficients
b(1:R) = 1;
a = 1;

y = filter(b,a,x);      % Filter signal N times
for ii = 2:N
    y = filter(b,a,y);    
end
end
