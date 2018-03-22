% Jeff Arata
% 1/17/18

function [ y ] = decimator( x, M )
% This function takes a signal x and decimates it by a rate of M, reducing
% its sampling rate by the same factor.
%
% Inputs:
% 
% x -       the input signal
% M -       the rate at which to decimate the signal 
% N -       the number of stages for the CIC filter; defaults to 2 
% 
% Outputs:
%
% y -       the output signal with a decreased sampling rate

if nargin < 3                   % Default number of stages to 2
   N = 2; 
end

% Filter method of decimation

y = CIC_filter(x, M, N);
y = downsampler(y, M);
y = y/(M^N);

% Ifft method of decimation
%{
X = fft(x);
len_X = length(X);
X_stuffing = zeros(1, M*len_X);
X_stuffing(1:floor(len_X/2)) = X(1:floor(len_X/2));
X_stuffing((M-1)*len_X+floor(len_X/2)+1:end) = X(floor(len_X/2)+1:end);
y = ifft(M*X_stuffing);
%}

end
