% Jeff Arata
% 1/17/18

function [ y ] = interpolator( x, R )
% This function interpolates a signal by upsampling by a rate of L and then
% filtering the upsampled signal with a lowpass filter.
%
% Inputs:
% 
% x -       the input signal
% R -       the interpolating/upsampling factor
% N -           the number of stages for the CIC filter; defaults to 2 
%
% Outputs:
%
% y -       the interpolated signal

if nargin < 3                   % Default number of stages to 2
   N = 2; 
end

% CIC Filtering Method of Interpolation

y = upsampler(x, R);
y = CIC_filter(y, R, N);
y = y/(R^(N-1));


% Interpolation via zero phase stuffing in frequency domain and then using
% the inverse FFT. Only works when original signal satisfies nyquist
% criterion.
%{
X = fft(x);
len_X = length(X);
X_stuffing = zeros(1, R*len_X);
X_stuffing(1:floor(len_X/2)) = X(1:floor(len_X/2));
X_stuffing((R-1)*len_X+floor(len_X/2)+1:end) = X(floor(len_X/2)+1:end);
y = ifft(R*X_stuffing);
%}

% Interpolation using lagrange interpolating polynomial coefficients as the
% FIR filter coefficients. Involves quite a bit of delay though
%{
y = upsampler(x, R);
Q = 50;                     % Use 25 samples for lagrange coeffs
b = lagrange_coeff(Q, R);
y = filter(b, 1, y);
%}

end
