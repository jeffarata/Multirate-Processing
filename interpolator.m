% Jeff Arata
% 1/17/18

function [ y ] = interpolator( x, R, F, filter_flag )
% This function interpolates a signal by upsampling by a rate of L and then
% filtering the upsampled signal with a lowpass filter.
%
% Inputs:
% 
% x -           the input signal
% R -           the interpolating/upsampling factor 
% N -           the number of stages for the CIC filter; defaults to 2 
% filter_flag - filters the signal if equal to 1
%
% Outputs:
%
% y -       the interpolated signal

if nargin < 3
    F = 3;
end
if nargin < 4
    filter_flag = 0;
end

y = x;
% CIC Filtering Method of Interpolation
y = upsampler(y, R);                % Add more samples
     
if ~filter_flag
    if F == 0
    else
        y = interp_filter_amp(y, R, F); % Filter and adjust amplitude   
    end
else
    if F > 2
        F = 2;
        y = interp_filter_amp(y, R, F); % Filter and adjust amplitude
    elseif (F == 2) || (F == 1)
        y = interp_filter_amp(y, R, F); % Filter and adjust amplitude
    else
    end
end


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
