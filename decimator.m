% Jeff Arata
% 1/17/18

function [ y ] = decimator( x, M, F, filter_flag )
% This function takes a signal x and decimates it by a rate of M, reducing
% its sampling rate by the same factor.
%
% Inputs:
% 
% x -           the input signal
% M -           the rate at which to decimate the signal 
% N -           the number of stages for the CIC filter; defaults to 2 
% filter_flag - filters the signal if equal to 0
% 
% Outputs:
%
% y -       the output signal with a decreased sampling rate

if nargin < 3
    F = 3;
end
if nargin < 4
    filter_flag = 1;
end

y = x;
% CIC Filtering Method of Decimation
if (filter_flag && (F>2))           % Filter here if R>M
    F = F-2;                        % Filter signal F-2 times only
    y = dec_filter_amp(y, M, F);    % Filter and adjust amplitude
end
y = downsampler(y, M);      


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
