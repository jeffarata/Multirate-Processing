% Jeff Arata
% 4/10/18

function [ h ] = sinc_filter( fc, N )
% This function implements a sinc filter and returns the impulse response
% of one with n samples.
% see https://tomroelandts.com/articles/how-to-create-a-configurable-filter-using-a-kaiser-window
% Also used the paper 'Design of Low pass FIR Filters using Kaiser Window 
% Function with variable parameter Beta (B)' by Rachna Arya and Shiva Jaiswal 
%
% Input:
% 
% fc -      filter's cutoff frequency as a fraction of half the sampling rate
% N -       the length of the filter
%
% Output:
% 
% h -       the filter coefficients

fc = fc/2;  % Adjusted for fraction of full sampling rate
% Filter coefficient calculation
h = (2*fc) * (sin(2*pi*fc*[-(N-1)/2:(N-1)/2]) ./ (2*pi*fc*[-(N-1)/2:(N-1)/2]));
h(ceil(N/2)) = 2*fc;    % Set middle coefficient correctly

end
