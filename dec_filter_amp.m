% Jeff Arata
% 7/11/18

function [ y ] = dec_filter_amp( x, M, F)
% This function takes in a signal x, and filters it F times using 
% decimation rate M to determine the coefficients of the CIC filter. It
% also scales the amplitude of the output signal based on the amount of
% filtering desired in the decimation.
%
% Input:
%
% x -       the input signal
% M -       the decimation rate
% F -       the number of filterings
%
% Output:
%
% y -       the filtered and scaled output signal

y = CIC_filter(x, M, F);   % Filter the signal F times
y = y/(M^(F));             % Adjust signal amplitude

end