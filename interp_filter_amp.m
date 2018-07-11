% Jeff Arata
% 7/10/18

function [ y ] = interp_filter_amp( x, R, F)
% This function takes in a signal x, and filters it F times using 
% interpolation rate R to determine the coefficients of the CIC filter. It
% also scales the amplitude of the output signal based on the amount of
% filtering desired in the interpolation.
%
% Input:
%
% x -       the input signal
% R -       the interpolation rate
% F -       the number of times to apply to CIC_filter
%
% Output:
% 
% y -       the filtered and scaled output signal

y = CIC_filter(x, R, F);        % Filter F times
y = y/(R^(F-1));                % Adjust amplitude of signal

end