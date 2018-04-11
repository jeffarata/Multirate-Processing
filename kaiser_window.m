% Jeff Arata
% 4/10/18

function [ K ] = kaiser_window( B, N )
% This function implements a Kaiser window of length N and with rolloff
% factor B.
% see https://en.wikipedia.org/wiki/Kaiser_window
% Also used the paper 'Design of Low pass FIR Filters using Kaiser Window 
% Function with variable parameter Beta (B)' by Rachna Arya and Shiva Jaiswal 
%
% Input:
%
% B -       rolloff factor
% N -       filter length
%
% Output:
% 
% K -       the values of the Kaiser window

% Calculate Kaiser Window values
K = besseli(0, B*sqrt(1-(2*[0:N-1]/(N-1)-1).^2)) ./ besseli(0, B);

end
