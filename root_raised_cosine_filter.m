% Jeff Arata
% 4/10/18

function [ h ] = root_raised_cosine_filter(  n, fc, fs, delta_f )
% This funciton implements a root raised cosine filter and outputs its
% coefficients.
%
% Input:
%
% n -       length of filter - n should be odd
% fc -      cutoff frequency - aka fn0 - fraction of sampling rate
% fs -      sampling frequency or symbol rate???
% R -       rolloff factor - between 0 and 1, smaller = faster roll off
% delta_f - maybe use this instead of R - transition bandwidth as fraction
% of sampling rate
%
% Output:
% 
% h -       filter coefficients



% can be found on wikipedia

% also:
% http://radio.feld.cvut.cz/matlab/toolbox/dspblks/digitalfirraisedcosinefil.html

% and:
% https://dspguru.com/dsp/reference/raised-cosine-and-root-raised-cosine-formulas/


ts = 1/fs;
%k = -ts*(n-1)/2:ts:ts*(n-1)/2; seems unneeded as the impulse response
%always includes k*ts in it;
k = -(n-1)/2:(n-1)/2;
t = k*ts;
R = 2*ts*delta_f;


h_num = (4*R*cos((1+R)*2*pi*t*fc) + (sin((1-R)*2*pi*t*fc) ./ (8*R*t*fc)));
h_den = pi*fs*sqrt(1/(2*fc))*((8*R*t*fc).^2-1);

h = h_num ./ h_den;




end

