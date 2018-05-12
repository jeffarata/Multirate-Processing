% Jeff Arata
% 4/10/18

function [ h ] = raised_cosine_filter( n, fc, fs, delta_f )
% This funciton implements a raised cosine filter and outputs its
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



% NOTE: SAMPLING RATE MUSt BE GREATER THAN TWICE SYMBOL RATE AT LEAST
% so maybe symbol rate is just the frequency of the signal? or the freq. of
% the max signal



ts = 1/fs;
%k = -ts*(n-1)/2:ts:ts*(n-1)/2; seems unneeded as the impulse response
%always includes k*ts in it;
k = -(n-1)/2:(n-1)/2;
t = k*ts;
R = 2*ts*delta_f;


h = (1/fs)*sinc(2*t*fc) .* cos(2*pi*R*t*fc) ./ (1-(4*R*t*fc).^2);

% may be unnecessary to declare this as it's just a property of the limit
% for t = 0;
%h(ceil(n/2)) = 1/fs; % ceil(n/2) should be middle term - essentially 0 on axis

%{
if sum(t == ts/(2*R)) == 1
    keyboard
    h(ts/(2*R) == t) = pi/(4*ts)*sinc(1/(2*R)); 
end
if sum(t == -ts/(2*R)) == 1
    keyboard
    h(-ts/(2*R) == t) = pi/(4*ts)*sinc(1/(2*R));
end
%}

if sum(t == 1/(4*R*fc)) == 1
    keyboard
    h(t == 1/(4*R*fc)) = (R/(2*fs))*sin(pi/(2*R));
end
if sum(t == -1/(4*R*fc)) == 1
    keyboard
    h(t == -1/(4*R*fc)) = (R/(2*fs))*sin(pi/(2*R));
end


% so R*fc == delta_f?????
% so 2*ts*delta_f*fc ==== delta_f????


% T - related to symbol rate on wikipedia should be the interpolation or
% decimation factor since we'll be using this for that. That will give us
% the correct cutoff frequency of pi/T

% I believe this filter and the root raised cosine filter wil haveto be
% windowed


% chec out:
% http://radio.feld.cvut.cz/matlab/toolbox/dspblks/digitalfirraisedcosinefil.html

% also:
% https://dspguru.com/dsp/reference/raised-cosine-and-root-raised-cosine-formulas/






end

