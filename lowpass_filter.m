% Jeff Arata
% 4/11/18

function [ h ] = lowpass_filter( fc, b, A )
% This function implemenets a lowpass filter using a Kaiser windowed sinc
% design as suggested in the following link:
% https://tomroelandts.com/articles/how-to-create-a-configurable-filter-using-a-kaiser-window
% along with some help from the paper 'Design of Low pass FIR Filters using
% Kaiser Window Function with variable parameter Beta (B)' by Rachna Arya
% and Shiva Jaiswal 
%
% Input:
% 
% fc -      the cutoff filter as a fraction of sampling rate
% b -       the transition bandwidth as a fraction of half the sampling rate
% A -       the desired attenuation of the passband in dB
% 
% Output:
% 
% h -       the filter coefficients satisfying the input parameters

b = b/2;    % Adjusted for fraction of full sampling rate

if A > 50                   % Set beta based on desired dB of attenuation
    beta = 0.1102*(A-8.7);  
elseif (A <= 50) && (A >= 21)
    beta = 0.5842*(A-21)^(0.4) + 0.07886*(A-21);
else
    beta = 0;    
end
% Filter length based on desired attenuation and transition bandwidth. N
% should always be odd for an odd filter length.
N = ceil((A-8)/(2.285*2*pi*b) + 1); 
if mod(N,2) == 0
    N = N + 1;
end

K = kaiser_window( beta, N );       % Finds Kaiser window
h_sinc = sinc_filter( fc, N );      % Gets sinc filter coefficients
h = h_sinc .* K;                    % Windows sinc filter with Kaiser window

end


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
