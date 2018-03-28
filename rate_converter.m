% Jeff Arata
% 3/13/18


% Another idea to try would be to go for implementing the interpolation and
% decimation in stages by doing it in prime factors of the desired values
% to see if that takes care of it. This may require filtering with each
% interp/decimation by a prime factor. Not sure though, will have to look
% into this.


function [ y ] = rate_converter( x, R, M, N )
% This function converts the sampling rate of a signal x changing it such
% that (if the sampling rate of x is 'fs_x'):
% fs_y = (R/M)*fs_x
% This is done by first interpolating x by a factor of R and then
% decimating the result by a factor of M;
%
% Inputs:
% 
% x -       the input signal
% R -       the rate to interpolate the signal x
% M -       the rate to decimate the resulting interpolated signal
% N -       the number of stages to for applying the CIC filter
%
% Outputs:
%
% y -       the sampling rate converted signal

filter_flag = R < M;    % If filter_flag is 1, then the filter is applied 
                        % in the decimator function. If it is 0, the
                        % filter is not applied in the decimator function.
if R == M               % If R and M are the same (no rate conversion), 
    filter_flag = [];   % then no filtering should take place and so y = x. 
end                     % There is always filtering in the interpolation 
                        % function, except the case of no rate conversion.

y = interpolator(x, R, N, filter_flag);
% try implementing the rate conversion in factors of the conversion rate
% for example, interpolating by a factor of 8 could be done by
% interpolating by 2, 3 separate times (8 = 2^3)
% Probably write a function that breaks a number down into its prime
% factors for this
y = decimator(y, M, N, filter_flag);

end


function [ y ] = interpolator( x, R, N, filter_flag )
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

if nargin < 3                   % Default number of stages to 2
   N = 2; 
end

y = upsampler(x, R);            % Add more samples

if isempty(filter_flag)
else % Filter here if R>M or R<M (everytime other than no rate conversion)
    y = CIC_filter(y, R, N);    % Filter signal
    y = y/(R^(N-1));            % Adjust amplitude of signal
end

end


function [ y ] = decimator( x, M, N, filter_flag )
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

if nargin < 3                   % Default number of stages to 2
   N = 2; 
end
% if filtering happens here, a second filtering means more delay (twice the
% delay???)
if filter_flag                  % Filter here if R<M
    y = CIC_filter(x, M, N);    % Filter the signal
    y = downsampler(y, M);      % Pick every Mth sample
    y = y/(M^N);                % Adjust signal amplitude
else
    y = downsampler(x, M);      
end

end


function [ y ] = CIC_filter( x, R, N )
% This function returns the filtered output of applying a CIC filter in N
% stages at a rate conversion of R
%
% Input:
%
% x -       the input signal
% R -       the rate conversion of the interpolator/decimator
% N -       the number of stages/number of times to filter the signal,
%           defaults to 2 (if it were 1, we would get zero order hold
%           interpolation i.e. a step function)
%
% Output:
%
% y -       CIC filtering output

if nargin < 3
   N = 2; 
end

b = zeros(1, 2*R);      % Create CIC filter coefficient vectors.
b(1:R) = 1;
a = 1;

y = filter(b,a,x);      % Filter the signal N times.
for ii = 2:N            
    y = filter(b,a,y);    
end

end


function [ y ] = upsampler( x, r, p )
% This function upsamples a signal by inserting r-1 zeros between the
% input samples, starting with the first sample. Optionally, it offsets the
% upsampling by n samples.
%
% Input:
% 
% x -       the input signal
% r -       the rate of upsampling, must be an integer
% p -       the offset of phase of the upsampled signal, defaults to 0
%           0 <= p <= r-1
%
% Output:
% 
% y -       the input signal upsampled by a rate of r and offset by n

if nargin < 2
    error('Not enough inputs. Expected at least 2 inputs.');
elseif mod(r,1) ~= 0
    error('Upsampling rate r must be an integer value.');
elseif nargin < 3
    p = 0;
elseif p > r-1
    error('Phase offset n must be less than r-1 samples.');
end

y = zeros(1, r*length(x));  % Initialize output to correct length
y(1+p:r:end) = x;           % Fill output with input spaced as desire by r
                            % and delayed by n
end


function [ y ] = downsampler( x, r, p )
% This function downsamples a signal by keeping every r-th data point of 
% the input samples, starting with the first sample. Optionally, it offsets
% the downsampling by n samples.
%
% Input:
% 
% x -       the input signal
% r -       the rate of downsampling, must be an integer
% p -       the offset of phase of the upsampled signal, defaults to 0
%           0 <= p <= r-1
%
% Output:
% 
% y -       the input signal downsampled by a rate of r and offset by n

if nargin < 2
    error('Not enough inputs. Expected at least 2 inputs.');
elseif mod(r,1) ~= 0
    error('Downsampling rate r must be an integer value.');
elseif nargin < 3
    p = 0;
elseif p > r-1
    error('Phase offset n must be less than r-1 samples.');
end

y = x(1+p:r:end);   % Index out every r-th sample, starting with the n+1-th
                    % sample.
end
