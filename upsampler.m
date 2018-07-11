% Jeff Arata
% 1/16/18

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
