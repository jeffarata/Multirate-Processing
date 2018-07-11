% Jeff Arata
% 1/16/18

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
