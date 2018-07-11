% Jeff Arata
% 3/13/18


% What if you gave the rate_converter function the option to put in the starting sampling
% rate and ending sampling rate instead of the interpolation and decimation
% rates. You would just have to calculate these rates from the input
% sampling rates. Would have issues with irrational combinations, so maybe
% throw an error? Or try to approximate something rational that's close? Or
% throw an error and say that we went for something close anyways (saying
% what the new rate was/new end sampling rate was)



% RAISED COSINE FILTERS - impulse response on wikipedia - an
% implementation (one kind of) nyquist filter
% ROOT RAISED COSINE FILTERS - impulse respone on wikipedia
% Nyquist filters - aka Mth band filters (more general version of haldband
% filters - there is an algorithm apparently for equiripple nyquist filters
% KAISER Window - on wikipedia, could be good for filtering, use the window
% and apply it the the impulse response in order to taper the impulse
% response before filtering
% It seems that Nyquist Filters are a general type of filter satisfying
% certain properties. Raised cosine filters satisfy these properties and so
% are nyquist filters. I believe root raised cosine filters also satisfy
% the nyquist filter conditions. An equiripple filter is something that
% comes out of adjusting an existing filter. In order to do that, look into
% the parks-mcclellan method and the remez exchange algorithm.


function [ y ] = rate_converter( x, R, M, N, F )
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
% N -       the number of stages of the rate conversion
% F -       the number of times to filter per stage
%
% Outputs:
%
% y -       the sampling rate converted signal
if nargin < 4       % Default to one stage for up/downsapling
    N = 1;
end
if nargin < 5       % Default to 3 filters per stage
    F = 3;
end

filter_flag = R > M;    % If filter_flag is 1, then 1 filter is applied in
                        % interpolation and (F-1) filters are applied in 
                        % decimation, per rate stage. If it is 0, all F 
                        % filters are applied in interpolation, per rate 
                        % stage                        
                        
if R == M               % If R and M are the same (no rate conversion), 
    filter_flag = [];   % then no filtering should take place and so y = x. 
end                     % There is always filtering in the interpolation 
                        % function, except the case of no rate conversion.

y = x;                  % Initialize output as the input signal x      
                        
if ~isempty(filter_flag)    % Does nothing if there is no rate conversion                     
    % Get the staging factors for the interpolation and decimation rates
    % Ensure that if the number of stages possible for either is lower than
    % expected, that both are factored to the same number of stages/factors.
    [sf_int, N_int] = rate_stager(R, N);  % N stages of upsampling to total to R.
    [sf_dec, N_dec] = rate_stager(M, N);  % N stages of downsampling to total to M.                      
    if N_int > N_dec
        [sf_int, N_int] = rate_stager(R, N_dec);
    elseif N_int < N_dec
        [sf_dec, N_dec] = rate_stager(M, N_int);
    end
    
    % Determines where filtering should occur at each stage. If one stage
    % has the same interpolation and decimation rate, then no filtering or 
    % up/downsampling will occur at that stage, for efficiency. 
    for ii = 1:N_int           
        if (sf_int(ii) ~= sf_dec(ii))        
            stage_filter_flag = sf_int(ii) > sf_dec(ii);
            y = interpolator(y, sf_int(ii), F, stage_filter_flag);
            y = decimator(y, sf_dec(ii), F, stage_filter_flag);
        else
        end
    end
end
end

function [ y ] = interpolator( x, R, F, filter_flag )
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


y = x;

y = upsampler(y, R);                % Add more samples
     
if ~filter_flag
    if F == 0
    else
        y = interp_filter_amp(y, R, F); % Filter and adjust amplitude   
    end
else
    if F > 2
        F = 2;
        y = interp_filter_amp(y, R, F); % Filter and adjust amplitude
    elseif (F == 2) || (F == 1)
        y = interp_filter_amp(y, R, F); % Filter and adjust amplitude
    else
    end
end
    
end


function [ y ] = interp_filter_amp( x, R, F)
% This function takes in a signal x, a interpolation rate R, and number of
% filterings F, and filters x using a CIC filter following by scaling its
% amplitude accordingly.
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


function [ y ] = decimator( x, M, F, filter_flag )
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

y = x;
if (filter_flag && (F>2))           % Filter here if R>M
    F = F-2;                        % Filter signal F-2 times only
    y = dec_filter_amp(y, M, F);    % Filter and adjust amplitude
end
y = downsampler(y, M);      

end


function [ y ] = dec_filter_amp( x, M, F)
% This function takes in a signal x, and filters it F times using 
% decimation rate M to determine the coefficients of the CIC filter. It
% also scales the amplitude of the output signal based on the amount of
% filtering.
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


function [ y ] = CIC_filter( x, R, F )
% This function returns the filtered output of applying a CIC filter in N
% stages at a rate conversion of R
%
% Input:
%
% x -       the input signal
% R -       the rate conversion of the interpolator/decimator
% F -       the number of times to filter the signal, defaults to 3 (if it
%           were 1, we would get zero order hold interpolation i.e. a step,
%           function, and 2 is not very smooth either)
%
% Output:
%
% y -       CIC filtering output

if nargin < 3
   F = 3; 
end

b = zeros(1, 2*R);      % Create CIC filter coefficient vectors.
b(1:R) = 1;
a = 1;

y = x;      
for ii = 1:F           % Filter the signal F times. 
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


function [ sf, Nout ] = rate_stager( R, f )
% This function takes in a number, breaks it into its prime factors, and
% condenses those to how ever many factors the user prefers. This will
% default to 3. More than the number of prime factors will print an error
% and instead give the most factors possible. These factors will be
% increasing.
%
% Inputs:
%
% R -       the rate of conversion to be broken into factors
% f -       the number of factors to break down R into
%
% Outputs:
% 
% sf -      the stage factors of R
% Nout -    the length of sf/the number of actual output stages

if nargin < 2
    f = 3;
end

pf = prime_factor(R);   % Get prime factors of R

if length(pf) > f       % When number of factors requested is less than in
    sf = pf(1:f);       % the number of prime factors
    
    % Multiply leftovers in pf to sf vector
    for ii = (f+1):length(pf)
        sf(mod(ii-1,f)+1) = sf(mod(ii-1,f)+1) * pf(ii);
    end             
    sf = sort(sf);      % Sort to be increasing order
    
elseif length(pf) == f
    sf = pf;    
else % when length(pf) < f
    sf = pf;
    warning(['Requested %1.f factor(s) but, the number %1.f has %1.f factor(s) ' ...
        'at most. Using %1.f factor(s).\n'], f, R, length(pf), length(pf)); 
end

Nout = length(sf);

end


function [ pf ] = prime_factor( x )
% This function finds the prime factors of x and outputs them as a vector
% pf.
%
% Inputs:
% 
% x -   the input number to find the prime factors of
%
% Outputs:
% 
% pf -  a vector of the prime factors of x

pf = [];        % initialize output

if x == 1       % if input is 1, output is 1
    pf = 1;
else
    while (mod(x, 2) == 0)  % Find all factors of 2 and add them to ouput
        pf = [pf, 2];
        x = x/2;
    end
    for ii = 3:2:sqrt(x)        % All other prime factors must be odd
        while (mod(x,ii) == 0)  % Find the other prime factors
            pf = [pf, ii];
            x = x/ii;
        end
    end
end

if x > 2                % If x is greater than 2 by this point, it must be
    pf = [pf x];        % be a prime input, so the output is just x.
end

end

