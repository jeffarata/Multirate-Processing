% Jeff Arata
% 4/3/18

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
