% Jeff Arata
% 3/21/18

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
