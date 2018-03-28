% Jeff Arata
% 2/28/18

function [ y ] = integrator_filter( x, R )
% This function implements a simplie filter integrator and applies it to x

b = 1;
a = [1, -1];

y = filter(b, a, x);

end