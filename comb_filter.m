% Jeff Arata
% 2/28/18

function [ y ] = comb_filter( x , R )
% This function implements a comb filter for use in interpolation where the
% filter coefficients are all 0 except the first and last being 1 and -1
% respectively. The number of coefficients is R+1

b = zeros(1, 2*R);
b(1) = 1;
b(R+1) = -1;
a = 1;

y = filter(b, a, x);

end