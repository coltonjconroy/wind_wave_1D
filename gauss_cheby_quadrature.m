function [xi,wi] = gauss_cheby_quadrature(n)

k  = 1:n;
xi = cos(pi*(2*k-1)/(2*n));
wi = ones(size(xi))*(pi/n);
%wi = wi./sqrt(1-xi.^2);
xi = sort(xi)';
wi = sort(wi)';

