function [x,dx] = geometric_distribution(XL,Xr,n,hx)
% n = number of elements
% XL = left coordinate
% Xr = right coordinate
if hx == 1
    gamma = 20/5;
elseif hx == 2
    gamma = 3/2;
elseif hx == 3
     gamma = 1680/1250;
elseif hx == 4
    gamma = 223/200;
elseif hx == 5
    gamma = 1.1; %gamma = 83/80;
elseif hx == 6
    gamma = 323/320;
end

dx = zeros(1,n);
x = zeros(1,n+1);

jj = 0;
for j = 0:n+1
    jj = jj + gamma^j;
end

xo = 50;


for j = 0:n;
    dx(j+1) = gamma^j*xo;
end

x(1) = XL;
for j = 2:n+1
    x(j) = x(j-1)+dx(j);
end
x(n+2) = Xr;