function xk = roots_cheby(n,a,b)

% k = 0:n;
% top = (2*n + 1 - 2*k)*pi;
% bot = 2*n + 2;
% t_k = cos(top/bot);
% %x_k = zeros(1,length(t_k)+2);
% x_k = (b-a)/2*t_k + ones(size(t_k))*(a+b)/2;
% %x_k(1) = a;
% %x_k(end) = b;
% dx_k(1) = (x_k(1)-a)+(x_k(2)-x_k(1))/2;
% j = 2:n;
% dx_k(j) = (x_k(j+1)-x_k(j))/2+(x_k(j)-x_k(j-1))/2;
% dx_k(n+1) = (b-x_k(n))+(x_k(n)-x_k(n-1))/2;

C1 = (b-a)/2;
C2 = (b+a)/2;
xk = zeros(n+1,1);

for k = 0:n
    top = (2*n+1 - 2*k)*pi;
    bot = 2*n+2;
    tk  = cos(top/bot);
    xk(k+1) = C1*tk + C2;
end