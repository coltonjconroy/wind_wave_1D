function [fp,alpha,U,E,nu,Hs] = hasselman_solns(x,t,Data,type)
%**************************************************************************
%
%                  two-parameter power law solution in
%          Hasselmann 1975 A simple parametric prediction model
%                             U = U(x,t)
%                     Written by Colton J. Conroy
%                               @ APAM
%                               9/9/15
%     
%**************************************************************************
Uo     = Data.Uo;
q      = Data.q;
g      = Data.g;
lambda = Data.lambda;
if strcmpi(type,'spatial')
    U     = Uo*(g*x./Uo^2).^q; % Wind field (spatially varying, frozen in time)
    A     = 2.84*(1+1.63*q)^(3/10);
    B     = 0.033*((1+1.50*q)/(1+1.63*q))^(1/2);
    C     = B*lambda;
    m     = (3/10)*(2*q-1);
    nu    = A*((g*x)./(Uo^2)).^m;  % Non-dimensional peak freq
    alpha = B*nu.^(2/3);      % Philips parameter
    epsil = C*nu.^(-10/3);    % Non-dimensional energy
    fp    = (nu.*g)./U;          % peak frequency
    E     = (epsil.*U.^4)./g^2;   % Energy density
    Hs    = 4.*E.^(1/2);
elseif strcmpi(type,'temporal')
    U     = Uo*(g*t/Uo)^q;
    A     = 16.8*(1+1.51*q)^(3/7);
    B     = 0.031*((1+1.33*q)/(1+1.51*q))^(1/2);
    C     = B*lambda;
    m     = (3/7)*(q-1);
    nu    = A*(g*t/Uo)^m; 
    alpha = B*nu^(2/3);      % Philips parameter
    epsil = C*nu^(-10/3);    % Non-dimensional energy
    fp    = nu*g/U;          % peak frequency
    E     = epsil*U^4/g^2;   % Energy density
    fp    = ones(size(x))*fp;
    alpha = ones(size(x))*alpha;
    E     = ones(size(x))*E;
end

