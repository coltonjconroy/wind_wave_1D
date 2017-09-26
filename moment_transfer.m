%--------------------------------
% Energy transfer between phases 
%--------------------------------
for j = 1:nelems.x
    m0s = [PHI.L2]*m0(:,j,irk+1);     m1s = [PHI.L2]*m1(:,j,irk+1);
    m0w = [PHI.L2]*Twind0(:,j,irk+1); m1w = [PHI.L2]*Twind1(:,j,irk+1);
    fms = m1s./m0s;                   fmw = m1w./m0w;
    uj  = [PHI.L2]*u10(:,j);
    %------------------
    % swell to wind-sea
    %------------------
    fcutoff = 0.13*g./(uj);
    ipt = find(fms > fcutoff);
    dm0 = zeros(size(m0s)); dm1 = zeros(size(m1s)); df = zeros(size(m0s));
    df(ipt)  = fms(ipt)./fmw(ipt); % variation in frequency
    dpt = find(df > 1); df(dpt) = 1.0;
    dm0(ipt) = df(ipt).*m0s(ipt);  % variation in moments
    dm1(ipt) = df(ipt).*m1s(ipt);
    m0s = m0s - dm0;               m1s = m1s - dm1;
    m0w = m0w + dm0;               m1w = m1w + dm1;
    %------------------
    % wind-sea to swell
    %------------------
%     fpm = g*(0.13)./uj;
%     kpt = find(fmw < fpm);
%     dw0 = zeros(size(m0w)); dw1 = zeros(size(m1w)); 
%     m0pm = 0.0045/g^2.*uj.^4;   % new wind energy
%     m1pm = m0pm.*fpm;           % 1st moment
%     dw0(kpt) = m0w(kpt) - m0pm(kpt); % variation in moments
%     dw1(kpt) = m1w(kpt) - m1pm(kpt); % variation in moments
%     m0w(kpt) = m0pm(kpt);        m1w(kpt) = m1pm(kpt);
%     m0s      = m0s + dw0;        m1s      = m1s + dw1; 
    %----------------
    % update solution
    %----------------
    m0(:,j,irk+1)     = L2*m0s;
    m1(:,j,irk+1)     = L2*m1s;
    Twind0(:,j,irk+1) = L2*m0w;
    Twind1(:,j,irk+1) = L2*m1w;
end
        