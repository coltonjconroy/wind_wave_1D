%**************************************************************************
%
%                      Wind-wave swell model
%                    (moment field equations) 
%                  Written by Colton J. Conroy
%                            @ APAM
%                           09/21/16
%     
%**************************************************************************
clear all; 
%--------------------------------------------------------------------------
% User Input 
%--------------------------------------------------------------------------
fig         = 'on';
vid         = 'off';
S           = 'on';
discret     = 'geo';
optimal     = 'on';
type        = 'spatial';  
basis       = 'Legendre'; % basis (Legendre polynomials or chebyshev?)
slimiting   = 'off';      % slope limiter
dtw         = 2*3600;     % equilibrium time
relax_time  = 4.0*3600;   % amount of time before swell flux is turned on
T           = 3.0*3600;  % total amount of time
wind        = 'off';
swell       = 'on';
wind_swell  = 'on';
windsea     = 'off';     
test        = 'temporal';  % spatial kernel or temporal kernel
wind_interp = 'elem';
windfig     = 'on';
swellfig    = 'off';     
transfer    = 'on';       % moment transfer for  wind-sea-swell
analytic    = 'off';
flux_type   = 'upwind';   % upwind or LLF       
C_plot      = 'off';      % plot of constants (Co & C1) for source term
if T == 0
    C_plot = 'off';
end
f_error    = 'on';        % 'on' = error in freq, 'off' = error in m1
%--------------------------------------------------------------------------
% Input Hasselmann problem data (power law test)
%--------------------------------------------------------------------------
DATA.geom1    = 10;      % can't use 0 b/c of singularirty in analytic soln
DATA.geom2    = 150e2;   % Right boundary location
DATA.Lx       = DATA.geom2 - DATA.geom1; % domain length
DATA.sigA     = 0.07;    % Avg. JONSWAP value
DATA.sigB     = 0.09;    % ``               "
DATA.sig      = (1/2)*(DATA.sigA + DATA.sigB);
DATA.gamma    = 3.30;    % Avg. JONSWAP value
DATA.K        = DATA.sig^2/(20*DATA.sig^2 + log(DATA.gamma));
DATA.lambda   = 1.58e-4; % avg JONSWAP shape parameter
DATA.lambda2  = 1.375e-4;
DATA.Uo       = 22.50;   % spatial wind parameters
DATA.q        = 0.11;    % "                     " 
DATA.g        = 9.81;    % acceleration due to gravity
DATA.A        = 2.84*(1+1.63*DATA.q)^(3/10); % spatial analytic soln
DATA.B        = 0.033*((1+1.50*DATA.q)/(1+1.63*DATA.q))^(1/2); % "   "
DATA.m        = (3/10)*(2*DATA.q-1);                           % "   "
lambda        = DATA.lambda;  % shape factor 
%--------------------------------------------------------------------------
% Input for DG method
%--------------------------------------------------------------------------
% spatial
%--------
p = 5; pd = 5; ndof_1d = p+1; L2_pts = 10;
%---------
% spectral
%---------
pf = 1; ndof_1df = pf+1; L2_ptsf = 10; Sf = 'on';
%--------------------------------------------------------------------------
% Source and flux functions
%--------------------------------------------------------------------------
Nu     = @(fp,U)  fp.*U/DATA.g;
Etilde = @(E,U)   E*DATA.g.^2./(U.^4);
Cg     = @(fp)    DATA.g./(4.*pi.*fp);
Te     = @(E)     (E./DATA.a).^(1/DATA.b);     
Tw     = @(f)     ((f)./DATA.c).^(1/DATA.d);
So     = @(x)     -(DATA.B.*DATA.lambda.*DATA.Uo.^5*(13.*DATA.m - ...
                  15.*DATA.q)*((DATA.g.*x)./DATA.Uo.^2).^(5.*DATA.q))./...
                  (12.*pi.*DATA.g.^2.*x.*(DATA.A.*((DATA.g.*x)./...
                  DATA.Uo.^2).^DATA.m).^(13/3));
S1     = @(x)     -(DATA.A.*DATA.B.*DATA.lambda.*DATA.Uo.^4.*(5.*DATA.m - 6.*...
                  DATA.q).*((DATA.g.*x)./DATA.Uo.^2).^(DATA.m + 4.*DATA.q))./...
                  (6.*pi.*DATA.g.*x.*(DATA.A.*((DATA.g.*x)./...
                  DATA.Uo.^2).^DATA.m).^(13/3));
Sot    = @(t,Uo,q,g) (3875*10^(4/21)*21^(2/3)*lambda*Uo^4*((133*q + 100)/(151*q + 100))^(1/2)*...
         (9*q + 5)*(151*q + 100)^(3/7)*((g*t)/Uo)^((31*q)/7 - 3/7))/...
         (87127488*g^2*t*((151*q + 100)^(3/7)*((g*t)/Uo)^((3*q)/7 - 3/7))^(13/3));
S1t    = @(t,Uo,q,g) (155*10^(1/3)*21^(2/3)*Uo^3*lambda*((133*q + 100)/(151*q + 100))^(1/2)*...
         (2*q + 1)*(151*q + 100)^(6/7)*((g*t)/Uo)^((27*q)/7 - 6/7))/...
         (592704*g*t*((151*q + 100)^(3/7)*((g*t)/Uo)^((3*q)/7 - 3/7))^(13/3));
E      = @(t,Uo,q,g) (31*Uo^4*lambda*(((133*q)/100 + 1)/((151*q)/100 + 1))^(1/2)*...
          ((g*t)/Uo)^(4*q))/(1000*g^2*((84*((151*q)/100 +...
          1)^(3/7)*((g*t)/Uo)^((3*q)/7 - 3/7))/5)^(10/3));    
Fm     = @(t,Uo,q,g) (42*10^(1/7)*g*(151*q + 100)^(3/7))/...
            (25*Uo*((g*t)/Uo)^((4*q)/7 + 3/7));     
Feq    = @(t,Uo,q,g) -(294*10^(1/7)*pi*g^2*t*(151*q + 100)^(3/7))/(25*Uo*(q...
          - 1)*((g*t)/Uo)^((4*q)/7 + 3/7));
Sw     = @(t,Uo,q,g) -(6*10^(1/7)*g*(4*q + 3)*(151*q + 100)^(3/7))/...
         (25*Uo*t*((g*t)/Uo)^((4*q)/7 + 3/7));    
Alpha  =  @(t,Uo,q,g) (31*10^(16/21)*21^(2/3)*((133*q + 100)/(151*q + ...
          100))^(1/2)*((151*q + 100)^(3/7)*((g*t)/Uo)^((3*q)/7 ...
          - 3/7))^(2/3))/25000;     
U10    = @(t,Uo,q,g) Uo.*(g.*t./Uo).^q;  

DcgDt  = @(t,Uo,q,g) (5.*10.^(6/7).*g.*((4.*q)./7 + 3/7).*((g.*t)./Uo).^((4.*q)./7 - ...
          4./7))./(336.*pi.*(151.*q + 100).^(3/7));
S0T    = @(t,Uo,q,g) -(19375*10^(4/21)*21^(2/3)*g*lambda*((133*q + 100)/...
          (151*q + 100))^(1/2)*((3*q)/7 - 3/7)*(151*q + 100)^(3/7)*...
          ((g*t)/Uo)^((3*q)/7 - 10/7))/(37340352*Uo*((151*q + ...
          100)^(3/7)*((g*t)/Uo)^((3*q)/7 - 3/7))^(13/3));
S1T    = @(t,Uo,q,g) -(155*10^(1/3)*21^(2/3)*lambda*((133*q + 100)/...
          (151*q + 100))^(1/2)*(151*q + 100)^(6/7)*(q - 1)*...
          ((g*t)/Uo)^((6*q)/7 - 6/7))/(592704*t*((151*q + 100)^(3/7)*...
          ((g*t)/Uo)^((3*q)/7 - 3/7))^(13/3));
S0nd   = @(t,Uo,q,g) -(19375*10^(4/21)*21^(2/3)*g*lambda*((133*q + 100)/(151*q + ...
           100))^(1/2)*((3*q)/7 - 3/7)*(151*q + 100)^(3/7)*((g*t)/...
           Uo)^((3*q)/7 - 10/7))/(37340352*Uo*((151*q + 100)^(3/7)*((g*t)...
           /Uo)^((3*q)/7 - 3/7))^(13/3));
       
S1nd   = @(t,Uo,q,g) -(155*10^(1/3)*21^(2/3)*lambda*((133*q + 100)/(151*q ...
           + 100))^(1/2)*(151*q + 100)^(6/7)*(q - 1)*((g*t)/Uo)^((6*q)/7 ...
           - 6/7))/(592704*t*((151*q + 100)^(3/7)*((g*t)/Uo)^((3*q)/7 - ...
           3/7))^(13/3));
       
%--------------------------------------------------------------------------
% Domain discretization
%--------------------------------------------------------------------------
% wind-sea
%--------
x0 = DATA.geom1; xN = DATA.geom2; 
hx = 1; 
dx = (xN-x0)/((4)*2^(hx-1)); dt = 1.00; 
X = x0:dx:xN;
if strcmpi(discret,'cheby')
    nelems.x = length(X)-1;
    sig = roots_cheby(nelems.x-1,x0,xN);
    X(1) = x0; X(2:nelems.x) = sig(1:nelems.x-1) ; X(end) = xN;
    XF = X;
elseif strcmpi(discret,'geo')
    nelems.x = length(X)-1;
    X = geometric_distribution(x0,xN,nelems.x-2,hx);
end
[nnodes,nelems,xNODE,~,xELEM] = XZ_1d_mesh(X,X);
%---------
% spectral
%---------
hf = 8; f0 = 0; fN = 10;
df = (fN-f0)/((4)*2^(hf-1)); f = f0:df:fN; 
[~,felems,~,~,fELEM] = XZ_1d_mesh(f,f);
%--------------------------------------------------------------------------
% Time stepping 
%--------------------------------------------------------------------------
if strcmpi(optimal,'on') % optimal time steppers 'on' or 'off'?
    if p > 0 && p < 2
        [alpha,beta,crk] = RKssp(p+1,p+2); nrk = length(alpha);
    elseif p >= 2
        [alpha,beta,crk] = RKssp(4,5);     nrk = length(alpha);
    else
        [alpha,beta,crk] = RKssp(p+1,p+1); nrk = length(alpha);
    end    
else
    [alpha,beta,crk] = SSPRK(p,p+1);       nrk = length(alpha); % number or stages
end   
%--------------------------------------------------------------------------
% Video setup 
%--------------------------------------------------------------------------
if strcmpi(vid,'on')
    writerObj = VideoWriter(['freq_int_1d_windsea_v6',num2str(p)]);
    writerObj.Quality = 100;
    open(writerObj);
end
%--------------------------------------------------------------------------
% Get LDG 1D matricies
%--------------------------------------------------------------------------
% spatial
%--------
if strcmpi(basis,'cheby')
    [L2,A,B,C,PHI,PSI] = LDG_matrices_1d(p,S,L2_pts);
    [L2d,Ad,Bd,Cd,~,~,DPHI] = CDG_matrices_1d(pd,S,L2_pts); % cheby for spectral anaylsis
else
    [L2,A,B,C,PHI,PSI] = LDG_matrices_1d(p,S,L2_pts);
    [L2d,Ad,Bd,Cd,~,~,DPHI] = LDG_matrices_1d(pd,S,L2_pts);
end
%---------
% spectral
%---------
[L2ptsf,L2wts] = gauss_jacobi_quad(L2_pts,0,0);
[L2f,~,~,~,PHIf,PSIf] = LDG_matrices_1d(pf,Sf,L2_ptsf);
%--------------------------------------------------------------------------
% Initialize variables
%--------------------------------------------------------------------------
m0     = zeros([p+1,nelems.x,nrk+1]); % swell phase
m1     = zeros([p+1,nelems.x,nrk+1]);
Fbar   = zeros([p+1,nelems.x]);
RHS_m0 = zeros([p+1,nelems.x,nrk]);   % rhs for swell phase
RHS_m1 = zeros([p+1,nelems.x,nrk]); 
u10    = zeros([p+1,nelems.x]);       % wind magnitude 
Twind0 = zeros([p+1,nelems.x,nrk+1]); % wind phase
Twind1 = zeros([p+1,nelems.x,nrk+1]); 
wind0  = zeros([p+1,nelems.x,nrk]);   % rhs for wind-sea
wind1  = zeros([p+1,nelems.x,nrk]);   
smult0 = zeros(L2_pts,nelems.x);      % source constants
smult1 = zeros(L2_pts,nelems.x);
%--------------------------------------------------------------------------
% Initial condition
%--------------------------------------------------------------------------
NT = floor(T/dt)-1; 
NTw = floor(dt/dtw); 
if strcmpi(wind_interp,'nodes')
    wind_spatial_1d
else
    wind_spatial_1dv2
end
fp = zeros(nnodes.x,1); ep = zeros(nnodes.x,1);
m0_j = zeros(p+1,1); m1_j = zeros(p+1,1); 
%tw = dtw*3.75;  % Uo = 10 q = 0.15
%tw = dtw/1.25; %Uo = 10 q = 0.05
%tw = dtw*1.7; %Uo = 10 q = 0.10
%tw = dtw*2.5; % Uo = 15 q = 0.10
%tw = dtw/1.4; %tw = dtw/1.785;
%tw = dtw/2.4;  % Uo = 4.3 q = 0.004
%tw  = dtw/2.235; % Uo = 8 q = 0.01
%tw  = dtw/2.365; % Uo = 8 q = 0.005
%tw  = dtw/1.38;   % Uo = 8 q = 0.05
%tw   = dtw/.65; % Uo = 8 q = 0.1
%tw = dtw/1.75; % Uo = 14.0 q = 0.02
tw = dtw*.4; % Uo = 22.50 q = 0.11 % elem
for j = 1:nelems.x
    if strcmpi(test,'spatial') || strcmpi(test,'temporal')
        x = [PSI.xa]*X(xELEM(j).nodes)';
        [fp,~,~,ep] = hasselman_solns(x,0,DATA,type);
        check = DATA.Uo;
        if check <= 3
            m0(:,j,1)  = L2*ones(size(ep))*1e-6;
            m1(:,j,1)  = L2*(ones(size(ep))*1e-6);
            if T < 1e-6 || relax_time < 1e-6
                Twind0(:,j,1)  = L2*ep;
                Twind1(:,j,1)  = L2*(ep.*fp);
            else
                Twind0(:,j,1)  = L2*ones(size(ep))*1e-6;
                Twind1(:,j,1)  = L2*(ones(size(ep))*1e-6.*fp);
            end
        else
            m0(:,j,1) = L2*ones(size(ep))*1e-7;
            m1(:,j,1)  = L2*(ones(size(ep))*1e-7.*fp);
            if T < 1e-6 || relax_time < 1e-6
                Twind0(:,j,1)  = L2*ep;
                Twind1(:,j,1)  = L2*(ep.*fp);
            else
                Twind0(:,j,1) = L2*ones(size(ep))*1e-7;
                Twind1(:,j,1) = L2*(ones(size(ep))*1e-6.*fp);
            end
        end
        if strcmpi(wind_interp,'elem')
            u10(:,j)= L2*(Uo(:,j).*(DATA.g*(tw)./Uo(:,j)).^(q(:,j)));
        else
            U0i  = U10(tw,Uo(1,j),q(1,j),DATA.g);
            U0i1 = U10(tw,Uo(1,j+1),q(1,j+1),DATA.g);
            U0j  = [PSI.xa]*[U0i; U0i1];
            u10(:,j) = L2*U0j;
        end
     elseif strcmpi(windsea,'on')
         x = [PSI.xa]*X(xELEM(j).nodes)';
         [fp,~,~,ep] = hasselman_solns(x,0,DATA,type);
         m0(:,j,1)   = L2*(ep.*1e-3);     Twind0(:,j,1) = m0(:,j,1);
         m1(:,j,1)   = L2*(fp.*ep.*1e-3); Twind1(:,j,1) = m1(:,j,1);
         u10(:,j)    = L2*(Uo(:,j).*(DATA.g*(tw)./Uo(:,j)).^(q(:,j)));
    else
        x = [PSI.xa]*X(xELEM(j).nodes)';
        [fp,~,~,ep] = hasselman_solns(x,0,DATA,type);
        m0(:,j,1)   = L2*(ep.*1e-3);
        m1(:,j,1)   = L2*(fp.*ep.*1e-3);
        u10(:,j)    = L2*(Uo(:,j).*(DATA.g*(tw)./Uo(:,j)).^(q(:,j)));
    end
  
end
%--------------------------------------------------------------------------
% Runge-Kutta time stepping
%--------------------------------------------------------------------------
 t = 0; m = 1; time = 0; plot_n = 0;
for n = 0:NT
    for irk = 1:nrk
        t = time + crk(irk)*dt;
        RHS_moment_field
        for i = 1:irk
            m0(:,:,irk+1) = m0(:,:,irk+1) + alpha(irk,i)*m0(:,:,i) ...
                + beta(irk,i)*dt*RHS_m0(:,:,i);
            m1(:,:,irk+1) = m1(:,:,irk+1) + alpha(irk,i)*m1(:,:,i) ...
                + beta(irk,i)*dt*RHS_m1(:,:,i);

            Twind0(:,:,irk+1) = Twind0(:,:,irk+1) + alpha(irk,i)*Twind0(:,:,i) ...
                + beta(irk,i)*dt*wind0(:,:,i);
            Twind1(:,:,irk+1) = Twind1(:,:,irk+1) + alpha(irk,i)*Twind1(:,:,i) ...
                + beta(irk,i)*dt*wind1(:,:,i);

        end
        if strcmpi(slimiting,'on') % slope limit the solution
            m0(:,:,irk+1) = slope_limiter(m0(:,:,irk),[PHI(1:2).b1],[PHI(1:2).b2],xELEM);
            m1(:,:,irk+1) = slope_limiter(m1(:,:,irk),[PHI(1:2).b1],[PHI(1:2).b2],xELEM);
        end
        if strcmpi(transfer,'on')
            moment_transfer;
        end
    end
    m0(:,:,1) = m0(:,:,nrk+1); m0(:,:,2:end) = 0; RHS_m0(:,:,:) = 0;
    m1(:,:,1) = m1(:,:,nrk+1); m1(:,:,2:end) = 0; RHS_m1(:,:,:) = 0;

    Twind0(:,:,1) = Twind0(:,:,nrk+1); Twind0(:,:,2:end) = 0; wind0(:,:,:) = 0;
    Twind1(:,:,1) = Twind1(:,:,nrk+1); Twind1(:,:,2:end) = 0; wind1(:,:,:) = 0;

    % Movie
    %----------------------------------------------------------------------
        if strcmpi(vid,'on')
            if n == plot_n
                if time >= relax_time
                    vid_1D_HISWA
                end
                plot_n = plot_n + floor(50/dt);
            end
        end
    time = time + dt;
end
if strcmpi(vid,'on')
    close(writerObj);
end
%--------------------------------------------------------------------------
% % Plot solution
%--------------------------------------------------------------------------
if strcmpi(fig,'on')
    if strcmpi(swellfig,'on')
        wind = 'off';
    else
        wind = 'on';
        swell = 'off';
        if relax_time >= T
            wind_swell = 'off';
        end
    end
    plot_1D_HISWA
end
%--------------------------------------------------------------------------
% L2 errors
%--------------------------------------------------------------------------
if relax_time > T % calc errors for wind-sea case (not wind-sea-swell)
    [x,w] = gauss_jacobi_quad(L2_pts,0,0);
    m0_L2_error = 0;
    m1_L2_error = 0;
    uw_L2_error = 0;
    for j = 1:nelems.x
        x  = [PSI.xa]*X(xELEM(j).nodes)';
        [fp,ALPHA,~,Em0,nu,hs]  = hasselman_solns(x,time,DATA,type);
        uw = DATA.Uo.*(DATA.g*x./DATA.Uo.^2).^(DATA.q);
        m0L2 = [PHI.L2]*Twind0(:,j,1); m1L2 = [PHI.L2]*Twind1(:,j,1); fbar = m1L2./m0L2;
        uwL2 = [PHI.L2]*u10(:,j);
        m0_L2_error = m0_L2_error + (1/xELEM(j).jacobian)*sum(w.*(Em0 - m0L2).^2);
        if strcmpi(f_error,'on')
            m1_L2_error = m1_L2_error + (1/xELEM(j).jacobian)*sum(w.*(fp - fbar).^2);
        else
            m1_L2_error = m1_L2_error + (1/xELEM(j).jacobian)*sum(w.*(fp.*Em0 - m1L2).^2);
        end
        uw_L2_error = uw_L2_error + (1/xELEM(j).jacobian)*sum(w.*(uw - uwL2).^2);
    end
    E_L2_error = sqrt(m0_L2_error)
    F_L2_error = sqrt(m1_L2_error)
    if strcmpi(windfig,'on')
        W_L2_error = sqrt(uw_L2_error)
    end
end