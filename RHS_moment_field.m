%--------------------------------------------------------------------------
% Loop over elements
%--------------------------------------------------------------------------
g = DATA.g;

if tw == 0
    tw = 5e-6;
end

for j = 1:nelems.x
    if strcmpi(test,'spatial')
        x = [PSI.xa]*X(xELEM(j).nodes)';             % x points for source winds
        [fxj,~,uw] = hasselman_solns(x,t,DATA,type); % wind evaluated at L2 pts
        S0j = So(x);
        S1j = S1(x);
    elseif strcmpi(test,'temporal')
        x = [PSI.xa]*X(xELEM(j).nodes)';             % x points for equilibrium
        [fxj,~,uw] = hasselman_solns(x,t,DATA,type); % kernel group freq
        cgj = g./(2*pi*fxj);
        kxj = (2*pi*fxj).^2./DATA.g;
        S0c = So(x); % to determine constants
        S1c = S1(x);
        if strcmpi(wind_interp,'elem')
            S0jt = zeros(size(fxj)); S1jt = zeros(size(fxj));
            for i = 1:L2_pts
                S0jt(i)   = Sot(tw,Uo(i,j),q(i,j),DATA.g);
                S1jt(i)   = S1t(tw,Uo(i,j),q(i,j),DATA.g);
            end
            S0j = cgj.*S0jt;
            S1j = (1/2*(cgj.*kxj./(2*pi*fxj)).^2).*S1jt;
            C0s  = S0c./S0j; C1s = S1c./S1j; % calculate constants
            S0j  = C0s.*S0j; S1j = C1s.*S1j; % multiply by constant
            if n == NT    % record constants
                if j == 1
                    Ck0 = zeros(L2_pts,nelems.x); Ck1 = zeros(L2_pts,nelems.x);
                end
                Ck0(:,j) = C0s; Ck1(:,j) = C1s;
            end
        elseif strcmpi(wind_interp,'nodes')
            S0i  = Sot(tw,Uo(1,j),q(1,j),DATA.g);
            S0i1 = Sot(tw,Uo(1,j+1),q(1,j+1),DATA.g);
            S0j  = [PSI.xa]*[S0i; S0i1];
            S0j  = cgj.*S0j;
            S1i  = S1t(tw,Uo(1,j),q(1,j),DATA.g);
            S1i1 = S1t(tw,Uo(1,j+1),q(1,j+1),DATA.g);
            S1j  = [PSI.xa]*[S1i; S1i1];
            S1j  = 1/2*(cgj.*kxj./(2*pi*fxj)).^2.*S1j;
            C0s  = S0c./S0j; C1s = S1c./S1j; % calculate constants
            S0j  = C0s.*S0j; S1j = C1s.*S1j; % multiply by constant
            if n == NT    % record constants
                if j == 1
                    Ck0 = zeros(L2_pts,nelems.x); Ck1 = zeros(L2_pts,nelems.x);
                end
                Ck0(:,j) = C0s; Ck1(:,j) = C1s;
            end
        end
    end
    m0j  = [PHI.elem]*m0(:,j,irk); w0j = [PHI.elem]*Twind0(:,j,irk);
    m1j  = [PHI.elem]*m1(:,j,irk); w1j = [PHI.elem]*Twind1(:,j,irk);
    % element flux (currently deep water)
    Jfb        = find(m0j < 1e-7); Ifb       = find(w0j < 1e-7);
    fbar       = m1j./m0j;         wbar      = w1j./w0j; 
    fbar(Jfb)  = 0;                wbar(Ifb) = 0;
    CGbar      = Cg(fbar);         cgwm      = Cg(wbar);
    CGbar(Jfb) = 0;                cgwm(Ifb) = 0;
    F0j = CGbar.*m0j; F1j = CGbar.*m1j;  % flux across elements
    W0j = cgwm.*w0j;  W1j = cgwm.*w1j;
    RHS_m0(:,j,irk) = xELEM(j).jacobian*A*(F0j);
    RHS_m1(:,j,irk) = xELEM(j).jacobian*A*(F1j);
    wind0(:,j,irk)  = xELEM(j).jacobian*A*(W0j);
    wind1(:,j,irk)  = xELEM(j).jacobian*A*(W1j);
    if strcmpi(wind,'on')
        RHS_m0(:,j,irk) = RHS_m0(:,j,irk) + C*S0j;
        RHS_m1(:,j,irk) = RHS_m1(:,j,irk) + C*S1j;
    elseif strcmpi(windsea,'on')
        RHS_m0(:,j,irk) = RHS_m0(:,j,irk) + C*((S0j));
        wind0(:,j,irk)  = wind0(:,j,irk)  + C*((S0j));
        RHS_m1(:,j,irk) = RHS_m1(:,j,irk) + C*((S1j));
        wind1(:,j,irk)  = wind1(:,j,irk)  + C*((S1j));
    elseif strcmpi(wind_swell,'on')
        if t > relax_time % transfer swell to wind-sea
            x   = [PSI.xa]*X(xELEM(j).nodes)';
            fp  = hasselman_solns(x,time,DATA,type);
            uj  = [PHI.L2]*u10(:,j);
            [B0,B1]=swell_wind(fp,0,0,uj);
            m0s = [PHI.L2]*m0(:,j,irk); m1s = [PHI.L2]*m1(:,j,irk);
            S0s = B0.*m0s;
            S1s = B1.*m0s;
            fc  = g*(0.13)./uj;
            fms = m1s./m0s; fipt = find(fms < fc);
            S0s(fipt) = 0.0; S1s(fipt) = 0.0;
        else  % no transer (b/c no swell present)
            fipt = [];
            S0s  = zeros(size(S0j)); S1s = S0s;
        end
        RHS_m0(:,j,irk) = RHS_m0(:,j,irk) + C*S0s;
        RHS_m1(:,j,irk) = RHS_m1(:,j,irk) + C*S1s;
        wind0(:,j,irk)  = wind0(:,j,irk)  + C*S0j; 
        wind1(:,j,irk)  = wind1(:,j,irk)  + C*S1j; 
    end
end
%--------------------------------------------------------------------------
% Loop over boundaries for parameters
%--------------------------------------------------------------------------
if strcmpi(wind,'on')
    [fL,~,~,m0L] = hasselman_solns(x0,tw,DATA,type);
    m1L = fL*m0L;
    CgL = Cg(fL);
    Fl = CgL*m0L;
elseif strcmpi(swell,'on')
    fL  = 0.05;
    Hs  = 1.00*cos(2*pi*fL*t/60);
    m0L = (1/4*Hs)^2;
    m1L = fL*m0L;
    CgL = Cg(fL);
    if t > relax_time
        Fl = CgL*m0L;
    else
        Fl = 0;
    end
    [fLw,~,~,m0Lw] = hasselman_solns(x0,tw,DATA,type);
    m1Lw = fLw*m0Lw;
    CgLw = Cg(fLw); Flw = CgLw*m0Lw;
elseif strcmpi(windsea,'on')
    [fLw,~,~,m0Lw] = hasselman_solns(x0,tw,DATA,type);
    m1Lw = fLw*m0Lw;
    CgLw = Cg(fLw); Flw = CgLw*m0Lw;
    fL  = 0.10;
    tl  = t - relax_time;
    Hs  = 0.250*cos(2*pi*fL*tl/60);
    m0L = (1/4*Hs)^2;
    m1L = fL*m0L;
    CgL = Cg(fL);
    if t > relax_time
        Fl = CgL*m0L + Flw;
    else
        Fl = 0;
    end
else
    Fl = 0; 
end
Fhat = Fl;
RHS_m0(:,1,irk) = RHS_m0(:,1,irk) + xELEM(1).jacobian*B.b1*Fhat; 
wind0(:,1,irk)  = wind0(:,1,irk)  + xELEM(1).jacobian*B.b1*(Flw);

if strcmpi(wind,'on')
    Fl = CgL*m1L;
elseif strcmpi(swell,'on')
    fL  = 0.05;
    tl  = t - relax_time;
    Hs  = 1.0*cos(2*pi*fL*tl/60);
    m0L = (1/4*Hs)^2;
    m1L = fL*m0L;
    CgL = Cg(fL);
    if t > relax_time
        Fl = CgL*m1L;
    else
        Fl = 0;
    end
elseif strcmpi(windsea,'on')
    if t > relax_time
        Fl = CgL*m1L + CgLw*m1Lw;
    else
        Fl = CgLw*m1Lw;
    end
else
    Fl = 0;
end
Fhat = Fl;
RHS_m1(:,1,irk) = RHS_m1(:,1,irk) + xELEM(1).jacobian*B.b1*Fhat; 
wind1(:,1,irk)  = wind1(:,1,irk)  + xELEM(1).jacobian*B.b1*(CgLw*m1Lw);
% Loop over nodes
for i = 2:nnodes.x-1
    jL  = xNODE(i).elems(1);          jR  = xNODE(i).elems(2); 
    m0L = [PHI.b2]*m0(:,jL,irk);     m0R  = [PHI.b1]*m0(:,jR,irk);
    m1L = [PHI.b2]*m1(:,jL,irk);     m1R  = [PHI.b1]*m1(:,jR,irk);
    w0L = [PHI.b2]*Twind0(:,jL,irk); w0R  = [PHI.b1]*Twind0(:,jR,irk);
    w1L = [PHI.b2]*Twind1(:,jL,irk); w1R  = [PHI.b1]*Twind1(:,jR,irk);
    
    ckL1 = isnan(m1L);    ckL3 = isnan(w1L);
    ckL2 = max(m0L);      ckL4 = max(w0L); 
    
    if max(ckL1) > 0; display('soln is NaN'); keyboard; end
    if ckL2 > 100; display('soln is blowing up'); keyboard; end
    if max(ckL3) > 0; display('soln is NaN'); keyboard; end
    if ckL4 > 100; display('soln is blowing up'); keyboard; end
    
    
    if m0L < 1e-7; fL  = 0; CgL = 0; else fL  = m1L/m0L; CgL = Cg(fL); end
    
    if m0R < 1e-7; fR  = 0; CgR = 0; else fR  = m1R/m0R; CgR = Cg(fR); end
    
    if w0L < 1e-7; wL  = 0; cgwL = 0; else wL  = w1L/w0L; cgwL = Cg(wL); end
    
    if w0R < 1e-7; wR  = 0; cgwR = 0; else wR  = w1R/w0R; cgwR = Cg(wR); end


    Fl  = CgL*m0L;                    Fr  = CgR*m0R;
    Flw = cgwL*w0L;                   Frw = cgwR*w0R;
    if strcmpi(flux_type,'LLF')
        Fhat = 1/2*((Fl+Fr) - max(abs(CgL),abs(CgR))*(m0R - m0L));
        What = 1/2*((Flw+Frw) - max(abs(cgwL),abs(cgwR))*(w0R - w0L));
    elseif strcmpi(flux_type,'upwind')
        Fhat = Fl; What = Flw;
    end
    RHS_m0(:,jL,irk) = RHS_m0(:,jL,irk) + xELEM(jL).jacobian*B.b2*Fhat;
    RHS_m0(:,jR,irk) = RHS_m0(:,jR,irk) + xELEM(jR).jacobian*B.b1*Fhat;
    wind0(:,jL,irk)  = wind0(:,jL,irk)  + xELEM(jL).jacobian*B.b2*What;
    wind0(:,jR,irk)  = wind0(:,jR,irk)  + xELEM(jR).jacobian*B.b1*What;
    
    Fl = CgL*m1L;                      Fr = CgR*m1R;
    Flw = cgwL*w1L;                   Frw = cgwR*w1R;
    if strcmpi(flux_type,'LLF')
        Fhat = 1/2*((Fl+Fr) - max(abs(CgL),abs(CgR))*(m1R - m1L));
        What = 1/2*((Flw+Frw) - max(abs(cgwL),abs(cgwR))*(w1R - w1L));
    elseif strcmpi(flux_type,'upwind')
        Fhat = Fl; What = Flw;
    end
    RHS_m1(:,jL,irk) = RHS_m1(:,jL,irk) + xELEM(jL).jacobian*B.b2*Fhat;
    RHS_m1(:,jR,irk) = RHS_m1(:,jR,irk) + xELEM(jR).jacobian*B.b1*Fhat;
    wind1(:,jL,irk)  = wind1(:,jL,irk)  + xELEM(jL).jacobian*B.b2*What;
    wind1(:,jR,irk)  = wind1(:,jR,irk)  + xELEM(jR).jacobian*B.b1*What;    
end
% Right bc
m0L = [PHI.b2]*m0(:,nelems.x,irk);     m1L = [PHI.b2]*m1(:,nelems.x,irk);
w0L = [PHI.b2]*Twind0(:,nelems.x,irk); w1L = [PHI.b2]*Twind1(:,nelems.x,irk);
fL  = m1L/m0L;  wL = w1L/w0L;

CgL = Cg(fL); cgwL = Cg(wL);
Fl  = CgL*m0L; Flw = cgwL*w0L; 

Fhat = Fl;  What = Flw;
RHS_m0(:,nelems.x,irk) = RHS_m0(:,nelems.x,irk) + xELEM(nelems.x).jacobian*B.b2*Fhat;
wind0(:,nelems.x,irk)  = wind0(:,nelems.x,irk)  + xELEM(nelems.x).jacobian*B.b2*What;

Fl = CgL*m1L; Flw = cgwL*w1L;
Fhat = Fl; What = Flw;
RHS_m1(:,nelems.x,irk) = RHS_m1(:,nelems.x,irk) + xELEM(nelems.x).jacobian*B.b2*Fhat;
wind1(:,nelems.x,irk)  = wind1(:,nelems.x,irk)  + xELEM(nelems.x).jacobian*B.b2*What;