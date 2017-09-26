%----------------------
% spatial distribution nodes
%----------------------
scaleF = 61/64;
Uxf = zeros(L2_pts,nelems.x);
Uxi = zeros(L2_pts,nelems.x);
for j = 1:nelems.x
    x  = [PSI.xa]*X(xELEM(j).nodes)';
    [~,~,Uxf(:,j)] = hasselman_solns(x,0,DATA,type);
   Uxi(:,j) = Uxf(:,j).*5e-3; 
end
%--------------------------------------------------------------------------
% Wind interpolation (Uo & q)
%--------------------------------------------------------------------------
options = optimoptions('fsolve');
options = optimoptions(options,'Display','off'); % turn off fsolve output.
Lw = 2; % number of observations
q  = zeros(L2_pts,nelems.x); Uo = zeros(L2_pts,nelems.x);
for i = 1:nelems.x
    for j = 1:L2_pts
        for k = 1:Lw-1
            Ul = Uxi(j,i);    Ur = Uxf(j,i);
            if Ul < 1e-7
                Ul = 5e-6;
            end
            if Ur < 1e-7
                Ur = 5e-6;
            end
            x1 = [0.1,.5*(Ul+Ur)];
            x  = fsolve(@(x) root2d(x,dtw,Ur,Ul),x1,options);
            q(j,i)  = x(1);
            Uo(j,i) = x(2);
        end
    end
end
