%----------------------
% spatial distribution nodes
%----------------------
Uxf = zeros(nnodes.x,1);
for j = 1:nnodes.x
    [~,~,Uxf(j)] = hasselman_solns(X(j),0,DATA,type);
end
Uxi = Uxf.*61/64;
%--------------------------------------------------------------------------
% Wind interpolation (Uo & q)
%--------------------------------------------------------------------------
options = optimoptions('fsolve');
options = optimoptions(options,'Display','off'); % turn off fsolve output.
Lw = 2; % number of observations
if strcmpi(wind_interp,'nodes')
     q  = zeros(Lw-1,nnodes.x); Uo = zeros(Lw-1,nnodes.x);
else
    qn  = zeros(nnodes.x);     Uon = zeros(nnodes.x);
end
for j = 1:nnodes.x
    for k = 1:Lw-1
        Ul = Uxi(j);    Ur = Uxf(j);
        if Ul < 1e-7
            Ul = 5e-6;
        end
        if Ur < 1e-7
            Ur = 5e-6;
        end
        x1 = [0.1,.5*(Ul+Ur)];
        x  = fsolve(@(x) root2d(x,dtw,Ur,Ul),x1,options);
        if strcmpi(wind_interp,'nodes')
            q(k,j)  = x(1);
            Uo(k,j) = x(2);
        else
             qn(j) = x(1);
            Uon(j) = x(2);
        end
    end
end