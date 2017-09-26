function u = slope_limiter(U,PHI1,PHI2,sELEM)

M = 0.0;
u = zeros(size(U));
nelems = length(U(1,:,1));

for j = 2:nelems-1
    dx   = 2/sELEM(j).jacobian;
    vh_L = PHI1*U(1:2,j);                        % Soln evaluated at left boundary
    vh_R = PHI2*U(1:2,j);                        % Soln evaluated at right boundary 
    a1   = vh_R-U(1,j);
    a2   = U(1,j)-U(1,j-1);
    a3   = U(1,j+1)-U(1,j); 
    M_R  = minmod_corrected(a1,a2,a3,M,dx);    % minmod due to Shu
    a1   = U(1,j)-vh_L;
    M_L  = minmod_corrected(a1,a2,a3,M,dx);
    uh_R = U(1,j) + M_R;                       % Slope limited soln at right boundary
    uh_L = U(1,j) - M_L;                       % Slope limited soln at left boundary
    check1 = abs(uh_R - vh_R);
    check2 = abs(uh_L - vh_L);
    
   % if abs(aof) < p+(1/2)
   %     u(:,j) = U(:,j);
   % else
        if check1 < 1e-6 && check2 < 1e-6        % No slope limiting required
            
            u(:,j) = U(:,j);
            
        else                                       % Limit soln
            
            u(1,j)     = U(1,j);
            a1         = U(2,j);
            u(2,j)     = minmod_corrected(a1,a2/(dx/2),a3/(dx/2),M,dx);
            u(3:end,j) = 0;
            
        end
   % end
    
end
u(:,1)      = U(:,1);
u(:,nelems) = U(:,nelems);
    

    