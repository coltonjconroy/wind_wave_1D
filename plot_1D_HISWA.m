figure
subplot(2,1,1)
eplot = zeros(size(nelems.x));
for j = 1:nelems.x
    if p == 0
        x = [PSI.c]*X(xELEM(j).nodes)';
        xplot = x;
        xexact = [PSI.xa]*X(xELEM(j).nodes)';
    else
        x = [PSI.xa]*X(xELEM(j).nodes)';
        xplot = [  X(j); x ; X(j+1)];
    end
    if strcmpi(wind,'on')
        if relax_time >= T
            if (p == 0)
                [fp,~,uw,Ee,nu]  = hasselman_solns(xexact,time,DATA,type);
                plot(xexact,Ee,'r:','LineWidth',5)
            else
                [fp,~,uw,Ee,nu]  = hasselman_solns(xplot,time,DATA,type);
                plot(xplot,Ee,'r:','LineWidth',5)
            end
        end
    end
    hold on
    if strcmpi(wind,'on')
        if p == 0
            uplot = [PHI.c]*Twind0(:,j,1);
            eplot(j) = uplot;
        else
            M0 = [PHI.b1]*Twind0(:,j,1); Mi = [PHI.L2]*Twind0(:,j,1); Mn = [PHI.b2]*Twind0(:,j,1);
            uplot = [M0; Mi; Mn];
        end
    else
        if p == 0
            uplot = [PHI.c]*m0(:,j,1);
            eplot(j) = uplot;
        else
            M0 = [PHI.b1]*m0(:,j,1); Mi = [PHI.L2]*m0(:,j,1); Mn = [PHI.b2]*m0(:,j,1);
            uplot = [M0; Mi; Mn];
        end
    end
    
    if p == 0
        plot(xplot,uplot,'bo',...
            'LineWidth',1,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',7)
    else
        plot([xplot(1),xplot(end)],[uplot(1),uplot(end)],'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',9)
        plot(xplot,uplot,'b-')
    end

end
xlabel(['fetch, $\chi$ (m)'],'FontSize',13,'FontWeight','demi','Interpreter','latex');
ylabel('$m_0$ ($m^2$)','FontSize',13,'FontWeight','demi','Interpreter','latex')%,'FontName','SansSerif'
if strcmpi(wind,'on')
    if strcmpi(wind_swell,'on')
        legend(['  DG, p = ',num2str(p)],'Location','NorthWest')
        title('Wind-sea-swell interacion test case','FontSize',13,'FontWeight','demi','Interpreter','latex')
    else
        legend('  Analytic Solution',['  DG, p = ',num2str(p)],'Location','NorthWest')
    end
else
    legend(['  DG, p = ',num2str(p)],'Location','NorthEast')
    title('Wind-sea-swell interacion test case','FontSize',13,'FontWeight','demi','Interpreter','latex')
end
subplot(2,1,2)
eplot = zeros(size(nelems.x)); fplot = zeros(size(nelems.x)); Uplot = zeros(size(nelems.x));
for j = 1:nelems.x
    if p == 0
        x = [PSI.c]*X(xELEM(j).nodes)';
        xplot = x;
        xexact = [PSI.xa]*X(xELEM(j).nodes)';
    else
        x = [PSI.xa]*X(xELEM(j).nodes)';
        xplot = [  X(j); x ; X(j+1)];
    end
    if strcmpi(wind,'on')
        if relax_time >= T
            if (p == 0)
                [fp,~,uw,Ee,nu]  = hasselman_solns(xexact,time,DATA,type);
                if strcmpi(f_error,'on')
                    plot(xexact,fp,'r:','LineWidth',5)
                else
                    plot(xexact,fp.*Ee,'r:','LineWidth',5)
                end
            else
                [fp,~,uw,Ee,nu]  = hasselman_solns(xplot,time,DATA,type);
                if strcmpi(f_error,'on')
                    plot(xplot,fp,'r:','LineWidth',5)
                else
                    plot(xplot,fp.*Ee,'r:','LineWidth',5)
                end
            end
        end
    end
    
    hold on

    if strcmpi(wind,'on')
        if p == 0
            Fplot = [PHI.c]*Twind1(:,j,1);
            Eplot = [PHI.c]*Twind0(:,j,1);
            eplot(j) = Eplot;
            fplot(j) = Fplot;
        else
            Fplot = [ [PHI.b1]*Twind1(:,j,1); [PHI.L2]*Twind1(:,j,1); [PHI.b2]*Twind1(:,j,1) ];
            Eplot = [ [PHI.b1]*Twind0(:,j,1); [PHI.L2]*Twind0(:,j,1); [PHI.b2]*Twind0(:,j,1) ];
        end
    else
        if p == 0
            Fplot = [PHI.c]*m1(:,j,1);
            Eplot = [PHI.c]*m0(:,j,1);
            eplot(j) = Eplot;
            fplot(j) = Fplot;
        else
            Fplot = [ [PHI.b1]*m1(:,j,1); [PHI.L2]*m1(:,j,1); [PHI.b2]*m1(:,j,1) ];
            Eplot = [ [PHI.b1]*m0(:,j,1); [PHI.L2]*m0(:,j,1); [PHI.b2]*m0(:,j,1) ];
        end
    end
    inan = find(Eplot < 1e-8);
    if strcmpi(f_error,'on')
        uplot = Fplot./Eplot;
    else
        uplot = Fplot;
    end
    uplot(inan) = 0;
    if p == 0
          plot(xplot,uplot,'bo',...
        'LineWidth',1,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',[.49 1 .63],...
        'MarkerSize',7)
    else
        plot([xplot(1),xplot(end)],[uplot(1),uplot(end)],'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',9)
        plot(xplot,uplot,'b-')
    end

end
xlabel(['fetch, $\chi$ (m)'],'FontSize',13,'FontWeight','demi','Interpreter','latex');
ylabel('frequency~(1/s)','FontSize',13,'FontWeight','demi','Interpreter','latex')
if strcmpi(wind,'on')
    if strcmpi(wind_swell,'on')
        legend(['  DG, p = ',num2str(p)],'Location','SouthEast')
    else
        legend('  Analytic Solution',['  DG, p = ',num2str(p)],'Location','NorthEast')
    end
else
    legend(['  DG, p = ',num2str(p)],'Location','NorthEast')
end

if strcmpi(windfig,'on')
    figure
    for j = 1:nelems.x
        
        x = [PSI.xa]*X(xELEM(j).nodes)';
        xplot = [  X(j); x ; X(j+1)];
        plot(x,DATA.Uo.*(DATA.g*x./DATA.Uo.^2).^(DATA.q),'r:','LineWidth',5)
        
        hold on
        uplot = [ [PHI.b1]*u10(:,j); ...
            [PHI.L2]*u10(:,j); ...
            [PHI.b2]*u10(:,j)];
        
        plot([xplot(1),xplot(end)],[uplot(1),uplot(end)],'bo',...
            'LineWidth',1,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',7)
        plot(xplot,uplot,'b-')
        
    end
    xlabel(['fetch, $\chi$ (m)'],'FontSize',13,'FontWeight','demi','Interpreter','latex');
    ylabel('wind speed (m/s)','FontSize',13,'FontWeight','demi','Interpreter','latex')
    legend('  fetch dependent',['  duration dependent, p = ',num2str(p)],'Location','NorthWest')
    title('Hasselman et. al. wind test case','FontSize',13,'FontWeight','demi','Interpreter','latex')
end
if strcmpi(C_plot,'on')
    figure
    for j = 1:nelems.x
        x = [PSI.xa]*X(xELEM(j).nodes)';
        plot(x,Ck0(:,j),'b')
        hold on
        plot(x,Ck1(:,j),'r')
    end
    xlabel(['fetch, $\chi$ (m)'],'FontSize',13,'FontWeight','demi','Interpreter','latex');
    ylabel('constants','FontSize',13,'FontWeight','demi','Interpreter','latex')
    legend('c_0','c_1','Location','Best')
    title('Source term constants for Hasselman et.~al.~wind test case','FontSize',13,'FontWeight','demi','Interpreter','latex')
end

