figure('Renderer','zbuffer');
set(gca,'NextPlot','replaceChildren'); hold on; 
set(gcf,'position',get(0,'screensize'))
set(gca,'FontSize',13,'FontWeight','demi')
box on

subplot(2,1,1)
for j = 1:nelems.x
    if p == 0
        xplot = [PSI.c]*X(xELEM(j).nodes)';
    else
        x = [PSI.xa]*X(xELEM(j).nodes)';
        xplot = [  X(j); x ; X(j+1)];
    end
    if strcmpi(analytic,'on')
        [fp,~,uw,Ee,nu]  = hasselman_solns(xplot,time,DATA,type);
        plot(xplot,Ee,'r:','LineWidth',5)
    end
    hold on
    if p == 0
        if strcmpi(wind,'on')
            uplot = [PHI.c]*Twind0(:,j,1);
        else
            uplot = [PHI.c]*m0(:,j,1);
        end
    else
        M0 = [PHI.b1]*m0(:,j,1); Mi = [PHI.L2]*m0(:,j,1); Mn = [PHI.b2]*m0(:,j,1);
        uplot = [M0; Mi; Mn];
    end
    
    if p == 0
        plot(xplot,uplot,'bo',...
            'LineWidth',1,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',7)
        plot(xplot,uplot,'b-')
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
ylabel('$m_0$ ($m^2$) ','FontSize',13,'FontWeight','demi','Interpreter','latex')%,'FontName','SansSerif'
legend(['  DG, p = ',num2str(p)],'Location','NorthWest')
title('Wind-swell interaction','FontSize',13,'FontWeight','demi','Interpreter','latex')
axis([x0, xN, -0.00001, 2.0])
%text(-0.90*umax,-4.5,['t = ',num2str(time/3600),' hrs'],'FontSize',13,'FontWeight','demi')
subplot(2,1,2)
for j = 1:nelems.x
    if p == 0
        xplot = [PSI.c]*X(xELEM(j).nodes)';
    else
        x = [PSI.xa]*X(xELEM(j).nodes)';
        xplot = [  X(j); x ; X(j+1)];
    end
    if strcmpi(analytic,'on')
        [fp,~,~,ep] = hasselman_solns(xplot,time,DATA,type);
        plot(xplot,fp.*ep,'r:','LineWidth',5)
    end
    
    hold on
    if p == 0
        if strcmpi(wind,'on')
            Eplot = [PHI.c]*Twind0(:,j,1);
            Fplot = [PHI.c]*Twind1(:,j,1);
        else
            Eplot = [PHI.c]*m0(:,j,1);
            Fplot = [PHI.c]*m1(:,j,1);
        end
        uplot = Fplot./Eplot;
    else
        Eplot = [ [PHI.b1]*m0(:,j,1); [PHI.L2]*m0(:,j,1); [PHI.b2]*m0(:,j,1) ];
        Fplot = [ [PHI.b1]*m1(:,j,1); [PHI.L2]*m1(:,j,1); [PHI.b2]*m1(:,j,1) ];
        uplot = Fplot;
    end
    
    if p == 0
        plot(xplot,uplot,'bo',...
            'LineWidth',1,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',7)
        plot(xplot,uplot,'b-')
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
ylabel('frequency ($1/s$)','FontSize',13,'FontWeight','demi','Interpreter','latex')
legend(['  DG, p = ',num2str(p)],'Location','NorthWest')
axis([x0, xN, -0.00001, 2.0])

frame = getframe(gcf);
writeVideo(writerObj,frame);
% plot_n = plot_n + floor(75/dt);
close