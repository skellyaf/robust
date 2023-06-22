function [out] = plotIterationHistory(xhist,simparams)
%plotIterationHistory Plots an initial guess and all iterations of a
%trajectory to get to the final
%   Detailed explanation goes here

pp=figure;
hold on
axis equal
grid on
for i = 1:size(xhist,3)
    

    [~, xi_t, it, it_s] = createStateHistory(xhist(:,:,i), simparams);
    if i == 1


        p1=plot3(xi_t(:,1), xi_t(:,2), xi_t(:,3),'--','Color',[1 0 0 .5],'linewidth',5,'DisplayName','Initial Guess');
    elseif i==size(xhist,3)


%         plotMultiSegTraj(xhist(:,:,i), xi_t, it_s, simparams);


    
        p2=plot3(xi_t(:,1), xi_t(:,2), xi_t(:,3),'Color',[0 0 1 1],'linewidth',4,'DisplayName','Converged Trajectory');
    elseif i == 2
        p3=plot3(xi_t(:,1), xi_t(:,2), xi_t(:,3),'Color',[0 0 0 0.1],'linewidth',2,'DisplayName','Intermediate Steps');
    else
        plot3(xi_t(:,1), xi_t(:,2), xi_t(:,3),'Color',[0 0 0 0.1],'linewidth',2)
    end



end


out = 0;
mu = simparams.mu;

options = simparams.options;

x_target = simparams.x_target;
x_init = simparams.x_init;

dynSys = simparams.dynSys;

if strcmp(dynSys,'2bp')
    [~, x_1] = ode113(@r2bp_de, [0,simparams.T0], x_init, options, mu);
    [~, x_2] = ode113(@r2bp_de, [0,simparams.T_target], x_target, options, mu);
elseif strcmp(dynSys,'cr3bp')
    [~, x_1] = ode113(@cr3bp_sFrame_nd_de, [0,simparams.T0], x_init, options, mu);
    [~, x_2] = ode113(@cr3bp_sFrame_nd_de, [0,simparams.T_target], x_target, options, mu);
end

p4=plot3(x_1(:,1), x_1(:,2), x_1(:,3),'--','Color',[1 0 1 0.7],'LineWidth',1,'DisplayName','Initial Orbit');
p5=plot3(x_2(:,1), x_2(:,2), x_2(:,3),'--','Color',[0 1 0 0.7],'LineWidth',1,'DisplayName','Target Orbit');




ip1=plot3(simparams.x_init(1),simparams.x_init(2),simparams.x_init(3),'.','Color','Green','MarkerSize',25,'DisplayName','Initial Position');
tp1=plot3(simparams.x_target(1),simparams.x_target(2),simparams.x_target(3),'.','Color','Red','MarkerSize',25,'DisplayName','Target Position');






if strcmp(simparams.dynSys,'2bp')
    xlabel('X (km)')
    ylabel('Y (km)')
    zlabel('Z (km)')
    pp.Position=[100, 100, 700, 700];
    fontsize(gcf,14,"points")
    [az,el]=view
%     view([-20, 15]);
    view([-227.5, 20]); % 
    % save

    pp2 = copyobj(pp,0);
%     legend([p1, p2, p3, p4, p5, ip1, tp1], 'Location','NorthEast')
    legend([p1, p2, p3, p4, p5], 'Location','NorthEast')
    



    view([-180, 0]);

%     view([0, 35]);
    


%     [xunit,yunit,zunit]=sphere;
%     m = surf(xunit*6378,yunit*6378,zunit*6378);
%     colormap('gray');


elseif strcmp(simparams.dynSys,'cr3bp')
%     view([-75, 10]);
%     view([80, 12]); %leo-nri view

    view([0 90]); % planar
    xlabel('X (nd)')
    ylabel('Y (nd)')
    zlabel('Z (nd)')



    moon_radius=1737/384400;
    % Plot location of moon
    [xunit,yunit,zunit]=sphere;
    m = surf(xunit*moon_radius + (1-simparams.mu),yunit*moon_radius,zunit*moon_radius);
    colormap('gray');
end


end