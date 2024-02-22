function plotMultiSegTraj(x, x_t, t_s, simparams, varargin)
%plotMultiSegTraj Plots a multi-segment trajectory
%   Doesn't open a figure or perform any of the hold commands
%   Just does the plotting when passed a multi-segment trajectory x
%   A different color for each segment, a dot for each leg intersection

options = simparams.options;
colorblind = simparams.colorblind;
% clf
hold on;
m = simparams.m;
n = simparams.n;
x_target = simparams.x_target;
x_init = simparams.x_init;
mu = simparams.mu;
x = reshape(x,simparams.m,simparams.n);

if nargin == 5
    tcm_idx = varargin{1};
end
    
% Plot the trajectory segments in different colors
for i = 1:n
    icolor = colorblind(mod(i,length(colorblind)-1)+1,:);
    r_plot = x_t(t_s==i, 1:3);    
    
   if size(r_plot,1) > 0
       plot3(r_plot(:,1), r_plot(:,2), r_plot(:,3),'LineWidth',3,'Color',icolor);
   end
       
end

% Plot the correction(s)
if exist('tcm_idx','var')
    for i = 1:length(tcm_idx)
        tcm_pos = x_t(tcm_idx(i),1:3);
        corrSeg = t_s(tcm_idx(i));
        seg_idxs = find(t_s==corrSeg);

        if tcm_idx(i) == seg_idxs(end)
            plot3(tcm_pos(1), tcm_pos(2), tcm_pos(3),'^','MarkerSize',12,'Color','Black','DisplayName','TCMr Location');
        else
            plot3(tcm_pos(1), tcm_pos(2), tcm_pos(3),'^','MarkerSize',12,'Color','Black','DisplayName','Correction Location');
        end
    end
end


%%% if the correction is co-located with a nominal maneuver
% if simparams.nom_dvctied && simparams.perform_correction
%     tcm_pos = x_t(t_corr==t,1:3);
%     cp1=plot3(tcm_pos(1), tcm_pos(2), tcm_pos(3),'o','MarkerSize',10,'Color','Black','DisplayName','Correction Location');
% %%%otherwise, if there is a nonzero tcm
% elseif tcm_3sigma > 0
%     tcm_pos = x_t(t_corr==t,1:3);
%     cp1=plot3(tcm_pos(1), tcm_pos(2), tcm_pos(3),'.','MarkerSize',25,'Color','Black','DisplayName','Correction Location');
% end




% Plot reference orbits
dynSys = simparams.dynSys;

if strcmp(dynSys,'2bp')
    [~, x_1] = ode113(@r2bp_de, [0,simparams.T0], x_init, options, mu);
    [~, x_2] = ode113(@r2bp_de, [0,simparams.T_target], x_target, options, mu);

    % Plot the initial orbit
    plot3(x_1(:,1), x_1(:,2), x_1(:,3),'Color','Black')
    % Plot the target orbit
    plot3(x_2(:,1), x_2(:,2), x_2(:,3),'Color','Black')
elseif strcmp(dynSys,'cr3bp')
    [~, x_1] = ode113(@cr3bp_sFrame_nd_de, [0,simparams.T0], x_init, options, mu);
    [~, x_2] = ode113(@cr3bp_sFrame_nd_de, [0,simparams.T_target], x_target, options, mu);

    % Plot the initial orbit
    plot3(x_1(:,1), x_1(:,2), x_1(:,3),'Color','Black')
    % Plot the target orbit
    plot3(x_2(:,1), x_2(:,2), x_2(:,3),'Color','Black')

elseif strcmp(dynSys,'br4bp_sb1')
    % Plot a circle for Earth orbit
    a4 = simparams.a4;
    r_b1 = [1-simparams.mub; 0; 0];
    viscircles(r_b1(1:2)', (1-mu)/a4, 'LineWidth',1,'Color','Black');

    % Plot a circle for lunar orbit
    viscircles(r_b1(1:2)', mu/a4, 'LineWidth',1,'Color','Black');

    


    % Plot a circular orbit at the initial earth position
    x_0 = x_t(1, :)';
    theta_0 = x_0(7);
    mub = simparams.mub;
    xe_0 = 1 - mub - 1/a4 * mu * cos(theta_0);
    ye_0 = -1/a4 * mu *sin(theta_0);
    re_0 = [xe_0; ye_0; 0];

    if existsAndTrue('fixed_initial_radius', simparams)
        viscircles(re_0(1:2)', simparams.fixed_initial_radius, 'LineWidth',1,'Color','Green');
    end    

    % Plot the initial earth position
    plot3(xe_0, ye_0, 0, '.', 'MarkerSize', 25, 'Color', 'Blue');

    % Plot circular orbit at the final moon position
    x_f = x_t(end, :)';
    theta_f = x_f(7);
    xm_f = 1 - mub + 1/a4 * (1-mu) * cos(theta_f);
    ym_f = 1/a4 * (1-mu) * sin(theta_f);
    rm_f = [xm_f; ym_f; 0];

    if existsAndTrue('fixed_final_radius', simparams)
        viscircles(rm_f(1:2)', simparams.fixed_final_radius, 'LineWidth',1,'Color','Red');
    end    

    % Plot the initial moon position
    plot3(xm_f, ym_f, 0, '.', 'MarkerSize', 25, 'Color', [.7 .7 .7]);
end


% Adjust the view
axis equal; grid on;

% view([-20, 15]);


if strcmp(dynSys,'2bp')
%     view([-20, 15]);
    xlabel('X (km)')
    ylabel('Y (km)')
    zlabel('Z (km)')
elseif strcmp(dynSys,'cr3bp')
%     view([-65, 15]); % nrho view
    view([80, 12]); %leo-nri view

%     view([0 90]); % planar

    % planar eed pef leo
%     xlim([-.75 1.1])
%     ylim([-.1 .9])

    
    xlabel('X (nd)')
    ylabel('Y (nd)')
    zlabel('Z (nd)')



    moon_radius=1737/384400;
    % Plot location of moon
    [xunit,yunit,zunit]=sphere;
    m = surf(xunit*moon_radius + (1-simparams.mu),yunit*moon_radius,zunit*moon_radius);
    colormap('gray');

elseif strcmp(dynSys,'br4bp_sb1')
    view([0 90])

end

% 
% annotation('textbox',[.85 .5 .1 .2], ...
%     'String',strcat('J:',num2str(round(J,4))),'EdgeColor','none')
% if exist('t_opt')
%     annotation('textbox',[.85 .45 .1 .2], ...
%     'String',strcat('t:',num2str(round(t_opt,2))),'EdgeColor','none')
% end






% Plot the initial position
if strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp')
    ip1=plot3(simparams.x_init(1),simparams.x_init(2),simparams.x_init(3),'.','Color','Green','MarkerSize',25,'DisplayName','Initial Position');

    % Plot the end of the trajectory (target after the final coast in most
    % instances)
    plot3(simparams.x_target(1),simparams.x_target(2),simparams.x_target(3),"square",'MarkerSize',20,'MarkerFaceColor','none','MarkerEdgeColor','Red')

end


% % Commenting out the red dot that gets plotted on the final TCM target
% if simparams.target_final_maneuver
%     finalManeuverSeg = simparams.maneuverSegments(end)-1;
%     r_finalMnvrSeg = x_t(t_s==finalManeuverSeg,1:3);
%     target = r_finalMnvrSeg(end,:);
% 
%     plot3(target(1),target(2),target(3),'.','Color','Red','MarkerSize',25,'DisplayName','Target Position');
% else
%     plot3(simparams.x_target(1),simparams.x_target(2),simparams.x_target(3),'.','Color','Red','MarkerSize',25,'DisplayName','Target Position');
% end



% Plot the nominal Delta V's
for i = 1:length(simparams.maneuverSegments)
    if simparams.maneuverSegments(i) == simparams.n + 1
        r_dv = x_t(end,:)';
    else
        r_dv = x(1:3,simparams.maneuverSegments(i));
    end
    plot3(r_dv(1),r_dv(2),r_dv(3),"pentagram",'MarkerSize',10,'MarkerFaceColor','Cyan','MarkerEdgeColor','Black');
end



if existsAndTrue('rdvz_flag',simparams)
    plot3(simparams.x0_target(1),simparams.x0_target(2),simparams.x0_target(3),"square",'MarkerSize',20,'MarkerFaceColor','none','MarkerEdgeColor','Green')
end





hold off;


end

