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

if nargin == 5
    tcm_idx = varargin{1};
end

% Calculate total impulsive delta V 
% [stm_i, x_i_f, x_t, stm_t, t, t_s] = createStateStmHistory(x, simparams);


    

for i = 1:n
    icolor = colorblind(mod(i,length(colorblind)-1)+1,:);
    r_plot = x_t(t_s==i, 1:3);    
    
   if size(r_plot,1) > 0
       plot3(r_plot(:,1), r_plot(:,2), r_plot(:,3),'LineWidth',3,'Color',icolor);
   end
       
end

% Plot the correction(s)

if exist('tcm_idx')
    for i = 1:length(tcm_idx)
        tcm_pos = x_t(tcm_idx(i),1:3);
        corrSeg = t_s(tcm_idx(i));
        seg_idxs = find(t_s==corrSeg);

        if tcm_idx(i) == seg_idxs(end)
            plot3(tcm_pos(1), tcm_pos(2), tcm_pos(3),'o','MarkerSize',10,'Color','Black','DisplayName','TCMr Location');
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
elseif strcmp(dynSys,'cr3bp')
    [~, x_1] = ode113(@cr3bp_sFrame_nd_de, [0,simparams.T0], x_init, options, mu);
    [~, x_2] = ode113(@cr3bp_sFrame_nd_de, [0,simparams.T_target], x_target, options, mu);
end

plot3(x_1(:,1), x_1(:,2), x_1(:,3),'Color','Black')
plot3(x_2(:,1), x_2(:,2), x_2(:,3),'Color','Black')

axis equal; grid on;

view([-20, 15]);


if strcmp(dynSys,'2bp')
%     view([-20, 15]);
    xlabel('X (km)')
    ylabel('Y (km)')
    zlabel('Z (km)')
elseif strcmp(dynSys,'cr3bp')
%     view([-75, 10]);
%     view([80, 12]); %leo-nri view

    view([0 90]); % planar
    xlabel('X (nd)')
    ylabel('Y (nd)')
    zlabel('Z (nd)')
end

% 
% annotation('textbox',[.85 .5 .1 .2], ...
%     'String',strcat('J:',num2str(round(J,4))),'EdgeColor','none')
% if exist('t_opt')
%     annotation('textbox',[.85 .45 .1 .2], ...
%     'String',strcat('t:',num2str(round(t_opt,2))),'EdgeColor','none')
% end







ip1=plot3(simparams.x_init(1),simparams.x_init(2),simparams.x_init(3),'.','Color','Green','MarkerSize',25,'DisplayName','Initial Position');

if simparams.target_final_maneuver
    finalManeuverSeg = simparams.maneuverSegments(end)-1;
    r_finalMnvrSeg = x_t(t_s==finalManeuverSeg,1:3);
    target = r_finalMnvrSeg(end,:);

    plot3(target(1),target(2),target(3),'.','Color','Red','MarkerSize',25,'DisplayName','Target Position');
else
    plot3(simparams.x_target(1),simparams.x_target(2),simparams.x_target(3),'.','Color','Red','MarkerSize',25,'DisplayName','Target Position');
end






% if exist('cp1')
%     legend([ip1 tp1 cp1],'Location','NorthEast')
% else
%     legend([ip1 tp1],'Location','NorthEast')
% end


hold off;


end

