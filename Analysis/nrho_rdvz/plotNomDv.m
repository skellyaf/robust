function plotNomDv(x, traj, simparams)


% Plot the nominal Delta V's
for i = 1:length(simparams.maneuverSegments)
    if simparams.maneuverSegments(i) == simparams.n + 1
        r_dv = traj.x_t(end,:)';
    else
        r_dv = x(1:3,simparams.maneuverSegments(i));
    end
    plot3(r_dv(1),r_dv(2),r_dv(3),"pentagram",'MarkerSize',10,'MarkerFaceColor','Cyan','MarkerEdgeColor','Black');
end

end