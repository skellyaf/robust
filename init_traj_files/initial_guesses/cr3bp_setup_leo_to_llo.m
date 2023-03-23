%% CR3BP leo to llo setup


% Dimensionless variable setup





%% Generating low lunar orbit state
a = moon.rad + 110; % Lunar altitude of 110km
e = .001; % Slightly eccentric
% Inclination of 27 degrees (27, 50, 76, & 86 enable extended LLO stays)
inc = (180 - 27) * pi/180;
% Right ascension - 180 degrees
W = 180 * pi/180;
% Arg of perigee - zero degrees
w = 0;
% True anomaly - 180 degrees;
nu = 180 * pi/180;

% 2 body COEs for LLO
[r_2b_llo,v_2b_llo] = orbel2rv(a,e,inc,W,w,nu,moon.mu);
r_2b_llo_nd = r_2b_llo / ndDist2km;
v_2b_llo_nd = v_2b_llo / ndVel2kms;

% Position of moon in nondimensional CR3BP is (1-mu)
% Velocity - all Y

r_3b_llo = r_2b_llo_nd + [1-mu; 0; 0];

X_target = [r_3b_llo; v_2b_llo_nd];

% Test propagate the LLO initial conditions for 2 hours
[~,x_llo_t] = stateProp(X_target, 2 / ndTime2hrs, mu,'cr3bp');


% figure
% plot3(x_llo_t(:,1),x_llo_t(:,2),x_llo_t(:,3))
% xlabel('X (nd)')
% ylabel('Y (nd)')
% zlabel('Z (nd)')
% title('Position in Earth-centered Earth-Moon rotating frame')
% axis equal
% grid on
% hold on
% plot3(X_target(1), X_target(2), X_target(3), '.','MarkerSize',25)

%% Generating LEO state

%% LEO Options
% leo_poincare_points = [];
% leo_initial_states = [];
% leo_integration_times = [];

a = earth.rad + 450;
e = .0001;
% inc = 5.14 * pi/180;
inc = .1 * pi/180;
W = 180 * pi/180;
w = 0;
nu = .3:.1:pi/2;
dv = 3:.1:3.5;


% for j = 1:length(nu)
%     
% 
%     [r_2b_leo,v_2b_leo] = orbel2rv(a, e, inc, W, w, nu(j), earth.mu);
%     r_2b_leo_nd = r_2b_leo / ndDist2km;
%     r_3b_leo = r_2b_leo_nd - [mu; 0; 0];
%     v_2b_leo_nd = v_2b_leo / ndVel2kms;
%     X_leo = [r_3b_leo; v_2b_leo_nd];
%     
% 
%     
%     
%     
%     
%     for i = 1:length(dv)
%         i_leo_v = X_leo(4:6) / norm(X_leo(4:6));
%         v_leo_nd = X_leo(4:6) + dv(i) * i_leo_v;
%         X_xfer = [X_leo(1:3); v_leo_nd];
%     
%         [~,x_xfer_t] = stateProp(X_xfer, 3.5 / ndTime2days, mu,'cr3bp');
%         plot3(x_xfer_t(:,1),x_xfer_t(:,2),x_xfer_t(:,3))
%     
%     
% 
%     
%     end
% 
% end


% Best solution corresponds to a little less than dv(1) and between nu(6) and nu(7)


% figure
% plot3(x_llo_t(:,1),x_llo_t(:,2),x_llo_t(:,3))
% xlabel('X (nd)')
% ylabel('Y (nd)')
% zlabel('Z (nd)')
% title('Position in Earth-centered Earth-Moon rotating frame')
% axis equal
% grid on
% hold on
% plot3(X_target(1), X_target(2), X_target(3), '.','MarkerSize',25)

% nu_sol = (nu(6) + nu(7))/2;
nu_sol = nu(7);
dv_sol = dv(1);


% Re-prop plot with the best guess
[r_2b_leo,v_2b_leo] = orbel2rv(a, e, inc, W, w, nu_sol, earth.mu);
r_2b_leo_nd = r_2b_leo / ndDist2km;
r_3b_leo = r_2b_leo_nd - [mu; 0; 0];
v_2b_leo_nd = v_2b_leo / ndVel2kms;
X_leo = [r_3b_leo; v_2b_leo_nd];

i_leo_v = X_leo(4:6) / norm(X_leo(4:6));
v_leo_nd = X_leo(4:6) + dv_sol * i_leo_v;
X_xfer = [X_leo(1:3); v_leo_nd];

[~,x_xfer_t, t] = stateProp(X_xfer, 3.5 / ndTime2days, mu,'cr3bp');
% plot3(x_xfer_t(:,1),x_xfer_t(:,2),x_xfer_t(:,3))
% hold on




%% Break guess up into n segments
X = [];
X = [x_xfer_t(1,:)'];
n = 10;

idx_sep = floor( length(t) / n );

iter = 1;
i = idx_sep;
X(7,1) = t(i);

while iter

    
    
    if i > length(t) - idx_sep
        X(7,end) = t(end) - t(i - idx_sep);
        iter = 0;
    else
        X(1:6,end+1) = x_xfer_t(i,:)';
        X(7,end) = t(i + idx_sep) - t(i);
    end
    i = i + idx_sep;   
    

end

% Remove last segment - goes too far
X = X(:,1:end-1);

% while iter
% 
%     X(1:6,end+1) = x_xfer_t(i,:)'
%     i = i + idx_sep;   
%     
%     if i > length(t)
%         X(7,end) = t(end) - t(i - idx_sep);
%         iter = 0;
%     else
%         X(7,end) = t(i) - t(i - idx_sep);
%     end
% 
% end

% figure
% 
% plotMultiSegTrajDtC(X(:),mu,'cr3bp')
% axis equal
% grid on
% hold on

% for i = 1:size(X,2)
%     [~,x_plot] = stateProp(X(1:6,i), X(7,i), mu,'cr3bp');
%     plot3(x_plot(:,1),x_plot(:,2),x_plot(:,3),'LineWidth',2)
% 
% 
% 
% end


% plot3(x_xfer_t(:,1),x_xfer_t(:,2),x_xfer_t(:,3),'--','LineWidth',2)




%% Add coasts in initial and target orbit

% backward prop LEO from the departure point

[X_leo_0] = stateProp(X_leo, -.25 / ndTime2hrs, mu,'cr3bp');

X = [ [X_leo_0; .25 / ndTime2hrs], X];

% Add LLO target state and time

X = [X, [X_target; .4 / ndTime2hrs]];
% Propagate forward to get actual target in LLO
X_target = stateProp(X_target, .4/ndTime2hrs,mu,'cr3bp');

