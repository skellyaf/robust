function [xfinal, stm, stt2, Q, Q2, xstmstt2QQ2_t, ti] = statestmsttQQ2Prop(xinitial, t, simparams, stm_initial, stt_initial, Q_initial, Q2_initial)

mu = simparams.mu;
nsv = simparams.nsv;
Qt = simparams.Qt;
options = simparams.options;
dynSys = simparams.dynSys;

if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') || strcmp(dynSys,'br4bp_sb1'))
    error('Error: Dynamical system identifier text not recognized!')
end


if nargin == 3
    stm_initial = eye(nsv);
    stt_initial = zeros(nsv,nsv,nsv);
    Q_initial = zeros(nsv,nsv);
    Q2_initial = zeros(nsv,nsv,nsv);
end
if nargin == 4
    stt_initial = zeros(nsv,nsv,nsv);
    Q_initial = zeros(nsv,nsv);
    Q2_initial = zeros(nsv,nsv,nsv);
end
if nargin == 5
    Q_initial = zeros(nsv,nsv);
    Q2_initial = zeros(nsv,nsv,nsv);
end
if nargin == 6
    Q2_initial = zeros(nsv,nsv,nsv);
end

if length(Q_initial(:))~=nsv^2
    assert(0);
end
xstmstt2QQ2 = [xinitial; stm_initial(:); stt_initial(:); Q_initial(:); Q2_initial(:)];
if t ~= 0
    if strcmp(dynSys,'2bp')
    % Two body numerical propagation of state and STM
        assert(0,'not set up yet!');
        [ti, xstmstt2QQ2_t] = ode113(@r2bp_stm_de, [0,t], xQQ2, options, mu);
    elseif strcmp(dynSys,'cr3bp')
    % Circular restricted three body numerical propagation of state and STM
        [ti, xstmstt2QQ2_t] = ode113(@statestmstt2QQ2dot, [0,t], xstmstt2QQ2, options, mu, Qt);
    elseif strcmp(dynSys,'br4bp_sb1')
        % Bicircular restricted four body problem, Sun-B1 frame
        [ti, xstmstt2QQ2_t] = ode113(@br4bp_sb1_state_stm_stt_QdQ_de, [0,t], xstmstt2QQ2, options, simparams.mub, mu, simparams.a4, simparams.theta_em_dot, Qt);

    end

    xfinal = reshape(xstmstt2QQ2_t(end,1:nsv),nsv,1);
    stm = reshape(xstmstt2QQ2_t(end,nsv+1:nsv+nsv*nsv),nsv,nsv);
    stt2 = reshape(xstmstt2QQ2_t(end,nsv+nsv*nsv+1:nsv+nsv^2+nsv^3),nsv,nsv,nsv);
    Q = reshape(xstmstt2QQ2_t(end,nsv+nsv^2+nsv^3+1:nsv+2*nsv^2+nsv^3),nsv,nsv);    
    
    Q2 = reshape(xstmstt2QQ2_t(end,nsv+2*nsv^2+nsv^3+1:nsv+2*nsv^2+2*nsv^3),nsv,nsv,nsv);
else
    Q = Q_initial;
    xfinal = xinitial;
    xstmstt2QQ2_t = xstmstt2QQ2;
end
end