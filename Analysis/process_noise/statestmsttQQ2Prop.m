function [xfinal, stm, stt2, Q, Q2, xstmstt2QQ2_t, ti] = statestmsttQQ2Prop(xinitial, t, simparams, stm_initial, stt_initial, Q_initial, Q2_initial)

mu = simparams.mu;
Qt = simparams.Qt;
options = simparams.options;
dynSys = simparams.dynSys;

if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') )
    error('Error: Dynamical system identifier text not recognized!')
end


if nargin == 3
    stm_initial = eye(6);
    stt_initial = zeros(6,6,6);
    Q_initial = zeros(6,6);
    Q2_initial = zeros(6,6,6);
end
if nargin == 4
    stt_initial = zeros(6,6,6);
    Q_initial = zeros(6,6);
    Q2_initial = zeros(6,6,6);
end
if nargin == 5
    Q_initial = zeros(6,6);
    Q2_initial = zeros(6,6,6);
end
if nargin == 6
    Q2_initial = zeros(6,6,6);
end

if length(Q_initial(:))~=36
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
    end

    xfinal = reshape(xstmstt2QQ2_t(end,1:6),6,1);
    stm = reshape(xstmstt2QQ2_t(end,7:42),6,6);
    stt2 = reshape(xstmstt2QQ2_t(end,43:258),6,6,6);
    Q = reshape(xstmstt2QQ2_t(end,259:294),6,6);    
    
    Q2 = reshape(xstmstt2QQ2_t(end,295:510),6,6,6);
else
    Q = Q_initial;
    xfinal = xinitial;
    xstmstt2QQ2_t = xstmstt2QQ2;
end
end