
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Mission Script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear;                                %%% Reset all allocated parameters %%%
    pack;                                 %%% Perform garbage collection     %%%
%%% clc;                                  %%% Clear the command window       %%%

    format long g;
    format compact;
    fprintf ('\n');
    mu = 3.9860050e+14;
    re = 6.3781370e+06;

%%% Define initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt =    999;
    ti =    3.03438601492081e+03;
    ri = [ -3.47795506927413e+06; -5.46034046309521e+06; -1.76857863720122e+06 ];
    vi = [  5.02029898803858e+03; -1.35535016661676e+03; -5.68801451441240e+03 ];

%%% Exercise Kepler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    crafti.t = ti;
    crafti.r = ri;
    crafti.v = vi;

    [ craftf, phik ] = kepler (-mu, dt, crafti);

%%% Exercise Lambert %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tf = craftf.t;
    rf = craftf.r;
    vf = craftf.v;
    uh = cross (ri, rf);
    uh = uh / norm (uh);

    [ vi, vf, phil ] = lambert (-mu, dt, ri, rf, uh);

%%% Compute Lambert transition matrix from phik %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k11 = phik(1:3,1:3); k12 = phik(1:3,4:6);
    k21 = phik(4:6,1:3); k22 = phik(4:6,4:6);

    k12inv = inv (k12);

    l11 = -k12inv * k11;
    l12 = +k12inv;
    l21 = -transpose (k12inv);
    l22 = +transpose (k22 * k12inv);

    phi = [ l11, l12; l21, l22 ];

    for i = 1 : 6,
        fprintf ('      %13.5e%13.5e%13.5e |%13.5e%13.5e%13.5e\n', ...
              phi(i,1), phi(i,2), phi(i,3), phi(i,4), phi(i,5), phi(i,6) );
        if ((i == 3) | (i == 6)),
            fprintf ('      -------------------------------------');
            fprintf ('---|---------------------------------------');
            fprintf ('\n');
        end;
    end;

%%% Exit gracefully %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf ('Done\n\n');

    return;













