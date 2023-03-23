
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Mission Script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear;                                %%% Reset all allocated parameters %%%
    pack;                                 %%% Perform garbage collection     %%%
%%% clc;                                  %%% Clear the command window       %%%

    global fput;

    true     = 1;
    false    = 0;
    degtorad = pi / 180;
    radtodeg = 180 / pi;
    mu  = 3.9860050e+14;
    re  = 6.3781370e+06;

    format long g;
    format compact;
    fprintf ('\n');
%%% fput = fopen ('main.out', 'w');

%%% Initialize craft and extrapolate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    craft.t =    3.03438601492081e+03;
    craft.r = [ -3.47795506927413e+06; -5.46034046309521e+06; -1.76857863720122e+06 ];
    craft.v = [  5.02029898803858e+03; -1.35535016661676e+03; -5.68801451441240e+03 ];
    craft.t =    3.03438601492081e+03;
    craft.r = [ 7000e3; 500e3; 100e3 ];
    craft.v = [ 500 ; 12000; 100 ];

    [ craftf, phi ] = kepler (-mu, 99999, craft);

%%% Exit gracefully %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf ('Done\n\n');

    return;













