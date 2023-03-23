
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Title     : Kepler                                                       %%%
%%% Purpose   : To propagate an inertial state vector through a conic        %%%
%%%             gravity field by a specified transfer time interval. This    %%%
%%%             routine also contains the option of computing the state      %%%
%%%             transition matrix corresponding to the transfer.             %%%
%%% Inputs    : mu      - Gravitational constant (neg -> phi)                %%%
%%%             dt      - (s)   Transfer time interval                       %%%
%%%             craft.t - (s)   Initial state time tag                       %%%
%%%             craft.r - (m)   Initial inertial position vector             %%%
%%%             craft.v - (m/s) Initial inertial velocity vector             %%%
%%% Outputs   : craft.t - (s)   Final state time tag                         %%%
%%%             craft.r - (m)   Final inertial position vector               %%%
%%%             craft.v - (m/s) Final inertial velocity vector               %%%
%%%             phi     - (var) Transfer state transition matrix (6x6)       %%%
%%% Comments  : The continued fraction used (n = 5) is defined as the        %%%
%%%             ratio of two hypergeometric functions as follows:            %%%
%%%                                                                          %%%
%%%                                     F(n,1,1+(n/2),q)                     %%%
%%%                 cf = G(n,0,n/2,q) = ----------------                     %%%
%%%                                       F(n,0,n/2,q)                       %%%
%%%                                                                          %%%
%%% Exception : Convergence failure                                          %%%
%%% Reference : Shepperd, S.W., "Universal Keplerian State Transition        %%%
%%%             Matrix," Celestial Mechanics, Vol. 35, 1985.                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ craft, phi ] = kepler (mu, dt, craft)

%%% Summarize Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    io = 0;
    ti = craft.t;
    ri = craft.r;
    vi = craft.v;

    if (io > 0),
        fprintf ('\n');
        fprintf ('==============\n');
        fprintf ('Conic : Kepler\n');
        fprintf ('==============\n\n');
        fprintf ('Input Parameters:\n');
        fprintf ('-----------------\n\n');
        fprintf ('    Transfer time interval   %19.11e\n', dt);
        fprintf ('    Gravitation constant     %19.11e\n', mu);
        fprintf ('    Initial state time tag   %19.11e\n', ti);
        fprintf ('    Initial position vector  %19.11e%19.11e%19.11e\n', ri(1), ri(2), ri(3));
        fprintf ('    Initial velocity vector  %19.11e%19.11e%19.11e\n', vi(1), vi(2), vi(3));
        fprintf ('\n');
    end;

    if (io > 1),
        fprintf ('                              I t');
        fprintf (' e r a t i o n   V a r i a b l e s\n');
        fprintf ('        -----------------------------------');
        fprintf ('-------------------------------------------\n');
        fprintf ('                 u                  t   ');
        fprintf ('               slope              t error\n');
    end;

%%% Initialization Iteration Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (mu > 0), sensitivity = 0; end;
    if (mu < 0), sensitivity = 1; mu = -mu; end;

    u         = 0;
    imax      = 20;
    umax      = +realmax;
    umin      = -realmax;
    tolerance = 1000 * eps;

    ww        = 0;
    orbits    = 0;
    tdesired  = dt;
    threshold = tolerance * abs (tdesired);

    r0        = norm (ri);
    n0        = dot (ri, vi);
    beta      = 2 * (mu / r0) - dot (vi, vi);

    if (beta ~= 0),
        umax = +1 / sqrt (abs (beta));
        umin = -1 / sqrt (abs (beta));
    end;

    if (beta > 0),
        orbits = beta * dt - 2 * n0;
        orbits = 1 + (orbits * sqrt (beta)) / (pi * mu);
        orbits = floor (orbits / 2);
        ww = 2 * pi * orbits / (beta * beta * sqrt (beta));
    end;

%%% Iterate Until Transfer Time Matches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1 : imax,

        q  = beta * u * u;
        q  = q / (1 + q);

        h0 = 1 - 2 * q;
        h1 = 2 * u * (1 - q);

        u0 = 2 * h0 * h0 - 1;
        u1 = 2 * h0 * h1;
        u2 = 2 * h1 * h1;
        uu = ww + 4 * h1 * u2 * u2 * cfgauss (5, q) / 15;
        u3 = beta * uu + u1 * u2 / 3;

        r1     = r0 * u0 + n0 * u1 + mu * u2;
        dt     = r0 * u1 + n0 * u2 + mu * u3;
        slope  = 4 * r1 / (1 + beta * u * u);
        terror = tdesired - dt;

        if (io > 1),
            fprintf ('%10d%19.11e%19.11e%19.11e%19.11e\n', i, u, dt, slope, terror);
        end;

        if (abs (terror) < threshold), break; end;
        if ((i > 1) & (u  ==  uold) ), break; end;
        if ((i > 1) & (dt == dtold) ), break; end;

        uold  = u;
        dtold = dt;
        ustep = terror / slope;

        if (ustep > 0),
            umin = u;
            u    = u + ustep;
            if (u > umax), u = (umin + umax) / 2; end;
        else
            umax = u;
            u    = u + ustep;
            if (u < umin), u = (umin + umax) / 2; end;
        end;

        if (i == imax),
            fprintf ('Kepler : Convergence failure\n\n');
        end;

    end;

%%% Final Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fm  = -(mu / r0) * u2;
    ggm = -(mu / r1) * u2;

    f   = 1 + fm;
    g   = r0 * u1 + n0 * u2;
    ff  = -mu * u1 / (r0 * r1);
    gg  = 1 + ggm;

    tf  = ti + dt;
    rf  = f  * ri + g  * vi;
    vf  = ff * ri + gg * vi;

    craft.t = tf;
    craft.r = rf;
    craft.v = vf;

%%% If Requested, Compute Transition Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    phi = zeros (6,6);

    if (sensitivity == 1),

        mt(1,1) = ff * (u0 / (r0 * r1) + 1 / (r0 * r0) + 1 / (r1 * r1));
        mt(1,2) = (ff * u1 + (ggm / r1)) / r1;
        mt(1,3) = ggm * u1 / r1;
        mt(2,1) = -(ff * u1 + (fm / r0)) / r0;
        mt(2,2) = -ff * u2;
        mt(2,3) = -ggm * u2;
        mt(3,1) = fm * u1 / r0;
        mt(3,2) = fm * u2;
        mt(3,3) = g * u2;

        a0      = mu / (r0 * r0 * r0);
        a1      = mu / (r1 * r1 * r1);
        uu      = g * u2 + 3 * mu * uu;

        mt(1,1) = mt(1,1) - a0 * a1 * uu;
        mt(1,3) = mt(1,3) - a1 * uu;
        mt(3,1) = mt(3,1) - a0 * uu;
        mt(3,3) = mt(3,3) - uu;

        mi      = [ ri, vi ];
        mf      = [ rf, vf ];

        phi11   =  mf * mt(2:3,1:2) * transpose (mi);
        phi12   =  mf * mt(2:3,2:3) * transpose (mi);
        phi21   = -mf * mt(1:2,1:2) * transpose (mi);
        phi22   = -mf * mt(1:2,2:3) * transpose (mi);

        for k = 1 : 3,
            phi11(k,k) = phi11(k,k) + f;
            phi12(k,k) = phi12(k,k) + g;
            phi21(k,k) = phi21(k,k) + ff;
            phi22(k,k) = phi22(k,k) + gg;
        end;

        phi = [ phi11, phi12; phi21, phi22 ];

    end;

%%% Summarize Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (io > 0),

        fprintf ('\n');
        fprintf ('Output Parameters:\n');
        fprintf ('------------------\n');
        fprintf ('\n');
        fprintf ('    Final state time tag     %19.11e\n', craft.t);
        fprintf ('    Final position vector    %19.11e%19.11e%19.11e\n', rf(1), rf(2), rf(3));
        fprintf ('    Final velocity vector    %19.11e%19.11e%19.11e\n', vf(1), vf(2), vf(3));
        fprintf ('\n');

        if (sensitivity == 0), return; end;

        for i = 1 : 6,
            fprintf ('      %13.5e%13.5e%13.5e |%13.5e%13.5e%13.5e\n', ...
                  phi(i,1), phi(i,2), phi(i,3), phi(i,4), phi(i,5), phi(i,6) );
            if ((i == 3) | (i == 6)),
                fprintf ('      -------------------------------------');
                fprintf ('---|---------------------------------------');
                fprintf ('\n');
            end;
        end;

        fprintf ('\n');

    end;

    return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Title     : CF Gauss                                                     %%%
%%% Purpose   : Compute the Gaussian continued fraction.                     %%%
%%% Inputs    : i - (na) Order of continued fraction                         %%%
%%%             q - (na) Argument of continued fraction                      %%%
%%% Outputs   : g - (na) Continued fraction result                           %%%
%%% Comments  : The continued fraction is defined as the ratio of two        %%%
%%%             hypergeometric functions as follows:                         %%%
%%%                                                                          %%%
%%%                                    F(i,1,1+(i/2),q)                      %%%
%%%                 g = G(i,0,i/2,q) = ----------------                      %%%
%%%                                      F(i,0,i/2,q)                        %%%
%%%                                                                          %%%
%%% Exception : Illegal argument                                             %%%
%%% Reference : Shepperd, S.W., "Universal Keplerian State Transition        %%%
%%%             Matrix," Celestial Mechanics, Vol. 35, 1985.                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = cfgauss (i, q)

    if (q > 0.5 + 16 * eps),
        fprintf ('CF Gauss : Illegal argument\n\n');
        keyboard;
    end;

    g = 1; r = 1; s = 1;
    n = 0; l = i - 2; d = i * (i - 2); k = 1 - 2 * i;

    while (1),
        k    = -k;
        l    = l + 2;
        d    = d + 4 * l;
        n    = n + (1 + k) * l;
        r    = d / (d - n * r * q);
        s    = (r - 1) * s;
        gold = g;
        g    = gold + s;
        if (g == gold), break; end;
    end;

    return;




