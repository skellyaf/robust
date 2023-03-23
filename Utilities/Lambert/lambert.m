
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Title     : Lambert                                                      %%%
%%% Purpose   : Given a conic gravity field, compute the required velocity   %%%
%%%             to transfer between initial and final inertial position      %%%
%%%             vectors in a specified transfer time interval. This routine  %%%
%%%             also contains the option of computing the sensitivity        %%%
%%%             matrix corresponding to the transfer.                        %%%
%%% Inputs    : mu  - Gravitational constant (neg -> phi)                    %%%
%%%             dt  - (s)   Transfer time interval                           %%%
%%%             ri  - (m)   Init. inertial position vector                   %%%
%%%             rf  - (m)   Final inertial position vector                   %%%
%%%             uh  - (na)  Unit angular momentum vector                     %%%
%%% Outputs   : vi  - (m/s) Init. inertial velocity vector (required)        %%%
%%%             vf  - (m/s) Final inertial velocity vector                   %%%
%%%             phi - (var) Transition (sensitivity) matrix (6x6)            %%%
%%% Comments  : 1) This routine does not have a multi-revolution capability  %%%
%%%             2) The unit angular momentum vector determines the plane,    %%%
%%%                and direction, of the transfer                            %%%
%%% Exception : Convergence failure                                          %%%
%%% Reference : Shepperd, S.W., unpublished notes.                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ vi, vf, phi ] = lambert (mu, dt, ri, rf, uh)

%%% Summarize Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    io = 0;

    if (io > 0),
        fprintf ('\n');
        fprintf ('===============\n');
        fprintf ('Conic : Lambert\n');
        fprintf ('===============\n');
        fprintf ('\n');
        fprintf ('Input Parameters:\n');
        fprintf ('-----------------\n');
        fprintf ('\n');
        fprintf ('    Transfer time interval   %19.11e\n', dt);
        fprintf ('    Gravitation constant     %19.11e\n', abs (mu));
        fprintf ('    Init. position vector    %19.11e%19.11e%19.11e\n', ri(1), ri(2), ri(3));
        fprintf ('    Final position vector    %19.11e%19.11e%19.11e\n', rf(1), rf(2), rf(3));
        fprintf ('    Unit angular momentum    %19.11e%19.11e%19.11e\n', uh(1), uh(2), uh(3));
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

    r0 = ri - dot (ri, uh) * uh;
    r1 = rf - dot (rf, uh) * uh;

    if (abs (dot (uh, uh) - 1) > 10000 * eps),
        vi = [ 0; 0; 0 ]; vf = [ 0; 0; 0 ]; phi = zeros (6,6);
        fprintf ('Lambert : Aborting : Normal vector not unitary\n');
        return;
    end;

    if (abs (dot (ri, uh)) > 0.00001),
        fprintf ('Lambert : Warning : Init. position not in plane\n');
    end;

    if (abs (dot (rf, uh)) > 0.00001),
        fprintf ('Lambert : Warning : Final position not in plane\n');
    end;

    imax      = 20;
    maxu      = realmax;
    tolerance = 1000 * eps;

    m0 = norm (r0);
    m1 = norm (r1);
    cc = norm (r1 - r0);
    ss = (m0 + m1 + cc) / 2;

    rc = ss / 2;
    vc = sqrt (mu / rc);
    wc = vc / rc;

    k1 = sqrt (abs ((ss - m0) * (ss - m1)));
    k  = dot  (cross (r0, r1), uh);

    if ((k1 * k1) > (ss * (ss - cc))),
        k  = k / (2 * k1 * ss);
        k2 = 1 - k * k;
    else
        k2 = cc / ss;
        if (k >= 0),
            k  = +sqrt (1 - k2);
        else
            k  = -sqrt (1 - k2);
        end;
    end;

    tdesired = wc * dt;

    if (tdesired > 4 * (1 - k * k * k) / 3),
        u = 0; umin = -1; umax = 1;
    else
        u = 1; umin =  1; umax = maxu;
    end;

%%% Iterate Until Transfer Time Matches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1 : imax,

        q = k * u;
        y = sqrt (q * q + k2);

        if (q <= 0),
            h1 = y - q;
        else
            h1 = k2 / (y + q);
        end;

        h0 = k + u * h1;
        q  = (1 - abs (h0)) / 2;
        uu = (16 / 15) * h1 ^ 5 * cfgauss (5, q);

        if (h0 < 0),
            pp = 2 * pi / sqrt (1 - u * u) ^ 5;
            uu = pp - uu;
        end;

        dt     = 4 * h1 * (k + h0 * h1 * h1 / 3) + (1 - u * u) * uu;
        slope  = 3 * u * uu - 4 * (h1 / y) * (k * k + (k * h0 + h1 * h1) * h1 * h1);
        terror = tdesired - dt;

        if (io > 1),
            fprintf ('%10d%19.11e%19.11e%19.11e%19.11e\n', i, u, dt, slope, terror);
        end;

        if (abs (terror) < wc),
            if (abs (terror) < tolerance * abs (tdesired )), break; end;
            if (abs (terror) < tolerance * abs (slope * u)), break; end;
        end;

        if ((i > 1) & (u  ==  uold)), break; end;
        if ((i > 1) & (dt == dtold)), break; end;

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
            fprintf ('Lambert : Convergence failure\n');
        end;

    end;

%%% Final Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h1 = h1 / vc;
    uu = uu / vc ^ 5;

    h  = k1 / h1;
    n0 = +(k * ss - m0 * h0) / h1;
    n1 = -(k * ss - m1 * h0) / h1;
    v0 = (n0 * r0 + h * cross (uh, r0)) / (m0 * m0);
    v1 = (n1 * r1 + h * cross (uh, r1)) / (m1 * m1);

    vi  = v0;
    vf  = v1;

%%% If Requested, Compute Transition Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    phi = zeros (6,6);

    if (sensitivity == 1),

        u2  = 2 * h1 * h1;
        u1  = 2 * h0 * h1;
        u0  = 2 * h0 * h0 - 1;

        fm  = -mu * (u2 / m0);
        ggm = -mu * (u2 / m1);

        f   = 1 + fm;
        g   = 2 * k * ss * h1;
        ff  = -mu * u1 / (m0 * m1);
        gg  = 1 + ggm;

        uu    =  g * u2 + 3 * mu * uu;
        omega = ff * uu - (fm + ggm) * u2;

        r0r1 = dot (r0, r1);
        r0v1 = dot (r0, v1);
        v0r1 = dot (v0, r1);
        v0v1 = dot (v0, v1);

        dd = g * g + omega * u2 * h * h ...
             + ((v0r1 - r0v1) * u2 - v0v1 * uu) * g;

        if ((g == 0) | (dd == 0)),

            fprintf ('Lambert : Singular sensitivity matrix\n');

        else

            m(1,1) = (omega * v0v1 - g * ff) * u2;
            m(1,2) = -(ggm / (m1 * m1)) * (((m0 + m1) ...
                     * (2 * m0 + n1 * u1) + g * n1) * u2 ...
                     - ((1 + u0) * n1 - 2 * n0) * uu);
            m(1,3) = -(omega * r0v1 + g * ggm) * u2;
            m(2,1) = -(fm / (m0 * m0)) * (((m0 + m1) ...
                     * (2 * m1 - n0 * u1) - g * n0) * u2 ...
                     + ((1 + u0) * n0 - 2 * n1) * uu);
            m(2,2) = (m0 * m1 * omega - g * g) * u2 + g * uu;
            m(2,3) = (m1 * u1 * omega - g * fm * u2);
            m(3,1) = -(omega * v0r1 - g * fm) * u2;
            m(3,2) = -(m0 * u1 * omega - g * ggm * u2);
            m(3,3) = (omega * r0r1 + g * g) * u2 - g * uu;

            m      = m / (g * dd);

            mi     = [ r0, v0 ];
            mf     = [ r1, v1 ];

            phi(1:3,1:3) = -mi * [ m(2,1), m(2,3); m(2,3), m(2,2) ] * transpose (mi);
            phi(1:3,4:6) = -mf * [ m(1,1), m(1,3); m(3,1), m(3,3) ] * transpose (mi);
            phi(4:6,1:3) = +mi * [ m(1,1), m(3,1); m(1,3), m(3,3) ] * transpose (mf);
            phi(4:6,4:6) = +mf * [ m(1,2), m(3,2); m(3,2), m(2,2) ] * transpose (mf);

            for k = 1 : 3,
                phi(k,k)     = phi(k,k)     - ( f / g);
                phi(k,k+3)   = phi(k,k+3)   + ( 1 / g);
                phi(k+3,k)   = phi(k+3,k)   - ( 1 / g);
                phi(k+3,k+3) = phi(k+3,k+3) + (gg / g);
            end;

        end;

    end;

%%% Summarize Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (io > 0),

        fprintf ('\n');
        fprintf ('Output Parameters:\n');
        fprintf ('------------------\n');
        fprintf ('\n');
        fprintf ('    Reqd. velocity vector    %19.11e%19.11e%19.11e\n', vi(1), vi(2), vi(3));
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




