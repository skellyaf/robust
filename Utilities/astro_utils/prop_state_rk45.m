%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnew = prop_state_rk45(diffeq,xold,dt);

h=dt;
hh = h/2;
h6 = h/6;

y = xold;
dydx = feval(diffeq,y);
yt = y + hh*dydx;

dyt = feval(diffeq,yt);
yt = y +hh*dyt; 

dym = feval(diffeq,yt);
yt = y +h*dym; 
dym = dyt + dym;

dyt = feval(diffeq,yt);
yout = y + h6*(dydx+dyt+2*dym);
xnew = yout; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%