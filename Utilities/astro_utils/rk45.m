%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnew = rk45(F,xold,t,dt)

y = xold;
x = t;
dydx = feval(F,y,x);

h=dt;
hh = h/2;
h6 = h/6;
xh = x + hh;

%keyboard
yt = y + hh*dydx;

dyt = feval(F,yt,xh);
yt = y +hh*dyt; 

dym = feval(F,yt,xh);
yt = y +h*dym; 
dym = dyt + dym;

dyt = feval(F,yt,x+h);
yout = y + h6*(dydx+dyt+2*dym);
xnew = yout; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%