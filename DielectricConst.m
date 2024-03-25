% Function for calculating the dielectric constant of a sample from KPFM
% scans. Caculations use an interative Newton Method to solve an expression
% derived for the extraction of the required parameters.
% Further detail can be found in my FYP - Peter Meehan

% x0 is the inital guess for x
% tol is the tolerance level for root convergence.
% G is the average value of the Lock-in amplitude for the sample.
% A0 is the free oscillation amplitude of the cantilever (near the set
% point)
% A is the steady state amplitude of the cantilever

% How to Use
% Type DielectricConst(x0,tol,G,A0,A,h,U,k) into the command window
% after running this code, inputting experimental values for G, A0, A, h, and U .
% Example: DielectricConst(0.1,1e-4,1.57E-3,16.48E-9,30.55E-9,100E-9,1)


function [E,iters] = DielectricConst(x0,tol,G,A0,A,h,U)

k = 2.8;
r = 25e-09;
o = 20;
y = r*(1-sind(o));
Q = 246.07;
a = A/y;


c=G*(2*A0*k*y)/((pi)*Q*(U^2)*(8.85e-12)*r*a);
f = @(x) (2/a^2)*(1/sqrt(1-(a^2)*x^2)-(1+x)/sqrt((1+x)^2-(a^2)*(x^2))) -c;
df = @(x) (2/a^2)*((a^2*x)/(1-a^2*x^2)^(3/2) + (1+x)*(-2*x*a^2 +2*(1+x))/(2*(-a^2*x^2+(1+x)^2)^(3/2)))-1/(sqrt(-a^2*x^2+(1+x)^2));

xg = x0; 
incr = 1;
iters = 0; 
m=1;

while (abs(incr) > tol)
    incr = -m*f(xg)/df(xg);
    xg = xg + incr;
    
    iters = iters + 1;
end

% Final value
x = xg;
X = abs(x);
p = y/X -100E-9;
E = h/p;


% Uncertainty Calculations
orad = o*(pi/180);

deltaA0 = 0.0005E-9;
deltak = 0.05;
deltar = 0.5E-9;
deltao = 0.009;
deltaQ = 0.005;
deltaU = 0.05;
deltaA = 0.0005E-9;
deltaG = 0.005E-3;
deltah = 0.005E-9;
deltaZ = 0.005E-9;

fractA0 = deltaA0/A0;
fractk = deltak/k;

deltay = sqrt((deltar)^2*(y/r)^2  +((r*cos(orad))*deltao)^2 );
fracty = deltay/y;
fractQ = deltaQ/Q;
fractU = deltaU/U;
fractr = deltar/r;
fractA = deltaA/A;
fractG = deltaG/G;



deltac = c*sqrt((fractA0)^2 + (fractk)^2 + 2*(fracty)^2 + (fractQ)^2 + (2*fractU)^2 + (fractr)^2 + (fractA)^2 + (fractG)^2);
deltax = deltac +1e-4;

fractx = deltax/x;
fracth = deltah/h;
deltap = deltaZ + (y/x)*sqrt((fracty)^2+(fractx)^2);
fractp = deltap/p;
deltaE = E*sqrt((fracth)^2+ (fractp)^2);

fprintf('The Dielectric Constant is %s +/- %d.\n',E,deltaE);
format short;