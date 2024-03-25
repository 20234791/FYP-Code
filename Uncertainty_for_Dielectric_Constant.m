fractA0 = deltaA0/A0;
fractk = deltak/k;

deltay = sqrt((deltar)^2*(y/r)^2  +(r*cos(orad))*deltao );
fracty = deltay/y;
fractQ = deltaQ/Q;
fractU = deltaU/U;
fractr = deltar/r;
fractA = detlaA/A;
fractG = deltaG/G;



deltac = c*sqrt((fractA0)^2 + (fractk)^2 + 2*(fracty)^2 + (fractQ)^2 + (2*fractU)^2 + (fractr)^2 + (fractA)^2 + (fractG)^2);
deltax = deltac +1e-4;

fractx = deltax/x;
fracth = deltah/h;
deltap = deltaZ + (y/x)*sqrt((fracty)^2+(fractx)^2);
fractp = deltap/p;
deltaE = E*sqrt((fracth)^2+ (fractp)^2);