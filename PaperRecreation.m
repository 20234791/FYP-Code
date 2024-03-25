%Calculation of Transfer Function and Verification of Method 
c=2.6;
Q=283.7;
u=(pi)*4E-7;
A = readmatrix('PaperPhase.txt');
B= readmatrix('PaperMagnetisation2.txt');
Anew =  A(1:615 , 1:615);
Bnew =  B(1:615 , 1:615);
%Fourier transform of the phase image
G=fftshift(fft2(Anew));

% Calculation of the wave vector k
x=  -307:1:307;
y=x;
[XX,YY]=meshgrid(x,y);
rr= (sqrt(XX.*XX+YY.*YY));
r = 2*pi./rr;
kk= (sqrt(XX.*XX+YY.*YY)).*(2*pi*5.11E-6)*(1/615);

h=130E-9;
z=60E-9;
% Implementation of equation 4.1.1
H=fftshift(fft2(Bnew).*(1-exp(-kk.*h)).*exp(-kk.*z).*0.5*500E3);
Hb=conj(H);

%Wiener deconvolution using the optimisation parameter from [11]
F=G.*(Hb./((H.*conj(H))+405.6));
 


TTF=conj(F)*(c*pi)/(180*Q*u*((10E-9)^2)*500E3);
T= ifftshift((ifft2(TTF)));

D = G.*(TTF./((F.*conj(F))+1E6));
HField = ifft2(ifftshift(D)).*(c*pi)/(360*Q*u*((10)^2)*20*500E3);


figure;
subplot(221);imagesc(abs(T));colorbar;title('TTF');
cbr = colorbar;
cbr.Label.String = 'A/m^2';
subplot(222);imagesc(real(HField));colorbar;title('Calculated H Field: Paper');
cbr = colorbar;
cbr.Label.String = 'kA/m';

E = readmatrix('PaperPhaseTest.txt');
Enew = E(200:500 , 200:500);
E1 = fftshift(fft2(Enew));
TTF300 = TTF(150:450 , 150:450);
F300 = F(150:450 , 150:450);
E2 = E1.*(TTF300./((F300.*conj(F300))+1E6));
MySample = ifft2(ifftshift(E2)).*(c*pi)/(360*Q*u*((10)^2)*20*500E3);
subplot(224);imagesc(real(MySample));colorbar;title('Test Sample: Paper')
cbr = colorbar;
cbr.Label.String = 'kA/m';



P = readmatrix('HorizontalLines_levelled_strayfield.txt');
P1 = fftshift(fft2(P));
TTF512 = TTF(51:562 , 51:562);
F512 = F(51:562 , 51:562);
P2= P1.*(TTF512./((F512.*conj(F512))+2.5E3));
MyScan = ifft2(ifftshift(P2)).*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
subplot(223);imagesc(real(MyScan));colorbar;title('My Sample: Horizontal Lines')
cbr = colorbar;
cbr.Label.String = 'A/m';


impixelinfo



%Uncertainty calculations

uG= ones(615).*0.005;
FuGR2= (real(G)).^2.*(uG.^2);
FuGS2= (imag(G)).^2.*(uG.^2);
FuHR2= kk.^2.*(real(H)).^2.*(0.5E-9)^2 + (kk.*exp(-kk.*h).*exp(-kk.*z))^2.*real(Bnew).^2.*(0.5E-9)^2;
FuHS2= kk.^2.*(imag(H)).^2.*(0.5E-9)^2 + (kk.*exp(-kk.*h).*exp(-kk.*z))^2.*imag(Bnew).^2.*(0.5E-9)^2;

RRG= real(H)./(real(H).^2+imag(H).^2+406.5);
RSG= imag(H)./(real(H).^2+imag(H).^2+406.5);
RRH= real(G)./(real(H).^2+imag(H).^2+406.5);
RSH= imag(G)./(real(H).^2+imag(H).^2+406.5);

FuFR2= RRG.^2.*FuGR2 + RSG.^2.*FuGS2 + RRH.^2.*FuHR2 + RSH.^2.*FuHS2 +405.6*100;
FuFS2= RSG.^2.*FuGR2 + RRG.^2.*FuGS2 + RSH.^2.*FuHR2 + RRH.^2.*FuHS2 +405.6*100;
FuF= sqrt(FuFR2 + FuFS2);


deltaQ = 0.005;
deltac = 0.05;
deltaPix = 0.025;

Pix1 = 10;
Pix2 = 50;
Pix3 = 5;
fractQ = deltaQ/Q;
fractc = deltac/c;
fractPix1 = deltaPix/Pix1;
fractPix2 = deltaPix/Pix2;
fractPix3 = deltaPix/Pix3;

FuTTF2 = abs(TTF).^2.*((fractQ)^2+(fractc)^2+(fractPix1)^2 + (0.5/500)^2) + FuF.^2.*(c*pi)/(180*Q*u*((10)^2)*500E3);
FuTTF= sqrt(FuTTF2);

FuF512 = FuF(51:562 , 51:562);
FuF256 = FuF(179:434 , 179:434);
FuF128 = FuF(243:370 , 243:370);


%Uncertainty Calculations for my Samples:

HszHLine = readmatrix("HorizontalLinesStrayField.txt");
uHsz2HLine = HszHLine.^2.*((fractQ)^2+(fractc)^2+(fractPix2)^2)*25 + FuF512.^2.*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
uHszHLine = sqrt(uHsz2HLine);

HszVLine = readmatrix("vertLines_strayfield.txt");
uHsz2VLine = HszVLine.^2.*((fractQ)^2+(fractc)^2+(fractPix2)^2)*25 + FuF512.^2.*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
uHszVLine = sqrt(uHsz2VLine);

%figure for lines stray fields
figure;
subplot(121);imagesc(HszHLine);colorbar;impixelinfo;title('Horizontal Line Structure')
subplot(122);imagesc(HszVLine);colorbar;impixelinfo;title('Vertical Line Structure')


%Uncertainty for Magnetised Sample
HszLine1 = readmatrix("line 1_strayfield.txt");
uHsz2Line1 = HszLine1.^2.*((fractQ)^2+(fractc)^2+(fractPix3)^2)*25 + FuF256.^2.*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
uHszLine1 = sqrt(uHsz2Line1);


HszLine3 = readmatrix("line 3_strayfield.txt");
uHsz2Line3 = HszLine3.^2.*((fractQ)^2+(fractc)^2+(fractPix3)^2)*25 + FuF256.^2.*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
uHszLine3 = sqrt(uHsz2Line3);


HszLine5 = readmatrix("line 5_strayfield.txt");
uHsz2Line5 = HszLine5.^2.*((fractQ)^2+(fractc)^2+(fractPix3)^2)*25 + FuF128.^2.*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
uHszLine5 = sqrt(uHsz2Line5);


HszLine7 = readmatrix("line 7_strayfield.txt");
uHsz2Line7 = HszLine7.^2.*((fractQ)^2+(fractc)^2+(fractPix3)^2)*25 + FuF256.^2.*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
uHszLine7 = sqrt(uHsz2Line7);


HszLine9 = readmatrix("line 9_strayfield.txt");
uHsz2Line9 = HszLine9.^2.*((fractQ)^2+(fractc)^2+(fractPix3)^2)*25 + FuF256.^2.*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
uHszLine9 = sqrt(uHsz2Line9);


HszLine11 = readmatrix("line 11_strayfield.txt");
uHsz2Line11 = HszLine11.^2.*((fractQ)^2+(fractc)^2+(fractPix3)^2)*25 + FuF256.^2.*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
uHszLine11 = sqrt(uHsz2Line11);


HszSpiral = readmatrix("Spiral_strayfield.txt");
uHsz2Spiral = HszSpiral.^2.*((fractQ)^2+(fractc)^2+(fractPix2)^2)*25 + FuF256.^2.*(c*pi)/(360*Q*u*((50)^2)*20*500E3);
uHszSpiral = sqrt(uHsz2Spiral);

%Plots for Inclusion in Report
figure;
subplot(121);imagesc(uHszHLine);colorbar;impixelinfo;title('Horizontal Line Structure')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(122);imagesc(uHszVLine);colorbar;impixelinfo;title('Vertical Line Structure')
cbr = colorbar;
cbr.Label.String = 'A/m';


figure;
subplot(121);imagesc(HszSpiral);colorbar;impixelinfo;title('Stray Field (A/m)')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(122);imagesc(uHszSpiral);colorbar;impixelinfo;title('Uncertainty (A/m)')
cbr = colorbar;
cbr.Label.String = 'A/m';



figure;
subplot(231);imagesc(HszLine1);colorbar;title('Line 1')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(232);imagesc(HszLine3);colorbar;title('Line 3')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(233);imagesc(HszLine5);colorbar;title('Line 5')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(234);imagesc(HszLine7);colorbar;title('Line 7')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(235);imagesc(HszLine9);colorbar;title('Line 9')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(236);imagesc(HszLine11);colorbar;title('Line 11')
cbr = colorbar;
cbr.Label.String = 'A/m';
impixelinfo


figure;
subplot(231);imagesc(uHszLine1);colorbar;title('Line 1')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(232);imagesc(uHszLine3);colorbar;title('Line 3')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(233);imagesc(uHszLine5);colorbar;title('Line 5')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(234);imagesc(uHszLine7);colorbar;title('Line 7')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(235);imagesc(uHszLine9);colorbar;title('Line 9')
cbr = colorbar;
cbr.Label.String = 'A/m';
subplot(236);imagesc(uHszLine11);colorbar;title('Line 11')
cbr = colorbar;
cbr.Label.String = 'A/m';
impixelinfo
