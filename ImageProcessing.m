% Matlab Transfer Fn
k=3.3;
Q=226;
u=(pi)*4E-7;
P= (360*u*(1E-6)^2)/(k*pi); 
A = readmatrix('TestSamplePhase.txt');
B= readmatrix('TestSampleField.txt');
G=fftshift(fft2(A));
H=fftshift(fft2(B)).*(2/500E6);
Hb=conj(H);

F=G.*(Hb./((H.*conj(H))+0.01));

TTF=conj(F)*(k*pi)/(180*Q*u*((5E-6)^2)*500E3);
T= ifftshift(ifft2(TTF));

C = H.*F;
Phase = ifft2(C);

D = G.*(TTF./((F.*conj(F))+1E6));
HField = ifft2(ifftshift(D));

figure(1);
subplot(221);imagesc(abs(T));colorbar;title('TTF');
subplot(222);imagesc(-abs(HField));colorbar;title('Calculated H Field');
subplot(223);imagesc(-abs(Phase));colorbar;title('Calculated Phase');

E = readmatrix('HorizontalLines_levelled_strayfield.txt');
E1 = fftshift(fft2(E));
E2 = E1.*(TTF./((F.*conj(F))+1E9));
MySample = ifft2(ifftshift(E2));
subplot(224);imagesc(-abs(MySample));colorbar;title('My Sample')


%Second Figure: My own samples using the TTF.

I = readmatrix('HorizontalLines_levelled_strayfield.txt');
J = readmatrix('HorizontalLinesStrayField.txt');
GFig2=fftshift(fft2(I));
HFig2=fftshift(fft2(J)).*(2/500E6);
HbFig2=conj(HFig2);


K = HFig2.*(TTF./((F.*conj(F))+1E9));
PhaseFig2 = ifft2(K);

L = GFig2.*(TTF./((F.*conj(F))+7.5E8));
HFieldFig2 = ifft2(ifftshift(L));


figure(2);
subplot(221);imagesc(abs(T));colorbar;title('TTF');
subplot(222);imagesc(real(HFieldFig2));colorbar;title('Calculated H Field');
subplot(223);imagesc(J);colorbar;title('Gwy H Field');


%Third Figure: Magnetised Sample
M = readmatrix('line 1_Phase.txt');
N = readmatrix('line 1_strayfield.txt');

GFig3=fftshift(fft2(M));
HFig3=fftshift(fft2(N)).*(2/500E6);
HbFig3=conj(HFig3);

TTF256 = TTF(128:383 , 128:383);
F256 = F(128:383 , 128:383);
T256= ifftshift(ifft2(TTF256));



O = GFig3.*(TTF256./((F256.*conj(F256))+7.5E8));
HFieldFig3 = ifft2(ifftshift(O));

figure(3);
subplot(221);imagesc(abs(T256));colorbar;title('TTF');
subplot(222);imagesc(-real(HFieldFig3));colorbar;title('Calculated H Field');
subplot(223);imagesc(N);colorbar;title('Gwy H Field');




impixelinfo
%% 
%  L-curve
k=2.6;
Q=226;
u=(pi)*4E-7;
P=((pi)*k)/(Q*180); 
A = readmatrix('MFM_Transfer_Fn.txt');
B= readmatrix('CalStrayField.txt');
G=fftshift(fft2(A));
H=fftshift(fft2(B))*2/500E3;
Hb=conj(H);
t=[1E-5,1E-4,1E-3,1E-2,0.1,1,10,100,1E3,1E4,1E5,1E6,1E7,1E8,1E9,1E10];
tv=t(:);
for i = 1:16
    F=G.*(Hb./((H.*conj(H))+t(1,i)));
    Y=norm(F);
    Y_vec(i) = Y;
    X=norm(G - F*H);
    X_vec(i) = X;
end
X_vec(:);
Y_vec(:);
norm(H,1);
logx = log(X_vec);
logy = log(Y_vec);

subplot(221);plot(logx)
subplot(224);plot(logy)
subplot(223);plot(logx,logy,'-o')
lcorner(logx,logy)