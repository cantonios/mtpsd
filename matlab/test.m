N=1000;
NFFT=N;
NW=2;

dt=0.001;
t=0:dt:(N-1)*dt;
Fs=1/dt;

x=sin(17.35*2*pi*t)+sin(28.35*2*pi*t);
sig=0.01;
x=x+sig*randn(1,N);

[S, f, Sc]=mtpsd(x,NW, floor(2*NW)-1,[],'unity');

f(f>=0.5) = f(f>=0.5)-1;
f = f/dt;
df=f(2)-f(1);

S = fftshift(S);
Sc = fftshift(Sc);
f = fftshift(f);

figure(1);
plot(f,10*log10(S),'-g',f,10*log10(Sc),'-b');
title('Spectrum of x');


% cross spectrum
y = 0.5*fftshift(x)+2;

[S, ~, Sc]=mtcpsd(x,y,NW);
S = fftshift(S);
Sc = fftshift(Sc);

figure(2);
plot(f,10*log10(abs(S)),'-g',f,10*log10(abs(Sc)),'-b');
title('Cross Spectrum of (x,y)');