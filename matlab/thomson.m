N=100;
NFFT=N;
NW=2;

dt=0.01;
t=0:dt:(N-1)*dt;
Fs=1/dt;

x=sin(30*2*pi*t)+1; %zeros(1,N);%+10*t;%sin(7.35*2*pi*t)+sin(8.35*2*pi*t);%+0.5*cos(23.23*2*pi*t)+3*sin(76.7*2*pi*t+0.2)+1.3;

sig=0;
x=x+sig*randn(1,N);

[S Sc ft Cf F]=mtpsdwf(x,NW, floor(2*NW)-1, NFFT, 'unity');
[S1 Pxxc]=pmtm(x,NW,NFFT,'twosided');   %rescale to proper units
S1=S1/Fs*2*pi;

f=0:NFFT;
f=f/NFFT/dt;
df=f(2)-f(1);

prange=1:NFFT;%(abs(f-60.13)<20);%(f > 350)&(f < 450);

figure(2);
plot(f(prange),10*log10(S(prange)),'-g',f(prange),10*log10(Sc(prange,:)),'-b', f(prange), 10*log10(S1),'r');
figure(1);
plot(F);


%"true" PSD:
% S=zeros(NFFT,1);
% S(1)=4;
% [mm jj]=min(abs(f-4));
% S(jj)=1^2/2;
% [mm jj]=min(abs(f-7.5/2));
% S(jj)=0.5^2/2;
% [mm jj]=min(abs(f-8));
% S(jj)=3^2/2;

%subplot(2,1,2);
%plot(f(prange),S1(prange),'r',f(prange),Sp(prange),'k',f(prange),Samt(prange),'b');
