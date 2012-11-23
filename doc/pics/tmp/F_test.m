## Copyright (C) 2010 Antonio Sánchez
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## fejer

## Author: Antonio Sánchez <antonio@bagheera>
## Created: 2010-10-19


margin=0.5;
width=4.5;
height=3;

fh=figure;
clf;
set(fh,'units','inches');
set(fh, 'position',[0 0 width+2*margin height+2*margin]);
fpos=get(fh,'position');
impos=[margin/fpos(3) margin/fpos(4) width/fpos(3) height/fpos(4)];
set(fh,'PaperPositionMode','manual');
papersize=[fpos(3) fpos(4)];
set(fh,'PaperPosition',[0 0 papersize]);
a=axes();
set(a,'FontSize',20);
set(a,'position',impos);

N=256;
Nlong=2048;
dt=1;
NW=2;
pf1=0.95;
pf2=0.99;

Fs=1;
t=0:N-1;

x = sin(2*pi*0.2*t) + cos(2*pi*0.4*t) + 0.1*randn(size(t));

[S1 f Sc F1]=mtpsd(x,NW,[],2*N,Fs,[],pf1,'eigen');
[S2 f Sc F2, a1, wk]=mtpsd(x,NW,[],2*N,Fs,[],pf1,'adapt');
[S3 f Sc F3]=mtpsd(x,NW,[],2*N,Fs,[],pf2,'eigen');
[S4 f Sc F4]=mtpsd(x,NW,[],2*N,Fs,[],pf2,'adapt');

f(f>=0.5)=f(f>=0.5)-1;
f=fftshift(f);
F1(:,1)=fftshift(F1(:,1));
F1(:,2)=fftshift(F1(:,2));
F2(:,1)=fftshift(F2(:,1));
F2(:,2)=fftshift(F2(:,2));
F3(:,1)=fftshift(F3(:,1));
F3(:,2)=fftshift(F3(:,2));
F4(:,1)=fftshift(F4(:,1));
F4(:,2)=fftshift(F4(:,2));

 b1=find(F1(:,1)>F1(:,2));
 b2=find(F2(:,1)>F2(:,2));
 b3=find(F3(:,1)>F3(:,2));
 b4=find(F4(:,1)>F4(:,2));

plot(f,log10([F1 F3(:,2)]), 'linewidth',2);
axis([-0.5 0.5 -2 3]);
xlabel('f');
ylabel('log10(F)');
filename=['Fstat1_',num2str(N)];
print('-depsc2','-r300',[filename,'.eps']);
system(['epspdf ',filename,'.eps']);

plot(f,log10([F2 F4(:,2)]), 'linewidth',2);
axis([-0.5 0.5 -2 3]);
xlabel('f');
ylabel('log10(F)');
filename=['Fstat2_',num2str(N)];
print('-depsc2','-r300',[filename,'.eps']);
system(['epspdf ',filename,'.eps']);
