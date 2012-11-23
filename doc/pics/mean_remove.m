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


N=128;
Nlong=2048;
dt=1;
NW=1.9;

Fs=1;
t=0:N-1;

x = 3 + sin(2*pi*0.016*t);
[S f]=mtpsd(x,NW,Nlong);
Sm=mtpsd(x,NW,Nlong,'mean');
f(f>=0.5)=f(f>=0.5)-1;
f=fftshift(f);
Sm=fftshift(Sm);
S=fftshift(S);

plot(f,10*log10([S Sm]), 'linewidth',2);
%hlegend=legend(['h_0';'h_1';'h_2']);
%set(hlegend,'FontSize',6);
axis([-0.51 0.51 -120 20]);
xlabel('f');
ylabel('dB');
filename=['mean_remove',num2str(N)];
print('-depsc2','-r300',[filename,'.eps']);
system(['epspdf ',filename,'.eps']);
%system(['pdfcrop ',filename,'.pdf ',filename,'.pdf']);