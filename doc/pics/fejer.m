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


margin=0.35;
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


N=64;
Nlong=1000;
dt=1;

f=linspace(-0.5,0.5,Nlong);
F=dt*sin(N*dt*pi*f).^2/N./sin(pi*f*dt).^2;

plot(f,10*log10(F), 'linewidth',2);
axis([-0.51 0.51 -40 20]);
xlabel('f');
ylabel('dB');
filename=['fejer',num2str(N)];
print('-depsc2','-r300',[filename,'.eps']);
system(['epspdf ',filename,'.eps']);
%system(['pdfcrop ',filename,'.pdf ',filename,'.pdf']);