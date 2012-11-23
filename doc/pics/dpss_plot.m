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
set(a,'FontSize',12);
set(a,'position',impos);

N=10^7;
NB=2^16;

NW=2.5;
R=[2 6];
plotrange=1:100:N;
Np=length(plotrange);

%[h l]=dpss(N, NW,R,'spline',NB,'trace');

plot(h(plotrange,:), 'linewidth',2);
legend(['h_1';'h_2 ';'h_3 ';'h_4 ';'h_5 ']);
axis([0 Np-1 -0.001 0.001]);
xlabel('t');
set(gca,'XTick',[0:(Np-1)/4:(Np-1)])
set(gca,'XTickLabel',['0'; '2.5x10^6'; '5x10^6'; '7.5x10^6'; '10^7']);

ylabel('h');
filename=['dpss',num2str(N)];
print('-depsc2','-r300',[filename,'.eps']);
system(['epspdf ',filename,'.eps']);
%system(['pdfcrop ',filename,'.pdf ',filename,'.pdf']);