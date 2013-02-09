function [S, f, Sc] = mtcpsd(x,y, NW, ntaper, NFFT, method, pconf)
%Computes the two-sided multi-taper power spectral density
% Inputs:
%   x,y the data
%   NW time-bandwidth product (2.5 is a good value)
%   ntaper number of tapers (2*NW-1 is good)
%   NFFT number of points to use for FFT
%   method equal or eigen for weighting
%   pconf confidence probability

% Outputs:
%   S the power spectral density
%   f the frequencies
%   Sc confidence intervals
    
    if (nargin < 5)
        NFFT = length(x);
    end
    if (nargin < 7)
        pconf = 0.98;
    end
    switch(nargin)
        case 1
            y = x;
            method='eigen';
            NW=4;
            ntaper=2*NW-1;
        case 2
            method='eigen';
            NW=4;
            ntaper=2*NW-1;
        case 3
            method='eigen';
            ntaper=2*NW-1;
        case 4
            method='eigen';
        case 5
            method='eigen';
        case 6
        case 7
        otherwise
            disp('Invalid number of arguments.');
            return;
    end
    
    N=length(x);        
    
    if (size(x,1)==1)
        x=transpose(x);     %don't want to change complex component
    end
    
    if (size(y,1)==1)
        y = transpose(y);
    end
    
    repk=ones(1,ntaper);            %used for duplicating vector in ntaper columns
    repn=ones(N,1);                 %used for duplicating vector in N rows
    reps=ones(NFFT,1);              %used for duplicating vector in NFFT rows
    
    [h,l]=dpss(N,NW,ntaper);        %Data taper
    
    f = linspace(0, 1, NFFT+1);
    f = f(1:end-1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Remove Mean to reduce bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mx=mean(x);                      %DC component
    my=mean(y);
    xx=x*repk;
    yy=y*repk;

    %Calculate weighted means
    leven=1:2:ntaper;               %zero out the mean for even DPSS to prevent leakage
    
    % remove mean
    wm=mx*repk;
    wm(leven)=sum(xx(:,leven).*h(:,leven),1)./sum(h(:,leven),1);
    xx=( xx - repn*wm );
    
    wm=my*repk;
    wm(leven)=sum(yy(:,leven).*h(:,leven),1)./sum(h(:,leven),1);
    yy =(yy-repn*wm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Calculate Eigenspectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Window series
    cwx=h.*xx;
    cwy=h.*yy;

    %Fourier, normalized
    Jkx=fft(cwx,NFFT,1);
    Jkx=Jkx/sqrt(N);             %normalized so independent of N
    Jkx(1,:)=mx;                 %replace mean
    
    Jky=fft(cwy,NFFT,1);
    Jky=Jky/sqrt(N);             %normalized so independent of N
    Jky(1,:)=my;                 %replace mean
    
    %Calculate Power
    Sk=conj(Jkx).*Jky;     

    
    if (strcmp(method,'equal')==1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Unit Weighting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        wk=reps*ones(1,ntaper)/ntaper;    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Weight by energy concentrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        wk=reps*(l'/sum(l));          
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Power Estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S=sum(wk.*Sk,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Confidence intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (nargout > 2)

        p=(1-pconf)/2;
        v=2./sum(wk.^2,2);  

        %correct degrees of freedom at zero and fn
    %     v(1)=v(1)/2;
    %     if (mod(NFFT,2)==0)
    %         v(NFFT/2+1)=v(NFFT/2+1)/2;
    %     end

        Ql=chisqinv(p,v);
        Qu=chisqinv(1-p,v);

        Sl=v.*S./Qu;
        Su=v.*S./Ql;

        Sc=[Sl Su];
        
    end
end
