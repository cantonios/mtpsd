function [ fb Cf ft] = linefreq2( S, Jk, H, pp )
%LINEPSD Uses the Harmonic F-Test to find significant line frequencies
%        ft is the list of estimated significant frequencies
%        fb is the list of significant frequencies rounded to nearest bin
%        Cf are the estimated fourier coefficients for fb
%
%        Jk is the transformed weighted data (i.e. Sk=|Jk|^2 )
%        wk are the weights
%        H  are the fourier of the DPSS
%        pp is the percentage point above which a frequency is deemed
%           significant in the F-test

    NFFT=size(Jk,1);
    K=size(Jk,2);
    reps=ones(NFFT,1);
    repk=ones(1,K);
    
    H0= reps*H(1,:);                    % Mean value of DPSSs
    
    H02 =sum(H0.*H0, 2);
    C= sum(Jk.*H0, 2) ./ H02;           %estimated fitting parameters
    Jke = (C*repk).*H0;                 %estimated fourier coefficient
    
    sige2=sum(abs(Jk-Jke).^2,2)/K;    %variance in Jke(f) estimator

    F=(K-1)*(abs(C).^2)./sige2/K.* H02;
    
    % NOTE #2: With adaptive weights, the degrees of freedom can get rather
    % small in regions that are dominated by bias.  This raises the px100%
    % percentage point. E.g. If a single eigenspectrum is favoured, wwsq->1, 
    % and percentage point of F -> infinity.  This shouldn't be a problem, 
    % however, since S(f) dominated by bias usually implies no significant 
    % energy at that frequency.  
    
    a=1-pp;
    F_dof=2*K - 2;  
    b=a.^(2./F_dof);
    
    Fu = F_dof.*(1-b) ./2./b;  % Percentage point threshold (f-independent)
    Fu=reps*Fu;
    
    ii=find(F>Fu);
    nii=length(ii);
    
    frange=zeros(nii,2);     %range of bump

    nline=0;
    jj=1;
    while (jj<nii)
        nline=nline+1;
        frange(nline,1)=max(1,ii(jj)-2);       
        while (ii(jj+1)-ii(jj)==1)
            jj=jj+1;
            if (jj+1>nii)
                break;
            end
        end
        frange(nline,2)=min(ii(jj)+2,NFFT);
        jj=jj+1;
    end
    if (jj==nii)
        nline=nline+1;
        frange(nline,1)=max(1,ii(jj)-2);       
        frange(nline,2)=min(ii(jj)+2,NFFT);
    end
    
    ft=zeros(nline,1);
    for ii=1:nline
        a=frange(ii,1):frange(ii,2);
        Sp=S(a);                        %estimate of power in small section
        ft(ii)=a*Sp/sum(Sp);            %centre of mass 
    end
    fb=round(ft);                               %rounded bins
    Cf=C(fb);


end

