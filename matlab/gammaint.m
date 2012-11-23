function [ out ] = gammaint( x, p )
% computes I(x,p)= 1/gamma(p) int ( exp(-t)*t^(p-1), t=0..x);
%                  Algorithm AS 73, G.P. Bhattacharjee (1970)

    if ( (x<p) || ( (x>=p)&&(x<=1) ) )
        out=gamma_int_series(x,p);
    else
        out=gamma_int_fraction(x,p);
    end
    
end

function gint=gamma_int_series( x, n)
    
    c = exp(-x)*x^n/gamma(n+1);
    itol= 1e-7/c;
    
    S=0; %Series sum
    term=1;     %Next term to add
    d=n+1;      %Next denominator term to include
    while (term>itol)
        S=S+term;
        term=term*x/d;
        d=d+1;
    end
    
    gint=c*S;   
end

function gint = gamma_int_fraction(x,n)
    
    oflow=1e20;        %prevent overflow
        
    c=exp(-x-log(gamma(n)))*x^n;  %multiplicative constant in expansion
    itol=1e-7/c;
    
    a=1.0-n;
    b=1.0+a+x;      
    pn=[1, x, x+1, x*b, 0, 0];

    term=0;
    
    gc=pn(3)/pn(4);  
    gl=gc+2*itol;    
    
    while ( abs(gc-gl) > itol)
        a=a+1.0;
        b=b+2.0;
        term=term+1.0;
        
        for ii=1:2
            pn(ii+4)=b*pn(ii+2)-a*term*pn(ii);
        end

        if (pn(6)~=0)
            gl=gc;
            gc=pn(5)/pn(6);
        end
    
        for ii=1:4
            pn(ii)=pn(ii+2);
        end

        if ( abs(pn(5)) >= oflow)
            for ii=1:4
                pn(ii)=pn(ii)/oflow;
            end
        end
        
    end
    
    gint=1.0-gc*c;
end