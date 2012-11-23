function [ x ] =chisqinv(p,v)
%chisqinv calculates the percentage point of a 
%       chi-squared distribution with v degrees
%       of freedom.  Based on three Fortran 77 routines:
%             Algorithm AS 91, D.J. Best and D.E. Roberts (1975)
%             Algorithm AS 70, R.E. Odeh, J.R. Evans (1974)
%             Algorithm AS 73, G.P. Bhattacharjee (1970)

    n=length(v);
    m=length(p);
    x=zeros([n m]);
    
    %can't vectorize because each element is iterative
    %although, maybe I can simplify if v is about the same
    %in several entries
    for ii=1:n
        for jj=1:m
            x(ii,jj)=pw_chisqinv(p(jj),v(ii));
        end
    end

end

function [ x ] =pw_chisqinv( p, v )
%pw_chisqinv calculates the percentage point of a 
%       chi-squared distribution with v degrees
%       of freedom in a pointwise fashion.
%       Algorithm AS 91, D.J. Best and D.E. Roberts (1975)
%             Algorithm AS 70, R.E. Odeh, J.R. Evans (1974)
    tol=1e-7;
    G=log(gamma(v/2));
    
    if (p<0.000002)
        x=0;
        return;
    elseif( (p>0.999998) || (v<=0) )
        x=NaN;
        return;
    end
        
    if ( v < -1.24*log(p))                    %P -> 0 (small z)
        z0 = (p*v*2^(v/2-1)*gamma(v/2))^(2/v);         
    elseif (v<=0.32)                          %v -> 0
        %estimate z0 using Newton-Raphson method
        z0=0.4;
        a=log(1-p);
        q=2*z0;        %forces loop to start
        while (abs(q/z0)-1>0.01)
            q=z0;
            p1=1+z0*(4.67+z0);
            p2=z0*(6.73+z0*(6.66+z0));
            t=-0.5+(4.67+2.0*z0)/p1-(6.73+z0*(13.32+3.0*z0))/p2;
            z0=z0-(1-exp(a+G+0.5*z0+log(2)*(v/2-1))*p2/p1)/t;
        end
    else
        xg=gaussinv(p);                      %General
        z0 = v*(xg*sqrt(2/9/v)+1-(2/9/v))^3;     
        if (z0 > 2.2*v+6)                     %P -> 1 (large z)
            z0 = -2*(log(1-p)-(v/2-1)*log(z0/2)+log(gamma(v/2)));
        end
    end
    
    %By now, z0 is defined, and ready for Taylor expansion
    
    if (z0 < tol)
        x=z0;       %small value of z ( --> p was small)
        return;
    end
    
    q=2*z0;
    while (abs(q/z0-1)>tol)
        q=z0;
        p1=z0/2;
        p2=p-gammaint(p1,v/2);
        t=p2*exp(v/2*log(2)+G+p1-(v/2-1)*log(z0));    
        c=(v/2-1);
        b=t/z0;
        a=0.5*t-b*c;
        s1=(210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
        s2=(420.0+a*(735.0+a*(966.0+a*(1141.0+a*1278.0)))) / 2520.0;
        s3=(210.0+a*(462.0+a*(707.0+932.0*a))) / 2520.0;
        s4=(252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a))) / 5040.0;
        s5=(84.0+264.0*a+c*(175.0+606.0*a)) /  2520.0;
        s6=(120.0+c*(346.0+127.0*c)) /  5040.0;
        z0=z0+t*(1.0+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
    end
    x=z0;
    
end
