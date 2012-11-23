function [ xp ] = gaussinv( p )
%gaussint(p)  given percentage point p, finds x s.t.
%             cdf_gauss(x)=p
%             Algorithm AS 70, R.E. Odeh, J.R. Evans (1974)

    alim=1e-20;
    p0=-0.322232431088;
    p1=-1;
    p2=-0.342242088547;
    p3=-0.0204231210245;
    p4=-0.453642210148e-4;
    q0=0.993484626060e-1;
    q1=0.588581570495;
    q2=0.531103462366;
    q3=0.103537752850;
    q4=0.38560700634e-2;
    
    ps=p;
    if (p>0.5)
        ps=1-p;
    end
    if (ps<alim)
        xp=NaN;
    elseif (ps==0.5)
        xp=0;
    else
       y=sqrt(log(1/(ps*ps)));
       xp=y+((((y*p4+p3)*y+p2)*y+p1)*y+p0)/((((y*q4+q3)*y+q2)*y+q1)*y+q0);
    end   
    if (p<0.5)
        xp=-xp;
    end

end

