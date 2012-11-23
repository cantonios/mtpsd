/*

Copyright (c) 2010 C. Antonio Sanchez

This file is part of MTPSD.

MTPSD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MTPSD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MTPSD.  If not, see <http://www.gnu.org/licenses/>.

*/

////////////////////////////////////////////////////////////////////
//  REFER TO appliedstats.h for information
//
//  Warning: Magic Numbers Everywhere!!!  Most of these are
//           coefficients in polynomial expansions (some fitted).
//           Please refer to the papers listed in appliedstats.h
//           for more information.
///////////////////////////////////////////////////////////////////

#include "applied_stats.h"
#include <cmath>

template <class T> T abs ( T a ) { return (0<a)?a:-a; }

// Inverts chi-squared CDF with v degrees of freedom
double chisquared_inv(double p, double v){

    const double p_bound=2e-6;  //limits on p
    double z0=0;                //centre for series expansion
    
    double p1,p2,t,q,a,b,c, s1,s2,s3,s4,s5,s6;    //temp variables
    double G=lgamma(v/2.0);
    double vo2=v/2.0;
    double lg2=log(2);
    
    if (p<p_bound)
        return 0;
    else if( ( p > 1.0-p_bound ) || ( v<= 0 ) )
        return -1; //error value
    
    if ( v < -1.24*log(p) ){    //P -> 0 (small z)
        z0 = pow( (p*vo2*exp(G+vo2*log(2.0)) ), (2.0/v));
    }
    else if ( v<=0.32 ){        //v -> 0
        //estimate z0 using Newton-Raphson method
        a = log(1.0-p);
        z0=0.4;
        q=2*z0;        //forces loop to start
        
        while (abs<double>(q/z0-1.0) > 0.01){  //small tolerance of 1% used since this is just an initial estimate
            q=z0;
            p1=1.0+z0*(4.67+z0);
            p2=z0*(6.73+z0*(6.66+z0));
            t=-0.5+(4.67+2.0*z0)/p1-(6.73+z0*(13.32+3.0*z0))/p2;
            z0=z0-(1.0-exp(a+G+0.5*z0+lg2*(vo2-1.0))*p2/p1)/t;
        }
    }
    else {
        a=gauss_inv(p);
        b=2.0/9.0/v;
        z0 = v * pow( ( a*sqrt(b)+1.0-b ), 3);
        
        if (z0 > 2.2*v+6.0)                 //P -> 1 (large z)
            z0 = -2.0*(log(1.0-p)-(vo2-1.0)*log(z0/2.0)+G);
    }
    
    //By now, z0 is defined, and ready for Taylor expansion
    if (z0 < AS_TOL)
        return z0;          //near-zero solution
    
    q=2*z0*(1+AS_TOL);      //forces loop to start
    c=(v/2.0-1.0);
    
    // Loops until estimate has been modified by a value less than AS_TOL
    while ( abs<double>(q/z0-1) > AS_TOL){
        q=z0;
        p1=z0/2.0;
        p2=p-gamma_int(p1,vo2);
        t=p2*exp(vo2*lg2+G+p1-c*log(z0));    
        b=t/z0;
        a=0.5*t-b*c;
        s1=(210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
        s2=(420.0+a*(735.0+a*(966.0+a*(1141.0+a*1278.0)))) / 2520.0;
        s3=(210.0+a*(462.0+a*(707.0+932.0*a))) / 2520.0;
        s4=(252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a))) / 5040.0;
        s5=(84.0+264.0*a+c*(175.0+606.0*a)) /  2520.0;
        s6=(120.0+c*(346.0+127.0*c)) /  5040.0;
        z0=z0+t*(1.0+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
    }
    
    return z0;
    
}

// Inverts chi-squared CDF for an array of DOFs
void chisquared_inv(double p, double *v, unsigned int n, double *x){
    
    for (unsigned int ii=0; ii<n; ii++)
        x[ii]=chisquared_inv(p,v[ii]);
    
}

void chisquared_inv(double *p, double *v, unsigned int np, unsigned int nv, double *x){
    
    for (unsigned int ii=0; ii<np; ii++){
        for (unsigned int jj=0; jj<nv; jj++){
            x[ii*nv+jj] = chisquared_inv( p[ii], v[jj] );
        }
    }
}

// Inverts standard Gaussian CDF
double gauss_inv(double p){
    
    const double p_limit=1e-10;  //lower probability limit
    
    double x=0, y=0;    //temp variables
    
    //constants in polynomial expansion
    const double p0=-0.322232431088;
    const double p1=-1.0;
    const double p2=-0.342242088547;
    const double p3=-0.204231210245e-1;
    const double p4=-0.453642210148e-4;
    const double q0= 0.993484626060e-1;
    const double q1= 0.588581570495;
    const double q2= 0.531103462366;
    const double q3= 0.103537752850;
    const double q4= 0.38560700634e-2;
    
    //only checks for p<1/2, then flips x if appropriate
    double ps=p;
    if (p>0.5)
        ps=1.0-p;

    if (ps<p_limit)
        ps=p_limit;
    else if (ps==0.5)
        x=0;
    else{
        y = sqrt(log(1.0/(ps*ps)));
        x =     ( ( ( ( y*p4 + p3 )*y + p2)*y + p1)*y + p0);
        x = x / ( ( ( ( y*q4 + q3 )*y + q2)*y + q1)*y + q0);
        x = x + y;
    }
    
    //flip x for lower probability
    if (p<0.5)
        x=-x;
    
    return x;
}

// Computes the incomplete gamma integral
double gamma_int(double x, double n){
    
    double out;
    if ( ( x < n ) || ( ( x >= n) && ( x <= 1 ) ) )
        out=gamma_int_series(x,n);
    else
        out=gamma_int_fraction(x,n);
    
    return out;
    
}

// Uses series expansion
double gamma_int_series(double x, double n){
    
    double c = exp(-x-lgamma(n+1.0))*pow(x,n); //multiplicative constant in expansion
    const double itol= AS_TOL/c;
    
    double S=0; //Series sum
    double term=1;     //Next term to add
    double d=n+1.0;      //Next denominator term to include
    while (term>itol){
        S=S+term;
        term=term*x/d;
        d=d+1;
    }
    
    return c*S;   
}

// Uses a continued fraction expansion
double gamma_int_fraction(double x, double n){
    
    const double oflow=1e20;        //prevent overflow
        
    double c=exp(-x-lgamma(n))*pow(x,n);  //multiplicative constant in expansion
    double itol=AS_TOL/c;
    
    //temporary variables
    double a=1.0-n;
    double b=1.0+a+x;      
    double pn[]={1.0, x, x+1.0, x*b, 0.0, 0.0};

    double term=0;
    
    double gc=pn[2]/pn[3];  //current gamma integral
    double gl=gc+2*itol;    //last gamma integral (initialized to force loop to start)
    
    while ( abs<double>(gc-gl) > itol){
        a=a+1.0;
        b=b+2.0;
        term=term+1.0;
        
        //calculate next two terms in expansion
        for (int ii= 0; ii<2; ii++)
            pn[ii+4]=b*pn[ii+2]-a*term*pn[ii];

        //update new estimate (as long as no division by zero)
        if (pn[5]!=0){
            gl=gc;
            gc=pn[4]/pn[5];
        }
    
        //shift terms
        for (int ii=0; ii < 4; ii++)
            pn[ii]=pn[ii+2];

        //if numbers are overflowing, reduce both numerator and denominator
        if ( abs<double>(pn[4]) >= oflow){
            for (int ii=0; ii<4; ii++)
                pn[ii]=pn[ii]/oflow;
        }
        
    }
    
    gc=1.0-gc*c;
    
    return gc;
}