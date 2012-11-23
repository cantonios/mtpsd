/*
 * 
 C opyright (c) 20*10 C. Antonio Sanchez
 
 This file is part of DPSS.
 
 DPSS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 DPSS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with DPSS.  If not, see <http://www.gnu.org/licenses/>.
 
 This component makes use of the FFTW3 package for computing FFTs
 <http://www.fftw.org>
 
 */  

#include "dpss.h"

#include <complex.h>
#ifndef complex
#define complex             // in msys, complex.h doesn't define "complex", needed by fftw
#endif
#include <fftw3.h>
#include <cmath>

// integrates   int( y(x), x=0..xfraction ), 
// where x assumed to be [0, (n-1)/n] (so dx=1/n)
// Approximates function with natural cubic splines
double spline_integral( double *y, uint_t n, double xfraction){
    
    double *c=new double[n];
    double *d=new double[n];
    
    //TRIDIAGONAL SOLVE to find splines
    //forward sweep
    c[0]=0.5; d[0]=0.5*3.0*n*(y[1]-y[0]);
    for ( uint_t ii=1; ii<(n-1); ii++ ){
        c[ii]= 1.0/(4.0-c[ii-1]);
        d[ii]= ( 3.0*n*(y[ii+1]-y[ii-1]) - d[ii-1] ) * c[ii];
    }
    c[n-1]=0;
    d[n-1]=( 3.0*n*( y[n-1]-y[n-2] ) - d[n-2])/(2.0-c[n-2]);
    
    //back substitution
    for (uint_t ii=n-1; ii>1; --ii){
        d[ii]=d[ii]-c[ii]*d[ii+1];
    }
    
    // Now we have coefficients of cubic splines, we should be able to integrate
    double xbound=xfraction*n;
    double integral=0;
    double dx, dy, aa, bb, cc, dd, step;
    for ( uint_t ii=0; ii<xbound; ii++){
        dy=y[ii+1]-y[ii];
        aa=y[ii];
        bb=d[ii];
        cc= 3*dy*n-2*d[ii]-d[ii+1];
        dd=-2*dy*n+d[ii]+d[ii+1];  
        
        if (ii < xbound-1)
            step=aa+bb/n/2.0+cc/n/3.0 +dd/n/4;
        else if (ii < xbound){
            dx=xbound-ii;
            step=aa*dx+bb/n/2.0*pow(dx,2)+cc/n/3.0*pow(dx,3) +dd/n/4.0*pow(dx,4);
        }
        integral+=step;
    }
    
    integral=integral/n;    //because of normalized x
    
    delete [] c;
    delete [] d;
    
    return integral;
}

double spline_integral_fraction( double *y, uint_t n, double xfraction){
    
    double *c=new double[n];
    double *d=new double[n];
    
    //TRIDIAGONAL SOLVE
    //forward sweep
    c[0]=0.5; d[0]=0.5*3.0*n*(y[1]-y[0]);
    for ( uint_t ii=1; ii<(n-1); ii++ ){
        c[ii]= 1.0/(4.0-c[ii-1]);
        d[ii]= ( 3.0*n*(y[ii+1]-y[ii-1]) - d[ii-1] ) * c[ii];
    }
    c[n-1]=0;
    d[n-1]=( 3.0*n*( y[n-1]-y[n-2] ) - d[n-2])/(2.0-c[n-2]);
    
    //back substitution
    for (uint_t ii=n-1; ii>1; --ii){
        d[ii]=d[ii]-c[ii]*d[ii+1];
    }
    
    // Now we have coefficients of cubic splines, we should be able to integrate
    double totalintegral=0, integral=0;
    double dx, dy,aa,bb,cc,dd;                     //coefficients of cubic
    double xbound=xfraction*(n);
    double step=0;
    for ( uint_t ii=0; ii<n-1; ii++){
        dy=y[ii+1]-y[ii];
        aa=y[ii];
        bb=d[ii];
        cc= 3*dy*n-2*d[ii]-d[ii+1];
        dd=-2*dy*n+d[ii]+d[ii+1];  
        step=aa+bb/n/2.0+cc/n/3.0 +dd/n/4.0;
        totalintegral+=step;
        
        if (ii < xbound-1)
            integral+=step;
        else if (ii < xbound){
            dx=xbound-ii;
            step=aa*dx+bb/n/2.0*pow(dx,2)+cc/n/3.0*pow(dx,3) +dd/n/4.0*pow(dx,4);
            integral+=step;
        }
    }
    
    delete [] c;
    delete [] d;
    
    return integral/totalintegral;
}

double linear_integral_fraction( double *y, uint_t n, double xfraction){
    
    // Now we have coefficients of cubic splines, we should be able to integrate
    double totalintegral=0, integral=0;
    double dx, yn;
    double xbound=xfraction*(n);
    double step=0;
    for ( uint_t ii=0; ii<n-1; ii++){
        step=(y[ii+1]+y[ii])/2;
        totalintegral+=step;
        
        if (ii < xbound-1)
            integral+=step;
        else if (ii < xbound){
            //linearly interpolate
            dx=xbound-ii;
            yn=(y[ii+1]-y[ii])*dx+y[ii];
            step=(yn+y[ii+1])/2*dx;
            integral+=step;
        }
    }
    
    return integral/totalintegral;
}

void compute_energy_concentrations(double *h, uint_t n, double nW, double *lambda, uint_t K){
    
    uint_t N=4*n;
    uint_t Nmid=floor(N/2);
    double W=nW/n;    
    double *S=new double[Nmid];
    
    fftw_complex *H;
    fftw_plan p;
    H = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
    
    //transform taper
    p = fftw_plan_dft_1d(N, H, H, FFTW_FORWARD,FFTW_ESTIMATE);
    
    for (uint_t ii=0; ii<K; ii++){
    
        for (uint_t jj=0; jj<n; jj++){
            H[jj]=h[jj+ii*n]+0*I;
        }
        for (uint_t jj=n; jj<N; jj++){
            H[jj]=0;
        }
        
        fftw_execute(p);
        
        for (uint_t jj=0; jj<Nmid; jj++){
            S[jj]=cabs(H[jj])*cabs(H[jj]);
        }
        
        //compute concentration
        lambda[ii]=spline_integral_fraction(S,Nmid,(W*N)/Nmid);
        //fix round-off
        if (lambda[ii]>1) lambda[ii]=1;
        else if (lambda[ii]<0 ) lambda[ii]=0;
        
    }
    
    //ensures strictly decreasing eigenvalues (biased upwards toward one)
    for (uint_t ii=K-1; ii>0; --ii){
        if (lambda[ii]<lambda[ii+1]) lambda[ii]=lambda[ii+1];
        if (lambda[ii] > 1.0 ) lambda[ii]=1.0;
    }
    
    //clean up 
    fftw_destroy_plan(p);
    fftw_free(H);
    fftw_cleanup();
    delete [] S;
    
    return;
    
}