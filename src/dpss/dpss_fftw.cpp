/*

Copyright (c) 2010 C. Antonio Sanchez

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

//computes the energy concentrations *lambda, of dpss sequences *h
void compute_energy_concentrations(double *h,uint_t n, uint_t K, double nW, double *lambda);


///////////////////////////////////////////////////////////////////////////////////////////
// CODE
///////////////////////////////////////////////////////////////////////////////////////////

void compute_energy_concentrations(double *h,uint_t n, uint_t K, double nW, double *lambda){
    
    //kernel and eigenvector
    double *ev=(double*) fftw_malloc(sizeof(double)*2*n);
    fftw_complex *KERN, *EV;
    
    fftw_plan p, pi;

    ev[0]=2.0*nW/n; ev[n]=0;
    for (uint_t ii=1;ii<n;ii++){
        ev[ii]=sin(2.0*PI*nW/n*ii)/(PI*ii);
        ev[2*n-ii]=ev[ii];
    }

    //for (uint_t ii=0; ii<2*n; ii++)
    //    printf(" %23.16e \n",ev[ii]);

    //transform kernel to FFT
    KERN = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(n+1));
    p = fftw_plan_dft_r2c_1d(2*n, ev, KERN, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    //for (uint_t ii=0; ii<n+1; ii++)
    //    printf(" %23.16e %23.16e i\n",creal(KERN[ii]), cimag(KERN[ii]));

    EV=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(n+1));
    
    p =fftw_plan_dft_r2c_1d(2*n, ev, EV, FFTW_ESTIMATE);
    pi=fftw_plan_dft_c2r_1d(2*n, EV, ev, FFTW_ESTIMATE);
    
    for (uint_t ii=0; ii<K; ii++){
    //transform eigenvector
        for (uint_t jj=0; jj<n; jj++){
            ev[jj]=h[ii*n+jj];
            ev[jj+n]=0;
        }
    
        //transform, point-wise multiply, then inverse transform
        fftw_execute(p);
        tmath::pw_mult<fftw_complex,fftw_complex>(EV,KERN,n+1,EV);
        fftw_execute(pi);
        tmath::scale<double>(ev,1.0/2.0/n,n,ev);  //rescale ifft
        
        lambda[ii]=tmath::dot_mult<double,double>(ev,&h[ii*n],n);
        
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
    fftw_destroy_plan(pi);
    fftw_free(KERN); fftw_free(EV); fftw_free(ev);
    fftw_cleanup();
    
}
