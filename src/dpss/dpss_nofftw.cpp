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

*/  


/////////////////////////////////////////////////////////////////////
//  NOTE: THIS FILE REPLACE dpss_fftw.cpp IF YOU DO NOT WISH TO USE
//        THE FFTW3 PACKAGE TO COMPUTE ENERGY CONCENTRATIONS.
//         
//        FFTW3 is a great package, and you should consider using it.
/////////////////////////////////////////////////////////////////////

#include "dpss.h"

// approximates the concentration of energy in bandwidth [-W W] of the dpss sequences in h
// h is in row-major format with n columns and nseq rows
// energies placed in *lambda
// Computes eigenvalue, which corresponds to energy concentration
void compute_energy_concentrations(double *h, uint_t n, uint_t K, double nW, double *lambda){
    
    //kernel
    double *kern= new double[n];
    
    kern[0]=2.0*nW/n;
    for (uint_t ii=1;ii<n;ii++){
        kern[ii]=sin(2.0*PI*nW/n*ii)/(PI*ii);
    }
    
    double v;
    uint_t ki;
    for (uint_t ii=0; ii<K; ii++){
        
        lambda[ii]=0;
        for (uint_t jj=0; jj<n; jj++){
            v=0;
            for (uint_t kk=0; kk<n; kk++){
                
                //absolute difference between unsigned ints
                if (jj>kk) ki=jj-kk;
                else ki=kk-jj;
                
                v+=kern[ki]*h[n*ii+kk];
            }
            lambda[ii]+=v*h[n*ii+jj];
        }
    }
    
    //fixes eigenvalues
    if (lambda[K-1] > 1.0 )  lambda[K-1]=1.0;
    else if (lambda[K-1]<0 ) lambda[K-1]=0;
    //ensures strictly decreasing eigenvalues (biased upwards toward one)
    for (uint_t ii=K-1; ii>0; --ii){
        if (lambda[ii]<lambda[ii+1]) 
            lambda[ii]=lambda[ii+1];
        if (lambda[ii] > 1.0 )  
            lambda[ii]=1.0;
        else if (lambda[ii]<0 ) 
            lambda[ii]=0;
    }
    
    delete [] kern;
    
    return;
    
}