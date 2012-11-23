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

/////////////////////////////////////////////////////////////////////////
// Basic math / vector math operations, including:
//       abs:       absolute value (real types only)
//       max:       maximum        (real types only)
//       min:       minimum        (real types only)
//       sum:       vector
//       pw_mult:   vector pointwise multiplication
//       dot_mult:  vector dot-product
//       scale:     vector scale by constant
//       mean:      vector mean
//       wmean:     vector weighted mean
//       var:       vector variance
//       mom2:      vector 2nd moment (non-centered)
////////////////////////////////////////////////////////////////////////

#ifndef MTPSD_TEMPLATE_MATH
#define MTPSD_TEMPLATE_MATH

#include <complex.h>

namespace tmath {

//absolute value function
template <class T> T abs ( T a ) { return (0<a)?a:-a; }
template <class T> T max ( T a, T b ) { return (b<a)?a:b; }
template <class T> T min ( T a, T b ) { return (b<a)?b:a; }

// basic sum algorithm
template <class T> T sum(const T *a, unsigned int n){
    T mysum=0;
    for (unsigned int ii=0; ii<n; ii++)
        mysum+=a[ii];
    
    return mysum;
}

// pointwise multiplication
template <class A, class B> void pw_mult(const A *a, const B *b, unsigned int n, A *ab ){
    
    for (unsigned int ii=0; ii<n; ii++)
        ab[ii]=a[ii]*b[ii];
    
}

//dot product
template <class A, class B> A dot_mult(const A *a, const B *b, unsigned int n){
    A dotsum=0;
    for (unsigned int ii=0; ii<n; ii++)
        dotsum+=a[ii]*b[ii];
    
    return dotsum;    
}

//scale by a double value
template <class T> void scale(const T *a, double b, unsigned int n, T *ab){
    for (unsigned int ii=0; ii<n; ii++)
        ab[ii]=a[ii]*b;    
}

//standard mean algorithm
template <class T> T mean(const T *a, unsigned int n){  
    return sum<T>(a,n)/n; 
} 
    
//weighted mean algorithm
template <class A, class B> A wmean(const A *a, const B *wts, unsigned int n){  
    return (A)(dot_mult<A,B>(a,wts,n) / sum<B>(wts,n)); 
}

template <class T> double var(const T *a, unsigned int n){
    double varsum=0;
    T m=mean<T>(a,n);
    
    double tmp;
    for (unsigned int ii=0; ii<n; ii++){
        tmp=cabs(a[ii]-m);
        varsum+= tmp*tmp;
    }
        
    return varsum/n; 
} 

// 2nd non-centered moment
template <class T> double mom2(const T *a, unsigned int n){
    
    double varsum=0;
    
    double tmp;
    for (unsigned int ii=0; ii<n; ii++){
        tmp=cabs(a[ii]);
        varsum+= tmp*tmp;
    }
    
    return varsum/n; 
} 

}

#else
#warning template_math already loaded
#endif
