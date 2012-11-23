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

///////////////////////////////////////////////////////////////////////////////////
//  dpss.cpp contains all functions required for calculating the discrete prolate
//           spheroidal sequences, apart from the LAPACK fortran routines:
//                     dlamch_, dstebz_, dstevr_, dstein_
//           and the routine to compute the energy concentrations of the DPSSs:
//                     compute_energy_concentrations
//           which can be found in either dpss_fftw.cpp, or dpss_nofftw.cpp
//
//           An interface class is also provided: class dpss{}
//
//////////////////////////////////////////////////////////////////////////////////

#if !defined(UINT_T)
#define UINT_T
typedef unsigned int uint_t;   //unsigned integer
#endif

#include "simple_error.h"
#include "template_math.h"
#include <cmath>

#define DPSS_VERSION 1.1
#define PI 3.141592653589793238463
#define MAX_N 131072-1          // My LAPACK fails with size larger than 2^16-1 (size of INTEGER*2)
                                // (but we use a splitting technique, so can get away with n=2^17-1)
                                // Feel free to change this is your LAPACK supports larger indices

//constants
enum ARGS { SEQ_LENGTH, TIME_HALF_BANDWIDTH, SEQL, SEQU, INTERP_METHOD, INTERP_BASE };
enum INTERP_TYPE { NONE, LINEAR, SPLINE };

//keeps track of calculation information
struct dpss_workspace{
    uint_t n;                    //length of sequence
    double nW;                   //time half-bandwidth
    uint_t seql;                 //lower dpss
    uint_t sequ;                 //upper dpss
    uint_t K;                    //number of sequences
    INTERP_TYPE interp_method;   //interpolation method (NONE if none)
    uint_t interp_base;          //base length of interpolation
    bool energy;                 //true if energy concentrations are to be computed
    
    // Sets default values
    dpss_workspace(): n(0), nW(1), seql(0), sequ(0), K(1), interp_method(NONE), interp_base(0), energy(true) { }
    dpss_workspace(uint_t n_, double nW_): n(n_), nW(nW_), seql(0), sequ(floor(2*nW)-2), K(floor(2*nW)-1), interp_method(NONE), interp_base(n_), energy(true) { }
};

// Adjusts settings in the dpss_workspace so that no information conflicts
// and the sequences and be calculated without throwing errors
void fix_workspace(dpss_workspace& ws);

class dpss{
public:
    dpss(uint_t n, double nW);      //constructor
    dpss(dpss_workspace work);

    ~dpss();                        //destructor
    void energize();                //computes energy concentrations, sets info.energy==true
    uint_t length();                //return n*K
    uint_t size(int dim);           //dim < 1 returns n, otherwise returns K
    void compute();                 //computes sequences (and energies if info.energy==true)
    double lambda(uint_t k);        //returns the kth energy concentration, seql<= k <= sequ
    double operator()(uint_t k, uint_t i);   //returns the ith value of sequence seql<= k <= sequ
    dpss_workspace getinfo();       //return information structure
    const double* ph();             //returns const pointer to filter h
    const double* pl();             //returns const pointer to energies l
	
protected:
    double *h;                      //filter coefficients in row-major format
    double *l;                      //energy concentrations
    dpss_workspace info;            //workspace information
};

// MAIN DPSS calculation routine
// Uses symmetric tridiagonal matrix formulation with even-odd splitting
// Fills h in row-major format with length n dpss sequences seql to sequ (indices starting at zero)
// and time-bandwidth product nW.
// NOTE: Does not check size of n.  If you would like to have automatic checking and interpolating,
//       use the dpss class.
void dpss_calc(uint_t n, double nW, int seql, int sequ, double *h);
// takes a function of the type eig_rrr eig_iit
void dpss_calc(uint_t n, double nW, int seql, int sequ, void (*eig_calc)(uint_t,double*,double*,uint_t,uint_t,double*, double*, uint_t), double *h);

//External function provided by either dpss_fftw.cpp or dpss_nofftw.cpp
extern void compute_energy_concentrations(double *h, uint_t n, uint_t K, double nW, double *l);

//Polarizes the sequences so even sequences have positive mean, and odd sequences have positive sum
// (n-2*i)*h[odd][i]
void polarize_dpss(double *h, uint_t n, uint_t kk);

//interpolation
void normalize_vec( double *h, uint_t n);                            //normalizes h
void linear_interp( double *y, uint_t n, uint_t nout, double *z);    //linearly interpolates y from n points to nout points, filling z
void spline_interp( double *y, uint_t n, uint_t nout, double *z);    //interpolates y using cubic splines

//calculate eigenvalues (interface to LAPACK routines
void eig_rrr(uint_t n, double *D, double *E, uint_t il, uint_t iu, double *eig_val, double *eig_vec, uint_t vec_length);
void eig_rrr(uint_t n, double *D, double *E, uint_t il, uint_t iu, double *eig_val, double *eig_vec);
void eig_iit(uint_t n, double *D, double *E, uint_t il, uint_t iu, double *eig_val, double *eig_vec, uint_t vec_length);
void eig_iit(uint_t n, double *D, double *E, uint_t il, uint_t iu, double *eig_val, double *eig_vec);

//External LAPACK functions
// Technically, dpss_calc only uses the dstebz_ and dstein_ routines, however, with one small change in code, 
// one can change to the dstevr.  dstebz/dstein uses bisection and inverse iteration, which takes longer, but is slightly 
// accurate.  dstevr uses relatively robust representations, which is faster, but slightly less accurate and 
// uses more memory.
extern "C" void dstevr_(char *jobz, char *range, uint_t *n, double *D, double *E, double *vl, double *vu, uint_t *il, uint_t *iu, double *abstol, uint_t *m, double *W, double *Z, uint_t *ldz, uint_t *ISUPPZ, double *WORK, uint_t *lwork, uint_t *IWORK, uint_t *liwork, uint_t *info );
extern "C" double dlamch_(char *type);
extern "C" void dstebz_(char *range, char *order, uint_t *n, double *vl, double *vu, uint_t *il, uint_t *iu, double *abstol, double *D, double *E, uint_t *m, uint_t *nsplit, double *W, uint_t *IBLOCK, uint_t *ISPLIT, double *WORK, uint_t *IWORK, uint_t *info );
extern "C" void dstein_(uint_t *n, double *D, double *E, uint_t *m, double *W, uint_t *IBLOCK, uint_t *ISPLIT, double *Z, uint_t *ldz, double *WORK, uint_t *IWORK, uint_t *IFAIL, uint_t *info );
